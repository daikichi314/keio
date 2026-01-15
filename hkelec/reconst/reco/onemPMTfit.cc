/**
 * @file onemPMTfit.cc
 * @brief 光源位置再構成クラスの実装ファイル
 *
 * TMinuitのFCN関数を用いて、観測データとモデルの残差(Chi2)を最小化します。
 *
 * [ロジックの概要]
 * 1. フィッティングパラメータ: 光源位置(x,y,z)、発光時刻(t0)、光量パラメータ(A,B)
 * 2. フィッティング関数(FCN):
 *   - 各PMTヒットに対して、期待される到達時間と電荷を計算
 *  - 計算には、選択された電荷モデル(FuncF/FuncG)とChi2計算モデルを使用
 *  - 観測値との残差を基にChi2を計算
 * 3. Minuitを用いてChi2を最小化し、最適なパラメータを探索
 * 4. フィッティング結果をFitResult構造体に格納して返す
 * 
 * [主な変更点 2025-01-08 v3]
 * 1. 電荷モデル FuncF / FuncG の分離を完全化。
 * それぞれのモデルで異なる距離依存係数(c0)と角度依存係数(epsilon params)を使用します。
 * 2. 時間フィットにおいても、選択した電荷モデルに対応する PMT半径(r_pmt) を使用。
 * これにより、中心位置(Z = 80.5 - r_pmt)がモデルごとに変化します。
 * 3. パラメータ B は常に 0 に固定されます。
 *
 * @author Gemini (Modified based on user request)
 */

#include "onemPMTfit.hh"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <regex>
#include <TMath.h> 

// グローバルポインタ
LightSourceFitter* gFitter = nullptr;

// =========================================================
// ファイル名パース関数の実装
// =========================================================
FilenameParams ParseFilename(const std::string& filename) {
    FilenameParams params;
    params.valid = false;
    
    // 正規表現で数値を抽出 (負の数、小数に対応)
    // パターン: x(数値)_y(数値)_z(数値)-(ラン番号)-(数値)dB
    std::regex pattern(R"(x([\d\.-]+)_y([\d\.-]+)_z([\d\.-]+)-(\d+)-([\d\.]+)dB)");
    std::smatch match;
    
    if (std::regex_search(filename, match, pattern)) {
        try {
            params.x = std::stod(match[1].str());
            params.y = std::stod(match[2].str());
            params.z = std::stod(match[3].str());
            params.db = std::stod(match[5].str());
            params.valid = true;
        } catch (const std::exception& e) {
            params.valid = false;
        }
    }
    
    return params;
}

// =========================================================
// パラメータ計算用関数
// f(q) = c0 * q^{-1/2} + c1 + c2 * q + c3 * q^2
// TW_PARAMS用の場合、上限値チェックが適用されます
// =========================================================
double CalcParametricValue(int ch, double charge, const double params[4][4]) {
    if (ch < 0 || ch >= 4) return 1.0; 
    
    double q = (charge > 1e-3) ? charge : 1e-3;

    double c0 = params[ch][0];
    double c1 = params[ch][1];
    double c2 = params[ch][2];
    double c3 = params[ch][3];

    double val = c0 / std::sqrt(q) + c1 + c2 * q + c3 * q * q;
    
    // TW_PARAMSの場合のみ上限値チェックを適用
    // (paramsのアドレスがTW_PARAMSと一致するか確認)
    if (params == TW_PARAMS) {
        if (val > TW_MAX_VALUES[ch]) {
            val = TW_MAX_VALUES[ch];
        }
    }
    
    return val;
}

double GetEMG_Sigma(int ch, double charge) { return 1.0; }
double GetEMG_Tau(int ch, double charge)   { return 1.0; }

// =========================================================
// ヘルパー関数: ベクトルのなす角(cos)を計算
// =========================================================
double CalculateCosAngle(double v1x, double v1y, double v1z, double v2x, double v2y, double v2z) {
    double dot = v1x * v2x + v1y * v2y + v1z * v2z;
    double mag1 = std::sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
    double mag2 = std::sqrt(v2x * v2x + v2y * v2y + v2z * v2z);
    
    if (mag1 == 0 || mag2 == 0) return -1.0; 
    
    double cos_a = dot / (mag1 * mag2);
    if (cos_a > 1.0) cos_a = 1.0;
    if (cos_a < -1.0) cos_a = -1.0;
    
    return cos_a;
}

// =========================================================
// 評価関数: EMGの負の対数尤度 (-2lnL)
// =========================================================
double CalcEMG_NLL(double t, double mu, double sigma, double tau) {
    if (tau <= 0 || sigma <= 0) return 1e9;
    double arg_erfc = (sigma/tau - (t - mu)/sigma) / std::sqrt(2.0);
    double term_exp = (sigma*sigma)/(2.0*tau*tau) - (t - mu)/tau;
    double val_erfc = std::erfc(arg_erfc);
    if (val_erfc <= 1e-15) val_erfc = 1e-15;
    double ln_f = -std::log(2.0 * tau) + term_exp + std::log(val_erfc);
    return -2.0 * ln_f;
}

// =========================================================
// Minuit用 目的関数 (Chi2計算)
// =========================================================
void fcn_wrapper(int& npar, double* gin, double& f, double* par, int iflag) {
    // フィッティングパラメータ
    double x = par[0];
    double y = par[1];
    double z = par[2]; // 光源位置
    double t0 = par[3]; // 発光時刻
    double A = par[4];  // 光量パラメータ
    // B = par[5] は使用せず (0固定)

    const auto& hits = gFitter->GetData();
    const auto& config = gFitter->GetConfig();

    // -------------------------------------------------------------
    // 1. モデルに基づく定数・パラメータ配列の選択
    // -------------------------------------------------------------
    double r_pmt_eff = 0.0;
    const double (*charge_c0_params)[4] = nullptr;   // 距離係数配列へのポインタ (便宜上配列として扱う)
    const double (*ang_params)[8] = nullptr; // 角度係数配列へのポインタ

    // モデルに応じた選択
    if (config.chargeModel == ChargeModelType::FuncG) {
        // --- FuncG ---
        r_pmt_eff = PMT_RADIUS_G; // 23.5 cm
        // キャストしてポインタを合わせる
        // 本来は構造体などでまとめるべきですが、既存コードの延長でポインタ操作で対応
        // CHARGE_RADIAL_PARAMS_FUNC_G は 1次元配列ですが、アドレスを渡してアクセスします
        // (ここでは個別アクセスするので単に配列を指すようにします)
        
        // 角度パラメータ
        ang_params = CHARGE_ANGULAR_PARAMS_FUNC_G;
        
    } else {
        // --- FuncF (Default) ---
        r_pmt_eff = PMT_RADIUS_F; // 28.5 cm
        ang_params = CHARGE_ANGULAR_PARAMS_FUNC_F;
    }

    // PMT球の中心Z座標 = 表面(80.5) - 半径
    double pmt_center_z = PMT_SURFACE_Z - r_pmt_eff;

    double chi2_total = 0.0;

    // -------------------------------------------------------------------
    // 2. 電荷 (Charge) に関する Chi2
    // -------------------------------------------------------------------
    if (config.chargeType != ChargeChi2Type::None) {
        for (const auto& hit : hits) {
            // ★追加: チャンネル番号の安全確認
            if (hit.ch < 0 || hit.ch >= 4) continue;
            
            // PMT球中心座標
            double pmt_cx = PMT_XY_POS[hit.ch][0];
            double pmt_cy = PMT_XY_POS[hit.ch][1];
            double pmt_cz = pmt_center_z;

            // 光源からPMT中心へのベクトル
            double vec_x = pmt_cx - x;
            double vec_y = pmt_cy - y;
            double vec_z = pmt_cz - z;
            double dist_center2 = vec_x*vec_x + vec_y*vec_y + vec_z*vec_z;
            double dist_center = std::sqrt(dist_center2);

            // 角度計算 cos(alpha)
            double cos_alpha = CalculateCosAngle(vec_x, vec_y, vec_z, 
                                                 hit.dir_x, hit.dir_y, hit.dir_z);
            
            // 角度依存項 epsilon (モデルごとに選択された配列を使用)
            const double* c_ang = ang_params[hit.ch];
            double epsilon = c_ang[0] + cos_alpha * (c_ang[1] + cos_alpha * (c_ang[2] + cos_alpha * (c_ang[3] + 
                             cos_alpha * (c_ang[4] + cos_alpha * (c_ang[5] + cos_alpha * (c_ang[6] + cos_alpha * c_ang[7]))))));
                            
            if (epsilon < 0) epsilon = 0.0;

            double mu = 0.0;
            double f_r = 0.0;

            // --- モデル式計算 ---
            if (config.chargeModel == ChargeModelType::FuncF) {
                // [FuncF] r_pmt=28.5
                // f(r) = c0 * (1 - sqrt(1 - (r_pmt/r)^2))
                double c0 = CHARGE_RADIAL_PARAMS_FUNC_F[hit.ch];

                // ルート内保護
                if (dist_center > r_pmt_eff + 0.001) {
                    double ratio = r_pmt_eff / dist_center;
                    f_r = c0 * (1.0 - std::sqrt(1.0 - ratio * ratio));
                } else {
                    f_r = c0; 
                }
                
                mu = A * f_r * epsilon;

            } else if (config.chargeModel == ChargeModelType::FuncG) {
                // [FuncG] r_pmt=23.5
                // f(r) = c0 / r^2
                double c0 = CHARGE_RADIAL_PARAMS_FUNC_G[hit.ch];
                
                // ゼロ除算防止
                if (dist_center2 < 1.0) dist_center2 = 1.0;
                
                f_r = c0 / dist_center2;
                mu = A * f_r * epsilon;
            }
            
            if (mu < 1e-9) mu = 1e-9; 

            // --- Chi2 加算 ---
            double n = hit.charge;
            if (config.chargeType == ChargeChi2Type::BakerCousins) {
                double term = 0.0;
                if (n > 1e-9) term = mu - n + n * std::log(n / mu);
                else term = mu; 
                chi2_total += 2.0 * term;
            } else {
                double sigma_q = 1.0; 
                chi2_total += std::pow(n - mu, 2) / (sigma_q * sigma_q);
            }
        }
    }

    // -------------------------------------------------------------------
    // 3. 時間 (Time) に関する Chi2
    // -------------------------------------------------------------------
    if (config.timeType != TimeChi2Type::None) {
        double goodness_sum = 0.0;
        
        for (const auto& hit : hits) {
            if (!hit.isHit) continue; 

            // ★追加: チャンネル番号の安全確認
            if (hit.ch < 0 || hit.ch >= 4) continue;

            // PMT球中心座標 (共通設定された pmt_center_z を使用)
            double pmt_cx = PMT_XY_POS[hit.ch][0];
            double pmt_cy = PMT_XY_POS[hit.ch][1];
            double pmt_cz = pmt_center_z;

            // 光源(x,y,z) と PMT球中心間の距離
            double dx = x - pmt_cx;
            double dy = y - pmt_cy;
            double dz = z - pmt_cz;
            double dist_center = std::sqrt(dx*dx + dy*dy + dz*dz);
            
            // 時間フィット用の飛行距離 = 中心距離 - 半径
            double dist_surface = dist_center - r_pmt_eff;
            
            // 物理的にあり得ない近距離の保護
            if (dist_surface < 0.1) dist_surface = 0.1;

            // 飛行時間
            double t_flight = dist_surface / C_LIGHT;
            
            // 期待時刻
            double tw_val = CalcParametricValue(hit.ch, hit.charge, TW_PARAMS);
            double t_corr_val = TIME_CORRECTION_VAL[hit.ch];
            
            double t_expected = t0 + t_flight + tw_val + t_corr_val;
            double t_obs = hit.time;

            // 時間分解能(Sigma)
            double sigma_t = CalcParametricValue(hit.ch, hit.charge, SIGMA_T_PARAMS);
            if (sigma_t < 0.1) sigma_t = 0.1; 

            // Chi2加算
            if (config.timeType == TimeChi2Type::Goodness) {
                double res2 = std::pow(t_obs - t_expected, 2);
                goodness_sum += std::exp( -res2 / (2.0 * sigma_t * sigma_t) );

            } else if (config.timeType == TimeChi2Type::EMG) {
                double sigma = sigma_t; 
                double tau = 1.0; 
                chi2_total += CalcEMG_NLL(t_obs, t_expected, sigma, tau);

            } else {
                // Gaussian
                chi2_total += std::pow(t_obs - t_expected, 2) / (sigma_t * sigma_t);
            }
        }

        if (config.timeType == TimeChi2Type::Goodness) {
            if (goodness_sum < 1e-9) goodness_sum = 1e-9;
            chi2_total += -2.0 * std::log(goodness_sum);
        }
    }

    f = chi2_total;
}

// LightSourceFitterクラスの実装
LightSourceFitter::LightSourceFitter() {
    fMinuit = new TMinuit(6);
    fMinuit->SetFCN(fcn_wrapper);
    fMinuit->SetPrintLevel(-1);  // 全ての出力を抑制（エラーは別途検出）
    gFitter = this;
}

LightSourceFitter::~LightSourceFitter() {
    delete fMinuit;
}

void LightSourceFitter::SetConfig(const FitConfig& config) {
    fConfig = config;
}

void LightSourceFitter::SetDataFilename(const std::string& filename) {
    fDataFilename = filename;
}

bool LightSourceFitter::FitEvent(const std::vector<PMTData>& eventHits, FitResult& res) {
    fCurrentHits = eventHits;
    InitializeParameters(eventHits);

    // 最小化実行
    double arglist[10];
    int ierflg = 0;
    arglist[0] = 100000;
    arglist[1] = 0.1;
    fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    // 結果取得
    double val, err;
    fMinuit->GetParameter(0, val, err); res.x = val; res.err_x = err;
    fMinuit->GetParameter(1, val, err); res.y = val; res.err_y = err;
    fMinuit->GetParameter(2, val, err); res.z = val; res.err_z = err;
    fMinuit->GetParameter(3, val, err); res.t = val; res.err_t = err;
    fMinuit->GetParameter(4, val, err); res.A = val;
    fMinuit->GetParameter(5, val, err); res.B = val;

    double fmin, fedm, errdef;
    int npari, nparx, istat;
    fMinuit->mnstat(fmin, fedm, errdef, npari, nparx, istat);
    res.chi2 = fmin;
    
    // NDF計算
    int nDataPoints = 0;
    if (fConfig.chargeType != ChargeChi2Type::None) nDataPoints += eventHits.size();
    if (fConfig.timeType != TimeChi2Type::None) {
        for(const auto& h : eventHits) if(h.isHit) nDataPoints++;
    }
    int nFreeParams = fMinuit->GetNumFreePars();
    res.ndf = nDataPoints - nFreeParams;
    res.status = istat;

    // // 収束失敗時にエラーメッセージを出力
    // if (istat != 3) {
    //     std::cerr << "WARNING: Fit did not converge for event. Status=" << istat 
    //               << " Chi2=" << fmin << " EDM=" << fedm << std::endl;
    // }

    // 変更前: 完全収束(3)のみ許可
    return (istat == 3);

    // 変更案1: 多少の近似や補正があっても、最小値が見つかっていればOKとする (推奨)
    // "Not pos-def" (2) や "Approximate" (1) も合格にする
    // return (istat >= 1);
}

void LightSourceFitter::InitializeParameters(const std::vector<PMTData>& hits) {
    // ファイル名から初期値を取得（設定されている場合）
    FilenameParams fnParams = ParseFilename(fDataFilename);
    
    double iniX, iniY, iniZ, iniA;
    
    if (fnParams.valid) {
        // ファイル名から取得した値を使用
        iniX = fnParams.x;
        iniY = fnParams.y;
        iniZ = fnParams.z;
        // A初期値: 15dBを基準 (15dBの時A=1.0)
        iniA = std::pow(10.0, (15.0 - fnParams.db) / 10.0);
    } else {
        // ファイル名がパースできない場合は重心計算にフォールバック
        double sumQ = 0, sumX = 0, sumY = 0, sumZ = 0;
        for (const auto& hit : hits) {
            if(hit.isHit) {
                sumQ += hit.charge;
                sumX += hit.x * hit.charge;
                sumY += hit.y * hit.charge;
                sumZ += hit.z * hit.charge;
            }
        }
        iniX = (sumQ > 0) ? sumX/sumQ : 0;
        iniY = (sumQ > 0) ? sumY/sumQ : 0;
        iniZ = 100.0;
        iniA = 1.0;
    }

    fMinuit->DefineParameter(0, "x", iniX, 1.0, -200, 200);
    fMinuit->DefineParameter(1, "y", iniY, 1.0, -200, 200);
    fMinuit->DefineParameter(2, "z", iniZ, 1.0, 0, 300);
    
    // 時間パラメータ
    if (fConfig.timeType == TimeChi2Type::None) {
        fMinuit->DefineParameter(3, "t", 0.0, 0.0, 0.0, 0.0);
        fMinuit->FixParameter(3);
    } else {
        fMinuit->DefineParameter(3, "t", 0.0, 0.1, -300, 300);
        fMinuit->Release(3);
    }

    // 電荷パラメータ
    if (fConfig.chargeType == ChargeChi2Type::None) {
        fMinuit->DefineParameter(4, "A", 0.0, 0.0, 0.0, 0.0);
        fMinuit->FixParameter(4);
        fMinuit->DefineParameter(5, "B", 0.0, 0.0, 0.0, 0.0);
        fMinuit->FixParameter(5);
    } else {
        fMinuit->DefineParameter(4, "A", iniA, 0.05, 0, 20);
        fMinuit->Release(4);
        
        // Bの固定: 今回のモデル(FuncF, FuncG)では常に0に固定
        fMinuit->DefineParameter(5, "B", 0.0, 0.0, 0.0, 0.0);
        fMinuit->FixParameter(5);
    }
}
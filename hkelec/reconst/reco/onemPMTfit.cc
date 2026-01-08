/**
 * @file onemPMTfit.cc
 * @brief 光源位置再構成クラスの実装ファイル
 *
 * TMinuitのFCN関数を用いて、観測データとモデルの残差(Chi2)を最小化します。
 *
 * [主な変更点 2025-01-08]
 * 1. 時間期待値計算に TimeWalk (TW) と TimeCorrection の項を追加・復帰させました。
 * 2. TWと時間分解能(Sigma)を電荷の関数として実装しました。
 * f(q) = c0*q^-0.5 + c1 + c2*q + c3*q^2
 * 3. 距離計算(dist)を、PMT半球中心からの距離から半径を引く方式に変更しました。
 * 4. 電荷モデルのR^2を分母に配置する修正を維持しています。
 *
 * @author Gemini (Modified based on user request)
 */

#include "onemPMTfit.hh"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <TMath.h> 

// グローバルポインタ
LightSourceFitter* gFitter = nullptr;

// =========================================================
// パラメータ計算用関数
// f(q) = c0 * q^{-1/2} + c1 + c2 * q + c3 * q^2
// =========================================================
double CalcParametricValue(int ch, double charge, const double params[4][4]) {
    if (ch < 0 || ch >= 4) return 1.0; // エラーガード
    
    // 電荷が0以下の場合の保護（平方根計算のため）
    // Unhitの場合はそもそも呼ばれない想定だが、念のため小さな正の数にする
    double q = (charge > 1e-3) ? charge : 1e-3;

    double c0 = params[ch][0];
    double c1 = params[ch][1];
    double c2 = params[ch][2];
    double c3 = params[ch][3];

    double val = c0 / std::sqrt(q) + c1 + c2 * q + c3 * q * q;
    return val;
}

// EMG用 (今回は使用しないが互換性のため残す)
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
    double B = par[5];  // バックグラウンド等

    const auto& hits = gFitter->GetData();
    const auto& config = gFitter->GetConfig();

    double chi2_total = 0.0;

    // -------------------------------------------------------------------
    // 1. 電荷 (Charge) に関する Chi2
    // -------------------------------------------------------------------
    if (config.chargeType != ChargeChi2Type::None) {
        for (const auto& hit : hits) {
            // 電荷モデルの距離計算は、従来通りPMT表面中心との距離を使用するか、
            // あるいは半球中心距離を使うか議論の余地がありますが、
            // 光量は「受光面への距離」に依存するため、ここではひとまず
            // ユーザー指定の半球モデル距離（最近接距離）を一貫して使用することにします。
            
            // PMT半球中心座標
            double pmt_cx = PMT_XY_POS[hit.ch][0];
            double pmt_cy = PMT_XY_POS[hit.ch][1];
            double pmt_cz = PMT_SPHERE_Z;

            double dx = x - pmt_cx;
            double dy = y - pmt_cy;
            double dz = z - pmt_cz;
            double dist_center = std::sqrt(dx*dx + dy*dy + dz*dz);
            
            // 半球表面までの距離 (光源から一番近い点)
            double dist_surface = dist_center - PMT_RADIUS;
            
            // 光源がPMT内部に入り込んだ場合の保護
            if (dist_surface < 1.0) dist_surface = 1.0; 

            double mu = 0.0;
            
            // 期待値計算: mu = (A * f(cos)) / R^2 + B
            if (config.chargeModel == ChargeModelType::Cosine) {
                // 入射角 cos(alpha) の計算
                // 光源からPMT中心に向かうベクトル vs PMTの向き(0,0,1)
                // ※ここでは簡易的に中心方向ベクトルを使用します
                // 厳密には「受光点での法線」など考慮が必要ですが、
                // 与えられた条件「光源とPMT中心の距離...」に合わせています。
                // ベクトル: 光源(x,y,z) -> PMT(hit.x, hit.y, hit.z) 
                // hit.x等は readData で表面中心が入っていますが、
                // 角度計算用ベクトルとしては (pmt_cx - x, ...) を用います
                
                double vec_s2p_x = pmt_cx - x;
                double vec_s2p_y = pmt_cy - y;
                double vec_s2p_z = pmt_cz - z;
                
                double cos_alpha = CalculateCosAngle(vec_s2p_x, vec_s2p_y, vec_s2p_z, 
                                                     hit.dir_x, hit.dir_y, hit.dir_z); // hit.dirは(0,0,1)
                
                // 角度依存係数 (以前のシグモイド関数)
                double ang_factor = 1.0 / (1.0 + std::exp(-6.0 * (cos_alpha - 1.0)));
                
                mu = (A * ang_factor) / (dist_surface * dist_surface) + B;

            } else {
                // 等方モデル
                mu = A / (dist_surface * dist_surface) + B;
            }
            
            if (mu < 1e-9) mu = 1e-9; 

            double n = hit.charge;

            if (config.chargeType == ChargeChi2Type::BakerCousins) {
                double term = 0.0;
                if (n > 1e-9) {
                    term = mu - n + n * std::log(n / mu);
                } else {
                    term = mu; 
                }
                chi2_total += 2.0 * term;
            } else {
                // Gaussian: sigmaは簡易的に1.0またはsqrt(n)など
                double sigma_q = 1.0; 
                chi2_total += std::pow(n - mu, 2) / (sigma_q * sigma_q);
            }
        }
    }

    // -------------------------------------------------------------------
    // 2. 時間 (Time) に関する Chi2
    // -------------------------------------------------------------------
    if (config.timeType != TimeChi2Type::None) {
        double goodness_sum = 0.0;
        int hit_count_time = 0;
        
        for (const auto& hit : hits) {
            if (!hit.isHit) continue; 

            // --- [変更点3] 距離計算の変更 ---
            // PMT半球中心座標 (Z=48.0)
            double pmt_cx = PMT_XY_POS[hit.ch][0];
            double pmt_cy = PMT_XY_POS[hit.ch][1];
            double pmt_cz = PMT_SPHERE_Z;

            // 光源(x,y,z) と PMT半球中心間の距離
            double dx = x - pmt_cx;
            double dy = y - pmt_cy;
            double dz = z - pmt_cz;
            double dist_center = std::sqrt(dx*dx + dy*dy + dz*dz);
            
            // 半径を引いて表面までの距離とする
            double dist_surface = dist_center - PMT_RADIUS;
            // 負の値にならないよう保護（光源が近すぎる場合）
            if (dist_surface < 0.1) dist_surface = 0.1;

            // 飛行時間
            double t_flight = dist_surface / C_LIGHT;
            
            // --- [変更点1] 期待時刻の計算 (TimeWalk, Correction追加) ---
            // t_expected = t0 + t_flight + TW(q) + Correction
            
            // TimeWalk補正値の計算
            double tw_val = CalcParametricValue(hit.ch, hit.charge, TW_PARAMS);
            
            // チャンネルごとの固定補正値 (ケーブル遅延など)
            double t_corr_val = TIME_CORRECTION_VAL[hit.ch];
            
            // 期待値
            double t_expected = t0 + t_flight + tw_val + t_corr_val;
            double t_obs = hit.time;

            // --- [変更点2] 時間分解能(Sigma)の計算 ---
            double sigma_t = CalcParametricValue(hit.ch, hit.charge, SIGMA_T_PARAMS);
            if (sigma_t < 0.1) sigma_t = 0.1; // 保護

            // Chi2加算
            if (config.timeType == TimeChi2Type::Goodness) {
                // Goodness (APFit like): sum( exp( -res^2 / 2sigma^2 ) )
                double res2 = std::pow(t_obs - t_expected, 2);
                goodness_sum += std::exp( -res2 / (2.0 * sigma_t * sigma_t) );
                hit_count_time++;

            } else if (config.timeType == TimeChi2Type::EMG) {
                // EMGは現在使わないとのことだが、ロジックは残す
                double sigma = sigma_t; // EMGのGaussian成分にsigma_tを流用
                double tau = 1.0;       // 仮
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

// ... (クラスの実装部分やInitializeParametersは前回と同様、変更なし) ...
// LightSourceFitterクラスの実装
LightSourceFitter::LightSourceFitter() {
    fMinuit = new TMinuit(6);
    fMinuit->SetFCN(fcn_wrapper);
    gFitter = this;
}

LightSourceFitter::~LightSourceFitter() {
    delete fMinuit;
}

void LightSourceFitter::SetConfig(const FitConfig& config) {
    fConfig = config;
}

bool LightSourceFitter::FitEvent(const std::vector<PMTData>& eventHits, FitResult& res) {
    fCurrentHits = eventHits;
    InitializeParameters(eventHits);

    // 最小化実行
    double arglist[10];
    int ierflg = 0;
    arglist[0] = 5000;
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

    return (istat == 3);
}

void LightSourceFitter::InitializeParameters(const std::vector<PMTData>& hits) {
    // 重心計算
    double sumQ = 0, sumX = 0, sumY = 0, sumZ = 0;
    for (const auto& hit : hits) {
        if(hit.isHit) {
            sumQ += hit.charge;
            sumX += hit.x * hit.charge;
            sumY += hit.y * hit.charge;
            sumZ += hit.z * hit.charge;
        }
    }
    double iniX = (sumQ > 0) ? sumX/sumQ : 0;
    double iniY = (sumQ > 0) ? sumY/sumQ : 0;
    double iniZ = 50.0;

    fMinuit->DefineParameter(0, "x", iniX, 1.0, -200, 200);
    fMinuit->DefineParameter(1, "y", iniY, 1.0, -200, 200);
    fMinuit->DefineParameter(2, "z", iniZ, 1.0, 0, 200);
    
    if (fConfig.timeType == TimeChi2Type::None) {
        fMinuit->DefineParameter(3, "t", 0.0, 0.0, 0.0, 0.0);
        fMinuit->FixParameter(3);
    } else {
        fMinuit->DefineParameter(3, "t", 0.0, 1.0, -100, 100);
        fMinuit->Release(3);
    }

    if (fConfig.chargeType == ChargeChi2Type::None) {
        fMinuit->DefineParameter(4, "A", 0.0, 0.0, 0.0, 0.0);
        fMinuit->FixParameter(4);
        fMinuit->DefineParameter(5, "B", 0.0, 0.0, 0.0, 0.0);
        fMinuit->FixParameter(5);
    } else {
        fMinuit->DefineParameter(4, "A", 10000.0, 100.0, 0, 1000000);
        fMinuit->Release(4);
        
        if (fConfig.chargeModel == ChargeModelType::ZeroIntercept) {
            fMinuit->DefineParameter(5, "B", 0.0, 0.0, 0.0, 0.0);
            fMinuit->FixParameter(5);
        } else {
            fMinuit->DefineParameter(5, "B", 0.0, 0.1, 0, 1000);
            fMinuit->Release(5);
        }
    }
}
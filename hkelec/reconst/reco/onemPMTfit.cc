/**
 * @file onemPMTfit.cc
 * @brief 光源位置再構成クラスの実装ファイル
 *
 * TMinuitのFCN関数を用いて、観測データとモデルの残差(Chi2)を最小化します。
 * * [重要修正]
 * - 電荷モデルに角度依存性 (Cosine dependence) を追加しました。
 * - 期待値計算において R^2 を正しく分母に配置しました。
 *
 * @author Gemini (Modified based on user request)
 * @date 2025-12-14
 */

#include "onemPMTfit.hh"
#include <iostream>
#include <cmath>
#include <algorithm> // std::max
#include <TMath.h> 

// グローバルポインタ
LightSourceFitter* gFitter = nullptr;

// =========================================================
// パラメータ取得関数 (現在は仮定数)
// =========================================================
double GetEMG_Sigma(int ch, double charge) { return 1.0; }
double GetEMG_Tau(int ch, double charge)   { return 1.0; }
double GetTimeResolution(int ch, double charge) { return 1.0; } 

// =========================================================
// ヘルパー関数: ベクトルのなす角(cos)を計算
// =========================================================
double CalculateCosAngle(double v1x, double v1y, double v1z, double v2x, double v2y, double v2z) {
    double dot = v1x * v2x + v1y * v2y + v1z * v2z;
    double mag1 = std::sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
    double mag2 = std::sqrt(v2x * v2x + v2y * v2y + v2z * v2z);
    
    if (mag1 == 0 || mag2 == 0) return -1.0; // Error case
    
    double cos_a = dot / (mag1 * mag2);
    if (cos_a > 1.0) cos_a = 1.0;
    if (cos_a < -1.0) cos_a = -1.0;
    
    return cos_a;
}

// =========================================================
// 評価関数: EMG
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
    double x = par[0];
    double y = par[1];
    double z = par[2];
    double t0 = par[3];
    double A = par[4];
    double B = par[5];

    const auto& hits = gFitter->GetData();
    const auto& config = gFitter->GetConfig();

    double chi2_total = 0.0;

    // -----------------------------------------------------
    // 1. 電荷(Charge)に関するChi2
    // -----------------------------------------------------
    if (config.chargeType != ChargeChi2Type::None) {
        for (const auto& hit : hits) {
            // 光源からPMTへのベクトル
            double dx = hit.x - x;
            double dy = hit.y - y;
            double dz = hit.z - z;
            double dist2 = dx*dx + dy*dy + dz*dz;
            double dist = std::sqrt(dist2);
            
            double mu = 0.0;
            
            // 距離が近すぎる場合の保護
            if (dist > 1e-3) {
                if (config.chargeModel == ChargeModelType::Cosine) {
                    // --- 新モデル (Cosine Dependent) ---
                    // 光の進行方向ベクトル (Source -> PMT)
                    // PMTの向きベクトル (hit.dir_x, ...)
                    // cos_alpha = (Light_vec . PMT_dir)
                    
                    // 光源からPMTに向かうベクトル (dx, dy, dz)
                    // PMTの受光面は通常、PMTの向きベクトルと対向する方向からの光を受けるときに最大感度
                    // PMTの向き(dir)が「受光面の法線ベクトル（外向き）」だとすると、
                    // 光のベクトル(Source->PMT)とのなす角は 180度(cos=-1) で正面入射となる定義の場合と、
                    // 単純にベクトルのなす角で計算する場合がある。
                    // 頂いたコードでは dx,dy,dz をそのまま使い、シグモイドの引数が (cos - 1) となっていたため、
                    // cos=1 (平行) のときに exp(0)=1 -> 1/(1+1)=0.5 となり、cosが小さくなると exp(正) で分母大 -> 値小
                    // つまり、ベクトル同士が「同じ向き」のときに感度が最大となるロジック。
                    // したがって、ここでは純粋に「光源->PMTベクトル」と「PMT向きベクトル」のcosを計算します。
                    
                    double cos_alpha = CalculateCosAngle(dx, dy, dz, hit.dir_x, hit.dir_y, hit.dir_z);
                    
                    // 角度依存係数 (Sigmoid function based on provided snippet)
                    // 1 / (1 + exp(-6 * (cos_alpha - 1)))
                    double ang_factor = 1.0 / (1.0 + std::exp(-6.0 * (cos_alpha - 1.0)));
                    
                    // 期待値 mu = (A * ang_factor) / R^2 + B
                    // ※ R^2 (dist2) は分母
                    mu = (A * ang_factor) / dist2 + B;

                } else {
                    // --- 従来モデル (Standard / ZeroIntercept) ---
                    // 等方発光と仮定: mu = A / R^2 + B
                    mu = A / dist2 + B;
                }
            }
            
            if (mu < 1e-9) mu = 1e-9; 

            double n = hit.charge;

            if (config.chargeType == ChargeChi2Type::BakerCousins) {
                // Baker-Cousins (Poisson)
                double term = 0.0;
                if (n > 1e-9) {
                    term = mu - n + n * std::log(n / mu);
                } else {
                    term = mu; 
                }
                chi2_total += 2.0 * term;
            } else {
                // Gaussian
                double sigma_q = 1.0; // 簡易値
                chi2_total += std::pow(n - mu, 2) / (sigma_q * sigma_q);
            }
        }
    }

    // -----------------------------------------------------
    // 2. 時間(Time)に関するChi2
    // -----------------------------------------------------
    if (config.timeType != TimeChi2Type::None) {
        double goodness_sum = 0.0;
        
        for (const auto& hit : hits) {
            if (!hit.isHit) continue; 

            // 光源からPMTまでの距離
            double dist = std::sqrt(std::pow(hit.x - x, 2) + 
                                    std::pow(hit.y - y, 2) + 
                                    std::pow(hit.z - z, 2));
            
            // 期待到達時間
            double t_flight = dist / C_LIGHT;
            double t_expected = t0 + t_flight;
            double t_obs = hit.time;

            if (config.timeType == TimeChi2Type::Goodness) {
                double sigma_t = GetTimeResolution(hit.ch, hit.charge);
                double res2 = std::pow(t_obs - t_expected, 2);
                goodness_sum += std::exp( -res2 / (2.0 * sigma_t * sigma_t) );

            } else if (config.timeType == TimeChi2Type::EMG) {
                double sigma = GetEMG_Sigma(hit.ch, hit.charge);
                double tau = GetEMG_Tau(hit.ch, hit.charge);
                chi2_total += CalcEMG_NLL(t_obs, t_expected, sigma, tau);

            } else {
                // Gaussian
                double sigma_t = 1.0; 
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

// ... (クラスの実装部分は変更なし、InitializeParametersでBの固定等はFitConfigに従う) ...
// ... (必要であれば InitializeParameters で B=0 固定のロジックを確認しますが、
//      ChargeModelType::Cosine の場合でも 3本ヒットなら B=0 にするべきです)

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
    // 重心計算等
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

        // Bの固定ルール:
        // ZeroInterceptモデルの場合 -> 強制固定
        // 3本ヒットオプションが有効で、かつユーザーが意図的にZeroInterceptを選んでいない場合でも、
        // 物理的にパラメータが決まらないためB=0推奨だが、
        // ここでは「ユーザー責任案B」に従い、モデル設定に忠実に従う（ZeroInterceptなら固定、Standard/Cosineなら解放）
        // ただし、CosineモデルでもBを固定したい場合は、main側でオプション制御が必要かもしれないが、
        // 今回は「CosineモデルはBを持つ(=解放)」として定義し、3本時は -m zeroB を使う運用とするか、
        // または CosineZeroB のようなモデルを作るか...
        // ひとまず Cosine モデルでは B を解放します（4本以上ヒットを想定）。
        // 3本ヒットでCosineモデルを使いたい場合、Bを固定するロジックが必要ですが、
        // 今回のコード構造では ChargeModelType::ZeroIntercept と分かれているため、
        // 「CosineモデルかつB固定」という組み合わせが現状のEnumでは表現しにくいです。
        
        // ★修正案:
        // chargeModel = Cosine の時も、もし3本ヒットならBを固定する、という自動制御を
        // 入れても良いですが、「ユーザー責任」の指示があったため、
        // 厳密には -m cosine_zeroB が必要になります。
        // ここでは簡易的に、Standard または Cosine の場合は B 解放、ZeroIntercept なら B 固定とします。
        
        if (fConfig.chargeModel == ChargeModelType::ZeroIntercept) {
            fMinuit->DefineParameter(5, "B", 0.0, 0.0, 0.0, 0.0);
            fMinuit->FixParameter(5);
        } else {
            fMinuit->DefineParameter(5, "B", 0.0, 0.1, 0, 1000);
            fMinuit->Release(5);
        }
    }
}
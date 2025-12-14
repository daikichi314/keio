/**
 * @file onemPMTfit.cc
 * @brief 光源位置再構成クラスの実装ファイル
 *
 * TMinuitのFCN関数を用いて、観測データとモデルの残差(Chi2または負の対数尤度)を最小化します。
 * * 主な機能:
 * - Baker-Cousins Chi2 の計算 (電荷用)
 * - EMG (Exponential Modified Gaussian) 分布の尤度計算 (時間用)
 * - 3本ヒット時のパラメータB固定処理
 * - Unhitチャンネルの時間フィット除外処理
 *
 * @author Gemini (Modified based on user request)
 * @date 2025-12-14
 */

#include "onemPMTfit.hh"
#include <iostream>
#include <cmath>
#include <TMath.h> // Erfc計算用

// グローバルポインタの実体定義
LightSourceFitter* gFitter = nullptr;

// =========================================================
// EMGパラメータ取得関数 (現在は仮実装)
// =========================================================
double GetEMG_Sigma(int ch, double charge) {
    // 本来は電荷Qやチャンネルごとの関数として実装
    return 1.0; // [ns]
}
double GetEMG_Tau(int ch, double charge) {
    // 本来は電荷Qやチャンネルごとの関数として実装
    return 1.0; // [ns]
}

// =========================================================
// EMG (Exponential Modified Gaussian) の負の対数尤度計算
// PDF: f(t) = (1/2tau) * exp( (sigma^2/2tau^2) - (t-mu)/tau ) * erfc(...)
// NLL: -ln(f(t)) を返す (Chi2として加算するため)
// =========================================================
double CalcEMG_NLL(double t, double mu, double sigma, double tau) {
    if (tau <= 0 || sigma <= 0) return 1e9; // 異常値回避

    // erfcの引数
    double arg_erfc = (sigma/tau - (t - mu)/sigma) / std::sqrt(2.0);
    // 指数部の引数 (符号に注意)
    double term_exp = (sigma*sigma)/(2.0*tau*tau) - (t - mu)/tau;
    
    // erfc計算
    double val_erfc = std::erfc(arg_erfc);
    if (val_erfc <= 1e-15) val_erfc = 1e-15; // log(0)回避

    // ln(f) = -ln(2tau) + term_exp + ln(erfc)
    double ln_f = -std::log(2.0 * tau) + term_exp + std::log(val_erfc);
    
    // Chi2と同様のスケールにするため -2倍して返す (-2lnL)
    return -2.0 * ln_f;
}

// =========================================================
// Minuit用 目的関数 (Chi2計算)
// =========================================================
void fcn_wrapper(int& npar, double* gin, double& f, double* par, int iflag) {
    // パラメータ展開
    double x = par[0];
    double y = par[1];
    double z = par[2];
    double t0 = par[3]; // 発光開始時刻
    double A = par[4];  // 光量係数 A
    double B = par[5];  // 背景光係数 B

    const auto& hits = gFitter->GetData();
    const auto& config = gFitter->GetConfig();

    double chi2_total = 0.0;

    // -----------------------------------------------------
    // 1. 電荷(Charge)に関するChi2
    // -----------------------------------------------------
    for (const auto& hit : hits) {
        // 幾何学的距離
        double dist = std::sqrt(std::pow(hit.x - x, 2) + 
                                std::pow(hit.y - y, 2) + 
                                std::pow(hit.z - z, 2));
        
        // 期待電荷モデル: mu = A / r^2 + B
        double mu = 0.0;
        if (dist > 1e-3) mu = A / (dist * dist) + B;
        if (mu < 1e-9) mu = 1e-9; // 0割・対数定義域エラー防止

        double n = hit.charge;

        // モデルによる切り替え
        if (config.chargeType == ChargeChi2Type::BakerCousins) {
            // Baker-Cousins Chi2 (Poisson like): 2 * (mu - n + n * ln(n / mu))
            double term = 0.0;
            if (n > 1e-9) {
                term = mu - n + n * std::log(n / mu);
            } else {
                // n=0 (Unhit) の極限: n*ln(n) -> 0 なので term = mu
                term = mu; 
            }
            chi2_total += 2.0 * term;
        } else {
            // Gaussian (従来): (n - mu)^2 / sigma^2
            // sigmaの扱いは要調整。ここでは簡易的に1.0またはmuの平方根などを想定
            double sigma_q = 1.0; 
            chi2_total += std::pow(n - mu, 2) / (sigma_q * sigma_q);
        }
    }

    // -----------------------------------------------------
    // 2. 時間(Time)に関するChi2
    // -----------------------------------------------------
    for (const auto& hit : hits) {
        // Unhitチャンネル (isHit == false) は時間情報を持たないので無視
        if (!hit.isHit) continue;

        double dist = std::sqrt(std::pow(hit.x - x, 2) + 
                                std::pow(hit.y - y, 2) + 
                                std::pow(hit.z - z, 2));
        
        double t_flight = dist / C_LIGHT;
        double t_expected = t0 + t_flight; // 期待到達時刻
        double t_obs = hit.time;           // 観測時刻

        // モデルによる切り替え
        if (config.timeType == TimeChi2Type::EMG) {
            double sigma = GetEMG_Sigma(hit.ch, hit.charge);
            double tau = GetEMG_Tau(hit.ch, hit.charge);
            // EMG分布のピーク/期待値の解釈については、
            // 「muは今と同じ」という指示に従い、Gaussian mean相当の位置にt_expectedを代入
            chi2_total += CalcEMG_NLL(t_obs, t_expected, sigma, tau);

        } else {
            // Gaussian (従来)
            double sigma_t = 1.0; // 簡易値
            chi2_total += std::pow(t_obs - t_expected, 2) / (sigma_t * sigma_t);
        }
    }

    f = chi2_total;
}

// =========================================================
// LightSourceFitter クラス実装
// =========================================================
LightSourceFitter::LightSourceFitter() {
    fMinuit = new TMinuit(6); // パラメータ数: x,y,z,t,A,B
    fMinuit->SetFCN(fcn_wrapper);
    gFitter = this; // グローバルポインタ設定
}

LightSourceFitter::~LightSourceFitter() {
    delete fMinuit;
}

void LightSourceFitter::SetConfig(const FitConfig& config) {
    fConfig = config;
}

bool LightSourceFitter::FitEvent(const std::vector<PMTData>& eventHits, FitResult& res) {
    fCurrentHits = eventHits;

    // パラメータ初期化 (重心計算など)
    InitializeParameters(eventHits);

    // -----------------------------------------------------
    // 3本ヒットの場合の特別処理: パラメータBを0に固定
    // -----------------------------------------------------
    int nHit = 0;
    for(const auto& hit : eventHits) {
        if(hit.isHit) nHit++;
    }

    if (nHit == 3) {
        // 3本ヒット時: B (index 5) を固定
        fMinuit->FixParameter(5);
        
        // 値を0に強制設定 (現在の値がずれている場合の対策)
        double cur_val, cur_err;
        fMinuit->GetParameter(5, cur_val, cur_err);
        if (std::abs(cur_val) > 1e-9) {
             fMinuit->DefineParameter(5, "B", 0.0, 0.0, 0.0, 0.0); 
        }
    } else {
        // 4本以上(通常): B を自由に動かす
        fMinuit->Release(5);
    }

    // 最小化実行 (MIGRAD法)
    double arglist[10];
    int ierflg = 0;
    arglist[0] = 5000; // Max calls
    arglist[1] = 0.1;  // Tolerance
    fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    // 結果取得
    double val, err;
    fMinuit->GetParameter(0, val, err); res.x = val; res.err_x = err;
    fMinuit->GetParameter(1, val, err); res.y = val; res.err_y = err;
    fMinuit->GetParameter(2, val, err); res.z = val; res.err_z = err;
    fMinuit->GetParameter(3, val, err); res.t = val; res.err_t = err;
    fMinuit->GetParameter(4, val, err); res.A = val;
    fMinuit->GetParameter(5, val, err); res.B = val;

    // 統計情報取得
    double fmin, fedm, errdef;
    int npari, nparx, istat;
    fMinuit->mnstat(fmin, fedm, errdef, npari, nparx, istat);
    
    res.chi2 = fmin;
    // 自由度 = データ点数(2*nHit: 電荷+時間) - パラメータ数
    // 3本時: パラメータ5個(B固定), 4本時: パラメータ6個
    // ※Unhitは電荷のみ(1点)寄与するが、時間(0点)は寄与しない等の厳密な自由度計算は
    // 解析の目的に応じて調整してください。ここでは簡易的に以下とします。
    // 有効ヒット数 nHit に対し、データ点は 2*nHit (Q,T)
    // Unhitがある場合、Qのみ増えるため、厳密には (2*nHit_real + 1*nUnhit) となる
    // ここでは単純化のため nHit (有効ヒット数) を基準に計算
    res.ndf = 2 * nHit - (nHit == 3 ? 5 : 6); 
    
    res.status = istat;

    return (istat == 3); // 3: 正常収束
}

void LightSourceFitter::InitializeParameters(const std::vector<PMTData>& hits) {
    // 電荷重心を計算して初期座標とする
    double sumQ = 0;
    double sumX = 0, sumY = 0, sumZ = 0;
    
    for (const auto& hit : hits) {
        if(hit.isHit) { // Hitのみ使う
            sumQ += hit.charge;
            sumX += hit.x * hit.charge;
            sumY += hit.y * hit.charge;
            sumZ += hit.z * hit.charge;
        }
    }
    
    double iniX = (sumQ > 0) ? sumX/sumQ : 0;
    double iniY = (sumQ > 0) ? sumY/sumQ : 0;
    double iniZ = 50.0; // Z方向の初期値は適当な正の値

    // パラメータ定義: ID, 名前, 初期値, ステップ, 最小, 最大
    fMinuit->DefineParameter(0, "x", iniX, 1.0, -200, 200);
    fMinuit->DefineParameter(1, "y", iniY, 1.0, -200, 200);
    fMinuit->DefineParameter(2, "z", iniZ, 1.0, 0, 200);
    fMinuit->DefineParameter(3, "t", 0.0, 1.0, -100, 100);
    fMinuit->DefineParameter(4, "A", 10000.0, 100.0, 0, 1000000);
    fMinuit->DefineParameter(5, "B", 0.0, 0.1, 0, 1000); 
}
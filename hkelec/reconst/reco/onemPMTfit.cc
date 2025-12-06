/*
 * id: onemPMTfit.cc
 * Place: /home/daiki/keio/hkelec/reconst/reco/
 * Author: Gemini 3 Pro
 * Last Edit: 2025-12-06
 *
 * 概要:
 * 光源フィッティングの実装ファイル
 * Minuit を用いた χ² 最小化によるフィット処理を実装します。
 * 光量モデル (距離減衰、角度依存) と時刻モデル (飛行時間) を定義し、
 * 観測データとの差異を計算して最適パラメータを推定します。
 *
 * 主な関数:
 * - CalculateExpectedCharge(): 予想光量計算
 * - CalculateChi2(): χ² 計算
 * - FitEvent(): イベント単位でのフィット実行
 */

#include "onemPMTfit.hh"
#include <iostream>
#include <cmath>
#include <TMath.h>

std::vector<PMTData> LightSourceFitter::g_hits;

LightSourceFitter::LightSourceFitter() {}
LightSourceFitter::~LightSourceFitter() {}

// 時刻分解能
double LightSourceFitter::GetSigmaTime(int ch, double charge) {
    return 1.0; 
}

// 光量モデル
double LightSourceFitter::CalculateExpectedCharge(const double* params, double distance, double cos_angle) {
    double A = params[4];
    double B = params[5];
    
    double dist2 = distance * distance;
    if (dist2 < 1.0) dist2 = 1.0;

    double ang_eff = 1.0 / (1.0 + std::exp(-6.0 * (cos_angle - 1.0)));
    
    return (A * ang_eff + B) / dist2;
}

double LightSourceFitter::CalculateChi2(double *par) {
    double chi2_total = 0.0;
    
    double t0 = par[0];
    double sx = par[1];
    double sy = par[2];
    double sz = par[3];

    for (const auto &hit : g_hits) {
        double dx = hit.x - sx;
        double dy = hit.y - sy;
        double dz = hit.z - sz;
        double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (dist < 1e-3) continue;

        // 1. 光量項
        double cos_theta = ((-dx)*hit.dir_x + (-dy)*hit.dir_y + (-dz)*hit.dir_z) / dist;
        double exp_Q = CalculateExpectedCharge(par, dist, cos_theta);
        double obs_Q = hit.charge;
        
        double sigma_Q2 = obs_Q;
        if (sigma_Q2 < 1.0) sigma_Q2 = 1.0;
        
        double chi2_Q = (obs_Q - exp_Q)*(obs_Q - exp_Q) / sigma_Q2;

        // 2. 時刻項
        double exp_T = (dist / C_LIGHT) + t0;
        double obs_T = hit.time - TIME_CORRECTION_VAL[hit.ch];
        double sigma_T = GetSigmaTime(hit.ch, hit.charge);

        double chi2_T = (obs_T - exp_T)*(obs_T - exp_T) / (sigma_T * sigma_T);

        chi2_total += chi2_Q + chi2_T;
    }
    return chi2_total;
}

void LightSourceFitter::FcnForMinuit(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag) {
    fval = CalculateChi2(par);
}

bool LightSourceFitter::FitEvent(const std::vector<PMTData> &hits, FitResult &result) {
    g_hits = hits;
    // 【修正】パラメータ数(6)に対してデータが少なすぎる場合はエラー
    if (hits.size() < 4) return false;

    TMinuit minuit(6);
    minuit.SetFCN(FcnForMinuit);
    minuit.SetPrintLevel(-1);

    double max_q = 0;
    for(const auto& h : hits) if(h.charge > max_q) max_q = h.charge;
    
    // パラメータ定義
    // 【修正】Zの初期値を100から170に変更 (光源位置の事前情報に近づける)
    minuit.DefineParameter(0, "t", 0.0, 0.1, -1000, 1000);
    minuit.DefineParameter(1, "x", 0.0, 1.0, -400, 400);
    minuit.DefineParameter(2, "y", 0.0, 1.0, -400, 400);
    minuit.DefineParameter(3, "z", 170.0, 1.0, -400, 400); 
    minuit.DefineParameter(4, "A", max_q * 10000, max_q*100, 0, 0); 
    minuit.DefineParameter(5, "B", 0.0, 1.0, 0, 0);

    int status = minuit.Migrad();

    double val, err;
    minuit.GetParameter(1, val, err); result.x = val; result.err_x = err;
    minuit.GetParameter(2, val, err); result.y = val; result.err_y = err;
    minuit.GetParameter(3, val, err); result.z = val; result.err_z = err;
    minuit.GetParameter(0, val, err); result.t = val; result.err_t = err;
    minuit.GetParameter(4, val, err); result.A = val;
    minuit.GetParameter(5, val, err); result.B = val;

    double fmin, fedm, errdef;
    int nvpar, nparx, istat;
    minuit.mnstat(fmin, fedm, errdef, nvpar, nparx, istat);
    result.chi2 = fmin;
    result.ndf = 2 * hits.size() - 6; 
    result.status = status;

    return (status == 0 || status == 4);
}
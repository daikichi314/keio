#include "onemPMTfit.hh"
#include <iostream>
#include <cmath>
#include <TMath.h>

// 静的メンバの実体定義
std::vector<PMTData> LightSourceFitter::g_hits;

// 時刻分解能を計算する関数
// 将来的にはPMTごとの関数形や電荷依存性をここに記述する
double GetSigmaTime(int ch, double charge) {
    // 例: return sqrt(pow(const_term, 2) + pow(stoch_term, 2)/charge);
    // 現在は固定値 1.0 ns を返す
    return 1.0;
}

LightSourceFitter::LightSourceFitter() {
}

LightSourceFitter::~LightSourceFitter() {
}

// ----------------------------------------------------------------------
// モデル関数: パラメータと幾何条件から予想光量を計算
// params: [t, x, y, z, A, B]
// ----------------------------------------------------------------------
double LightSourceFitter::CalculateExpectedCharge(const double* params, double distance, double cos_angle) {
    // シグモイドモデルを使用
    double A = params[4];
    double B = params[5];
    
    // distanceの2乗で減衰 (立体角)
    double distance2 = distance * distance;
    if (distance2 < 1.0) distance2 = 1.0; 

    // 角度依存性項 (係数は仮)
    double angular_acceptance = (1.0 / (1.0 + std::exp(-6.0 * (cos_angle - 1.0))));
    
    // モデル式
    double model = (A * angular_acceptance + B) / distance2;
    
    if (model < 0) model = 0;
    return model;
}

// ----------------------------------------------------------------------
// Chi2計算関数
// ----------------------------------------------------------------------
double LightSourceFitter::CalculateChi2(double *par) {
    double chi2_total = 0.0;
    
    double t0 = par[0];
    double src_x = par[1];
    double src_y = par[2];
    double src_z = par[3];

    for (const auto &hit : g_hits) {
        // --- 幾何計算 ---
        double dx = hit.x - src_x;
        double dy = hit.y - src_y;
        double dz = hit.z - src_z;
        double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
        
        if (dist < 1e-3) continue;

        // 入射角のコサイン計算 (PMT法線とPMT->光源ベクトルのなす角)
        double vec_pmt_to_src_x = -dx;
        double vec_pmt_to_src_y = -dy;
        double vec_pmt_to_src_z = -dz;
        double cos_theta = (vec_pmt_to_src_x*hit.dir_x + vec_pmt_to_src_y*hit.dir_y + vec_pmt_to_src_z*hit.dir_z) / dist;

        // --- 1. 光量のChi2 ---
        double expected_L = CalculateExpectedCharge(par, dist, cos_theta);
        double obs_L = hit.charge;
        
        // 誤差評価 (観測光量を分散の代用とする簡易法)
        double sigma_L2 = obs_L; 
        if (sigma_L2 < 1.0) sigma_L2 = 1.0; 

        double chi2_charge = (obs_L - expected_L) * (obs_L - expected_L) / sigma_L2;

        // --- 2. 時刻のChi2 ---
        double flight_time = dist / C_LIGHT;
        // time_diffの期待値 = 飛行時間 + 発光時刻t0
        double expected_time = flight_time + t0;
        
        // PMTごとの補正値を取得
        double time_correction = TIME_CORRECTION_VAL[hit.ch];
        
        // 観測値 (測定時刻 - 補正値)
        double obs_time = hit.time - time_correction;
        
        // 分解能を取得 (電荷依存性を考慮可能)
        double sigma_t = GetSigmaTime(hit.ch, hit.charge);
        
        double chi2_time = (obs_time - expected_time) * (obs_time - expected_time) / (sigma_t * sigma_t);

        // --- 合計 ---
        chi2_total += chi2_charge + chi2_time;
    }

    return chi2_total;
}

void LightSourceFitter::FcnForMinuit(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag) {
    fval = CalculateChi2(par);
}

bool LightSourceFitter::FitEvent(const std::vector<PMTData> &hits, FitResult &result) {
    g_hits = hits; 

    if (hits.size() < 1) return false;

    TMinuit minuit(6);
    minuit.SetFCN(FcnForMinuit);
    minuit.SetPrintLevel(-1);

    // 初期値の設定
    double init_x = 0;
    double init_y = 0;
    double init_z = 100;
    double init_t = 0;
    
    double max_q = 0;
    for(const auto& h : hits) if(h.charge > max_q) max_q = h.charge;
    double init_A = max_q * 10000; 
    double init_B = 0;

    minuit.DefineParameter(0, "t", init_t, 0.1, -1000, 1000); 
    minuit.DefineParameter(1, "x", init_x, 1.0, -400, 400);   
    minuit.DefineParameter(2, "y", init_y, 1.0, -400, 400);
    minuit.DefineParameter(3, "z", init_z, 1.0, -400, 400);
    minuit.DefineParameter(4, "A", init_A, init_A*0.1, 0, 0); 
    minuit.DefineParameter(5, "B", init_B, 1.0, 0, 0);

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
    result.ndf = 2 * hits.size() - minuit.GetNumFreePars();
    result.status = status;

    return (status == 0 || status == 4); 
}
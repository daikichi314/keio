#ifndef ONEMPMTFIT_HH
#define ONEMPMTFIT_HH

#include "fittinginput.hh"
#include <vector>
#include <TMinuit.h>
#include <cmath>

// /////////////////////////////////////////////////////////
// 追加: フィッティングモードの定義
// /////////////////////////////////////////////////////////
enum FitMode {
    kChi2_Gaussian,      // [Default] 従来のカイ二乗法 (ガウス分布仮定)
    kLikelihood_Poisson, // ポアソン尤度 (Unhit考慮 + Baker-Cousins Chi2)
    kLikelihood_EMG,     // EMG時間項 + ポアソン電荷 (fiTQun風)
    kCustom_Template     // 将来のためのカスタムテンプレート
};
// /////////////////////////////////////////////////////////

class LightSourceFitter {
public:
    LightSourceFitter();
    ~LightSourceFitter();

    // イベントを受け取ってフィットを実行するメイン関数
    bool FitEvent(const std::vector<PMTData> &hits, FitResult &result);
    
    // /////////////////////////////////////////////////////////
    // 追加: フィッティングモードを切り替える静的関数
    // 例: LightSourceFitter::SetFitMode(kLikelihood_EMG);
    // /////////////////////////////////////////////////////////
    static void SetFitMode(FitMode mode) { g_fitMode = mode; }

    // TMinuitが呼び出す関数 (staticである必要があります)
    static void FcnForMinuit(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);

private:
    // TMinuitからアクセスするための静的メンバ
    static std::vector<PMTData> g_hits;
    static FitMode g_fitMode; // 現在のモード

    // ヘルパー関数群
    static double GetSigmaTime(int ch, double charge);
    static double CalculateExpectedCharge(const double* params, double distance, double cos_angle);
    
    // カイ二乗(または負の対数尤度)を計算するコア関数
    static double CalculateChi2(double *par);
    
    // /////////////////////////////////////////////////////////
    // 追加: EMG分布の負の対数尤度 (-2lnL) を計算する関数
    // /////////////////////////////////////////////////////////
    static double CalculateNLL_EMG(double t, double mu, double sigma, double tau);
};

#endif
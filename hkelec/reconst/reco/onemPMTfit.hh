/*
 * id: onemPMTfit.hh
 * Place: /home/daiki/keio/hkelec/reconst/reco/
 * Author: Gemini 3 Pro
 * Last Edit: 2025-12-06
 *
 * 概要:
 * 光源フィッティングクラスのヘッダーファイル
 * LightSourceFitter クラスの定義と Minuit 用コールバック関数の宣言を含みます。
 */

#ifndef ONEMPMTFIT_HH
#define ONEMPMTFIT_HH // インクルードガード

#include "fittinginput.hh"
#include <vector>
#include <TMinuit.h> // ROOTの最小化ライブラリ

// フィッティングを行うクラス
class LightSourceFitter {
public:
    LightSourceFitter();  // コンストラクタ
    ~LightSourceFitter(); // デストラクタ

    // イベントごとのフィットを実行する関数
    // hits: 入力データのリスト, result: 結果格納用構造体
    bool FitEvent(const std::vector<PMTData> &hits, FitResult &result);

    // Minuitから呼び出される静的関数 (目的関数FCN)
    // 引数の形式はTMinuitの仕様で決まっている
    static void FcnForMinuit(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);

private:
    static std::vector<PMTData> g_hits; // Minuitに渡すためにデータを保持する静的変数

    // 時刻分解能を取得する関数
    static double GetSigmaTime(int ch, double charge);

    // 予想される光量を計算するモデル関数
    static double CalculateExpectedCharge(const double* params, double distance, double cos_angle);
    
    // カイ二乗値を計算するヘルパー関数
    static double CalculateChi2(double *par);
};

#endif // ONEMPMTFIT_HH
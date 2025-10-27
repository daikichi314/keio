/*
 * id: tts_fitter.h
 * Place: ~/hkelec/DiscreteSoftware/Analysis/macro/fit_results/
 * Last Edit: 2025-10-18 Gemini
 *
 * 概要: 時間分布ヒストグラムをフィットするための関数群を定義したヘッダーファイル。
 * 2段階フィットとフィット範囲の安全機構を実装。
 * AsymGaus, EMG, ExpGausの関数を選択可能。
 * コンパイル不要 (ヘッダーファイル)
 */
#include "TF1.h"
#include "TH1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMath.h"
#include <algorithm> // for std::max, std::min

// 1. フィット結果を格納するための構造体
struct TTSFitResult {
    double tts = 0;
    double sigma = 0;
    double fwhm = 0;
    double peak = 0;
    double chi2 = 0;
    int ndf = 0;
    double tau = 0; 
};

// 2. フィット用の関数定義
Double_t fitf_asymgaus(Double_t *x, Double_t *par) {
    Double_t ret = 0;
    if (x[0] < par[1]) { ret = par[0] * TMath::Gaus(x[0], par[1], par[2]); }
    else { ret = par[0] * TMath::Gaus(x[0], par[1], par[3]); }
    return ret;
}
Double_t fitf_emg(Double_t *x, Double_t *par) {
    Double_t den = par[3] > 0 ? par[3] : 1e-9;
    Double_t arg = (par[2]*par[2]/2. + par[3]*(par[1]-x[0]))/den;
    Double_t ret = par[0] * TMath::Exp(arg) * TMath::Erfc( (par[1]-x[0])/sqrt(2.)/par[2] + par[2]/sqrt(2.)/den );
    return ret;
}
Double_t fitf_expgaus(Double_t *x, Double_t *par) {
    Double_t den = par[3] > 0 ? par[3] : 1e-9;
    Double_t arg = par[2]*par[2]/2./den/den - (x[0]-par[1])/den;
    Double_t ret = par[0] * TMath::Exp(arg) * ( 1-TMath::Erf( (par[2]/sqrt(2.)/den - (x[0]-par[1])/sqrt(2.)/par[2]) ) );
    return ret;
}

// 3. 時間分布フィットを実行するメイン関数
TTSFitResult perform_tts_fit(TH1* hist) {
    TTSFitResult result;
    if (!hist || hist->GetEntries() < 100) return result;

    // --- ステップ1: プリフィットのための準備 ---
    double hist_min_x = hist->GetXaxis()->GetXmin();
    double hist_max_x = hist->GetXaxis()->GetXmax();
    int bmax = hist->GetMaximumBin();
    double xmax = hist->GetXaxis()->GetBinCenter(bmax);
    double xrms = hist->GetRMS();
    if (xrms == 0) return result;

    // --- ステップ2: 広い範囲でガウス関数によるプリフィットを実行 ---
    double pre_fit_min = std::max(hist_min_x, xmax - 5 * xrms);
    double pre_fit_max = std::min(hist_max_x, xmax + 5 * xrms);
    TF1* f_prefit = new TF1("f_prefit", "gaus", pre_fit_min, pre_fit_max);
    TFitResultPtr pre_fit_res = hist->Fit(f_prefit, "QNRS");

    if (!pre_fit_res->IsValid()) {
        delete f_prefit;
        return result;
    }

    // --- ステップ3: 最終フィットの範囲と初期値を決定 ---
    double refined_mean = pre_fit_res->Parameter(1);
    double refined_sigma = pre_fit_res->Parameter(2);
    if (refined_sigma == 0) { delete f_prefit; return result; }
    
    double final_fit_min = std::max(hist_min_x, refined_mean - 3 * refined_sigma);
    double final_fit_max = std::min(hist_max_x, refined_mean + 3 * refined_sigma);

    // --- ステップ4: 最終的なフィット関数を選択して実行 ---
    Bool_t IsAsymGaus = kFALSE;
    Bool_t IsEMG      = kTRUE;
    Bool_t IsExpGaus  = kFALSE;
    
    TF1* fitFunc = nullptr;
    if (IsAsymGaus) {
        fitFunc = new TF1("fitFunc", fitf_asymgaus, final_fit_min, final_fit_max, 4);
        fitFunc->SetParameters(pre_fit_res->Parameter(0), refined_mean, refined_sigma, refined_sigma);
    } else if (IsEMG) {
        fitFunc = new TF1("fitFunc", fitf_emg, final_fit_min, final_fit_max, 4);
        fitFunc->SetParameters(pre_fit_res->Parameter(0), refined_mean, refined_sigma, refined_sigma);
        fitFunc->SetParLimits(2, 0, 1e10);// 10/20:エラーがあるとするとこの辺が怪しい?
        fitFunc->SetParLimits(3, 0, 1e10);
    } else if (IsExpGaus) {
        fitFunc = new TF1("fitFunc", fitf_expgaus, final_fit_min, final_fit_max, 4);
        fitFunc->SetParameters(pre_fit_res->Parameter(0), refined_mean, refined_sigma, refined_sigma);
        fitFunc->SetParLimits(2, 0, 1e10);
        fitFunc->SetParLimits(3, 0, 1e10);
    }
    
    if (!fitFunc) { delete f_prefit; return result; }
    
    TFitResultPtr fit_res = hist->Fit(fitFunc, "SQR");

    // --- ステップ5: 結果を抽出 ---
    if (fit_res.Get() && fit_res->IsValid()) {
        result.peak = fitFunc->GetMaximumX(final_fit_min, final_fit_max);
        result.sigma = fit_res->Parameter(2);
        result.chi2 = fit_res->Chi2();
        result.ndf = fit_res->Ndf();
        
        double half_max = fitFunc->GetMaximum(final_fit_min, final_fit_max) / 2.0;
        double x1 = fitFunc->GetX(half_max, final_fit_min, result.peak);
        // ★★★ ここが修正点です: fit_max -> final_fit_max ★★★
        double x2 = fitFunc->GetX(half_max, result.peak, final_fit_max);
        result.fwhm = x2 - x1;

        if (IsAsymGaus) {
            result.tts = (fit_res->Parameter(2) > fit_res->Parameter(3)) ? fit_res->Parameter(2) : fit_res->Parameter(3);
        } else {
            result.tau = fit_res->Parameter(3);
            result.tts = sqrt(pow(result.sigma, 2) + pow(result.tau, 2));
        }
        hist->GetListOfFunctions()->Add(fitFunc->Clone());
    }

    delete f_prefit;
    return result;
}
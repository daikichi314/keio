/*
 * id: peakfinder.C
 * Place: ~/hkelec/DiscreteSoftware/Analysis/macro/fit_results/
 * Last Edit: 2025-10-29 Gemini
 *
 * 概要: 電荷はシンプルなピーク検出、時間は高精度なTTSフィットを行う。
 * (ttshistofit.C のロジックを統合し、tts_fitter.h を不要にしたバージョン)
 *
 * コンパイル:
 * g++ peakfinder.C -o peakfinder $(root-config --cflags --glibs)
 */

// 1. ヘッダーファイルのインクルード
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
#include <algorithm>
// (tts_fitter.h は不要になりました)

// --- 2. グローバル設定 (ttshistofit.C より) ---
Bool_t IsAsymGaus = kFALSE;
Bool_t IsEMG = kTRUE;
Bool_t IsExpGaus = kFALSE;


// --- 3. ユーティリティ関数 (ttshistofit.C より) ---

// 3a. EMG (Exponentially Modified Gaussian) 関数の定義
// par[0]: #mu (ガウス中心), par[1]: #gamma (振幅スケール), par[2]: #sigma (ガウス幅), par[3]: #lambda (1/tau)
Double_t EMG(Double_t *x, Double_t *par)
{
    if (par[2] == 0) return 0;
    if (par[3] == 0) return 0;
    return 0.5*par[3]*exp(0.5*par[3]*(2*par[0]+par[3]*par[2]*par[2]-2.*x[0]))
           *TMath::Erfc((par[0]+par[3]*par[2]*par[2]-x[0])/(sqrt(2.)*par[2]))*par[1];
}

// 3b. FWHM (半値全幅) を計算する関数
Double_t GetFWHM(TF1 *f)
{
    if (!f) return 0;
    double peak_pos = f->GetMaximumX();
    double half_max = f->GetMaximum() * 0.5;
    double x_min = f->GetXaxis()->GetXmin();
    double x_max = f->GetXaxis()->GetXmax();
    double x1 = f->GetX(half_max, x_min, peak_pos);
    double x2 = f->GetX(half_max, peak_pos, x_max);
    return x2 - x1;
}

// 3c. EMG関数のピーク位置を解析的に計算する関数
Double_t GetPeak(TF1 *f)
{
    if (!f) return 0;
    
    // EMG関数のパラメータを取得
    double mu = f->GetParameter(0);     // ガウス中心
    double sigma = f->GetParameter(2);  // ガウス幅
    double lambda = f->GetParameter(3); // 1/tau
    
    // EMG関数の最大値の位置を解析的に計算
    // peak = μ + σ²λ
    return mu + sigma * sigma * lambda;
}

// --- 4. ファイル名から電圧を抽出する関数 (変更なし) ---
double get_voltage_from_filename(const std::string& filename) {
    std::regex re("(\\d+)V");
    std::smatch match;
    if (std::regex_search(filename, match, re) && match.size() > 1) {
        return std::stod(match.str(1));
    }
    return -1.0;
}

// --- 5. 電荷情報のピークを探す関数 (変更なし) ---
// (gausfit.C とは異なり、シンプルな最大ビン検出を使用)
void find_charge_peaks(TString input_filename) {
    auto infile = TFile::Open(input_filename, "READ");
    if (!infile || infile->IsZombie()) return;

    TString output_txt_filename = input_filename;
    output_txt_filename.ReplaceAll("_eventhist.root", "_peak.txt");
    std::ofstream outfile(output_txt_filename.Data());
    outfile << "# ch,type,voltage,peak_pos" << std::endl;
    
    double voltage = get_voltage_from_filename(input_filename.Data());
    std::vector<std::string> hist_types = {"hgain", "lgain", "tot"};
    
    for (int ch = 0; ch < 12; ++ch) {
        for (const auto& type : hist_types) {
            TString hist_name = Form("h_%s_ch%d", type.c_str(), ch);
            auto hist = infile->Get<TH1D>(hist_name);
            if (!hist || hist->GetEntries() == 0) continue;
            int max_bin = hist->GetMaximumBin();
            double peak_pos = hist->GetXaxis()->GetBinCenter(max_bin);
            outfile << ch << "," << type << "," << voltage << "," << peak_pos << "," << input_filename.Data() << std::endl;
        }
    }
    std::cout << "電荷ピーク検出完了 -> " << output_txt_filename << std::endl;
    outfile.close();
    infile->Close();
}


// --- 6. ★★★ 時間フィット関数 (gausfit.C と同様のロジックに置換) ★★★ ---
void fit_time(TString input_filename, bool save_pdf) {
    auto infile = TFile::Open(input_filename, "READ");
    if (!infile || infile->IsZombie()) return;

    TString output_txt_filename = input_filename;
    output_txt_filename.ReplaceAll("_eventhist.root", "_timefit.txt");
    std::ofstream outfile(output_txt_filename.Data());
    
    // 6a. 出力ヘッダー (gausfit.C / ttshistofit.C と同一)
    outfile << "# ch,type,voltage,tts(sigma),sigma,fwhm(calc),peak(calc),tau(1/lambda),chi2_ndf" << std::endl;

    double voltage = get_voltage_from_filename(input_filename.Data());
    
    if (save_pdf) {
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(1);
    }

    // 6b. time_diff のみフィット対象とする (tdc_diff を除外)
    std::vector<std::string> hist_types = {"time_diff"};

    for (int ch = 0; ch < 12; ++ch) {
        for (const auto& type : hist_types) {
            TString hist_name = Form("h_%s_ch%d", type.c_str(), ch);
            auto hist = infile->Get<TH1D>(hist_name);
            
            // 6c. エントリ数チェック (ttshistofit.C に倣い 100 に)
            if (!hist || hist->GetEntries() < 100) continue;

            // 6d. ttshistofit.C のフィッティングロジック
            double hist_min = hist->GetXaxis()->GetXmin();
            double hist_max = hist->GetXaxis()->GetXmax();

            // ステップA: ガウス関数によるプレフィット
            TF1 *fgaus = new TF1("fgaus", "gaus", hist_min, hist_max);
            fgaus->SetLineColor(kCyan);
            fgaus->SetLineWidth(1);
            fgaus->SetParameter(1, hist->GetBinCenter(hist->GetMaximumBin()));
            fgaus->SetParameter(2, hist->GetRMS());
            hist->Fit(fgaus, "QN", "", hist_min, hist_max);
            
            double pre_amp = fgaus->GetParameter(0);
            double pre_mean = fgaus->GetParameter(1);
            double pre_sigma = TMath::Abs(fgaus->GetParameter(2));
            if (pre_sigma == 0) { delete fgaus; continue; }

            // ステップB: EMG関数による本フィット
            TF1 *emg = nullptr;
            TFitResultPtr fit_result = nullptr;

            if (IsEMG) {
                emg = new TF1("emg", EMG, hist_min, hist_max, 4);
                emg->SetLineColor(kRed);
                emg->SetLineStyle(2);
                emg->SetNpx(2000);
                
                emg->SetParName(0, "#mu");
                emg->SetParName(1, "#gamma");
                emg->SetParName(2, "#sigma");
                emg->SetParName(3, "#lambda");
                
                emg->SetParameter(0, pre_mean);
                emg->SetParameter(1, pre_amp * 10.0);
                emg->SetParameter(2, pre_sigma * 0.7);
                emg->SetParameter(3, (pre_sigma > 1e-9) ? (1. / pre_sigma) : 1.0);

                emg->SetParLimits(1, 1, 1e9);
                emg->SetParLimits(2, 0.01, 100);
                emg->SetParLimits(3, 0.001, 500);

                fit_result = hist->Fit(emg, "SQR", "", hist_min, hist_max);
            }
            
            // 6e. 結果の抽出とテキストファイルへの出力
            if (emg && fit_result.Get() && fit_result->IsValid() && fit_result->Ndf() > 0) {
                double sigma = emg->GetParameter(2);
                double lambda = emg->GetParameter(3);
                double tau = (lambda > 1e-9) ? (1.0 / lambda) : 0.0;
                double tts = sigma;
                double fwhm = GetFWHM(emg);
                double peak = GetPeak(emg);
                double chi2_ndf = fit_result->Chi2() / fit_result->Ndf();

                outfile << ch << "," << type << "," << voltage << ","
                        << tts << "," << sigma << "," << fwhm << ","
                        << peak << "," << tau << "," << chi2_ndf << std::endl;
            }

            // 6f. PDFの保存
            if (save_pdf && emg && fit_result.Get() && fit_result->IsValid()) {
                TCanvas* canvas = new TCanvas("c", "c", 800, 600);
                hist->GetXaxis()->SetRangeUser(hist->GetBinCenter(hist->GetMaximumBin())-15, hist->GetBinCenter(hist->GetMaximumBin())+20);
                hist->Draw();
                emg->Draw("same");
                
                TString pdf_name = input_filename;
                pdf_name.ReplaceAll("_eventhist.root", Form("_%s_fit.pdf", hist_name.Data()));
                canvas->SaveAs(pdf_name);
                delete canvas;
            }

            // 6g. メモリ解放
            delete fgaus;
            if (emg) delete emg;
        }
    }
    std::cout << "時間フィット完了 -> " << output_txt_filename << std::endl;
    outfile.close();
    infile->Close();
}

// --- 7. main関数 (変更なし) ---
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "使い方: " << argv[0] << " <input.root> [--fit-charge | --fit-time | --fit-all] [--no-pdf]" << std::endl;
        return 1;
    }

    std::string fit_mode = "--fit-charge"; // デフォルト
    bool save_pdf = true;
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--fit-charge" || arg == "--fit-time" || arg == "--fit-all") {
            fit_mode = arg;
        }
        if (arg == "--no-pdf") {
            save_pdf = false;
        }
    }
    
    if (fit_mode == "--fit-charge" || fit_mode == "--fit-all") {
        find_charge_peaks(argv[1]);
    }
    if (fit_mode == "--fit-time" || fit_mode == "--fit-all") {
        fit_time(argv[1], save_pdf);
    }

    return 0;
}
/*
 * id: meanfinder.C
 * Place: ~/hkelec/DiscreteSoftware/Analysis/macro/fit_results/
 * Last Edit: 2025-10-29 Gemini
 *
 * 概要: 電荷はヒストグラムの平均値 (GetMean()) で算出し、
 * 時間は高精度なEMGフィット (ttshistofit.Cベース) を行う。
 *
 * 電荷計算ロジック:
 * 1. hgain と lgain の両方の平均値を取得。
 * 2. hgain の最終ビンが飽和しているか判定
 * (最終ビン > (最終-1)ビン * 5.0)
 * 3. 飽和していたら lgain を、していなければ hgain を採用。
 *
 * コンパイル:
 * g++ meanfinder.C -o meanfinder $(root-config --cflags --glibs)
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

// --- 2. グローバル設定 (ttshistofit.C より) ---
Bool_t IsAsymGaus = kFALSE;
Bool_t IsEMG = kTRUE;
Bool_t IsExpGaus = kFALSE;


// --- 3. ユーティリティ関数 (ttshistofit.C ベース) ---

// 3a. EMG (Exponentially Modified Gaussian) 関数の定義
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

// 3c. ピーク位置を計算する関数
Double_t GetPeak(TF1 *f)
{
    if (!f) return 0;
    return f->GetMaximumX();
}

// --- 4. 電圧取得関数 (変更なし) ---
double get_voltage_from_filename(const std::string& filename) {
    std::regex re("(\\d+)V");
    std::smatch match;
    if (std::regex_search(filename, match, re) && match.size() > 1) {
        return std::stod(match.str(1));
    }
    return -1.0;
}

// --- 5. ★★★ 電荷計算関数 (新規ロジック) ★★★ ---
// (fit_charge から calculate_charge_mean に変更)
void calculate_charge_mean(TString input_filename) {
    auto infile = TFile::Open(input_filename, "READ");
    if (!infile || infile->IsZombie()) return;

    TString output_txt_filename = input_filename;
    // 1. 出力ファイル名を _mean.txt に変更
    output_txt_filename.ReplaceAll("_eventhist.root", "_mean.txt");
    std::ofstream outfile(output_txt_filename.Data());
    
    // 2. 出力ヘッダーを変更 (peak, sigma ではなく mean, mean_err)
    outfile << "# ch,type,voltage,mean,mean_err,rms" << std::endl;
    
    double voltage = get_voltage_from_filename(input_filename.Data());
    
    // tot も含めて計算
    std::vector<std::string> hist_types = {"hgain", "lgain", "tot"};
    
    for (int ch = 0; ch < 12; ++ch) {
        
        // 3. hgain/lgain の平均値とRMSを先に取得
        double hgain_mean = 0, hgain_mean_err = 0, hgain_rms = 0;
        double lgain_mean = 0, lgain_mean_err = 0, lgain_rms = 0;
        bool hgain_saturated = false;

        auto hist_hgain = infile->Get<TH1D>(Form("h_hgain_ch%d", ch));
        if (hist_hgain && hist_hgain->GetEntries() > 0) {
            hgain_mean = hist_hgain->GetMean();
            hgain_mean_err = hist_hgain->GetMeanError();
            hgain_rms = hist_hgain->GetRMS();
            
            // 4. ★★★ hgain 飽和判定ロジック ★★★
            // （hgainのヒストグラムの一番右端のビンの値が右端から2番目のビンの値の5倍以上）
            int n_bins = hist_hgain->GetNbinsX();
            if (n_bins >= 2) {
                double last_bin_content = hist_hgain->GetBinContent(n_bins);
                double second_last_bin_content = hist_hgain->GetBinContent(n_bins - 1);
                
                // 5. 0除算を避け、5倍の条件を確認
                if (second_last_bin_content > 1e-9 && (last_bin_content > second_last_bin_content * 5.0)) {
                    hgain_saturated = true;
                }
            }
        }

        auto hist_lgain = infile->Get<TH1D>(Form("h_lgain_ch%d", ch));
        if (hist_lgain && hist_lgain->GetEntries() > 0) {
            lgain_mean = hist_lgain->GetMean();
            lgain_mean_err = hist_lgain->GetMeanError();
            lgain_rms = hist_lgain->GetRMS();
        }

        // 6. 飽和判定に基づき、"hgain" または "lgain" として出力
        if (hgain_saturated) {
            // lgain を採用
            outfile << ch << ",lgain," << voltage << "," << lgain_mean << "," << lgain_mean_err << "," << lgain_rms << std::endl;
        } else {
            // hgain を採用
            outfile << ch << ",hgain," << voltage << "," << hgain_mean << "," << hgain_mean_err << "," << hgain_rms << std::endl;
        }

        // 7. tot は個別に処理
        auto hist_tot = infile->Get<TH1D>(Form("h_tot_ch%d", ch));
        if (hist_tot && hist_tot->GetEntries() > 0) {
            outfile << ch << ",tot," << voltage << "," 
                    << hist_tot->GetMean() << "," 
                    << hist_tot->GetMeanError() << "," 
                    << hist_tot->GetRMS() << std::endl;
        }
    }
    std::cout << "電荷 平均値計算完了 -> " << output_txt_filename << std::endl;
    outfile.close();
    infile->Close();
}


// --- 6. 時間フィット関数 (gausfit.C と同一) ---
void fit_time(TString input_filename, bool save_pdf) {
    auto infile = TFile::Open(input_filename, "READ");
    if (!infile || infile->IsZombie()) return;

    TString output_txt_filename = input_filename;
    output_txt_filename.ReplaceAll("_eventhist.root", "_timefit.txt");
    std::ofstream outfile(output_txt_filename.Data());
    outfile << "# ch,type,voltage,tts(sigma),sigma,fwhm(calc),peak(calc),tau(1/lambda),chi2_ndf" << std::endl;

    double voltage = get_voltage_from_filename(input_filename.Data());
    
    if (save_pdf) {
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(1);
    }

    std::vector<std::string> hist_types = {"time_diff"};

    for (int ch = 0; ch < 12; ++ch) {
        for (const auto& type : hist_types) {
            TString hist_name = Form("h_%s_ch%d", type.c_str(), ch);
            auto hist = infile->Get<TH1D>(hist_name);
            
            if (!hist || hist->GetEntries() < 100) continue;

            double hist_min = hist->GetXaxis()->GetXmin();
            double hist_max = hist->GetXaxis()->GetXmax();

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

            delete fgaus;
            if (emg) delete emg;
        }
    }
    std::cout << "Time fit completed. -> " << output_txt_filename << std::endl;
    outfile.close();
    infile->Close();
}

// --- 7. main関数 (オプション解析) ---
int main(int argc, char* argv[]) {
    if (argc < 2) {
        // 1. 使い方を meanfinder に修正
        std::cerr << "使い方: " << argv[0] << " <input.root> [--fit-charge | --fit-time | --fit-all] [--no-pdf]" << std::endl;
        return 1;
    }

    // 2. デフォルトを --fit-charge に設定
    std::string fit_mode = "--fit-charge";
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
    
    // 3. --fit-charge または --fit-all の場合
    if (fit_mode == "--fit-charge" || fit_mode == "--fit-all") {
        // 4. 電荷計算関数 (calculate_charge_mean) を呼び出す
        calculate_charge_mean(argv[1]);
    }
    
    // 5. --fit-time または --fit-all の場合
    if (fit_mode == "--fit-time" || fit_mode == "--fit-all") {
        // 6. 時間フィット関数 (fit_time) を呼び出す
        //    (電荷計算とPDF保存オプションを共有)
        fit_time(argv[1], save_pdf);
    }

    return 0;
}
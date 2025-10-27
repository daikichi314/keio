/*
 * id: gausfit.C
 * Place: ~/hkelec/DiscreteSoftware/Analysis/macro/fit_results/
 * Last Edit: 2025-10-16 Gemini
 * (修正: 2025-10-27 Gemini (tdc_diff のフィットをスキップ) )
 *
 * 概要: 信号データ(..._eventhist.root)を読み込み、電荷または時間のヒストグラムをフィットする。
 * オプションで解析対象を選択可能。
 * コンパイル可能
 */
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
#include <algorithm>
// 1. 新しいヘッダーをインクルード
#include "tts_fitter.h" 

// --- (get_voltage_from_filename関数は変更なし) ---
double get_voltage_from_filename(const std::string& filename) {
    std::regex re("(\\d+)V");
    std::smatch match;
    if (std::regex_search(filename, match, re) && match.size() > 1) {
        return std::stod(match.str(1));
    }
    return -1.0;
}

// 2. 電荷フィットを行う関数 (以前のprocess_file)
void fit_charge(TString input_filename, bool save_pdf) {
    auto infile = TFile::Open(input_filename, "READ");
    if (!infile || infile->IsZombie()) return;

    TString output_txt_filename = input_filename;
    output_txt_filename.ReplaceAll("_eventhist.root", "_gausfit.txt");
    std::ofstream outfile(output_txt_filename.Data());
    outfile << "# ch,type,voltage,peak,peak_err,sigma,sigma_err,chi2_ndf,rough_sigma" << std::endl;
    
    double voltage = get_voltage_from_filename(input_filename.Data());
    gStyle->SetOptFit(1111);

    std::vector<std::string> hist_types = {"hgain", "lgain", "tot"};
    for (int ch = 0; ch < 12; ++ch) {
        for (const auto& type : hist_types) {
            // ... (gausfitのロジックは変更なし)
            TString hist_name = Form("h_%s_ch%d", type.c_str(), ch);
            auto hist = infile->Get<TH1D>(hist_name);
            if (!hist || hist->GetEntries() < 200) continue;
            int max_bin = hist->GetMaximumBin();
            double rough_peak_pos = hist->GetXaxis()->GetBinCenter(max_bin);
            double rough_sigma = hist->GetStdDev();
            if (rough_sigma == 0) continue;
            TF1* f_prefit = new TF1("f_prefit", "gaus", std::max(hist->GetXaxis()->GetXmin(), rough_peak_pos - 5*rough_sigma), std::min(hist->GetXaxis()->GetXmax(), rough_peak_pos + 5*rough_sigma));
            TFitResultPtr pre_fit_res = hist->Fit(f_prefit, "QNRS");
            if(!pre_fit_res->IsValid()) {delete f_prefit; continue;}
            double r_mean = pre_fit_res->Parameter(1), r_sigma = pre_fit_res->Parameter(2);
            if(r_sigma==0) {delete f_prefit; continue;}
            TF1* f_final = new TF1("f_final", "gaus", std::max(hist->GetXaxis()->GetXmin(), r_mean - 2*r_sigma), std::min(hist->GetXaxis()->GetXmax(), r_mean + 2*r_sigma));
            TFitResultPtr fit_result = hist->Fit(f_final, "SQR");
            if (fit_result.Get() && fit_result->IsValid() && fit_result->Ndf() > 0) {
                 outfile << ch << "," << type << "," << voltage << "," << fit_result->Parameter(1) << "," << fit_result->ParError(1) << "," << std::abs(fit_result->Parameter(2)) << "," << fit_result->ParError(2) << "," << fit_result->Chi2() / fit_result->Ndf() << "," << rough_sigma << std::endl;
            }
            if (save_pdf && fit_result.Get() && fit_result->IsValid()) {
                TCanvas* canvas = new TCanvas("c", "c", 800, 600); hist->Draw(); f_final->Draw("same");
                TString pdf_name = input_filename; pdf_name.ReplaceAll("_eventhist.root", Form("_%s_fit.pdf", hist_name.Data())); canvas->SaveAs(pdf_name); delete canvas;
            }
            delete f_prefit; delete f_final;
        }
    }
    std::cout << "Charge fit completed. -> " << output_txt_filename << std::endl;
    outfile.close();
    infile->Close();
}

// 3. 時間フィットを行う関数 (新規追加)
void fit_time(TString input_filename, bool save_pdf) {
    auto infile = TFile::Open(input_filename, "READ");
    if (!infile || infile->IsZombie()) return;

    TString output_txt_filename = input_filename;
    output_txt_filename.ReplaceAll("_eventhist.root", "_timefit.txt");
    std::ofstream outfile(output_txt_filename.Data());
    outfile << "# ch,type,voltage,tts,sigma,fwhm,peak,tau,chi2_ndf" << std::endl;

    double voltage = get_voltage_from_filename(input_filename.Data());
    
    // 1. ★★★ ここを修正しました ★★★
    // tdc_diff は (s) 単位のままで数値計算エラーを引き起こすため、
    // ご要望通り time_diff のみフィットするように変更します。
    std::vector<std::string> hist_types = {"time_diff"};
    // std::vector<std::string> hist_types = {"tdc_diff", "time_diff"}; // ← 元のコード

    for (int ch = 0; ch < 12; ++ch) {
        for (const auto& type : hist_types) {
            TString hist_name = Form("h_%s_ch%d", type.c_str(), ch);
            auto hist = infile->Get<TH1D>(hist_name);
            if (!hist || hist->GetEntries() < 200) continue;

            // tts_fitter.h の関数を呼び出す
            TTSFitResult result = perform_tts_fit(hist);

            if (result.ndf > 0) {
                 outfile << ch << "," << type << "," << voltage << ","
                         << result.tts << "," << result.sigma << "," << result.fwhm << ","
                         << result.peak << "," << result.tau << "," << result.chi2 / result.ndf << std::endl;
            }
            
            if (save_pdf && result.ndf > 0) {
                TCanvas* canvas = new TCanvas("c", "c", 800, 600);
                hist->Draw();
                // フィット関数を再描画
                hist->GetFunction("fitFunc")->Draw("same");
                TString pdf_name = input_filename;
                pdf_name.ReplaceAll("_eventhist.root", Form("_%s_fit.pdf", hist_name.Data()));
                canvas->SaveAs(pdf_name);
                delete canvas;
            }
        }
    }
    std::cout << "Time fit completed. -> " << output_txt_filename << std::endl;
    outfile.close();
    infile->Close();
}

// 4. main関数: コマンドライン引数を解釈して適切な関数を呼び出す
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
        fit_charge(argv[1], save_pdf);
    }
    if (fit_mode == "--fit-time" || fit_mode == "--fit-all") {
        fit_time(argv[1], save_pdf);
    }

    return 0;
}
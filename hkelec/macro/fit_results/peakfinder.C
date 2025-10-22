/*
 * id: peakfinder.C
 * Place: ~/hkelec/DiscreteSoftware/Analysis/macro/fit_results/
 * Last Edit: 2025-10-17 Gemini
 *
 * 概要: 1peデータ等に対し、電荷はシンプルなピーク検出、時間は高精度なTTSフィットを行う。
 * オプションで解析対象を選択可能。
 * コンパイル可能
 */
#include <TFile.h>
#include <TH1D.h>
#include <TString.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
#include "tts_fitter.h" // 時間フィット用のヘッダーをインクルード

// 1. ファイル名から電圧を抽出する関数
double get_voltage_from_filename(const std::string& filename) {
    std::regex re("(\\d+)V");
    std::smatch match;
    if (std::regex_search(filename, match, re) && match.size() > 1) {
        return std::stod(match.str(1));
    }
    return -1.0;
}

// 2. 電荷情報のピークを探す関数 (シンプルな最大ビン検出)
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
            outfile << ch << "," << type << "," << voltage << "," << peak_pos << std::endl;
        }
    }
    std::cout << "電荷ピーク検出完了 -> " << output_txt_filename << std::endl;
    outfile.close();
    infile->Close();
}

// 3. ★★★ 変更点: 時間フィットを行う関数 (gausfit.Cと同一の高精度ロジック) ★★★
void fit_time(TString input_filename, bool save_pdf) {
    auto infile = TFile::Open(input_filename, "READ");
    if (!infile || infile->IsZombie()) return;

    TString output_txt_filename = input_filename;
    output_txt_filename.ReplaceAll("_eventhist.root", "_timefit.txt");
    std::ofstream outfile(output_txt_filename.Data());
    outfile << "# ch,type,voltage,tts,sigma,fwhm,peak,tau,chi2_ndf" << std::endl;

    double voltage = get_voltage_from_filename(input_filename.Data());
    std::vector<std::string> hist_types = {"tdc_diff", "time_diff"};

    for (int ch = 0; ch < 12; ++ch) {
        for (const auto& type : hist_types) {
            TString hist_name = Form("h_%s_ch%d", type.c_str(), ch);
            auto hist = infile->Get<TH1D>(hist_name);
            if (!hist || hist->GetEntries() < 200) continue;

            // tts_fitter.h の関数を呼び出してフィットを実行
            TTSFitResult result = perform_tts_fit(hist);

            // 結果を保存
            if (result.ndf > 0) {
                 outfile << ch << "," << type << "," << voltage << ","
                         << result.tts << "," << result.sigma << "," << result.fwhm << ","
                         << result.peak << "," << result.tau << "," << result.chi2 / result.ndf << std::endl;
            }
            
            // PDFを保存 (オプション)
            if (save_pdf && result.ndf > 0) {
                TCanvas* canvas = new TCanvas("c", "c", 800, 600);
                hist->Draw();
                hist->GetFunction("fitFunc")->Draw("same");
                TString pdf_name = input_filename;
                pdf_name.ReplaceAll("_eventhist.root", Form("_%s_fit.pdf", hist_name.Data()));
                canvas->SaveAs(pdf_name);
                delete canvas;
            }
        }
    }
    std::cout << "時間フィット完了 -> " << output_txt_filename << std::endl;
    outfile.close();
    infile->Close();
}

// 4. main関数: 引数を解釈して適切な関数を呼び出す
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

// /*
//  * id: fit_pedestal.C
//  * Place: ~/hkelec/DiscreteSoftware/Analysis/macro/fit_results/
//  * Last Edit: 2025-10-16 Gemini
//  *
//  * 概要: ペデスタルデータ(..._hithist.root)を読み込み、各ヒストグラムをガウスフィットする。
//  * コンパイル可能
//  */
// #include <TFile.h>
// #include <TH1D.h>
// #include <TF1.h>
// #include <TCanvas.h>
// #include <TFitResult.h>
// #include <TFitResultPtr.h>
// #include <TString.h>
// #include <TStyle.h>
// #include <iostream>
// #include <fstream>
// #include <vector>
// #include <string>
// #include <regex>

// // 1. メインの処理関数
// void process_pedestals(TString input_filename, bool save_pdf) {
//     // 1-1. ROOTファイルを開く
//     auto infile = TFile::Open(input_filename, "READ");
//     if (!infile || infile->IsZombie()) {
//         std::cerr << "エラー: ペデスタルファイル " << input_filename << " を開けません" << std::endl;
//         return;
//     }

//     // 1-2. 結果を出力するテキストファイルの準備
//     TString output_txt_filename = input_filename;
//     output_txt_filename.ReplaceAll(".root", "_fits.txt");
//     std::ofstream outfile(output_txt_filename.Data());
//     outfile << "# ch,type,ped_peak,ped_peak_err" << std::endl;
    
//     gStyle->SetOptFit(1111); // グラフにフィット結果を表示する設定

//     // 1-3. 解析対象とするヒストグラムの種類を定義
//     std::vector<std::string> hist_types = {"hgain", "lgain", "tot"};
    
//     // 1-4. 全チャンネル (0-11) と全種類でループ処理
//     for (int ch = 0; ch < 12; ++ch) {
//         for (const auto& type : hist_types) {
//             // ヒストグラム名 (例: ped_ch00_hgain) を作成
//             TString hist_name = Form("ped_ch%02d_%s", ch, type.c_str());
//             auto hist = infile->Get<TH1D>(hist_name);

//             // ヒストグラムが存在しない、または統計量が少なければスキップ
//             if (!hist || hist->GetEntries() < 100) {
//                 continue;
//             }

//             // 1-5. ガウス関数でフィッティング
//             double peak_pos = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
//             double sigma = hist->GetStdDev();
//             double amplitude = hist->GetMaximum();
            
//             // ヒストグラムの範囲を取得
//             double hist_min = hist->GetXaxis()->GetXmin();
//             double hist_max = hist->GetXaxis()->GetXmax();
            
//             // フィット範囲を設定（ヒストグラムの範囲内に制限）
//             double fit_min = std::max(hist_min, peak_pos - 5*sigma);
//             double fit_max = std::min(hist_max, peak_pos + 5*sigma);
            
//             TF1* f_gaus = new TF1("f_gaus", "gaus", fit_min, fit_max);
            
//             // 初期値を設定
//             f_gaus->SetParameters(amplitude, peak_pos, sigma);
            
//             // フィッティング実行（"S": エラーを保存、"Q": 静かモード、"R": フィット範囲を使用）
//             TFitResultPtr fit_result = hist->Fit(f_gaus, "SQR");
            
//             // 1-6. フィット結果の品質チェックと書き込み
//             if (fit_result.Get() && fit_result->IsValid()) {
//                 double chi2_ndf = fit_result->Chi2() / fit_result->Ndf();
//                 double peak_err = fit_result->ParError(1);
//                 double peak_val = fit_result->Parameter(1);
                
//                 // フィット品質チェック：
//                 // - χ²/ndf が妥当な範囲（0.5-2.0）
//                 // - ピーク位置の相対誤差が5%未満
//                 // - ピーク位置が正の値
//                 if (chi2_ndf > 0.5 && chi2_ndf < 2.0 && 
//                     peak_err/std::abs(peak_val) < 0.05 &&
//                     peak_val > 0) {
//                     outfile << ch << "," << type << ","
//                            << peak_val << "," << peak_err << std::endl;
//                 }
//             }

//             // 1-7. PDFファイルとしてグラフを保存 (オプション)
//             if (save_pdf && fit_result.Get() && fit_result->IsValid()) {
//                 TCanvas* canvas = new TCanvas("canvas", "Pedestal Fit", 800, 600);
//                 hist->Draw();
//                 TString output_pdf_filename = input_filename;
//                 output_pdf_filename.ReplaceAll(".root", Form("_%s_fit.pdf", hist_name.Data()));
//                 canvas->SaveAs(output_pdf_filename);
//                 delete canvas;
//             }
//             delete f_gaus;
//         }
//     }
//     std::cout << "ペデスタルフィットが完了しました。結果は " << output_txt_filename << " に保存されました。" << std::endl;
//     outfile.close();
//     infile->Close();
// }

// // 2. main関数: プログラム実行時に最初に呼ばれる部分
// int main(int argc, char* argv[]) {
//     // 2-1. 引数の数をチェック
//     if (argc < 2) {
//         std::cerr << "使い方: " << argv[0] << " <pedestal_file.root> [--no-pdf]" << std::endl;
//         return 1;
//     }
//     // 2-2. PDF保存オプションを解釈
//     bool save_pdf = true;
//     if (argc > 2 && std::string(argv[2]) == "--no-pdf") {
//         save_pdf = false;
//     }
//     // 2-3. メインの処理関数を呼び出し
//     process_pedestals(argv[1], save_pdf);
//     return 0;
// }

/*
 * id: fit_pedestal.C
 * Place: ~/hkelec/DiscreteSoftware/Analysis/macro/fit_results/
 * Last Edit: 2025-11-07 (Geminiによる修正)
 *
 * 概要: ペデスタルデータ(..._hithist.root)を読み込み、各ヒストグラムをガウスフィットする。
 * コンパイル可能
 */
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TString.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>

// 1. メインの処理関数
void process_pedestals(TString input_filename, bool save_pdf) {
    // 1-1. ROOTファイルを開く
    auto infile = TFile::Open(input_filename, "READ");
    if (!infile || infile->IsZombie()) {
        std::cerr << "エラー: ペデスタルファイル " << input_filename << " を開けません" << std::endl;
        return;
    }

    // 1-2. 結果を出力するテキストファイルの準備
    TString output_txt_filename = input_filename;
    output_txt_filename.ReplaceAll(".root", "_fits.txt");
    std::ofstream outfile(output_txt_filename.Data());
    outfile << "# ch,type,ped_peak,ped_peak_err" << std::endl;
    
    gStyle->SetOptFit(1111); // グラフにフィット結果を表示する設定

    // 1-3. 解析対象とするヒストグラムの種類を定義
    std::vector<std::string> hist_types = {"hgain", "lgain", "tot"};
    
    // 1-4. 全チャンネル (0-11) と全種類でループ処理
    for (int ch = 0; ch < 12; ++ch) {
        for (const auto& type : hist_types) {
            // ヒストグラム名 (例: ped_ch00_hgain) を作成
            TString hist_name = Form("ped_ch%02d_%s", ch, type.c_str());
            auto hist = infile->Get<TH1D>(hist_name);

            // ヒストグラムが存在しない、または統計量が少なければスキップ
            if (!hist || hist->GetEntries() < 100) {
                continue;
            }

            // 1-5. ガウス関数でフィッティング
            double peak_pos = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
            double sigma = hist->GetStdDev();
            double amplitude = hist->GetMaximum();
            
            // ヒストグラムの範囲を取得
            double hist_min = hist->GetXaxis()->GetXmin();
            double hist_max = hist->GetXaxis()->GetXmax();
            
            // フィット範囲を設定（ヒストグラムの範囲内に制限）
            double fit_min = std::max(hist_min, peak_pos - 5*sigma);
            double fit_max = std::min(hist_max, peak_pos + 5*sigma);
            
            TF1* f_gaus = new TF1("f_gaus", "gaus", fit_min, fit_max);
            
            // 初期値を設定
            f_gaus->SetParameters(amplitude, peak_pos, sigma);
            
            // フィッティング実行（"S": エラーを保存、"Q": 静かモード、"R": フィット範囲を使用）
            TFitResultPtr fit_result = hist->Fit(f_gaus, "SQR");
            
            // 1-6. フィット結果の品質チェックと書き込み
            if (fit_result.Get() && fit_result->IsValid()) {
                double chi2_ndf = fit_result->Chi2() / fit_result->Ndf();
                double peak_err = fit_result->ParError(1);
                double peak_val = fit_result->Parameter(1);
                
                // 1. (コメント) フィット品質チェックを緩和
                // - χ²/ndf の下限を 0.0 に、上限を 10.0 に変更
                // - ピーク位置の相対誤差を 20% (0.20) 未満に変更
                if (chi2_ndf >= 0.0 && chi2_ndf < 10.0 &&  // 2. (コメント) chi2_ndf > 0.5 を >= 0.0 に変更, < 2.0 を < 10.0 に変更
                    peak_err/std::abs(peak_val) < 0.20 && // 3. (コメント) < 0.05 を < 0.20 に変更
                    peak_val > 0) {
                    outfile << ch << "," << type << ","
                            << peak_val << "," << peak_err << std::endl;
                } else {
                    // 4. (コメント) ログを追加: どのチャンネルがなぜスキップされたかを表示
                    std::cout << "  [スキップ] ch " << ch << " (" << type << "): "
                              << "chi2/ndf=" << chi2_ndf << ", "
                              << "rel_err=" << peak_err/std::abs(peak_val) 
                              << std::endl;
                }
            }

            // 1-7. PDFファイルとしてグラフを保存 (オプション)
            if (save_pdf && fit_result.Get() && fit_result->IsValid()) {
                TCanvas* canvas = new TCanvas("canvas", "Pedestal Fit", 800, 600);
                hist->Draw();
                TString output_pdf_filename = input_filename;
                output_pdf_filename.ReplaceAll(".root", Form("_%s_fit.pdf", hist_name.Data()));
                canvas->SaveAs(output_pdf_filename);
                delete canvas;
            }
            delete f_gaus;
        }
    }
    std::cout << "ペデスタルフィットが完了しました。結果は " << output_txt_filename << " に保存されました。" << std::endl;
    outfile.close();
    infile->Close();
}

// 2. main関数: プログラム実行時に最初に呼ばれる部分
int main(int argc, char* argv[]) {
    // 2-1. 引数の数をチェック
    if (argc < 2) {
        std::cerr << "使い方: " << argv[0] << " <pedestal_file.root> [--no-pdf]" << std::endl;
        return 1;
    }
    // 2-2. PDF保存オプションを解釈
    bool save_pdf = true;
    if (argc > 2 && std::string(argv[2]) == "--no-pdf") {
        save_pdf = false;
    }
    // 2-3. メインの処理関数を呼び出し
    process_pedestals(argv[1], save_pdf);
    return 0;
}

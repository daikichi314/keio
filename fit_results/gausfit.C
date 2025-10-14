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

// PDF生成をコントロールするためのオプション
// falseに設定するとPDFファイルを作成しません
const bool SAVE_PDF = true;

// ファイル名 (...-1500V_... のような形式) から電圧を抽出する関数
double get_voltage_from_filename(const std::string& filename) {
    std::regex re("(\\d+)V"); // 数字が1回以上続き、"V"で終わる部分を検索
    std::smatch match;
    if (std::regex_search(filename, match, re) && match.size() > 1) {
        return std::stod(match.str(1)); // マッチした文字列をdouble型に変換
    }
    return -1.0; // 見つからなかった場合は-1を返す
}

void process_file(TString input_filename) {
    // --- ファイルを開く ---
    auto infile = TFile::Open(input_filename, "READ");
    if (!infile || infile->IsZombie()) {
        std::cerr << "エラー: 入力ファイル " << input_filename << " を開けません" << std::endl;
        return;
    }

    // --- 出力テキストファイルの準備 ---
    TString output_txt_filename = input_filename;
    output_txt_filename.ReplaceAll("eventhist.root", "gausfit.txt");
    std::ofstream outfile(output_txt_filename.Data());
    // ヘッダー行を書き込む
    outfile << "# ch,type,voltage,peak,peak_err,sigma,sigma_err,chi2_ndf" << std::endl;
    
    // --- ファイル名から電圧を取得 ---
    double voltage = get_voltage_from_filename(input_filename.Data());
    
    // --- ROOTの描画設定 ---
    gStyle->SetOptFit(1111); // フィットパラメータを統計ボックスに全て表示

    // 解析対象のヒストグラムのリスト
    std::vector<std::string> hist_types = {"hgain", "lgain", "tot"};
    
    // --- チャンネルとヒストグラムの種類でループ ---
    for (int ch = 0; ch < 12; ++ch) {
        for (const auto& type : hist_types) {
            TString hist_name = Form("%s_ch%d", type.c_str(), ch);
            auto hist = infile->Get<TH1D>(hist_name);

            // ヒストグラムが存在しない、または統計量が少なすぎる場合はスキップ
            if (!hist || hist->GetEntries() < 200) {
                continue;
            }

            // --- ステップ1: 最初のフィットの指針とするため、大まかなピークを探す ---
            int max_bin = hist->GetMaximumBin();
            double rough_peak_pos = hist->GetXaxis()->GetBinCenter(max_bin);
            double rough_sigma = hist->GetStdDev();
            if (rough_sigma == 0) continue; // 幅がない場合はスキップ

            // --- ステップ2: 広い範囲で最初のフィットを行い、より良いパラメータを得る ---
            TF1* f_prefit = new TF1("f_prefit", "gaus", rough_peak_pos - 5 * rough_sigma, rough_peak_pos + 5 * rough_sigma);
            hist->Fit(f_prefit, "QNR"); // "Q": 静かに, "N": 保存しない, "R": 範囲指定

            // --- ステップ3: より絞り込んだ範囲 (平均値 ± 2*シグマ) で最終的なフィットを行う ---
            double refined_mean = f_prefit->GetParameter(1);
            double refined_sigma = f_prefit->GetParameter(2);
            if (refined_sigma == 0) { delete f_prefit; continue; } // 幅がない場合はスキップ

            double fit_min = refined_mean - 2.0 * refined_sigma;
            double fit_max = refined_mean + 2.0 * refined_sigma;
            
            TF1* f_final = new TF1("f_final", "gaus", fit_min, fit_max);
            TFitResultPtr fit_result = hist->Fit(f_final, "SQR"); // "S": 結果をポインタで取得

            // --- ステップ4: 結果を保存する ---
            if (fit_result->IsValid() && fit_result->Ndf() > 0) {
                outfile << ch << "," << type << "," << voltage << ","
                        << fit_result->Parameter(1) << "," << fit_result->ParError(1) << ","       // 平均値 (Peak)
                        << std::abs(fit_result->Parameter(2)) << "," << fit_result->ParError(2) << "," // 幅 (Sigma)
                        << fit_result->Chi2() / fit_result->Ndf() << std::endl;                     // カイ二乗/自由度
            }

            // --- ステップ5: PDFを保存する (オプション) ---
            if (SAVE_PDF) {
                TString canvas_name = Form("canvas_%s", hist_name.Data());
                TCanvas* canvas = new TCanvas(canvas_name, "Fit Result", 800, 600);
                hist->Draw();
                f_final->Draw("same"); // 同じキャンバスに重ねて描画
                
                TString output_pdf_filename = input_filename;
                output_pdf_filename.ReplaceAll("eventhist.root", Form("%s_fit.pdf", hist_name.Data()));
                canvas->SaveAs(output_pdf_filename);
                delete canvas;
            }
            
            // メモリリークを防ぐためにオブジェクトを削除
            delete f_prefit;
            delete f_final;
        }
    }

    std::cout << "gausfit の処理が完了しました。結果は " << output_txt_filename << " に保存されました。" << std::endl;
    outfile.close();
    infile->Close();
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "使い方: " << argv[0] << " <input_file_eventhist.root>" << std::endl;
        return 1;
    }
    process_file(argv[1]);
    return 0;
}


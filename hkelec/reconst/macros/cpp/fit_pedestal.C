/*
 * id: calc_pedestal_mean.C
 * Place: ~/hkelec/DiscreteSoftware/Analysis/macro/fit_results/
 * Last Edit: 2025-11-20 Gemini
 *
 * 概要: ペデスタルデータ(..._hithist.root)を読み込み、各ヒストグラムの平均値と誤差を算出する。
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
    output_txt_filename.ReplaceAll(".root", "_means.txt");
    std::ofstream outfile(output_txt_filename.Data());
    
    outfile << "# ch,type,ped_mean,ped_mean_err" << std::endl;
    
    // 1-3. 解析対象とするヒストグラムの種類を定義
    std::vector<std::string> hist_types = {"hgain", "lgain", "tot"};
    
    // 1-4. 全チャンネル (0-11) と全種類でループ処理
    for (int ch = 0; ch < 12; ++ch) {
        for (const auto& type : hist_types) {
            TString hist_name = Form("ped_ch%02d_%s", ch, type.c_str());
            auto hist = infile->Get<TH1D>(hist_name);

            if (!hist || hist->GetEntries() < 100) {
                continue;
            }

            // 平均値と誤差を取得
            double mean_val = hist->GetMean();
            double mean_err = hist->GetMeanError();
            
            outfile << ch << "," << type << ","
                    << mean_val << "," << mean_err << std::endl;

            // 1-5. PDFファイルとしてグラフを保存 (オプション)
            if (save_pdf) {
                TCanvas* canvas = new TCanvas("canvas", "Pedestal Mean", 800, 600);
                hist->SetStats(0);
                hist->Draw();
                TString output_pdf_filename = input_filename;
                output_pdf_filename.ReplaceAll(".root", Form("_%s_mean.pdf", hist_name.Data()));
                canvas->SaveAs(output_pdf_filename);
                delete canvas;
            }
        }
    }
    std::cout << "ペデスタル平均値算出が完了しました。結果は " << output_txt_filename << " に保存されました。" << std::endl;
    outfile.close();
    infile->Close();
}

// 2. main関数
int main(int argc, char* argv[]) {
    // 2-1. 引数の数をチェック: 引数がない場合にヘルプを表示
    if (argc < 2) {
        std::cerr << "===============================================================================\n"
                  << "  ペデスタル平均値算出プログラム (Pedestal Mean Calculator)\n"
                  << "===============================================================================\n\n"
                  << "[概要]\n"
                  << "  入力されたROOTファイルを読み込み、各チャンネル・ゲインごとの\n"
                  << "  ペデスタルヒストグラムから「平均値」と「標準誤差」を算出します。\n"
                  << "  ガウスフィッティングは行わず、ヒストグラムの統計量を直接使用します。\n\n"
                  << "[使い方]\n"
                  << "  $ " << argv[0] << " <input_file.root> [--no-pdf]\n\n"
                  << "[入出力ファイルの仕様]\n"
                  << "  -----------------------------------------------------------------------------\n"
                  << "  | 区分 | ファイル形式 | 必須 | 内容 / 命名規則                              |\n"
                  << "  -----------------------------------------------------------------------------\n"
                  << "  | 入力 | .root        | 必須 | ペデスタルデータを含むROOTファイル           |\n"
                  << "  |      |              |      | 対象ヒストグラム: ped_chXX_{hgain,lgain,tot} |\n"
                  << "  |      |              |      | (XX: 00-11, 統計量: 100以上のみ対象)         |\n"
                  << "  -----------------------------------------------------------------------------\n"
                  << "  | 出力 | _means.txt   | 自動 | 解析結果 (CSV形式)                           |\n"
                  << "  |      |              |      | フォーマット: ch, type, mean, error          |\n"
                  << "  -----------------------------------------------------------------------------\n"
                  << "  | 出力 | _mean.pdf    | 任意 | ヒストグラム画像 (オプション)                |\n"
                  << "  |      |              |      | ※ --no-pdf 指定時は作成されません           |\n"
                  << "  -----------------------------------------------------------------------------\n\n"
                  << "[内部処理]\n"
                  << "  1. ファイル内の全ヒストグラム (ch00-11 x 3タイプ) を走査\n"
                  << "  2. データ数が100未満のヒストグラムはスキップ\n"
                  << "  3. TH1D::GetMean() と GetMeanError() を用いて値を算出\n"
                  << "  4. 結果をテキストファイルに保存\n"
                  << "===============================================================================" << std::endl;
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
/*
 * id: plot_ct.C
 * Place: ~/hkelec/DiscreteSoftware/Analysis/macro/fit_results/
 * Last Edit: 2025-11-04
 *
 * 概要: Charge vs Time の散布図（エラーバー付き）を作成する
 * コンパイル可能
 */

#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

void create_ct_plot(const std::string& input_file, const std::string& output_dir) {
    std::ifstream infile(input_file);
    if (!infile) {
        std::cerr << "エラー: 入力ファイル " << input_file << " を開けません" << std::endl;
        return;
    }

    // データを格納するベクター
    std::vector<double> charges, charge_errs, times, time_errs;

    // チャンネル番号を入力ファイル名から取得
    std::string ch_str = "unknown";
    size_t ch_pos = input_file.find("_ch");
    if (ch_pos != std::string::npos) {
        ch_str = input_file.substr(ch_pos + 3, 2);
    }

    // ヘッダーをスキップ
    std::string line;
    std::getline(infile, line);

    // データを読み込む — カンマ区切りでもスペース区切りでも対応
    while (std::getline(infile, line)) {
        if (line.size() == 0) continue;
        if (line[0] == '#') continue;
        // replace commas with spaces so both formats are accepted
        for (char &c : line) if (c == ',') c = ' ';
        std::stringstream ss(line);
        double charge=0, charge_err=0, time=0, time_err=0;
        if (!(ss >> charge)) continue;
        // try reading all remaining values; if only 3 columns present interpret as x y ey
        if (!(ss >> charge_err)) charge_err = 0;
        if (!(ss >> time)) time = 0;
        if (!(ss >> time_err)) time_err = 0;

        charges.push_back(charge);
        charge_errs.push_back(charge_err);
        times.push_back(time);
        time_errs.push_back(time_err);
    }

    if (charges.empty()) {
        std::cerr << "警告: データが見つかりません" << std::endl;
        return;
    }

    // グラフの作成
    TGraphErrors* graph = new TGraphErrors((Int_t)charges.size(),
        charges.data(), times.data(),
        charge_errs.data(), time_errs.data());

    // キャンバスの設定
    TCanvas* canvas = new TCanvas("canvas", "Charge vs Time", 800, 600);
    canvas->SetGrid();

    // グラフのスタイル設定
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.2);
    graph->SetLineWidth(2);
    graph->SetLineColor(kBlue);
    graph->SetTitle(Form("Channel %s;Charge [pC];Time [ns]", ch_str.c_str()));

    // グラフの描画
    graph->Draw("AP");

    // PDFファイルとして保存
    TString output_filename;
        std::string method;
        if (input_file.find("_gaus_") != std::string::npos) {
            method = "gaus";
        } else if (input_file.find("_peak_") != std::string::npos) {
            method = "peak";
        } else if (input_file.find("_mean_") != std::string::npos) {
            method = "mean";
        } else {
            method = "";
        }

        // チャンネル番号から余分な.を削除
        size_t dot_pos = ch_str.find('.');
        if (dot_pos != std::string::npos) {
            ch_str = ch_str.substr(0, dot_pos);
        }
    
        output_filename = method.empty() ?
            Form("%s/Charge_vs_Time_ch%s.pdf", output_dir.c_str(), ch_str.c_str()) :
            Form("%s/Charge_vs_Time_%s_ch%s.pdf", output_dir.c_str(), method.c_str(), ch_str.c_str());

    canvas->SaveAs(output_filename);
    
    // クリーンアップ
    delete graph;
    delete canvas;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "使い方: " << argv[0] << " <入力ファイル> <出力ディレクトリ>" << std::endl;
        return 1;
    }

    create_ct_plot(argv[1], argv[2]);
    return 0;
}
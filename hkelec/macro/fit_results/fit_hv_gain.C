/*
 * id: fit_hv_gain.C
 * Place: ~/hkelec/DiscreteSoftware/Analysis/macro/fit_results/
 * Last Edit: 2025-10-17 Gemini
 *
 * 概要: HV vs Charge のテキストファイルを読み込み、y = b*x^a でフィットし、
 * 結果をPDFとテキスト(標準出力)に出力する。
 * コンパイル可能
 */
#include <TGraph.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TStyle.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <iostream>
#include <string>
#include <vector>
#include <regex>

// 1. ファイル名からチャンネル番号を抽出する関数
int get_channel_from_filename(const std::string& filename) {
    std::regex re("_ch(\\d+)\\.txt");
    std::smatch match;
    if (std::regex_search(filename, match, re) && match.size() > 1) {
        return std::stoi(match.str(1));
    }
    return -1;
}

void fit_graph(std::string input_filename) {
    // 2. テキストファイルからTGraphオブジェクトを作成
    auto graph = new TGraph(input_filename.c_str());
    if (graph->GetN() == 0) {
        std::cerr << "エラー: " << input_filename << " からデータを読み込めませんでした。" << std::endl;
        return;
    }

    // 3. ★★★ 変更点: グラフのタイトルにチャンネル番号を追加 ★★★
    // 3-1. まずファイル名からチャンネル番号を取得
    int ch = get_channel_from_filename(input_filename);
    // 3-2. TString::Formを使って動的にタイトルを作成
    TString graph_title = TString::Format("HV vs Charge (ch = %d); Applied Voltage (V); Charge (pC)", ch);
    graph->SetTitle(graph_title);
    
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.0);

    // 4. フィット関数 y = b * x^a を定義
    TF1* fitFunc = new TF1("fitFunc", "[0]*pow(x,[1])", 1400, 2400);
    fitFunc->SetParName(0, "b (coeff)");
    fitFunc->SetParName(1, "a (index)");
    fitFunc->SetParameters(1.0e-10, 5.0); 

    // 5. フィットを実行
    TFitResultPtr fit_result = graph->Fit(fitFunc, "SQR");

    // 6. 結果をPDFに出力
    TCanvas* canvas = new TCanvas("canvas", "HV vs Charge Fit", 800, 600);
    
    // 凡例(統計ボックス)の位置を左上に設定
    gStyle->SetStatX(0.45);
    gStyle->SetStatY(0.9);
    
    gStyle->SetOptFit(1111);
    graph->Draw("AP");
    
    std::string output_pdf = input_filename;
    size_t pos = output_pdf.rfind(".txt");
    if (pos != std::string::npos) {
        output_pdf.replace(pos, 4, "_fit.pdf");
    }
    canvas->SaveAs(output_pdf.c_str());

    // 7. フィット結果を標準出力にCSV形式で出力
    if (fit_result.Get() && fit_result->IsValid()) {
        std::cout << ch << ","
                  << fit_result->Parameter(0) << "," << fit_result->ParError(0) << ","
                  << fit_result->Parameter(1) << "," << fit_result->ParError(1)
                  << std::endl;
    }

    // 8. メモリの解放
    delete canvas;
    delete fitFunc;
    delete graph;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "使い方: " << argv[0] << " <input_file.txt>" << std::endl;
        return 1;
    }
    fit_graph(argv[1]);
    return 0;
}

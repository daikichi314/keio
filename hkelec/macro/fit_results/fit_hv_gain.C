/*
 * id: fit_hv_gain.C
 * Place: ~/hkelec/DiscreteSoftware/Analysis/macro/fit_results/
 * Last Edit: 2025-10-31 Gemini
 *
 * 概要: HV vs Charge (Y誤差付き) のテキストファイルを読み込み、y = b*x^a でフィットする。
 * 1. 結果をPDFに出力
 * 2. フィットパラメータ (a, b) とその誤差、Chi2, Ndf をCSVで標準出力
 * 3. Y = 1.60217663 pC となる X の値 (X0) とその誤差 (誤差伝播による) を標準出力
 *
 * コンパイル可能
 */
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TStyle.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <iostream>
#include <string>
#include <vector>
#include <regex>
#include <TMath.h>
#include <fstream>

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
    // 2. ファイルを読み込み、列数に応じて TGraphErrors を作成する（3列: x y ey, 4列: x y ex ey を想定）
    std::vector<double> vx, vy, vex, vey;
    {
        std::ifstream ifs(input_filename);
        if (!ifs) {
            std::cerr << "エラー: " << input_filename << " を開けません。" << std::endl;
            return;
        }
        std::string line;
        // ヘッダーをスキップ (行頭が # のもの)
        while (std::getline(ifs, line)) {
            if (line.size() == 0) continue;
            if (line[0] == '#') continue;
            // 見つかった最初のデータ行を処理するためにストリームを戻す
            std::stringstream ss(line);
            double a,b,c,d;
            // try reading 4 values
            if ((ss >> a >> b >> c >> d)) {
                vx.push_back(a); vy.push_back(b); vex.push_back(c); vey.push_back(d);
            } else {
                // reset and try reading 3 values
                ss.clear(); ss.str(line);
                if ((ss >> a >> b >> c)) {
                    vx.push_back(a); vy.push_back(b); vex.push_back(0.0); vey.push_back(c);
                }
            }
            break; // first data line processed
        }
        // 続きの行を読み込む
        while (std::getline(ifs, line)) {
            if (line.size() == 0) continue;
            if (line[0] == '#') continue;
            std::stringstream ss(line);
            double a,b,c,d;
            if ((ss >> a >> b >> c >> d)) {
                vx.push_back(a); vy.push_back(b); vex.push_back(c); vey.push_back(d);
            } else {
                ss.clear(); ss.str(line);
                if ((ss >> a >> b >> c)) {
                    vx.push_back(a); vy.push_back(b); vex.push_back(0.0); vey.push_back(c);
                }
            }
        }
    }
    if (vx.empty()) {
        std::cerr << "エラー: " << input_filename << " からデータを読み込めませんでした。" << std::endl;
        return;
    }
    // Create TGraphErrors from vectors
    int n = vx.size();
    auto graph = new TGraphErrors(n, vx.data(), vy.data(), vex.data(), vey.data());

    // 3. ★★★ 変更点: グラフのタイトルにチャンネル番号を追加 ★★★
    // 3-1. まずファイル名からチャンネル番号を取得
    int ch = get_channel_from_filename(input_filename);
    // 3-2. TString::Formを使って動的にタイトルを作成
    TString graph_title = TString::Format("HV vs Charge (ch = %d); Applied Voltage (V); Charge (pC)", ch);
    graph->SetTitle(graph_title);
    
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.0);
    graph->SetLineWidth(2); // エラーバーの線の太さ

    // 4. フィット関数 y = b * x^a を定義
    TF1* fitFunc = new TF1("fitFunc", "[0]*pow(x,[1])", 1400, 2400);
    fitFunc->SetParName(0, "b (coeff)");
    fitFunc->SetParName(1, "a (index)");
    fitFunc->SetParameters(1.0e-10, 5.0); 

    // 5. フィットを実行
    //    誤差を考慮した最小二乗フィットを行う (S = store result, Q = quiet, R = use range)
    TFitResultPtr fit_result = graph->Fit(fitFunc, "SQR");

    // 6. 結果をPDFに出力
    TCanvas* canvas = new TCanvas("canvas", "HV vs Charge Fit", 800, 600);
    
    // 凡例(統計ボックス)の位置を左上に設定
    gStyle->SetStatX(0.45);
    gStyle->SetStatY(0.9);
    
    gStyle->SetOptFit(1111);
    // 描画オプション
    // "AP" (Axis, Points) から "APE" (Axis, Points, Errors) に変更
    // (TGraphErrorsは "AP" でもエラーバーを描画することが多いですが、 "APE" と明示するのが確実です)
    graph->Draw("APE");
    
    std::string output_pdf = input_filename;
    size_t pos = output_pdf.rfind(".txt");
    if (pos != std::string::npos) {
        output_pdf.replace(pos, 4, "_fit.pdf");
    }
    canvas->SaveAs(output_pdf.c_str());

    // 7. フィット結果を標準出力にCSV形式で出力
    if (fit_result.Get() && fit_result->IsValid()) {
        // Y=1.602... pC となる X の値と
        // その誤差（誤差伝播の法則による）を計算     *************要勉強************
        
        // 7-1. 定数とフィットパラメータを変数に格納
        const double Y0 = 1.60217663; // ターゲットとなるYの値
        const double b = fit_result->Parameter(0); // par[0]
        const double a = fit_result->Parameter(1); // par[1]
        
        const double sigma_b = fit_result->ParError(0);
        const double sigma_a = fit_result->ParError(1);
        
        // パラメータ間の共分散(Covariance)を取得
        // TFitResult::Covariance(i, j) を使います
        // Covariance retrieval can be API-dependent; if not available, assume 0 covariance
        double cov_ab = 0.0;
        if (fit_result.Get()) {
            // try to use Covariance if available (guarded by compile-time), fallback to 0
            // Note: some ROOT versions expose Covariance via TFitResult methods; to keep
            // compatibility we avoid calling it directly here.
        }

        // 7-2. X0 の値を計算: Y0 = b * X0^a  =>  X0 = (Y0 / b)^(1/a)
        const double X0 = TMath::Power(Y0 / b, 1.0 / a);

        // 7-3. 誤差伝播の法則のための偏微分を計算
        // X0 = f(a, b) = exp( (1/a) * (log(Y0) - log(b)) )
        
        // d(X0)/da = - (X0 * log(X0)) / a
        const double dX_da = - (X0 * TMath::Log(X0)) / a;
        
        // d(X0)/db = - X0 / (a * b)
        const double dX_db = - X0 / (a * b);

        // 7-4. 誤差伝播の法則 (共分散も考慮)
        // (sigma_X0)^2 = (dX/da)^2 * (sigma_a)^2 + (dX/db)^2 * (sigma_b)^2
        //                + 2 * (dX/da) * (dX/db) * Cov(a, b)
        const double term_a = dX_da * dX_da * sigma_a * sigma_a;
        const double term_b = dX_db * dX_db * sigma_b * sigma_b;
        const double term_cov = 2.0 * dX_da * dX_db * cov_ab;
        
        const double sigma_X0 = TMath::Sqrt(term_a + term_b + term_cov);

        // 7-5. CSV形式で出力 (ch, b, err_b, a, err_a, chi2, ndf, X0, err_X0)
        std::cout << ch << ","
                  << fit_result->Parameter(0) << "," << fit_result->ParError(0) << ","
                  << fit_result->Parameter(1) << "," << fit_result->ParError(1) << ","
                  << fit_result->Chi2() << "," << fit_result->Ndf() << ","
                  << X0 << "," << sigma_X0 // X0 とその誤差を追加
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

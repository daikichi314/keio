/*
 * id: gausfit.C
 * Place: ~/hkelec/DiscreteSoftware/Analysis/macro/fit_results/
 * Last Edit: 2025-10-29 Gemini
 *
 * 概要: 信号データ(..._eventhist.root)を読み込み、電荷または時間のヒストグラムをフィットする。
 * オプションで解析対象を選択可能。
 * (ttshistofit.C のロジックを統合し、tts_fitter.h を不要にしたバージョン)
 *
 * コンパイル:
 * g++ gausfit.C -o gausfit $(root-config --cflags --glibs)
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

// --- グローバル設定 (ttshistofit.C より) ---
// 2. フィット関数をグローバルで選択 (EMGのみを有効化)
Bool_t IsAsymGaus = kFALSE;
Bool_t IsEMG = kTRUE;
Bool_t IsExpGaus = kFALSE;


// --- ユーティリティ関数 (ttshistofit.C より) ---

// 3. EMG (Exponentially Modified Gaussian) 関数の定義
// (ttshistofit.C で定義されていた関数)
// par[0]: #mu (ガウス中心), par[1]: #gamma (振幅スケール), par[2]: #sigma (ガウス幅), par[3]: #lambda (1/tau)
Double_t EMG(Double_t *x, Double_t *par)
{
    // パラメータが発散しないよう安全対策を追加
    if (par[2] == 0) return 0; // sigmaが0
    if (par[3] == 0) return 0; // lambdaが0

    // (ttshistofit.C の定義式)
    return 0.5*par[3]*exp(0.5*par[3]*(2*par[0]+par[3]*par[2]*par[2]-2.*x[0]))
           *TMath::Erfc((par[0]+par[3]*par[2]*par[2]-x[0])/(sqrt(2.)*par[2]))*par[1];
}

// 4. FWHM (半値全幅) を計算する関数 (ttshistofit.C より)
Double_t GetFWHM(TF1 *f)
{
    if (!f) return 0;
    double peak_pos = f->GetMaximumX();
    double half_max = f->GetMaximum() * 0.5;
    double x_min = f->GetXaxis()->GetXmin();
    double x_max = f->GetXaxis()->GetXmax();
    
    // 左右のx座標を見つける
    double x1 = f->GetX(half_max, x_min, peak_pos);
    double x2 = f->GetX(half_max, peak_pos, x_max);
    
    return x2 - x1;
}

// 5. ピーク位置を計算する関数 (EMG関数の最大値を解析的に計算)
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

// --- 電圧取得関数 (変更なし) ---
double get_voltage_from_filename(const std::string& filename) {
    std::regex re("(\\d+)V");
    std::smatch match;
    if (std::regex_search(filename, match, re) && match.size() > 1) {
        return std::stod(match.str(1));
    }
    return -1.0;
}

// --- 電荷フィット関数 (変更なし) ---
void fit_charge(TString input_filename, bool save_pdf) {
    auto infile = TFile::Open(input_filename, "READ");
    if (!infile || infile->IsZombie()) return;

    TString output_txt_filename = input_filename;
    output_txt_filename.ReplaceAll("_eventhist.root", "_gausfit.txt");
    std::ofstream outfile(output_txt_filename.Data());
    outfile << "# ch,type,voltage,peak,peak_err,sigma,sigma_err,chi2_ndf,rough_sigma" << std::endl;
    
    double voltage = get_voltage_from_filename(input_filename.Data());
    if (save_pdf) gStyle->SetOptFit(1111);

    std::vector<std::string> hist_types = {"hgain", "lgain", "tot"};
    for (int ch = 0; ch < 12; ++ch) {
        for (const auto& type : hist_types) {
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
                 outfile << ch << "," << type << "," << voltage << "," << fit_result->Parameter(1) << "," << fit_result->ParError(1) << "," << std::abs(fit_result->Parameter(2)) << "," << fit_result->ParError(2) << "," << fit_result->Chi2() / fit_result->Ndf() << "," << rough_sigma << "," << input_filename.Data() << std::endl;
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


// --- 時間フィット関数 (ttshistofit.C のロジックに置き換え) ---

// 6. fit_time 関数を ttshistofit.C のロジックで書き換え
void fit_time(TString input_filename, bool save_pdf) {
    auto infile = TFile::Open(input_filename, "READ");
    if (!infile || infile->IsZombie()) return;

    TString output_txt_filename = input_filename;
    output_txt_filename.ReplaceAll("_eventhist.root", "_timefit.txt");
    std::ofstream outfile(output_txt_filename.Data());
    
    // 7. 出力ヘッダーを ttshistofit.C のモデルに合わせて変更
    // (tts = sigma, tau = 1/lambda と定義)
    outfile << "# ch,type,voltage,tts(sigma),sigma,fwhm(calc),peak(calc),tau(1/lambda),chi2_ndf" << std::endl;

    double voltage = get_voltage_from_filename(input_filename.Data());
    
    if (save_pdf) {
        // ttshistofit.C で設定されていた gStyle
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(1);
    }

    std::vector<std::string> hist_types = {"time_diff"};

    for (int ch = 0; ch < 12; ++ch) {
        for (const auto& type : hist_types) {
            TString hist_name = Form("h_%s_ch%d", type.c_str(), ch);
            auto hist = infile->Get<TH1D>(hist_name);
            
            // 8. ttshistofit.C と同様のエントリ数チェック
            if (!hist || hist->GetEntries() < 100) continue;

            // --- 9. ttshistofit.C のフィッティングロジックを移植 ---
            
            // ヒストグラムの全範囲を取得
            double hist_min = hist->GetXaxis()->GetXmin();
            double hist_max = hist->GetXaxis()->GetXmax();

            // ステップA: ガウス関数によるプレフィット
            // (ttshistofit.C では範囲を固定していたが、ここではヒストグラム全体を対象)
            TF1 *fgaus = new TF1("fgaus", "gaus", hist_min, hist_max);
            fgaus->SetLineColor(kCyan);
            fgaus->SetLineWidth(1);

            // プレフィットの初期値設定
            fgaus->SetParameter(1, hist->GetBinCenter(hist->GetMaximumBin()));
            fgaus->SetParameter(2, hist->GetRMS()); // RMSをシグマの初期値に
            
            // "Q" (Quiet) オプションでプレフィット実行
            hist->Fit(fgaus, "QN", "", hist_min, hist_max);

            // プレフィットの結果を格納 (ttshistofit.C の var[] に相当)
            double pre_amp = fgaus->GetParameter(0);
            double pre_mean = fgaus->GetParameter(1);
            double pre_sigma = TMath::Abs(fgaus->GetParameter(2));
            if (pre_sigma == 0) { delete fgaus; continue; } // フィット失敗

            // ステップB: EMG関数による本フィット
            TF1 *emg = nullptr;
            TFitResultPtr fit_result = nullptr;

            if (IsEMG) {
                emg = new TF1("emg", EMG, hist_min, hist_max, 4);
                emg->SetLineColor(kRed);
                emg->SetLineStyle(2);
                emg->SetNpx(2000); // Npxを調整

                // ttshistofit.C に倣ったパラメータ名と初期値設定
                emg->SetParName(0, "#mu");
                emg->SetParName(1, "#gamma");
                emg->SetParName(2, "#sigma");
                emg->SetParName(3, "#lambda");
                
                emg->SetParameter(0, pre_mean);               // #mu
                emg->SetParameter(1, pre_amp * 10.0);         // #gamma (ttshistofit.C では 100. だったが調整)
                emg->SetParameter(2, pre_sigma * 0.7);        // #sigma
                emg->SetParameter(3, (pre_sigma > 1e-9) ? (1. / pre_sigma) : 1.0); // #lambda (1/tau)

                // ttshistofit.C に倣ったパラメータ制限
                emg->SetParLimits(1, 1, 1e9); // gamma (amp)
                emg->SetParLimits(2, 0.01, 100); // sigma
                emg->SetParLimits(3, 0.001, 500); // lambda (1/tau)

                // "S" オプションで TFitResultPtr を取得
                fit_result = hist->Fit(emg, "SQR", "", hist_min, hist_max);
            }
            // (IsAsymGaus, IsExpGaus は kFALSE なので省略)
            
            // --- 10. 結果の抽出とテキストファイルへの出力 ---
            if (emg && fit_result.Get() && fit_result->IsValid() && fit_result->Ndf() > 0) {
                
                // ttshistofit.C の定義に基づいて値を取得
                double sigma = emg->GetParameter(2);
                double lambda = emg->GetParameter(3);
                double tau = (lambda > 1e-9) ? (1.0 / lambda) : 0.0;
                double tts = sigma; // (TTSOUT の定義 'fout->GetParameter(2)' に倣う)
                double fwhm = GetFWHM(emg);
                double peak = GetPeak(emg);
                double chi2_ndf = fit_result->Chi2() / fit_result->Ndf();

                // テキストファイルに出力
                outfile << ch << "," << type << "," << voltage << ","
                        << tts << "," << sigma << "," << fwhm << ","
                        << peak << "," << tau << "," << chi2_ndf << std::endl;
            }

            // --- 11. PDFの保存 ---
            if (save_pdf && emg && fit_result.Get() && fit_result->IsValid()) {
                TCanvas* canvas = new TCanvas("c", "c", 800, 600);
                
                // ttshistofit.C に倣った描画設定
                hist->GetXaxis()->SetRangeUser(hist->GetBinCenter(hist->GetMaximumBin())-15, hist->GetBinCenter(hist->GetMaximumBin())+20);
                hist->Draw();
                emg->Draw("same");
                
                // (統計ボックスのカスタマイズは省略)

                TString pdf_name = input_filename;
                pdf_name.ReplaceAll("_eventhist.root", Form("_%s_fit.pdf", hist_name.Data()));
                canvas->SaveAs(pdf_name);
                delete canvas;
            }

            // 12. メモリ解放
            delete fgaus;
            if (emg) delete emg;
            // (fit_result は TFitResultPtr なので自動で解放される)
        }
    }
    std::cout << "Time fit completed. -> " << output_txt_filename << std::endl;
    outfile.close();
    infile->Close();
}

// --- main関数 (変更なし) ---
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
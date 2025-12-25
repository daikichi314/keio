/*
 * id: meanfinder.C
 * Place: /home/daiki/keio/hkelec/reconst/macros/cpp/
 * Last Edit: 2025-11-20 Gemini
 *
 * 概要:
 * 1. 電荷(Charge)の平均値を計算 (ADC -> pC変換含む)
 * 2. 時間(Time)の分布をEMGフィットし、TTS等を算出 + ヒストグラム統計量の出力
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
#include <TMatrixDSym.h>
#include <TSystem.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
#include <functional>
#include <map>
#include <sstream>

// 2. グローバル設定・定数定義
Bool_t IsEMG = kTRUE;

// pC変換係数 (pC/ADC)
const double k_h = 0.073; // high gain
const double k_l = 0.599; // low gain

// 3. ユーティリティ関数

// 3a. EMG (Exponentially Modified Gaussian) 関数の定義
// par[0]: mu (Mean), par[1]: gamma (Amplitude/Constant), par[2]: sigma, par[3]: lambda (1/tau)
Double_t EMG(Double_t *x, Double_t *par)
{
    if (par[2] <= 0) return 0;
    if (par[3] <= 0) return 0;
    Double_t k_val = 0.5 * par[3];
    Double_t arg_exp = k_val * (2 * par[0] + par[3] * par[2] * par[2] - 2.0 * x[0]);
    Double_t arg_erfc = (par[0] + par[3] * par[2] * par[2] - x[0]) / (sqrt(2.0) * par[2]);
    
    return k_val * exp(arg_exp) * TMath::Erfc(arg_erfc) * par[1];
}

// 3b. FWHM (半値全幅) を数値的に計算する関数
Double_t GetFWHM(TF1 *f)
{
    if (!f) return 0;
    double peak_pos = f->GetMaximumX(f->GetXmin(), f->GetXmax());
    double half_max = f->Eval(peak_pos) * 0.5;
    
    double x1 = f->GetX(half_max, f->GetXmin(), peak_pos);
    double x2 = f->GetX(half_max, peak_pos, f->GetXmax());
    
    return x2 - x1;
}

// 3c. EMG関数のピーク位置を数値的に取得する関数
Double_t GetPeak(TF1 *f)
{
    if (!f) return 0;
    // 数値計算で関数の最大値を与えるX座標（ピーク位置）を探す
    return f->GetMaximumX(f->GetXmin(), f->GetXmax());
}

// 3d. 誤差伝播用の数値微分関数
Double_t GetDerivedError(TF1 *f, const TMatrixDSym &cov, std::function<Double_t(TF1*)> func)
{
    int nPar = f->GetNpar();
    std::vector<Double_t> params(nPar);
    for(int i=0; i<nPar; ++i) params[i] = f->GetParameter(i);

    Double_t epsilon = 1e-4; 
    std::vector<Double_t> grad(nPar);

    for (int i = 0; i < nPar; ++i) {
        Double_t original_val = params[i];
        
        f->SetParameter(i, original_val + epsilon);
        Double_t val_plus = func(f);
        
        f->SetParameter(i, original_val - epsilon);
        Double_t val_minus = func(f);
        
        grad[i] = (val_plus - val_minus) / (2.0 * epsilon);
        f->SetParameter(i, original_val); 
    }

    Double_t variance = 0.0;
    for (int i = 0; i < nPar; ++i) {
        for (int j = 0; j < nPar; ++j) {
            variance += grad[i] * grad[j] * cov(i, j);
        }
    }

    return (variance > 0) ? TMath::Sqrt(variance) : 0.0;
}

// 4. 電圧取得関数
double get_voltage_from_filename(const std::string& filename) {
    std::regex re("(\\d+)V");
    std::smatch match;
    if (std::regex_search(filename, match, re) && match.size() > 1) {
        return std::stod(match.str(1));
    }
    return -1.0;
}

// 5. 飽和判定関数 (修正版: 最終データビンと直前データビンの比率で判定)
bool check_saturation(const std::string& root_file, int ch, const std::string& type) {
    TFile* file = TFile::Open(root_file.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "警告: ROOTファイル " << root_file << " を開けません" << std::endl;
        return false;
    }

    // ヒストグラム名の探索
    std::vector<TString> names_to_try;
    names_to_try.push_back(Form("h_%s_ch%d", type.c_str(), ch));
    names_to_try.push_back(Form("h_%s_ch%02d", type.c_str(), ch));
    names_to_try.push_back(Form("%s_ch%d", type.c_str(), ch));
    names_to_try.push_back(Form("%s_ch%02d", type.c_str(), ch));
    names_to_try.push_back(Form("ch%02d_%s", ch, type.c_str()));

    TH1* hist = nullptr;
    for (const auto& hn : names_to_try) {
        hist = (TH1*)file->Get(hn);
        if (hist) { break; }
    }
    
    if (!hist) {
        file->Close();
        return false;
    }

    // --- 修正ロジック ---
    
    // 1. 中身が入っている「一番右のビン」（実質的な最終ビン）を探す
    int last_filled_bin = hist->FindLastBinAbove(0);
    
    // データが全くない場合は false
    if (last_filled_bin < 0) {
        file->Close();
        return false;
    }

    double last_content = hist->GetBinContent(last_filled_bin);

    // 2. その手前にある「中身が入っているビン」（右から2番目のビン）を探す
    //    (間に空のビンがあってもスキップして、データがあるところと比較する)
    int second_last_filled_bin = -1;
    for (int i = last_filled_bin - 1; i >= 1; --i) {
        if (hist->GetBinContent(i) > 0) {
            second_last_filled_bin = i;
            break;
        }
    }

    bool is_saturated = false;

    // 3. 比較判定
    if (second_last_filled_bin > 0) {
        double second_last_content = hist->GetBinContent(second_last_filled_bin);
        
        // 条件: 最後のビンのカウント数が、一つ手前のビンの5倍より大きい場合
        if (last_content > second_last_content * 5.0) {
            is_saturated = true;
        }
    }
    // else: データが1ビンしか入っていない場合などは、比較対象がないため false (飽和とみなさない) とする

    file->Close();
    return is_saturated;
}

// 6. ペデスタルデータの読み込み構造体と関数
struct PedestalInfo {
    double mean;
    double err;
};

std::map<std::pair<int, std::string>, PedestalInfo> load_pedestal_file(TString input_root_filename) {
    std::map<std::pair<int, std::string>, PedestalInfo> peds;
    
    TString dir = gSystem->DirName(input_root_filename);
    TString ped_filename = dir + "/hkelec_pedestal_hithist_means.txt";

    std::ifstream infile(ped_filename.Data());
    if (!infile.is_open()) {
        std::cerr << "警告: ペデスタルファイルが見つかりません: " << ped_filename << std::endl;
        return peds;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        std::string segment;
        std::vector<std::string> cols;
        while(std::getline(ss, segment, ',')) {
            cols.push_back(segment);
        }
        
        if (cols.size() >= 4) {
            int ch = std::stoi(cols[0]);
            std::string type = cols[1];
            double mean = std::stod(cols[2]);
            double err = std::stod(cols[3]);
            peds[{ch, type}] = {mean, err};
        }
    }
    return peds;
}

// 7. 電荷平均計算関数
void calculate_charge_mean(TString input_filename) {
    auto infile = TFile::Open(input_filename, "READ");
    if (!infile || infile->IsZombie()) return;

    auto peds = load_pedestal_file(input_filename);

    TString output_txt_filename = input_filename;
    output_txt_filename.ReplaceAll("_eventhist.root", "_mean.txt");
    std::ofstream outfile(output_txt_filename.Data());
    
    outfile << "# ch,type,mean,mean_err,rms,root_file" << std::endl;
    
    std::vector<std::string> base_types = {"hgain", "lgain", "tot"};
    
    for (int ch = 0; ch < 12; ++ch) {
        for (const auto& type : base_types) {
            TString hist_name = Form("h_%s_ch%d", type.c_str(), ch);
            auto hist = infile->Get<TH1D>(hist_name);
            
            if (hist && hist->GetEntries() > 0) {
                outfile << ch << "," << type << ","
                        << hist->GetMean() << "," << hist->GetMeanError() << "," 
                        << hist->GetRMS() << ","
                        << input_filename.Data() << std::endl;
            }
        }

        bool is_saturated = check_saturation(input_filename.Data(), ch, "hgain");
        
        std::string pc_type;
        double adc_mean = 0;
        double adc_err = 0;
        double adc_rms = 0;
        double ped_mean = 0;
        double ped_err = 0;
        double k_val = 0;

        if (is_saturated) {
            pc_type = "pc_by_l";
            TString h_name = Form("h_lgain_ch%d", ch);
            auto h_l = infile->Get<TH1D>(h_name);
            if (h_l) {
                adc_mean = h_l->GetMean();
                adc_err = h_l->GetMeanError();
                adc_rms = h_l->GetRMS();
            }
            
            if (peds.count({ch, "lgain"})) {
                ped_mean = peds[{ch, "lgain"}].mean;
                ped_err = peds[{ch, "lgain"}].err;
            }
            k_val = k_l;
        } 
        else {
            pc_type = "pc_by_h";
            TString h_name = Form("h_hgain_ch%d", ch);
            auto h_h = infile->Get<TH1D>(h_name);
            if (h_h) {
                adc_mean = h_h->GetMean();
                adc_err = h_h->GetMeanError();
                adc_rms = h_h->GetRMS();
            }

            if (peds.count({ch, "hgain"})) {
                ped_mean = peds[{ch, "hgain"}].mean;
                ped_err = peds[{ch, "hgain"}].err;
            }
            k_val = k_h;
        }

        if (k_val > 0 && adc_mean != 0) {
            double pc_mean = (adc_mean - ped_mean) * k_val;
            double pc_err = sqrt(pow(adc_err, 2) + pow(ped_err, 2)) * k_val;
            double pc_rms = adc_rms * k_val;

            outfile << ch << "," << pc_type << ","
                    << pc_mean << "," << pc_err << ","
                    << pc_rms << ","
                    << input_filename.Data() << std::endl;
        }
    }
    
    std::cout << "Charge mean calc completed -> " << output_txt_filename << std::endl;
    outfile.close();
    infile->Close();
}

// 8. 時間フィット関数 (アップデート版)
void fit_time(TString input_filename, bool save_pdf) {
    auto infile = TFile::Open(input_filename, "READ");
    if (!infile || infile->IsZombie()) return;

    TString output_txt_filename = input_filename;
    output_txt_filename.ReplaceAll("_eventhist.root", "_timefit.txt");
    std::ofstream outfile(output_txt_filename.Data());
    
    // ヘッダー (変更なし)
    outfile << "# ch,peak,peak_err,tts(fwhm),mu,gamma,sigma,lambda,tts_err,chi2,ndf,mean,mean_err,rms,rms_err" << std::endl;

    if (save_pdf) {
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(1);
    }

    std::vector<std::string> hist_types = {"time_diff"};

    for (int ch = 0; ch < 12; ++ch) {
        for (const auto& type : hist_types) {
            TString hist_name = Form("h_%s_ch%d", type.c_str(), ch);
            auto hist = infile->Get<TH1D>(hist_name);
            
            // データ数が少なすぎる場合は完全にスキップ（Mean計算も信頼できないため）
            if (!hist || hist->GetEntries() < 10) continue;

            // 1. ヒストグラム統計量の取得 (必ず取得)
            double h_mean = hist->GetMean();
            double h_mean_err = hist->GetMeanError();
            double h_rms = hist->GetRMS();
            double h_rms_err = hist->GetRMSError();

            // フィットパラメータの初期化 (失敗時用の値: -9999)
            double peak = -9999, peak_err = 0;
            double fwhm = -9999, fwhm_err = 0;
            double mu = -9999, gamma = -9999, sigma = -9999, lambda = -9999;
            double chi2 = -1;
            int ndf = -1;
            bool fit_success = false;

            // フィット処理 (データ数が十分ある場合のみ試行)
            if (hist->GetEntries() >= 100) {
                double hist_min = hist->GetXaxis()->GetXmin();
                double hist_max = hist->GetXaxis()->GetXmax();

                TF1 *fgaus = new TF1("fgaus", "gaus", hist_min, hist_max);
                fgaus->SetParameter(1, hist->GetBinCenter(hist->GetMaximumBin()));
                fgaus->SetParameter(2, hist->GetRMS());
                hist->Fit(fgaus, "QN0", "", hist_min, hist_max);
                
                double pre_amp = fgaus->GetParameter(0);
                double pre_mean = fgaus->GetParameter(1);
                double pre_sigma = TMath::Abs(fgaus->GetParameter(2));
                delete fgaus;

                if (pre_sigma > 0) {
                    TF1 *emg = new TF1("emg", EMG, hist_min, hist_max, 4);
                    emg->SetLineColor(kRed);
                    emg->SetNpx(2000);
                    // パラメータ設定 (省略せず記述)
                    emg->SetParName(0, "#mu"); emg->SetParName(1, "#gamma");
                    emg->SetParName(2, "#sigma"); emg->SetParName(3, "#lambda");
                    emg->SetParameter(0, pre_mean);
                    emg->SetParameter(1, pre_amp * 10.0);
                    emg->SetParameter(2, pre_sigma * 0.7);
                    emg->SetParameter(3, (pre_sigma > 1e-9) ? (1.0 / pre_sigma) : 1.0);
                    emg->SetParLimits(2, 0.01, 100);
                    emg->SetParLimits(3, 0.001, 1000);

                    TFitResultPtr fit_result = hist->Fit(emg, "SQR", "", hist_min, hist_max);
                    
                    if (fit_result.Get() && fit_result->IsValid() && fit_result->Ndf() > 0) {
                        fit_success = true;
                        TMatrixDSym cov = fit_result->GetCovarianceMatrix();
                        
                        mu = emg->GetParameter(0);
                        gamma = emg->GetParameter(1);
                        sigma = emg->GetParameter(2);
                        lambda = emg->GetParameter(3);
                        chi2 = fit_result->Chi2();
                        ndf = fit_result->Ndf();
                        
                        peak = GetPeak(emg);
                        fwhm = GetFWHM(emg);
                        peak_err = GetDerivedError(emg, cov, GetPeak);
                        fwhm_err = GetDerivedError(emg, cov, GetFWHM);

                        if (save_pdf) {
                            TCanvas* canvas = new TCanvas("c", "c", 800, 600);
                            hist->GetXaxis()->SetRangeUser(peak - 15, peak + 20);
                            hist->Draw();
                            emg->Draw("same");
                            TString pdf_name = input_filename;
                            pdf_name.ReplaceAll("_eventhist.root", Form("_%s_fit.pdf", hist_name.Data()));
                            canvas->SaveAs(pdf_name);
                            delete canvas;
                        }
                    }
                    delete emg;
                }
            }

            // 出力 (フィットの成否に関わらず出力する)
            outfile << ch << "," 
                    << peak << "," << peak_err << "," 
                    << fwhm << "," 
                    << mu << "," << gamma << "," << sigma << "," << lambda << "," 
                    << fwhm_err << "," 
                    << chi2 << "," << ndf << ","
                    << h_mean << "," << h_mean_err << "," 
                    << h_rms << "," << h_rms_err << std::endl;
        }
    }
    std::cout << "Time fit completed -> " << output_txt_filename << std::endl;
    outfile.close();
    infile->Close();
}

// 9. main関数
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "===============================================================================\n"
                  << "  MeanFinder & TimeFitter - 電荷平均計算および時間分解能解析プログラム\n"
                  << "===============================================================================\n\n"
                  << "[概要]\n"
                  << "  入力されたイベントヒストグラムROOTファイルを解析し、以下の処理を行います。\n"
                  << "  1. 電荷(Charge): ADC平均値の算出、ペデスタル減算、pCへの単位変換\n"
                  << "  2. 時間(Time)  : 時間分布のEMGフィット、TTS(FWHM)、ピーク位置の算出、基本統計量の出力\n\n"
                  << "[使い方]\n"
                  << "  $ " << argv[0] << " <input_file.root> [オプション]\n\n"
                  << "[オプション]\n"
                  << "  --fit-charge : 電荷の計算のみ実行 (デフォルト)\n"
                  << "  --fit-time   : 時間フィットのみ実行\n"
                  << "  --fit-all    : 両方を実行\n"
                  << "  --no-pdf     : PDF画像を出力しない (デフォルトは出力する)\n\n"
                  << "[入出力ファイルの仕様]\n"
                  << "  -----------------------------------------------------------------------------\n"
                  << "  | 区分 | ファイル形式     | 必須 | 内容 / 命名規則                            |\n"
                  << "  -----------------------------------------------------------------------------\n"
                  << "  | 入力 | .root            | 必須 | イベントデータ (h_hgain_chXX 等を含む)     |\n"
                  << "  | 入力 | ...means.txt     | 自動 | ペデスタル情報 (hkelec_pedestal...means.txt)|\n"
                  << "  |      |                  |      | ※入力ROOTと同じディレクトリから自動探索   |\n"
                  << "  -----------------------------------------------------------------------------\n"
                  << "  | 出力 | _mean.txt        | 自動 | 電荷計算結果 (CSV形式)                     |\n"
                  << "  | 出力 | _timefit.txt     | 自動 | 時間フィット結果 (CSV形式)                 |\n"
                  << "  |      |                  |      | ※EMGパラメータに加え、Mean, RMSも出力     |\n"
                  << "  | 出力 | _fit.pdf         | 任意 | フィット結果のプロット画像                 |\n"
                  << "  -----------------------------------------------------------------------------\n\n"
                  << "===============================================================================" << std::endl;
        return 1;
    }

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
    
    if (fit_mode == "--fit-charge" || fit_mode == "--fit-all") {
        calculate_charge_mean(argv[1]);
    }
    
    if (fit_mode == "--fit-time" || fit_mode == "--fit-all") {
        fit_time(argv[1], save_pdf);
    }

    return 0;
}
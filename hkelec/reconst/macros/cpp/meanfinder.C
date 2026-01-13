/*
 * id: meanfinder.C
 * Place: /home/daiki/keio/hkelec/reconst/macros/cpp/
 * Last Edit: 2025-12-25 Gemini
 *
 * 概要:
 * 1. 電荷(Charge)の平均値を計算 (ADC -> pC変換含む)
 * 2. 時間(Time)の分布に対し、以下の2種類のフィットを行う
 * - ガウスフィット (範囲: Mean ± 3*RMS)
 * - EMGフィット (全範囲, Hits>100の場合のみ)
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

// 5. 飽和判定関数 (修正版: ヒストグラムの末尾形状から判定)
bool check_saturation(const std::string& root_file, int ch, const std::string& type) {
    TFile* file = TFile::Open(root_file.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "警告: ROOTファイル " << root_file << " を開けません" << std::endl;
        return false;
    }

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

    // 1. 中身が入っている「一番右のビン」（実質的な最終ビン）を探す
    int last_filled_bin = hist->FindLastBinAbove(0);
    
    // データが全くない場合は false
    if (last_filled_bin < 0) {
        file->Close();
        return false;
    }

    double last_content = hist->GetBinContent(last_filled_bin);

    // 2. その手前にある「中身が入っているビン」（右から2番目のビン）を探す
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

// 8. 時間フィット関数 (アップデート版: EMG + Gaussian Fit)
void fit_time(TString input_filename, bool save_pdf) {
    auto infile = TFile::Open(input_filename, "READ");
    if (!infile || infile->IsZombie()) return;

    TString output_txt_filename = input_filename;
    output_txt_filename.ReplaceAll("_eventhist.root", "_timefit.txt");
    std::ofstream outfile(output_txt_filename.Data());
    
    // ヘッダー: 
    // 0-10: EMG結果 (peak, tts, mu, gamma, sigma, lambda等)
    // 11-14: ヒストグラム統計量 (mean, rms等)
    // 15-22: ガウスフィット結果 (amp, mu, sigma, chi2等)
    outfile << "# ch,peak,peak_err,tts(fwhm),mu,gamma,sigma,lambda,tts_err,chi2,ndf," 
            << "mean,mean_err,rms,rms_err,"
            << "g_amp,g_amp_err,g_mu,g_mu_err,g_sigma,g_sigma_err,g_chi2,g_ndf" << std::endl;

    if (save_pdf) {
        gStyle->SetOptStat(1111);
        gStyle->SetOptFit(0); // 複数のフィット線を描くため、統計ボックスは一旦消す
    }

    std::vector<std::string> hist_types = {"time_diff"};

    for (int ch = 0; ch < 12; ++ch) {
        for (const auto& type : hist_types) {
            TString hist_name = Form("h_%s_ch%d", type.c_str(), ch);
            auto hist = infile->Get<TH1D>(hist_name);
            
            // データ数が少なすぎる場合はスキップ
            if (!hist || hist->GetEntries() < 10) continue;

            // 1. ヒストグラム統計量の取得
            double h_mean = hist->GetMean();
            double h_mean_err = hist->GetMeanError();
            double h_rms = hist->GetRMS();
            double h_rms_err = hist->GetRMSError();

            // 変数の初期化 (失敗時用の値: -9999)
            // -- EMG用 --
            double emg_peak = -9999, emg_peak_err = 0;
            double emg_fwhm = -9999, emg_fwhm_err = 0;
            double emg_mu = -9999, emg_gamma = -9999, emg_sigma = -9999, emg_lambda = -9999;
            double emg_chi2 = -1;
            int emg_ndf = -1;

            // -- Gaussian用 --
            double g_amp = -9999, g_amp_err = 0;
            double g_mu = -9999, g_mu_err = 0;
            double g_sigma = -9999, g_sigma_err = 0;
            double g_chi2 = -1;
            int g_ndf = -1;

            TF1 *fgaus = nullptr;
            TF1 *femg = nullptr;

            // データがある程度ある場合のみフィットを実行
            if (hist->GetEntries() >= 50) {

                // --- 2. ガウスフィット (範囲: Mean ± 3*RMS) ---
                double fit_min = h_mean - 3.0 * h_rms;
                double fit_max = h_mean + 3.0 * h_rms;
                
                // ヒストグラムの範囲内に収める
                if (fit_min < hist->GetXaxis()->GetXmin()) fit_min = hist->GetXaxis()->GetXmin();
                if (fit_max > hist->GetXaxis()->GetXmax()) fit_max = hist->GetXaxis()->GetXmax();

                fgaus = new TF1("fgaus", "gaus", fit_min, fit_max);
                fgaus->SetLineColor(kBlue);
                fgaus->SetLineWidth(2);
                fgaus->SetParameter(1, h_mean); // 初期値: ヒストグラムMean
                fgaus->SetParameter(2, h_rms);  // 初期値: ヒストグラムRMS
                
                TFitResultPtr g_res = hist->Fit(fgaus, "SQR", "", fit_min, fit_max);///////////////
                
                if (g_res.Get() && g_res->IsValid()) {
                    g_amp = fgaus->GetParameter(0); g_amp_err = fgaus->GetParError(0);
                    g_mu = fgaus->GetParameter(1);  g_mu_err = fgaus->GetParError(1);
                    g_sigma = fgaus->GetParameter(2); g_sigma_err = fgaus->GetParError(2);
                    g_chi2 = fgaus->GetChisquare();
                    g_ndf = fgaus->GetNDF();
                }

                // --- 3. EMGフィット (全範囲) ---
                double full_min = hist->GetXaxis()->GetXmin();
                double full_max = hist->GetXaxis()->GetXmax();
                
                femg = new TF1("femg", EMG, full_min, full_max, 4);
                femg->SetLineColor(kRed);
                femg->SetLineWidth(2);
                femg->SetParName(0, "#mu");
                femg->SetParName(1, "#gamma");
                femg->SetParName(2, "#sigma");
                femg->SetParName(3, "#lambda");

                // 初期値: ガウスが成功していればそれを利用、だめならHist統計量
                double init_mu = (g_mu > -9000) ? g_mu : h_mean;
                double init_sigma = (g_sigma > 0) ? fabs(g_sigma) : h_rms;
                double init_amp = (g_amp > 0) ? g_amp : hist->GetMaximum();

                femg->SetParameter(0, init_mu);
                femg->SetParameter(1, init_amp * 10.0); // Gammaのスケール調整
                femg->SetParameter(2, init_sigma * 0.7);
                femg->SetParameter(3, (init_sigma > 1e-9) ? (1.0 / init_sigma) : 1.0);
                
                femg->SetParLimits(2, 0.01, 100);
                femg->SetParLimits(3, 0.001, 1000);

                // "+" オプションでフィット関数リストに追加 (ガウスを消さないため)
                TFitResultPtr e_res = hist->Fit(femg, "SQR0+", "", full_min, full_max);/////////////////////
                // TFitResultPtr e_res = hist->Fit(femg, "SQR+", "", full_min, full_max);

                if (e_res.Get() && e_res->IsValid() && e_res->Ndf() > 0) {
                    TMatrixDSym cov = e_res->GetCovarianceMatrix();

                    emg_mu = femg->GetParameter(0);
                    emg_gamma = femg->GetParameter(1);
                    emg_sigma = femg->GetParameter(2);
                    emg_lambda = femg->GetParameter(3);
                    
                    emg_chi2 = femg->GetChisquare();
                    emg_ndf = femg->GetNDF();
                    
                    emg_peak = GetPeak(femg);
                    emg_fwhm = GetFWHM(femg);

                    emg_peak_err = GetDerivedError(femg, cov, GetPeak);
                    emg_fwhm_err = GetDerivedError(femg, cov, GetFWHM);
                }
            }

            // 出力 (各フィットの成否に関わらず、変数の値を出力)
            outfile << ch << "," 
                    << emg_peak << "," << emg_peak_err << "," 
                    << emg_fwhm << "," 
                    << emg_mu << "," << emg_gamma << "," << emg_sigma << "," << emg_lambda << "," 
                    << emg_fwhm_err << "," 
                    << emg_chi2 << "," << emg_ndf << ","
                    << h_mean << "," << h_mean_err << "," 
                    << h_rms << "," << h_rms_err << ","
                    << g_amp << "," << g_amp_err << ","
                    << g_mu << "," << g_mu_err << ","
                    << g_sigma << "," << g_sigma_err << ","
                    << g_chi2 << "," << g_ndf
                    << std::endl;
            
            // PDF保存
            if (save_pdf) {
                TCanvas* canvas = new TCanvas("c", "c", 800, 600);
                
                // 描画範囲を調整 (ガウスまたはEMGのピーク周辺)
                double center = (g_mu > -9000) ? g_mu : (emg_peak > -9000 ? emg_peak : h_mean);
                hist->GetXaxis()->SetRangeUser(center - 20, center + 25);
                
                hist->Draw();
                if (fgaus) fgaus->Draw("same"); // 青
                //if (femg) femg->Draw("same");   // 赤///////////////
                
                TString pdf_name = input_filename;
                pdf_name.ReplaceAll("_eventhist.root", Form("_%s_fit.pdf", hist_name.Data()));
                canvas->SaveAs(pdf_name);
                delete canvas;
            }
            if (fgaus) delete fgaus;
            if (femg) delete femg;
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
                  << "  2. 時間(Time)  : 時間分布のガウスフィット、EMGフィット、TTS(FWHM)、ピーク位置の算出\n\n"
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
                  << "  |      |                  |      | ※EMG, Hist統計量, ガウスパラメータを出力  |\n"
                  << "  | 出力 | _fit.pdf         | 任意 | フィット結果のプロット画像                 |\n"
                  << "  -----------------------------------------------------------------------------\n\n"
                  << "[内部処理の詳細]\n"
                  << "  1. 電荷計算 (--fit-charge)\n"
                  << "     - hgain/lgain/tot の平均値、誤差、RMSを算出\n"
                  << "     - hgainの飽和判定を行い、pC計算時に hgain/lgain を自動選択\n"
                  << "     - pC = (ADC_mean - Pedestal_mean) * k (k=0.073[Hi] or 0.599[Lo])\n\n"
                  << "  2. 時間フィット (--fit-time)\n"
                  << "     - ガウスフィット: 範囲 [Mean - 3*RMS, Mean + 3*RMS], 初期値 Mean/RMS\n"
                  << "     - EMGフィット   : 全範囲, ガウス結果を初期値に利用\n"
                  << "     - 誤差伝播      : 共分散行列を用いた数値微分により Peak, FWHM の誤差を算出\n"
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
/*
 * id: plot_summary.C
 * Place: /home/daiki/keio/hkelec/reconst/macros/cpp/
 * Last Edit: 2025-12-25 Gemini
 *
 * 概要: 指定ディレクトリ内の _mean.txt (電荷) と _timefit.txt (時間) を集計し、
 * Charge vs 各種パラメータのグラフを作成する。
 * * * 作成されるグラフ:
 * - Hist統計量: Mean, RMS
 * - Fitパラメータ(EMG) : Peak, TTS(FWHM), Mu, Sigma, Gamma, Tau(1/lambda)
 * - Fitパラメータ(Gaus): GausAmp, GausMu, GausSigma
 * * * グラフは7次多項式でフィットし、結果をPDFとテキストに出力する。
 * * フィット関数の最小値(min_val)とその時の電荷(at_charge)も算出して出力する。
 * * 描画範囲はデータに合わせて動的に決定する。
 *
 * コンパイル:
 * g++ plot_summary.C -o plot_summary $(root-config --cflags --glibs)
 */

#include <TSystem.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TAxis.h>
#include <TString.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMatrixDSym.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <algorithm>
#include <functional>

// 誤差伝播: 派生量 (例: 最小値) の不確かさを数値微分で推定
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

// 2. データ構造体の定義
struct ChargeData {
    double val;
    double err;
    bool valid;
};

struct TimeData {
    // フィット由来パラメータ (EMG)
    double peak_val, peak_err;
    double tts_val, tts_err;
    double mu, gamma, sigma, lambda;
    bool fit_valid; // EMGフィットが成功しているか

    // ヒストグラム統計量 (常に有効)
    double h_mean, h_mean_err;
    double h_rms, h_rms_err;
    bool hist_valid;

    // ガウスフィットパラメータ
    double g_amp, g_amp_err;
    double g_mu, g_mu_err;
    double g_sigma, g_sigma_err;
    bool g_valid; // ガウスフィットが成功しているか
};

// 3. ヘルプ表示関数
void print_usage(const char* prog_name) {
    std::cerr << "===============================================================================\n"
              << "  Summary Plotter - Charge vs TimeParams グラフ作成ツール\n"
              << "===============================================================================\n"
              << "  [使い方] $ " << prog_name << " <target_dir> [--no-pdf]\n"
              << "===============================================================================" << std::endl;
}

// 4. ファイル読み込みヘルパー
std::string get_root_key(std::string filename) {
    size_t last_slash = filename.find_last_of("/");
    if (last_slash != std::string::npos) filename = filename.substr(last_slash + 1);
    
    std::vector<std::string> suffixes = {"_mean.txt", "_timefit.txt", "_eventhist.root"};
    for (const auto& suf : suffixes) {
        size_t pos = filename.find(suf);
        if (pos != std::string::npos) {
            filename = filename.substr(0, pos);
            break;
        }
    }
    return filename;
}

// 5. メイン処理関数
void process_directory(TString target_dir, bool save_pdf) {
    std::map<int, std::map<std::string, ChargeData>> charge_map;
    std::map<int, std::map<std::string, TimeData>> time_map;

    void* dirp = gSystem->OpenDirectory(target_dir);
    if (!dirp) {
        std::cerr << "エラー: ディレクトリを開けません: " << target_dir << std::endl;
        return;
    }

    const char* entry;
    // 5-1. ファイル走査とデータ読み込み
    while ((entry = gSystem->GetDirEntry(dirp))) {
        TString filename = entry;
        TString fullpath = target_dir + "/" + filename;

        // --- Chargeデータの読み込み (_mean.txt) ---
        if (filename.EndsWith("_mean.txt")) {
            std::ifstream infile(fullpath.Data());
            std::string line;
            while (std::getline(infile, line)) {
                if (line.empty() || line[0] == '#') continue;
                std::stringstream ss(line);
                std::string segment;
                std::vector<std::string> cols;
                while(std::getline(ss, segment, ',')) cols.push_back(segment);

                if (cols.size() >= 6) {
                    int ch = std::stoi(cols[0]);
                    std::string type = cols[1];
                    if (type == "pc_by_h" || type == "pc_by_l") {
                        double val = std::stod(cols[2]);
                        double err = std::stod(cols[3]);
                        std::string root_file = cols[5];
                        std::string key = get_root_key(root_file);
                        charge_map[ch][key] = {val, err, true};
                    }
                }
            }
        }
        // --- Timeデータの読み込み (_timefit.txt) ---
        else if (filename.EndsWith("_timefit.txt")) {
            std::ifstream infile(fullpath.Data());
            std::string line;
            while (std::getline(infile, line)) {
                if (line.empty() || line[0] == '#') continue;
                std::stringstream ss(line);
                std::string segment;
                std::vector<std::string> cols;
                while(std::getline(ss, segment, ',')) cols.push_back(segment);

                // format: 
                // 0-10: EMG, 11-14: Hist, 15-22: Gaus
                if (cols.size() >= 23) {
                    int ch = std::stoi(cols[0]);
                    std::string key = get_root_key(filename.Data());

                    TimeData td;
                    
                    // --- EMG Parameters ---
                    double check_emg = std::stod(cols[1]); // peak
                    if (check_emg > -9000) { 
                        td.peak_val = std::stod(cols[1]); td.peak_err = std::stod(cols[2]);
                        td.tts_val  = std::stod(cols[3]); td.tts_err  = std::stod(cols[8]);
                        td.mu = std::stod(cols[4]); td.gamma = std::stod(cols[5]);
                        td.sigma = std::stod(cols[6]); td.lambda = std::stod(cols[7]);
                        td.fit_valid = true;
                    } else {
                        td.fit_valid = false;
                    }

                    // --- Hist Stats ---
                    td.h_mean = std::stod(cols[11]); td.h_mean_err = std::stod(cols[12]);
                    td.h_rms  = std::stod(cols[13]); td.h_rms_err  = std::stod(cols[14]);
                    td.hist_valid = true;

                    // --- Gaussian Parameters ---
                    double check_gaus = std::stod(cols[17]); // g_mu
                    if (check_gaus > -9000) {
                        td.g_amp   = std::stod(cols[15]); td.g_amp_err   = std::stod(cols[16]);
                        td.g_mu    = std::stod(cols[17]); td.g_mu_err    = std::stod(cols[18]);
                        td.g_sigma = std::stod(cols[19]); td.g_sigma_err = std::stod(cols[20]);
                        td.g_valid = true;
                    } else {
                        td.g_valid = false;
                    }

                    time_map[ch][key] = td;
                }
            }
        }
    }
    gSystem->FreeDirectory(dirp);

    // 5-2. 出力ファイルの準備
    TString csv_path = target_dir + "/summary_all_data.csv";
    std::ofstream csv_outfile(csv_path.Data());
    csv_outfile << "ch,key,charge,charge_err,"
                << "h_mean,h_mean_err,h_rms,h_rms_err,"
                << "peak,peak_err,tts,tts_err,"
                << "g_amp,g_amp_err,g_mu,g_mu_err,g_sigma,g_sigma_err," // 追加
                << "mu,gamma,sigma,tau,fit_valid,g_valid" << std::endl;

    TString out_txt_path = target_dir + "/fit_results_summary.txt";
    std::ofstream outfile(out_txt_path.Data());
    // ヘッダーに min_val, min_err, at_charge を追加（パラメータは p0*x^{-1/2} + p1 + p2*x + p3*x^2）
    outfile << "# ch,graph_type,p0,p0_err,p1,p1_err,p2,p2_err,p3,p3_err,chi2,ndf,min_val,min_err,at_charge" << std::endl;

    // 5-3. チャンネルごとの処理ループ
    for (int ch = 0; ch < 12; ++ch) {
        if (charge_map.count(ch) == 0 || time_map.count(ch) == 0) continue;

        struct GraphSet {
            std::vector<double> x, ex, y, ey;
        };
        std::map<std::string, GraphSet> graphs;
        // グラフタイプのリスト (Gaussian parameters を追加)
        std::vector<std::string> graph_types = {
            "Mean", "RMS", "Peak", "TTS", "Mu", "Sigma", "Gamma", "Tau",
            "GausAmp", "GausMu", "GausSigma"
        };

        for (auto const& [key, c_data] : charge_map[ch]) {
            if (time_map[ch].count(key)) {
                TimeData t = time_map[ch][key];
                
                double tau = (t.fit_valid && t.lambda > 1e-9) ? (1.0 / t.lambda) : -9999;
                
                // CSV出力
                csv_outfile << ch << "," << key << ","
                            << c_data.val << "," << c_data.err << ","
                            << t.h_mean << "," << t.h_mean_err << ","
                            << t.h_rms << "," << t.h_rms_err << ","
                            << (t.fit_valid ? t.peak_val : -9999) << "," << (t.fit_valid ? t.peak_err : 0) << ","
                            << (t.fit_valid ? t.tts_val : -9999) << "," << (t.fit_valid ? t.tts_err : 0) << ","
                            << (t.g_valid ? t.g_amp : -9999) << "," << (t.g_valid ? t.g_amp_err : 0) << ","
                            << (t.g_valid ? t.g_mu : -9999) << "," << (t.g_valid ? t.g_mu_err : 0) << ","
                            << (t.g_valid ? t.g_sigma : -9999) << "," << (t.g_valid ? t.g_sigma_err : 0) << ","
                            << (t.fit_valid ? t.mu : -9999) << "," 
                            << (t.fit_valid ? t.gamma : -9999) << "," 
                            << (t.fit_valid ? t.sigma : -9999) << "," 
                            << tau << ","
                            << t.fit_valid << "," << t.g_valid << std::endl;

                // グラフデータ蓄積
                // 1. ヒストグラム統計量 (常にプロット)
                if (t.hist_valid) {
                    graphs["Mean"].x.push_back(c_data.val); graphs["Mean"].ex.push_back(c_data.err);
                    graphs["Mean"].y.push_back(t.h_mean);   graphs["Mean"].ey.push_back(t.h_mean_err);

                    graphs["RMS"].x.push_back(c_data.val);  graphs["RMS"].ex.push_back(c_data.err);
                    graphs["RMS"].y.push_back(t.h_rms);     graphs["RMS"].ey.push_back(t.h_rms_err);
                }

                // 2. EMGフィットパラメータ (成功時のみプロット)
                if (t.fit_valid) {
                    graphs["Peak"].x.push_back(c_data.val); graphs["Peak"].ex.push_back(c_data.err);
                    graphs["Peak"].y.push_back(t.peak_val); graphs["Peak"].ey.push_back(t.peak_err);

                    graphs["TTS"].x.push_back(c_data.val);  graphs["TTS"].ex.push_back(c_data.err);
                    graphs["TTS"].y.push_back(t.tts_val);   graphs["TTS"].ey.push_back(t.tts_err);

                    graphs["Mu"].x.push_back(c_data.val);   graphs["Mu"].ex.push_back(c_data.err);
                    graphs["Mu"].y.push_back(t.mu);         graphs["Mu"].ey.push_back(0);

                    graphs["Sigma"].x.push_back(c_data.val); graphs["Sigma"].ex.push_back(c_data.err);
                    graphs["Sigma"].y.push_back(t.sigma);    graphs["Sigma"].ey.push_back(0);

                    graphs["Gamma"].x.push_back(c_data.val); graphs["Gamma"].ex.push_back(c_data.err);
                    graphs["Gamma"].y.push_back(t.gamma);    graphs["Gamma"].ey.push_back(0);

                    graphs["Tau"].x.push_back(c_data.val);   graphs["Tau"].ex.push_back(c_data.err);
                    graphs["Tau"].y.push_back(tau);          graphs["Tau"].ey.push_back(0);
                }

                // 3. Gaussianフィットパラメータ (成功時のみプロット)
                if (t.g_valid) {
                    graphs["GausAmp"].x.push_back(c_data.val);   graphs["GausAmp"].ex.push_back(c_data.err);
                    graphs["GausAmp"].y.push_back(t.g_amp);      graphs["GausAmp"].ey.push_back(t.g_amp_err);

                    graphs["GausMu"].x.push_back(c_data.val);    graphs["GausMu"].ex.push_back(c_data.err);
                    graphs["GausMu"].y.push_back(t.g_mu);        graphs["GausMu"].ey.push_back(t.g_mu_err);

                    graphs["GausSigma"].x.push_back(c_data.val); graphs["GausSigma"].ex.push_back(c_data.err);
                    graphs["GausSigma"].y.push_back(t.g_sigma);  graphs["GausSigma"].ey.push_back(t.g_sigma_err);
                }
            }
        }

        // グラフ描画とフィッティング
        for (const auto& type : graph_types) {
            auto& g = graphs[type];
            if (g.x.empty()) continue;

            // 1. データの範囲を動的に取得
            double x_min_data = *std::min_element(g.x.begin(), g.x.end());
            double x_max_data = *std::max_element(g.x.begin(), g.x.end());

            // 2. 描画・フィット範囲の決定 (マージン10%程度)
            double range_min = (x_min_data < 0) ? x_min_data * 1.1 : 0.0;
            double range_max = (x_max_data > 0) ? x_max_data * 1.1 : 100.0;
            // x^{-1/2} を含むため 0 を跨がないように下限を微小正数へシフト
            if (range_min <= 0) range_min = 1e-6;

            // 3. TGraphErrorsの作成
            TGraphErrors* gr = new TGraphErrors(g.x.size(), g.x.data(), g.y.data(), g.ex.data(), g.ey.data());
            
            // 軸ラベル設定
            std::string y_unit = "[ns]";
            if (type == "Gamma" || type.find("Amp") != std::string::npos) y_unit = "[arb. units]";
            if (type == "Mean" || type == "Peak" || type == "GausMu") y_unit = "[ns (abs)]";
            gr->SetTitle(Form("Ch%d %s;Charge [pC];%s %s", ch, type.c_str(), type.c_str(), y_unit.c_str()));
            gr->SetMarkerStyle(20);
            gr->SetMarkerSize(0.8);

            // モデル: p0*x^{-1/2} + p1 + p2*x + p3*x^2
            // 3パラメータ版に戻す場合は下記をコメントアウト解除
            // TF1* f_model = new TF1("f_model", "[0]*pow(x,-0.5) + [1] + [2]*x", range_min, range_max);
            TF1* f_model = new TF1("f_model", "[0]*pow(x,-0.5) + [1] + [2]*x + [3]*x*x", range_min, range_max);
            f_model->SetLineColor(kRed);
            
            // 1回目のフィット (結果を取得)
            TFitResultPtr r1 = gr->Fit(f_model, "QS", "", range_min, range_max);
            
            // 1回目の結果を初期値として2回目のフィット
            double p0_init = f_model->GetParameter(0);
            double p1_init = f_model->GetParameter(1);
            double p2_init = f_model->GetParameter(2);
            double p3_init = f_model->GetParameter(3);
            f_model->SetParameters(p0_init, p1_init, p2_init, p3_init);
            TFitResultPtr r2 = gr->Fit(f_model, "S", "APE", range_min, range_max);

            // 最小値の算出 (fit範囲内で)
            double min_val = f_model->GetMinimum(range_min, range_max);
            double at_charge = f_model->GetMinimumX(range_min, range_max);

            // 最小値の誤差 (パラメータ共分散による誤差伝播)
            TMatrixDSym cov = r2->GetCovarianceMatrix();
            auto min_value_func = [&](TF1* func) -> Double_t {
                double xm = func->GetMinimumX(range_min, range_max);
                return func->Eval(xm);
            };
            double min_err = GetDerivedError(f_model, cov, min_value_func);

            // パラメータ出力
            outfile << ch << "," << type;
            for(int i=0; i<4; ++i) {
                outfile << "," << f_model->GetParameter(i) << "," << f_model->GetParError(i);
            }
            // 最後の列に最小値情報を追加 (min_val, min_err, at_charge)
            outfile << "," << f_model->GetChisquare() << "," << f_model->GetNDF() 
                    << "," << min_val << "," << min_err << "," << at_charge << std::endl;

            // PDF出力
            if (save_pdf) {
                TCanvas* c = new TCanvas("c", "c", 800, 600);
                c->SetGrid();
                gr->GetXaxis()->SetLimits(range_min, range_max);
                gr->Draw("APE"); // エラーバー付きで描画
                f_model->Draw("same"); // フィット曲線を上書き
                
                TString pdf_name = Form("%s/Charge_vs_%s_ch%02d.pdf", target_dir.Data(), type.c_str(), ch);
                c->SaveAs(pdf_name);
                delete c;
            }
            delete f_model;
            delete gr;
        }
    }
    
    std::cout << "Processing completed." << std::endl;
    std::cout << " - CSV Data  : " << csv_path << std::endl;
    std::cout << " - Results   : " << out_txt_path << std::endl;
    
    csv_outfile.close();
    outfile.close();
}

// 6. main関数
int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    TString target_dir = argv[1];
    bool save_pdf = true;

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--no-pdf") save_pdf = false;
    }

    process_directory(target_dir, save_pdf);
    return 0;
}
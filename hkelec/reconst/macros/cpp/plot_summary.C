/*
 * id: plot_summary.C
 * Place: /home/daiki/keio/hkelec/reconst/macros/cpp/
 * Last Edit: 2025-11-21 Gemini
 *
 * 概要: 指定ディレクトリ内の _mean.txt (電荷) と _timefit.txt (時間) を集計し、
 * Charge vs 各種パラメータのグラフを作成する。
 * * 作成されるグラフ:
 * - Hist統計量: Mean, RMS
 * - Fitパラメータ: Peak, TTS(FWHM), Mu, Sigma, Gamma, Tau(1/lambda)
 * * グラフは7次多項式でフィットし、結果をPDFとテキストに出力する。
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
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <algorithm>

// 2. データ構造体の定義
struct ChargeData {
    double val;
    double err;
    bool valid;
};

struct TimeData {
    // フィット由来パラメータ
    double peak_val, peak_err;
    double tts_val, tts_err;
    double mu, gamma, sigma, lambda;
    bool fit_valid; // フィットが成功しているか

    // ヒストグラム統計量 (常に有効)
    double h_mean, h_mean_err;
    double h_rms, h_rms_err;
    bool hist_valid;
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
                // 0:ch, 1:peak, 2:peak_err, 3:tts, 4:mu, 5:gamma, 6:sigma, 7:lambda, 
                // 8:tts_err, 9:chi2, 10:ndf, 11:mean, 12:mean_err, 13:rms, 14:rms_err
                if (cols.size() >= 15) {
                    int ch = std::stoi(cols[0]);
                    std::string key = get_root_key(filename.Data());

                    TimeData td;
                    // ヒストグラム統計量 (常に有効と仮定)
                    td.h_mean = std::stod(cols[11]);
                    td.h_mean_err = std::stod(cols[12]);
                    td.h_rms = std::stod(cols[13]);
                    td.h_rms_err = std::stod(cols[14]);
                    td.hist_valid = true;

                    // フィットパラメータ (失敗時は -9999)
                    double check_val = std::stod(cols[1]); // peak
                    if (check_val > -9000) { 
                        td.peak_val = std::stod(cols[1]);
                        td.peak_err = std::stod(cols[2]);
                        td.tts_val = std::stod(cols[3]);
                        td.tts_err = std::stod(cols[8]);
                        td.mu = std::stod(cols[4]);
                        td.gamma = std::stod(cols[5]);
                        td.sigma = std::stod(cols[6]);
                        td.lambda = std::stod(cols[7]);
                        td.fit_valid = true;
                    } else {
                        td.fit_valid = false;
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
                << "mu,gamma,sigma,tau,fit_valid" << std::endl;

    TString out_txt_path = target_dir + "/fit_results_summary.txt";
    std::ofstream outfile(out_txt_path.Data());
    // ヘッダーに min_val, at_charge を追加
    outfile << "# ch,graph_type,p0,p0_err,p1,p1_err,p2,p2_err,p3,p3_err,p4,p4_err,p5,p5_err,p6,p6_err,p7,p7_err,chi2,ndf,min_val,at_charge" << std::endl;

    // 5-3. チャンネルごとの処理ループ
    for (int ch = 0; ch < 12; ++ch) {
        if (charge_map.count(ch) == 0 || time_map.count(ch) == 0) continue;

        struct GraphSet {
            std::vector<double> x, ex, y, ey;
        };
        std::map<std::string, GraphSet> graphs;
        std::vector<std::string> graph_types = {"Mean", "RMS", "Peak", "TTS", "Mu", "Sigma", "Gamma", "Tau"};

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
                            << (t.fit_valid ? t.mu : -9999) << "," 
                            << (t.fit_valid ? t.gamma : -9999) << "," 
                            << (t.fit_valid ? t.sigma : -9999) << "," 
                            << tau << ","
                            << t.fit_valid << std::endl;

                // グラフデータ蓄積
                if (t.hist_valid) {
                    graphs["Mean"].x.push_back(c_data.val); graphs["Mean"].ex.push_back(c_data.err);
                    graphs["Mean"].y.push_back(t.h_mean);   graphs["Mean"].ey.push_back(t.h_mean_err);

                    graphs["RMS"].x.push_back(c_data.val);  graphs["RMS"].ex.push_back(c_data.err);
                    graphs["RMS"].y.push_back(t.h_rms);     graphs["RMS"].ey.push_back(t.h_rms_err);
                }

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
            // 下限: データが負なら拡張、正なら基本的に0から (ただしデータが0から遠い場合は要調整だが、現状は0-start基準)
            double range_min = (x_min_data < 0) ? x_min_data * 1.1 : 0.0;
            double range_max = (x_max_data > 0) ? x_max_data * 1.1 : 100.0; // データがない場合等の保護

            TGraphErrors* gr = new TGraphErrors(g.x.size(), g.x.data(), g.y.data(), g.ex.data(), g.ey.data());
            
            // 軸ラベル設定
            std::string y_unit = "[ns]";
            if (type == "Gamma") y_unit = "[arb. units]";
            if (type == "Mean" || type == "Peak") y_unit = "[ns (abs)]";
            gr->SetTitle(Form("Ch%d %s;Charge [pC];%s %s", ch, type.c_str(), type.c_str(), y_unit.c_str()));
            gr->SetMarkerStyle(20);
            gr->SetMarkerSize(0.8);

            // 7次関数フィッティング (動的範囲)
            TF1* f7 = new TF1("f7", "pol7", range_min, range_max); 
            gr->Fit(f7, "Q", "", range_min, range_max);

            // 最小値の算出 (fit範囲内で)
            double min_val = f7->GetMinimum(range_min, range_max);
            double at_charge = f7->GetMinimumX(range_min, range_max);

            // パラメータ出力
            outfile << ch << "," << type;
            for(int i=0; i<8; ++i) {
                outfile << "," << f7->GetParameter(i) << "," << f7->GetParError(i);
            }
            // 最後の列に最小値情報を追加
            outfile << "," << f7->GetChisquare() << "," << f7->GetNDF() 
                    << "," << min_val << "," << at_charge << std::endl;

            // PDF出力
            if (save_pdf) {
                TCanvas* c = new TCanvas("c", "c", 800, 600);
                c->SetGrid();
                // 軸範囲も設定するとより見やすい
                gr->GetXaxis()->SetLimits(range_min, range_max);
                gr->Draw("AP");
                
                TString pdf_name = Form("%s/Charge_vs_%s_ch%02d.pdf", target_dir.Data(), type.c_str(), ch);
                c->SaveAs(pdf_name);
                delete c;
            }
            delete f7;
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
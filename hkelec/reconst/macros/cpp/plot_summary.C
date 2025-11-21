/*
 * id: plot_summary.C
 * Place: /home/daiki/keio/hkelec/reconst/macros/cpp/
 * Last Edit: 2025-11-21 Gemini
 *
 * 概要: 指定ディレクトリ内の _mean.txt (電荷) と _timefit.txt (時間) を集計し、
 * Charge vs Time および Charge vs TTS のグラフを作成する。
 * グラフは7次多項式でフィットし、結果をPDFとテキストに出力する。
 * ★追加: 全データをまとめたCSV (summary_all_data.csv) も出力する。
 *
 * コンパイル:
 * g++ plot_summary.C -o plot_summary $(root-config --cflags --glibs)
 */

// 1. ヘッダーファイルのインクルード
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TAxis.h>
#include <TString.h>
#include <TStyle.h>
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
    double time_val;
    double time_err;
    double tts_val;
    double tts_err;
    bool valid;
};

// 3. ヘルプ表示関数
void print_usage(const char* prog_name) {
    std::cerr << "===============================================================================\n"
              << "  Summary Plotter - Charge vs Time/TTS グラフ作成ツール\n"
              << "===============================================================================\n\n"
              << "[概要]\n"
              << "  指定されたディレクトリ内の解析結果ファイル (_mean.txt, _timefit.txt) を読み込み、\n"
              << "  チャンネルごとに以下の処理を行います。\n"
              << "  1. 全データを結合したCSV (summary_all_data.csv) を出力\n"
              << "  2. Charge vs Time および Charge vs TTS のグラフを作成\n"
              << "  3. 7次多項式 (pol7) でフィッティング\n"
              << "  4. 結果をPDF画像とテキストファイルに出力\n\n"
              << "[使い方]\n"
              << "  $ " << prog_name << " <target_dir> [--no-pdf]\n\n"
              << "[入出力ファイルの仕様]\n"
              << "  -----------------------------------------------------------------------------\n"
              << "  | 区分 | ファイル形式        | 必須 | 内容                                    |\n"
              << "  -----------------------------------------------------------------------------\n"
              << "  | 入力 | *_mean.txt          | 必須 | meanfinderの出力 (電荷情報)             |\n"
              << "  | 入力 | *_timefit.txt       | 必須 | meanfinderの出力 (時間情報)             |\n"
              << "  -----------------------------------------------------------------------------\n"
              << "  | 出力 | summary_all_data.csv| 自動 | 全データのまとめCSV (Python解析用)      |\n"
              << "  | 出力 | Charge_vs_*.pdf     | 自動 | フィット結果のプロット (chごとに生成)   |\n"
              << "  | 出力 | fit_results.txt     | 自動 | フィットパラメータのまとめ              |\n"
              << "  -----------------------------------------------------------------------------\n"
              << "===============================================================================" << std::endl;
}

// 4. ファイル読み込みヘルパー
// ファイル名からROOTファイル部分を特定してキーにする
std::string get_root_key(std::string filename) {
    // 例: /path/to/run01_eventhist_mean.txt -> run01
    // ディレクトリセパレータ削除
    size_t last_slash = filename.find_last_of("/");
    if (last_slash != std::string::npos) filename = filename.substr(last_slash + 1);
    
    // 拡張子削除 (_mean.txt, _timefit.txt, _eventhist.root など)
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
    // データ格納用マップ: ch -> (key -> data)
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

                // format: ch,type,mean,mean_err,rms,root_file
                if (cols.size() >= 6) {
                    int ch = std::stoi(cols[0]);
                    std::string type = cols[1];
                    // pCデータのみを使用
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

                // format: ch,peak,peak_err,tts(fwhm),mu,gamma,sigma,lambda,tts_err,chi2,ndf
                if (cols.size() >= 9) {
                    int ch = std::stoi(cols[0]);
                    // ファイル名からキーを生成
                    std::string key = get_root_key(filename.Data());

                    double peak = std::stod(cols[1]);
                    double peak_err = std::stod(cols[2]);
                    double tts = std::stod(cols[3]);
                    // tts_err is col 8
                    double tts_err = std::stod(cols[8]);

                    time_map[ch][key] = {peak, peak_err, tts, tts_err, true};
                }
            }
        }
    }
    gSystem->FreeDirectory(dirp);

    // 5-2. 出力ファイルの準備 (CSV, Summary Text)
    TString csv_path = target_dir + "/summary_all_data.csv";
    std::ofstream csv_outfile(csv_path.Data());
    // ヘッダー出力 (keyはファイル識別に有用なため追加しました)
    csv_outfile << "ch,key,charge,charge_err,time,time_err,tts,tts_err" << std::endl;

    TString out_txt_path = target_dir + "/fit_results_summary.txt";
    std::ofstream outfile(out_txt_path.Data());
    outfile << "# ch,graph_type,p0,p0_err,p1,p1_err,...,p7,p7_err,chi2,ndf" << std::endl;

    // 5-3. データ結合・出力・グラフ作成ループ
    for (int ch = 0; ch < 12; ++ch) {
        if (charge_map.count(ch) == 0 || time_map.count(ch) == 0) continue;

        // データポイントのペアリング
        std::vector<double> x_q, ex_q, y_t, ey_t, y_tts, ey_tts;
        
        for (auto const& [key, c_data] : charge_map[ch]) {
            // Chargeデータに対応するTimeデータが存在するか確認
            if (time_map[ch].count(key)) {
                TimeData t_data = time_map[ch][key];
                
                // ベクターに追加 (グラフ用)
                x_q.push_back(c_data.val);
                ex_q.push_back(c_data.err);
                
                y_t.push_back(t_data.time_val);
                ey_t.push_back(t_data.time_err);
                
                y_tts.push_back(t_data.tts_val);
                ey_tts.push_back(t_data.tts_err);

                // CSVファイルへ書き出し (Python解析用)
                csv_outfile << ch << "," << key << ","
                            << c_data.val << "," << c_data.err << ","
                            << t_data.time_val << "," << t_data.time_err << ","
                            << t_data.tts_val << "," << t_data.tts_err << std::endl;
            }
        }

        if (x_q.empty()) continue;

        // グラフ作成とフィッティングを行うヘルパーラムダ式
        auto process_graph = [&](const std::string& type_name, std::vector<double>& y, std::vector<double>& ey) {
            TGraphErrors* gr = new TGraphErrors(x_q.size(), x_q.data(), y.data(), ex_q.data(), ey.data());
            gr->SetTitle(Form("Channel %d %s;Charge [pC];%s [ns]", ch, type_name.c_str(), type_name.c_str()));
            gr->SetMarkerStyle(20);
            gr->SetMarkerSize(1.0);

            // 7次関数フィッティング
            TF1* f7 = new TF1("f7", "pol7", 0, 2000); 
            gr->Fit(f7, "Q");

            // フィットパラメータ出力
            outfile << ch << "," << type_name;
            for(int i=0; i<8; ++i) {
                outfile << "," << f7->GetParameter(i) << "," << f7->GetParError(i);
            }
            outfile << "," << f7->GetChisquare() << "," << f7->GetNDF() << std::endl;

            // PDF出力
            if (save_pdf) {
                TCanvas* c = new TCanvas("c", "c", 800, 600);
                c->SetGrid();
                gr->Draw("AP");
                TString pdf_name = Form("%s/Charge_vs_%s_ch%02d.pdf", target_dir.Data(), type_name.c_str(), ch);
                c->SaveAs(pdf_name);
                delete c;
            }
            delete f7;
            delete gr;
        };

        process_graph("Time", y_t, ey_t);
        process_graph("TTS", y_tts, ey_tts);
    }
    
    std::cout << "Summary processing completed." << std::endl;
    std::cout << " - CSV Data  : " << csv_path << std::endl;
    std::cout << " - Fit Resuls: " << out_txt_path << std::endl;
    std::cout << " - Plots     : " << target_dir << "/*.pdf" << std::endl;
    
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
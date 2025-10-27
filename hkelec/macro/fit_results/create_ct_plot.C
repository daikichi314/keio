/*
 * id: create_ct_plot.C
 * Place: ~/hkelec/DiscreteSoftware/Analysis/macro/fit_results/
 * Last Edit: 2025-10-18 Gemini
 *
 * 概要: 電荷フィットのサマリーと時間フィットのサマリーを読み込み、
 * チャンネルごとに Charge vs Time のグラフ用データを作成する。
 * コンパイル可能
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <algorithm>
#include <regex>

// --- 構造体の定義 ---
// 電荷フィットの結果を保持
struct ChargeResult {
    double peak = -1;
    double rough_sigma = -1;
    bool found = false;
};
// 時間フィットの結果を保持
struct TimeResult {
    double peak = -1;
    bool found = false;
};
// ペデスタルの値を保持
struct Pedestal {
    double hgain_peak = 0.0;
    double lgain_peak = 0.0;
};

// --- サチュレーション判定関数 (select_gain.Cと同一) ---
bool check_saturation(double peak, double sigma, double rough_sigma) {
    bool is_width_suspicious = false;
    if (rough_sigma > 1e-6 && (sigma / rough_sigma < 0.3)) {
        is_width_suspicious = true;
    }
    bool is_peak_at_max = peak > 4150.0;
    return is_width_suspicious || is_peak_at_max;
}


void create_plots(const std::string& charge_summary_file, const std::string& time_summary_file, const std::string& pedestal_file, const std::string& output_dir) {
    // 1. 各種入力ファイルを読み込み、データをmapに格納する
    std::map<std::string, Pedestal> pedestals;
    std::map<std::string, ChargeResult> charge_data;
    std::map<std::string, TimeResult> time_data;
    std::string line, key;

    // 1-1. ペデスタルファイル
    std::ifstream ped_fs(pedestal_file);
    if(ped_fs){
        std::getline(ped_fs, line); // header
        while (std::getline(ped_fs, line)) {
            std::stringstream ss(line); std::string item;
            int ch; std::string type; double peak;
            std::getline(ss, item, ','); ch = std::stoi(item);
            std::getline(ss, item, ','); type = item;
            std::getline(ss, item, ','); peak = std::stod(item);
            if(type == "hgain") pedestals["ch"+std::to_string(ch)].hgain_peak = peak;
            if(type == "lgain") pedestals["ch"+std::to_string(ch)].lgain_peak = peak;
        }
    }

    // 1-2. 電荷フィットサマリー
    std::ifstream charge_fs(charge_summary_file);
    std::getline(charge_fs, line); // header
    while (std::getline(charge_fs, line)) {
        std::stringstream ss(line); std::string item;
        int ch; std::string type; double voltage;
        std::getline(ss, item, ','); ch = std::stoi(item);
        std::getline(ss, item, ','); type = item;
        std::getline(ss, item, ','); voltage = std::stod(item);
        key = type + "_ch" + std::to_string(ch) + "_v" + std::to_string((int)voltage);
        charge_data[key].found = true;
        std::getline(ss, item, ','); charge_data[key].peak = std::stod(item);
        std::getline(ss, item, ','); // peak_err
        std::getline(ss, item, ','); // sigma
        std::getline(ss, item, ','); // sigma_err
        std::getline(ss, item, ','); // chi2_ndf
        std::getline(ss, item, ','); charge_data[key].rough_sigma = std::stod(item);
    }
    
    // 1-3. 時間フィットサマリー
    std::ifstream time_fs(time_summary_file);
    std::getline(time_fs, line); // header
    while (std::getline(time_fs, line)) {
        std::stringstream ss(line); std::string item;
        int ch; std::string type; double voltage;
        // time_diffの結果のみを使う
        std::getline(ss, item, ','); ch = std::stoi(item);
        std::getline(ss, item, ','); type = item;
        if(type != "time_diff") continue;
        std::getline(ss, item, ','); voltage = std::stod(item);
        key = "ch" + std::to_string(ch) + "_v" + std::to_string((int)voltage);
        time_data[key].found = true;
        std::getline(ss, item, ','); // tts
        std::getline(ss, item, ','); // sigma
        std::getline(ss, item, ','); // fwhm
        std::getline(ss, item, ','); time_data[key].peak = std::stod(item);
    }

    // 2. チャンネルごとにループして、Charge vs Time のデータを作成
    for (int ch = 0; ch < 12; ++ch) {
        std::string output_filename = output_dir + "/Charge_vs_Time_ch" + std::to_string(ch) + ".txt";
        std::ofstream outfile(output_filename);
        outfile << "# Charge(pC), Time_peak(ns)" << std::endl;

        // 電圧のリストを重複なく取得
        std::vector<double> voltages;
        for(auto const& [key_str, val] : charge_data) {
             if (key_str.find("ch"+std::to_string(ch)) != std::string::npos) {
                std::smatch match;
                if(std::regex_search(key_str, match, std::regex("_v(\\d+)"))) {
                    voltages.push_back(std::stod(match.str(1)));
                }
             }
        }
        std::sort(voltages.begin(), voltages.end());
        voltages.erase(std::unique(voltages.begin(), voltages.end()), voltages.end());

        // 各電圧について処理
        for(double volt : voltages) {
            std::string v_str = std::to_string((int)volt);
            ChargeResult hgain_res = charge_data["hgain_ch"+std::to_string(ch)+"_v"+v_str];
            ChargeResult lgain_res = charge_data["lgain_ch"+std::to_string(ch)+"_v"+v_str];
            TimeResult time_res = time_data["ch"+std::to_string(ch)+"_v"+v_str];

            // 時間データがなければスキップ
            if (!time_res.found) continue;

            // select_gain と同じロジックで電荷を計算
            const double k_hgain = 0.073;
            const double k_lgain = 0.599;
            double hgain_ped = pedestals["ch"+std::to_string(ch)].hgain_peak;
            double lgain_ped = pedestals["ch"+std::to_string(ch)].lgain_peak;
            double selected_charge = -1.0;

            if (hgain_res.found) {
                bool is_sat = check_saturation(hgain_res.peak, 0, hgain_res.rough_sigma); // sigmaは使わないので0
                 if (!is_sat) {
                    selected_charge = (hgain_res.peak - hgain_ped) * k_hgain;
                } else if (lgain_res.found) {
                    selected_charge = (lgain_res.peak - lgain_ped) * k_lgain;
                }
            } else if (lgain_res.found) {
                selected_charge = (lgain_res.peak - lgain_ped) * k_lgain;
            }

            if (selected_charge > 0) {
                outfile << selected_charge << " " << time_res.peak << std::endl;
            }
        }
        outfile.close();
        if(voltages.size() > 0) {
            std::cout << "Charge vs Time グラフ用ファイルを作成しました: " << output_filename << std::endl;
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "使い方: " << argv[0] << " <charge_summary.txt> <time_summary.txt> <pedestal_fits.txt> <output_dir>" << std::endl;
        return 1;
    }
    create_plots(argv[1], argv[2], argv[3], argv[4]);
    return 0;
}

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
#include <filesystem>
#include <cmath>

// --- 構造体の定義 ---
// 電荷フィットの結果を保持
struct ChargeResult {
    double peak = -1;
    double peak_err = 0;
    double rough_sigma = -1;
    bool found = false;
};
// 時間フィットの結果を保持
struct TimeResult {
    double peak = -1;
    double peak_err = 0;
    bool found = false;
};
// ペデスタルの値を保持
struct Pedestal {
    double hgain_peak = 0.0;
    double lgain_peak = 0.0;
    double hgain_peak_err = 0.0;
    double lgain_peak_err = 0.0;
};

// ADC->pC conversion constants (put here so they're easy to change)
const double k_hgain = 0.073; // pC/ADC for high gain
const double k_lgain = 0.599; // pC/ADC for low gain

// --- サチュレーション判定関数 (select_gain.Cと同一) ---
bool check_saturation(double peak, double sigma, double rough_sigma) {
    bool is_width_suspicious = false;
    if (rough_sigma > 1e-6 && (sigma / rough_sigma < 0.3)) {
        is_width_suspicious = true;
    }
    bool is_peak_at_max = peak > 4150.0;
    return is_width_suspicious || is_peak_at_max;
}


void create_plots(const std::string& charge_summary_file, const std::string& time_summary_file, const std::string& pedestal_file, const std::string& output_dir, const std::string& method = "") {
    // 1. 各種入力ファイルを読み込み、データをmapに格納する
    std::map<std::string, Pedestal> pedestals;
    std::map<std::string, ChargeResult> charge_data;
    // selected (already chosen gain) per-channel HV vs Charge files
    std::map<std::string, ChargeResult> selected_charge_data;
    std::map<std::string, TimeResult> time_data;
    std::string line, key;

    // 1-1. ペデスタルファイル
    std::ifstream ped_fs(pedestal_file);
    if(ped_fs){
        std::getline(ped_fs, line); // header
        while (std::getline(ped_fs, line)) {
            std::stringstream ss(line); std::string item;
            int ch; std::string type; double peak; double peak_err = 0.0;
            std::getline(ss, item, ','); ch = std::stoi(item);
            std::getline(ss, item, ','); type = item;
            std::getline(ss, item, ','); peak = std::stod(item);
            // try to read peak_err if present
            if (std::getline(ss, item, ',')) {
                try { peak_err = std::stod(item); } catch(...) { peak_err = 0.0; }
            }
            if(type == "hgain") {
                pedestals["ch"+std::to_string(ch)].hgain_peak = peak;
                pedestals["ch"+std::to_string(ch)].hgain_peak_err = peak_err;
            }
            if(type == "lgain") {
                pedestals["ch"+std::to_string(ch)].lgain_peak = peak;
                pedestals["ch"+std::to_string(ch)].lgain_peak_err = peak_err;
            }
        }
    }

    // 1-2. 電荷フィットサマリー
    // Charge summary input: support multiple CSV styles. Common formats encountered:
    // 1) ch,type,voltage,mean,mean_err,...  (e.g. from meanfinder with explicit type)
    // 2) ch,voltage,charge_mean,charge_mean_err  (simpler summary with no "type" field)
    // We try to parse any comma-delimited data lines and normalize to keys like
    //   <type>_ch<ch>_v<voltage>
    std::ifstream charge_fs(charge_summary_file);
    bool parsed_csv = false;
    if (charge_fs) {
        // read header (if present)
        std::getline(charge_fs, line);
        while (std::getline(charge_fs, line)) {
            if (line.size() == 0) continue;
            if (line[0] == '#') continue;
            // if the line contains commas, treat as CSV and attempt to parse
            if (line.find(',') == std::string::npos) continue;
            parsed_csv = true;
            // split fields
            std::vector<std::string> fields;
            std::stringstream ss(line); std::string item;
            while (std::getline(ss, item, ',')) fields.push_back(item);
            // Handle format with and without explicit 'type'
            try {
                if (fields.size() >= 4) {
                    // If second field is numeric (voltage), assume format: ch,voltage,mean,...
                    bool second_is_number = std::regex_match(fields[1], std::regex("^[0-9]+$"));
                    // If the second field is non-numeric, only accept it when it is a known charge type
                    bool second_is_charge_type = false;
                    if (!second_is_number) {
                        std::string t = fields[1];
                        if (t == "hgain" || t == "lgain" || t == "tot" || t == "mean") second_is_charge_type = true;
                    }
                    int ch = std::stoi(fields[0]);
                    if (second_is_number) {
                        int voltage = std::stoi(fields[1]);
                        std::string type = "tot"; // no explicit type available
                        key = type + "_ch" + std::to_string(ch) + "_v" + std::to_string(voltage);
                        charge_data[key].found = true;
                        // mean value is fields[2]
                        try { charge_data[key].peak = std::stod(fields[2]); } catch(...) { charge_data[key].peak = -1; }
                        if (fields.size() > 3) { try { charge_data[key].peak_err = std::stod(fields[3]); } catch(...) { charge_data[key].peak_err = 0; } }
                        // rough_sigma may be in later columns
                        for (size_t i = 3; i < fields.size(); ++i) { try { double v = std::stod(fields[i]); charge_data[key].rough_sigma = v; break; } catch(...) { } }
                    } else if (second_is_charge_type) {
                        // format: ch,type,voltage,mean,...
                        std::string type = fields[1];
                        int voltage = std::stoi(fields[2]);
                        key = type + "_ch" + std::to_string(ch) + "_v" + std::to_string(voltage);
                        charge_data[key].found = true;
                        try { charge_data[key].peak = std::stod(fields[3]); } catch(...) { charge_data[key].peak = -1; }
                        if (fields.size() > 4) { try { charge_data[key].peak_err = std::stod(fields[4]); } catch(...) { charge_data[key].peak_err = 0; } }
                        for (size_t i = 4; i < fields.size(); ++i) { try { double v = std::stod(fields[i]); charge_data[key].rough_sigma = v; break; } catch(...) { } }
                    } else {
                        // unrecognized CSV line (e.g. timefit block), skip
                    }
                }
            } catch(...) {
                // skip malformed lines
            }
        }
    }

    // Also try to find per-channel HV_vs_ChargeSelected files (they may exist alongside summaries)
    namespace fs = std::filesystem;
    try {
        for (auto &p : fs::directory_iterator(output_dir)) {
            std::string fname = p.path().filename().string();
            std::smatch m;
            std::regex re("HV_vs_ChargeSelected_.*_ch(\\d+)\\.txt");
            if (std::regex_search(fname, m, re)) {
                int ch = std::stoi(m.str(1));
                std::ifstream hvfs(p.path());
                if (!hvfs) continue;
                std::string hdr;
                std::getline(hvfs, hdr); // header
                std::string l;
                while (std::getline(hvfs, l)) {
                    if (l.size() == 0) continue;
                    if (l[0] == '#') continue;
                    std::stringstream ss(l);
                    double hv, charge, charge_err;
                    if (!(ss >> hv)) continue;
                    if (!(ss >> charge)) continue;
                    if (!(ss >> charge_err)) charge_err = 0;
                    std::string k = "selected_ch" + std::to_string(ch) + "_v" + std::to_string((int)hv);
                    selected_charge_data[k].found = true;
                    selected_charge_data[k].peak = charge;
                    selected_charge_data[k].peak_err = charge_err;
                }
            }
        }
    } catch (...) {
        // directory iteration may fail on old systems; ignore and continue
    }
    
    // 1-3. 時間フィットサマリー
    // Expected CSV header (examples vary):
    // ch,type,voltage,tts,sigma,fwhm,peak,peak_err,tau,chi2_ndf
    // We'll explicitly pick the 'peak' column if present (typically the 7th field, index 6).
    std::ifstream time_fs(time_summary_file);
    if (time_fs) {
        std::getline(time_fs, line); // header
        while (std::getline(time_fs, line)) {
            if (line.size() == 0) continue;
            if (line[0] == '#') continue;
            // split CSV fields
            std::vector<std::string> fields;
            std::stringstream ss(line); std::string item;
            while (std::getline(ss, item, ',')) fields.push_back(item);
            if (fields.size() < 4) continue;
            int ch = 0; double voltage = 0; std::string type="";
            try {
                ch = std::stoi(fields[0]);
                type = fields[1];
                voltage = std::stod(fields[2]);
            } catch(...) { continue; }
            if (type != "time_diff") continue;
            key = "ch" + std::to_string(ch) + "_v" + std::to_string((int)voltage);
            time_data[key].found = true;
            // prefer the explicit 'peak' column if available (index 6), otherwise fall back to last numeric
            double peak_val = -1;
            if (fields.size() > 6) {
                try { peak_val = std::stod(fields[6]); } catch(...) { peak_val = -1; }
            }
            if (peak_val < 0) {
                // fallback: find last numeric field in the row
                for (int i = (int)fields.size()-1; i >= 3; --i) {
                    try { peak_val = std::stod(fields[i]); break; } catch(...) { }
                }
            }
            if (peak_val >= 0) time_data[key].peak = peak_val;
            // peak_err not always available; leave as 0
        }
    }

    // 2. チャンネルごとにループして、Charge vs Time のデータを作成
    for (int ch = 0; ch < 12; ++ch) {
        // If a method string is provided, include it in the filename so files are distinguishable
        std::string output_filename;
        if (!method.empty()) {
            output_filename = output_dir + "/Charge_vs_Time_" + method + "_ch" + std::to_string(ch) + ".txt";
        } else {
            output_filename = output_dir + "/Charge_vs_Time_ch" + std::to_string(ch) + ".txt";
        }
        std::ofstream outfile(output_filename);
        outfile << "# Charge(pC), Time_peak(ns)" << std::endl;

        // 電圧のリストを重複なく取得 (charge_data と selected_charge_data の両方を参照)
        std::vector<double> voltages;
        for(auto const& [key_str, val] : charge_data) {
             if (key_str.find("ch"+std::to_string(ch)) != std::string::npos) {
                std::smatch match;
                if(std::regex_search(key_str, match, std::regex("_v(\\d+)"))) {
                    voltages.push_back(std::stod(match.str(1)));
                }
             }
        }
        for(auto const& [key_str, val] : selected_charge_data) {
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
            ChargeResult sel_res = selected_charge_data["selected_ch"+std::to_string(ch)+"_v"+v_str];
            TimeResult time_res = time_data["ch"+std::to_string(ch)+"_v"+v_str];

            // 時間データがなければスキップ
            if (!time_res.found) continue;

            // select_gain と同じロジックで電荷を計算
            double hgain_ped = pedestals["ch"+std::to_string(ch)].hgain_peak;
            double lgain_ped = pedestals["ch"+std::to_string(ch)].lgain_peak;
            double hgain_ped_err = pedestals["ch"+std::to_string(ch)].hgain_peak_err;
            double lgain_ped_err = pedestals["ch"+std::to_string(ch)].lgain_peak_err;
            double selected_charge = -1.0;
            double selected_charge_err = 0.0;

            // If a 'selected' HV_vs_Charge file exists for this (ch,volt), prefer it
            if (sel_res.found) {
                selected_charge = sel_res.peak;
                selected_charge_err = sel_res.peak_err;
            } else {
                if (hgain_res.found) {
                    bool is_sat = check_saturation(hgain_res.peak, 0, hgain_res.rough_sigma); // sigmaは使わないので0
                    if (!is_sat) {
                        selected_charge = (hgain_res.peak - hgain_ped) * k_hgain;
                        // propagate peak_err and pedestal_err
                        selected_charge_err = k_hgain * std::sqrt(
                            hgain_res.peak_err * hgain_res.peak_err +
                            hgain_ped_err * hgain_ped_err
                        );
                    } else if (lgain_res.found) {
                        selected_charge = (lgain_res.peak - lgain_ped) * k_lgain;
                        selected_charge_err = k_lgain * std::sqrt(
                            lgain_res.peak_err * lgain_res.peak_err +
                            lgain_ped_err * lgain_ped_err
                        );
                    }
                } else if (lgain_res.found) {
                    selected_charge = (lgain_res.peak - lgain_ped) * k_lgain;
                    selected_charge_err = k_lgain * std::sqrt(
                        lgain_res.peak_err * lgain_res.peak_err +
                        lgain_ped_err * lgain_ped_err
                    );
                }
            }

            if (selected_charge > 0) {
                // write: charge, charge_err, time_peak, time_peak_err
                outfile << selected_charge << " " << selected_charge_err << " " << time_res.peak << " " << time_res.peak_err << std::endl;
            }
        }
        outfile.close();
        if(voltages.size() > 0) {
            std::cout << "Charge vs Time グラフ用ファイルを作成しました: " << output_filename << std::endl;
        }
    }
}

int main(int argc, char* argv[]) {
    // accept either 4 args (old behaviour) or 5 args (with method)
    if (argc != 5 && argc != 6) {
        std::cerr << "使い方: " << argv[0] << " <charge_summary.txt> <time_summary.txt> <pedestal_fits.txt> <output_dir> [method]" << std::endl;
        return 1;
    }
    std::string method = "";
    if (argc == 6) method = argv[5];
    create_plots(argv[1], argv[2], argv[3], argv[4], method);
    return 0;
}

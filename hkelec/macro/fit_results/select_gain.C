/*
 * id: select_gain.C
 * Place: ~/hkelec/DiscreteSoftware/Analysis/macro/fit_results/
 * Last Edit: 2025-10-16 Gemini
 *
 * 概要: gausfitのサマリーとペデスタルの結果を読み込み、サチュレーション判定を行い、
 * 最適なADC値を選択して電荷[pC]に変換し、グラフ用データを作成する。
 * コンパイル可能
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <algorithm>

// 1. ADC->pC 変換係数 (k)
// 将来、値が変わる場合はこの部分を編集してください
const double k_hgain = 0.073; // pC/ADC for high gain
const double k_lgain = 0.599; // pC/ADC for low gain

// 2. 解析結果とペデスタル値を保持する構造体
struct FitResult {
    int ch;
    std::string type;
    double voltage;
    double peak;
    double sigma;
    double rough_sigma;
};
struct Pedestal {
    double peak = 0.0;
};

// 3. サチュレーション判定を行う関数   // 改善の余地あり
bool check_saturation(const FitResult& result) {
    bool is_width_suspicious = false;
    if (result.rough_sigma > 1e-6) {
        if (result.sigma / result.rough_sigma < 0.3) is_width_suspicious = true;
    }
    bool is_peak_at_max = result.peak > 4150.0;
    return is_width_suspicious || is_peak_at_max;
}

void process_summary(const std::string& summary_filename, const std::string& pedestal_filename, const std::string& output_dir) {
    // 4. ペデスタルファイルを読み込む
    std::map<std::string, Pedestal> pedestals;
    std::ifstream ped_file(pedestal_filename);
    if (!ped_file) {
        std::cerr << "警告: ペデスタルファイル " << pedestal_filename << " を開けません。ペデスタル=0として処理を続行します。" << std::endl;
    } else {
        std::string line;
        std::getline(ped_file, line);
        while (std::getline(ped_file, line)) {
            std::stringstream ss(line);
            std::string item, key;
            int ch; std::string type; double peak;
            std::getline(ss, item, ','); ch = std::stoi(item);
            std::getline(ss, item, ','); type = item;
            std::getline(ss, item, ','); peak = std::stod(item);
            key = type + "_" + std::to_string(ch);
            pedestals[key].peak = peak;
        }
        ped_file.close();
    }

    // 5. サマリーファイルを読み込む
    std::ifstream infile(summary_filename);
    if (!infile) {
        std::cerr << "エラー: サマリーファイル " << summary_filename << " を開けません。" << std::endl;
        return;
    }
    std::map<int, std::vector<FitResult>> data_by_ch;
    std::string line;
    std::getline(infile, line);
    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string item;
        FitResult res;
        std::getline(ss, item, ','); res.ch = std::stoi(item);
        std::getline(ss, item, ','); res.type = item;
        std::getline(ss, item, ','); res.voltage = std::stod(item);
        std::getline(ss, item, ','); res.peak = std::stod(item);
        std::getline(ss, item, ',');
        std::getline(ss, item, ','); res.sigma = std::stod(item);
        std::getline(ss, item, ',');
        std::getline(ss, item, ',');
        std::getline(ss, item, ','); res.rough_sigma = std::stod(item);
        data_by_ch[res.ch].push_back(res);
    }
    infile.close();

    // 6. チャンネルごとに処理
    for (auto const& [ch, results] : data_by_ch) {
        std::string output_filename = output_dir + "/HV_vs_Charge_ch" + std::to_string(ch) + ".txt";
        std::ofstream outfile(output_filename);
        outfile << "# HV(V), Charge(pC), source(hgain=1_lgain=0)" << std::endl;

        std::map<double, std::vector<FitResult>> data_by_voltage;
        for(const auto& res : results) data_by_voltage[res.voltage].push_back(res);

        for(auto const& [volt, res_pair] : data_by_voltage){
            FitResult hgain_res, lgain_res;
            bool hgain_found = false, lgain_found = false;
            for(const auto& res : res_pair){
                if(res.type == "hgain") { hgain_res = res; hgain_found = true; }
                if(res.type == "lgain") { lgain_res = res; lgain_found = true; }
            }
            
            // 7. ADCから電荷[pC]への変換ロジック
            double hgain_ped = pedestals["hgain_" + std::to_string(ch)].peak;
            double lgain_ped = pedestals["lgain_" + std::to_string(ch)].peak;
            double selected_charge = -1.0;
            int source_flag = -1;

            if (hgain_found) {
                if (!check_saturation(hgain_res)) {
                    selected_charge = (hgain_res.peak - hgain_ped) * k_hgain;
                    source_flag = 1;
                } else if (lgain_found) {
                    selected_charge = (lgain_res.peak - lgain_ped) * k_lgain;
                    source_flag = 0;
                }
            } else if (lgain_found) {
                selected_charge = (lgain_res.peak - lgain_ped) * k_lgain;
                source_flag = 0;
            }
            
            if(source_flag != -1) {
                outfile << volt << " " << selected_charge << " " << source_flag << std::endl;
            }
        }
        outfile.close();
        std::cout << "チャンネル " << ch << " のグラフ用ファイルを作成しました: " << output_filename << std::endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "使い方: " << argv[0] << " <summary_file.txt> <pedestal_file.txt> <output_dir>" << std::endl;
        return 1;
    }
    process_summary(argv[1], argv[2], argv[3]);
    return 0;
}


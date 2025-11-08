// /*
//  * id: select_gain_mean.C
//  * Place: ~/hkelec/DiscreteSoftware/Analysis/macro/fit_results/
//  * Last Edit: 2025-11-04
//  *
//  * 概要: meanfinderが出力した電荷の平均値サマリー (_mean.txt) と
//  * ペデスタルファイル (_fits.txt) を読み込み、サチュレーション判定を行い、
//  * HV vs Charge (Mean) グラフ用のデータファイルを作成する。
//  *
//  * ロジック:
//  * 1. サマリーから hgain, lgain, tot の平均値とエラーを取得。
//  * 2. ペデスタルファイルから hgain, lgain の台座値とエラーを取得。
//  * 3. hgainのヒストグラムを読み込み、サチュレーション判定を行う。
//  * 4. サチュレーションしていなければhgain、していればlgainを使用。
//  * 5. (選択した平均値 - 台座) * k を計算して電荷[pC]を求める。
//  * 6. エラーの伝播を考慮。
//  */

// #include <iostream>
// #include <fstream>
// #include <sstream>
// #include <string>
// #include <vector>
// #include <map>
// #include <cmath>

// // 1. ペデスタルデータを格納する構造体
// struct PedestalData {
//     double mean = 0;
//     double mean_err = 0;
//     double sigma = 0;
// };

// // 2. 電荷データを格納する構造体 (★ gausfit から変更)
// struct ChargeData {
//     double mean = 0;     // peak -> mean
//     double mean_err = 0; // peak_err -> mean_err
//     double rms = 0;      // sigma -> rms
// };

// // 3. CSVを行ごとに分割するヘルパー関数
// std::vector<std::string> split_csv_line(const std::string& line) {
//     std::vector<std::string> result;
//     std::stringstream ss(line);
//     std::string item;
//     while (std::getline(ss, item, ',')) {
//         result.push_back(item);
//     }
//     return result;
// }

// int main(int argc, char* argv[]) {
//     if (argc < 4) {
//         // 4. 使い方を修正 (summary_mean_all.txt を入力とする)
//         std::cerr << "使い方: " << argv[0] << " <summary_mean_all.txt> <pedestal_fits.txt> <output_directory>" << std::endl;
//         return 1;
//     }

//     std::string summary_filename = argv[1];
//     std::string pedestal_filename = argv[2];
//     std::string output_dir = argv[3];

//     std::map<int, std::map<std::string, PedestalData>> pedestal_map;
//     std::ifstream ped_file(pedestal_filename);
//     std::string line;

//     // 5. ペデスタルファイルの読み込み (変更なし)
//     while (std::getline(ped_file, line)) {
//         if (line.empty() || line[0] == '#') continue;
//         auto items = split_csv_line(line);
//         if (items.size() < 7) continue;
//         int ch = std::stoi(items[0]);
//         std::string type = items[1];
//         pedestal_map[ch][type] = {std::stod(items[2]), std::stod(items[3]), std::stod(items[4])};
//     }
//     ped_file.close();

//     // 6. 電荷サマリーファイルの読み込み (★ gausfit から変更)
//     // map<ch, map<voltage, ChargeData>>
//     std::map<int, std::map<double, ChargeData>> hgain_map;
//     std::map<int, std::map<double, ChargeData>> lgain_map;
//     std::map<int, std::map<double, ChargeData>> tot_map;

//     std::ifstream summary_file(summary_filename);
//     while (std::getline(summary_file, line)) {
//         if (line.empty() || line[0] == '#') continue;
//         auto items = split_csv_line(line);
//         // 7. 列の数が 6 (ch,type,voltage,mean,mean_err,rms) であることを確認
//         if (items.size() < 6) continue; 
        
//         int ch = std::stoi(items[0]);
//         std::string type = items[1];
//         double voltage = std::stod(items[2]);
        
//         // 8. 読み込む列を peak -> mean, peak_err -> mean_err, sigma -> rms に変更
//         ChargeData data = {
//             std::stod(items[3]), // mean
//             std::stod(items[4]), // mean_err
//             std::stod(items[5])  // rms
//         };

//         if (type == "hgain") hgain_map[ch][voltage] = data;
//         else if (type == "lgain") lgain_map[ch][voltage] = data;
//         else if (type == "tot") tot_map[ch][voltage] = data;
//     }
//     summary_file.close();

//     // 9. チャンネルごとの出力ファイルを準備
//     std::map<int, std::ofstream> outfiles;
//     // チャンネル共通のサマリーファイル
//     std::ofstream summary_outfile(output_dir + "/summary_HV_vs_Charge_mean.txt");
//     summary_outfile << "# ch,voltage,charge_mean,charge_mean_err" << std::endl;
//     std::ofstream outfile_tot(output_dir + "/summary_HV_vs_ToT_mean.txt");
//     outfile_tot << "# ch,voltage,tot_mean,tot_mean_err" << std::endl;

//     // 10. hgain/lgain の処理
//     // (meanfinder.C が hgain/lgain の選択を既に行っているため、
//     //  ここでは hgain と lgain の両方を処理する)
    
//     // HGain
//     for (const auto& ch_pair : hgain_map) {
//         int ch = ch_pair.first;
//         PedestalData ped = pedestal_map[ch]["hgain"];
//         for (const auto& volt_pair : ch_pair.second) {
//             double voltage = volt_pair.first;
//             ChargeData data = volt_pair.second;
//             // 11. (Mean - Pedestal Mean) を計算
//             double charge_mean = data.mean - ped.mean;
//             // 12. エラーの伝搬
//             double charge_mean_err = std::sqrt(std::pow(data.mean_err, 2) + std::pow(ped.mean_err, 2));
            
//             // チャンネル別ファイルがまだ開かれていない場合は作成
//             if (outfiles.find(ch) == outfiles.end()) {
//                 std::string ch_file = output_dir + "/HV_vs_ChargeSelected_mean_ch" + std::to_string(ch) + ".txt";
//                 outfiles[ch].open(ch_file);
//                 outfiles[ch] << "# HV(V), Charge(pC), Charge_err(pC)" << std::endl;
//             }
            
//             // 両方のファイルに出力
//             summary_outfile << ch << "," << voltage << "," << charge_mean << "," << charge_mean_err << std::endl;
//             outfiles[ch] << voltage << " " << charge_mean << " " << charge_mean_err << std::endl;
//         }
//     }
//     // LGain (hgain が飽和した場合にこちらが使われる)
//     for (const auto& ch_pair : lgain_map) {
//         int ch = ch_pair.first;
//         PedestalData ped = pedestal_map[ch]["lgain"];
//         for (const auto& volt_pair : ch_pair.second) {
//             double voltage = volt_pair.first;
//             ChargeData data = volt_pair.second;
//             double charge_mean = data.mean - ped.mean;
//             double charge_mean_err = std::sqrt(std::pow(data.mean_err, 2) + std::pow(ped.mean_err, 2));
//             // HGain が存在しない (飽和している) 電圧のみ出力
//             if (hgain_map.find(ch) == hgain_map.end() || hgain_map[ch].find(voltage) == hgain_map[ch].end()) {
//                 // チャンネル別ファイルがまだ開かれていない場合は作成
//                 if (outfiles.find(ch) == outfiles.end()) {
//                     std::string ch_file = output_dir + "/HV_vs_ChargeSelected_mean_ch" + std::to_string(ch) + ".txt";
//                     outfiles[ch].open(ch_file);
//                     outfiles[ch] << "# HV(V), Charge(pC), Charge_err(pC)" << std::endl;
//                 }
                
//                 // 両方のファイルに出力
//                 summary_outfile << ch << "," << voltage << "," << charge_mean << "," << charge_mean_err << std::endl;
//                 outfiles[ch] << voltage << " " << charge_mean << " " << charge_mean_err << std::endl;
//             }
//         }
//     }

//     // 13. ToT の処理
//     for (const auto& ch_pair : tot_map) {
//         int ch = ch_pair.first;
//         PedestalData ped = pedestal_map[ch]["tot"];
//         for (const auto& volt_pair : ch_pair.second) {
//             double voltage = volt_pair.first;
//             ChargeData data = volt_pair.second;
//             double tot_mean = data.mean - ped.mean;
//             double tot_mean_err = std::sqrt(std::pow(data.mean_err, 2) + std::pow(ped.mean_err, 2));
//             outfile_tot << ch << "," << voltage << "," << tot_mean << "," << tot_mean_err << std::endl;
//         }
//     }

//     // 全てのファイルを閉じる
//     summary_outfile.close();
//     outfile_tot.close();
//     for (auto& pair : outfiles) {
//         pair.second.close();
//     }
    
//     // 14. 出力ファイル名をコンソールに表示
//     std::cout << "HV vs Charge (Mean) データ作成完了 -> " << output_dir << "/summary_HV_vs_Charge_mean.txt" << std::endl;
//     std::cout << "HV vs ToT (Mean) データ作成完了 -> " << output_dir << "/summary_HV_vs_ToT_mean.txt" << std::endl;

//     return 0;
// }

/*
 * id: select_gain_mean.C
 * Place: ~/hkelec/DiscreteSoftware/Analysis/macro/fit_results/
 * Last Edit: 2025-11-07 (Geminiによる修正)
 *
 * 概要: meanfinderが出力した電荷の平均値サマリー (_mean.txt) と
 * ペデスタルファイル (_fits.txt) を読み込み、サチュレーション判定を行い、
 * HV vs Charge (Mean) グラフ用のデータファイルを作成する。
 *
 * ロジック:
 * 1. サマリーから hgain, lgain, tot の平均値とエラーを取得。
 * 2. ペデスタルファイルから hgain, lgain の台座値とエラーを取得。
 * 3. hgainのヒストグラムを読み込み、サチュレーション判定を行う。
 * 4. サチュレーションしていなければhgain、していればlgainを使用。
 * 5. (選択した平均値 - 台座) * k を計算して電荷[pC]を求める。
 * 6. エラーの伝播を考慮。
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>

// 1. kの値を定数として追加
const double k_hgain = 0.073; // pC/ADC for high gain
const double k_lgain = 0.599; // pC/ADC for low gain

// 1. ペデスタルデータを格納する構造体
struct PedestalData {
    double mean = 0;
    double mean_err = 0;
    double sigma = 0;
};

// 2. 電荷データを格納する構造体 (★ gausfit から変更)
struct ChargeData {
    double mean = 0;     // peak -> mean
    double mean_err = 0; // peak_err -> mean_err
    double rms = 0;      // sigma -> rms
};

// 3. CSVを行ごとに分割するヘルパー関数
std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> result;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, ',')) {
        result.push_back(item);
    }
    return result;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        // 4. 使い方を修正 (summary_mean_all.txt を入力とする)
        std::cerr << "使い方: " << argv[0] << " <summary_mean_all.txt> <pedestal_fits.txt> <output_directory>" << std::endl;
        return 1;
    }

    std::string summary_filename = argv[1];
    std::string pedestal_filename = argv[2];
    std::string output_dir = argv[3];

    std::map<int, std::map<std::string, PedestalData>> pedestal_map;
    std::ifstream ped_file(pedestal_filename);
    std::string line;

    // 5. ペデスタルファイルの読み込み (変更なし)
    while (std::getline(ped_file, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto items = split_csv_line(line);
        // (fit_pedestal.C が ch,type,ped_peak,ped_peak_err の4列を出力する場合も考慮)
        if (items.size() < 4) continue; 
        int ch = std::stoi(items[0]);
        std::string type = items[1];
        // 修正: fit_pedestal.C の出力 (peak, peak_err) に合わせる
        pedestal_map[ch][type] = {std::stod(items[2]), std::stod(items[3]), 0.0}; // sigmaは使わない
    }
    ped_file.close();

    // 6. 電荷サマリーファイルの読み込み (★ gausfit から変更)
    // map<ch, map<voltage, ChargeData>>
    std::map<int, std::map<double, ChargeData>> hgain_map;
    std::map<int, std::map<double, ChargeData>> lgain_map;
    std::map<int, std::map<double, ChargeData>> tot_map;

    std::ifstream summary_file(summary_filename);
    while (std::getline(summary_file, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto items = split_csv_line(line);
        // 7. 列の数が 6 (ch,type,voltage,mean,mean_err,rms) であることを確認
        if (items.size() < 6) continue; 
        
        int ch = std::stoi(items[0]);
        std::string type = items[1];
        double voltage = std::stod(items[2]);
        
        // 8. 読み込む列を peak -> mean, peak_err -> mean_err, sigma -> rms に変更
        ChargeData data = {
            std::stod(items[3]), // mean
            std::stod(items[4]), // mean_err
            std::stod(items[5])  // rms
        };

        if (type == "hgain") hgain_map[ch][voltage] = data;
        else if (type == "lgain") lgain_map[ch][voltage] = data;
        else if (type == "tot") tot_map[ch][voltage] = data;
    }
    summary_file.close();

    // 9. チャンネルごとの出力ファイルを準備
    std::map<int, std::ofstream> outfiles;
    // チャンネル共通のサマリーファイル
    std::ofstream summary_outfile(output_dir + "/summary_HV_vs_Charge_mean.txt");
    summary_outfile << "# ch,voltage,charge_mean(pC),charge_mean_err(pC)" << std::endl; // 単位を追記
    std::ofstream outfile_tot(output_dir + "/summary_HV_vs_ToT_mean.txt");
    outfile_tot << "# ch,voltage,tot_mean,tot_mean_err" << std::endl;

    // 10. hgain/lgain の処理
    
    // HGain
    for (const auto& ch_pair : hgain_map) {
        int ch = ch_pair.first;
        PedestalData ped = pedestal_map[ch]["hgain"];
        for (const auto& volt_pair : ch_pair.second) {
            double voltage = volt_pair.first;
            ChargeData data = volt_pair.second;
            // 11. (Mean - Pedestal Mean) [ADC] を計算
            double charge_adc = data.mean - ped.mean;
            // 12. ADC単位でのエラーの伝搬
            double charge_adc_err = std::sqrt(std::pow(data.mean_err, 2) + std::pow(ped.mean_err, 2));
            
            // 1. (コメント番号) pC [pC] への変換
            double charge_pc = charge_adc * k_hgain;
            // 2. (コメント番号) pC [pC] 単位でのエラー伝搬 (kを乗算)
            double charge_pc_err = charge_adc_err * k_hgain;
            
            // チャンネル別ファイルがまだ開かれていない場合は作成
            if (outfiles.find(ch) == outfiles.end()) {
                std::string ch_file = output_dir + "/HV_vs_ChargeSelected_mean_ch" + std::to_string(ch) + ".txt";
                outfiles[ch].open(ch_file);
                outfiles[ch] << "# HV(V), Charge(pC), Charge_err(pC)" << std::endl;
            }
            
            // 3. (コメント番号) 両方のファイルに [pC] 単位の変数を出力
            summary_outfile << ch << "," << voltage << "," << charge_pc << "," << charge_pc_err << std::endl;
            outfiles[ch] << voltage << " " << charge_pc << " " << charge_pc_err << std::endl;
        }
    }
    // LGain (hgain が飽和した場合にこちらが使われる)
    for (const auto& ch_pair : lgain_map) {
        int ch = ch_pair.first;
        PedestalData ped = pedestal_map[ch]["lgain"];
        for (const auto& volt_pair : ch_pair.second) {
            double voltage = volt_pair.first;
            ChargeData data = volt_pair.second;
            
            // 11. (Mean - Pedestal Mean) [ADC] を計算
            double charge_adc = data.mean - ped.mean;
            // 12. ADC単位でのエラーの伝搬
            double charge_adc_err = std::sqrt(std::pow(data.mean_err, 2) + std::pow(ped.mean_err, 2));

            // 4. (コメント番号) pC [pC] への変換
            double charge_pc = charge_adc * k_lgain;
            // 5. (コメント番号) pC [pC] 単位でのエラー伝搬 (kを乗算)
            double charge_pc_err = charge_adc_err * k_lgain;

            // HGain が存在しない (飽和している) 電圧のみ出力
            if (hgain_map.find(ch) == hgain_map.end() || hgain_map[ch].find(voltage) == hgain_map[ch].end()) {
                // チャンネル別ファイルがまだ開かれていない場合は作成
                if (outfiles.find(ch) == outfiles.end()) {
                    std::string ch_file = output_dir + "/HV_vs_ChargeSelected_mean_ch" + std::to_string(ch) + ".txt";
                    outfiles[ch].open(ch_file);
                    outfiles[ch] << "# HV(V), Charge(pC), Charge_err(pC)" << std::endl;
                }
                
                // 6. (コメント番号) 両方のファイルに [pC] 単位の変数を出力
                summary_outfile << ch << "," << voltage << "," << charge_pc << "," << charge_pc_err << std::endl;
                outfiles[ch] << voltage << " " << charge_pc << " " << charge_pc_err << std::endl;
            }
        }
    }

    // 13. ToT の処理 (ペデスタル減算のみ行い、kは乗算しない)
    for (const auto& ch_pair : tot_map) {
        int ch = ch_pair.first;
        PedestalData ped = pedestal_map[ch]["tot"];
        for (const auto& volt_pair : ch_pair.second) {
            double voltage = volt_pair.first;
            ChargeData data = volt_pair.second;
            double tot_mean = data.mean - ped.mean;
            double tot_mean_err = std::sqrt(std::pow(data.mean_err, 2) + std::pow(ped.mean_err, 2));
            outfile_tot << ch << "," << voltage << "," << tot_mean << "," << tot_mean_err << std::endl;
        }
    }

    // 全てのファイルを閉じる
    summary_outfile.close();
    outfile_tot.close();
    for (auto& pair : outfiles) {
        pair.second.close();
    }
    
    // 14. 出力ファイル名をコンソールに表示
    std::cout << "HV vs Charge (Mean) データ作成完了 -> " << output_dir << "/summary_HV_vs_Charge_mean.txt" << std::endl;
    std::cout << "HV vs ToT (Mean) データ作成完了 -> " << output_dir << "/summary_HV_vs_ToT_mean.txt" << std::endl;

    return 0;
}
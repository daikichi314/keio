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
#include <TFile.h>
#include <TH1.h>
#include <TString.h>
#include <cmath>

// 1. ADC->pC 変換係数 (k)
const double k_hgain = 0.073; // pC/ADC for high gain
const double k_lgain = 0.599; // pC/ADC for low gain

// 2. 解析結果とペデスタル値を保持する構造体
struct FitResult {
    int ch;
    std::string type;
    double voltage;
    double peak;
    double peak_err;
    double sigma;
    double sigma_err;
    std::string hist_filename;  // ヒストグラムファイル名を保持
};

struct Pedestal {
    double peak = 0.0;
    double peak_err = 0.0;
};

// 3. サチュレーション判定を行う関数（新基準）
bool check_saturation(const std::string& root_file, int ch, const std::string& type) {
    TFile* file = TFile::Open(root_file.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "警告: ROOTファイル " << root_file << " を開けません" << std::endl;
        return false;
    }

    // try several common histogram name patterns used in different macros
    std::vector<TString> names_to_try;
    names_to_try.push_back(Form("h_%s_ch%d", type.c_str(), ch));
    names_to_try.push_back(Form("h_%s_ch%02d", type.c_str(), ch));
    names_to_try.push_back(Form("%s_ch%d", type.c_str(), ch));
    names_to_try.push_back(Form("%s_ch%02d", type.c_str(), ch));
    names_to_try.push_back(Form("ch%02d_%s", ch, type.c_str()));

    TH1* hist = nullptr;
    TString used_name = "";
    for (const auto& hn : names_to_try) {
        hist = (TH1*)file->Get(hn);
        if (hist) { used_name = hn; break; }
    }
    if (!hist) {
        std::cerr << "警告: ヒストグラム (any of common patterns) for type=" << type << " ch=" << ch << " が見つかりません" << std::endl;
        file->Close();
        return false;
    }

    // 最後のビンとその前のビンの内容を比較
    int last_bin = hist->GetNbinsX();
    double last_bin_content = hist->GetBinContent(last_bin);
    
    // 右から順に最初の非空ビンを探す
    int reference_bin = last_bin - 1;
    while(reference_bin > 0 && hist->GetBinContent(reference_bin) == 0) {
        reference_bin--;
    }
    
    if(reference_bin <= 0) {
        file->Close();
        return false;  // 有効なビンが見つからない場合
    }
    
    double reference_content = hist->GetBinContent(reference_bin);
    bool is_saturated = (last_bin_content > reference_content * 5.0);
    
    file->Close();
    return is_saturated;
}

void process_summary(const char* summary_file, const char* pedestal_file, const char* output_dir, const char* method) {
    // 4. ペデスタルファイルを読み込む
    std::map<std::string, Pedestal> pedestals;
    std::ifstream ped_file(pedestal_file);
    if (!ped_file) {
        std::cerr << "警告: ペデスタルファイル " << pedestal_file << " を開けません。ペデスタル=0として処理を続行します。" << std::endl;
    } else {
        std::string line;
        std::getline(ped_file, line);
        while (std::getline(ped_file, line)) {
            std::stringstream ss(line);
            std::string item, key;
            int ch; std::string type; double peak, peak_err;
            std::getline(ss, item, ','); ch = std::stoi(item);
            std::getline(ss, item, ','); type = item;
            std::getline(ss, item, ','); peak = std::stod(item);
            std::getline(ss, item, ','); peak_err = std::stod(item);
            key = type + "_" + std::to_string(ch);
            pedestals[key].peak = peak;
            pedestals[key].peak_err = peak_err;
        }
        ped_file.close();
    }

    // 5. サマリーファイルを読み込む
    std::ifstream infile(summary_file);
    if (!infile) {
        std::cerr << "エラー: サマリーファイル " << summary_file << " を開けません。" << std::endl;
        return;
    }
    std::map<int, std::vector<FitResult>> data_by_ch;
    std::string line;
    std::getline(infile, line); // skip header
    while (std::getline(infile, line)) {
        // split by comma into tokens
        std::vector<std::string> tok;
        std::stringstream lss(line);
        std::string f;
        while (std::getline(lss, f, ',')) tok.push_back(f);

        // tolerate empty/comment lines
        if (tok.size() < 4) continue;

        FitResult res;
        // common fields: ch,type,voltage
        res.ch = std::stoi(tok[0]);
        res.type = tok[1];
        res.voltage = std::stod(tok[2]);

        // Two expected formats:
        // gausfit-style: ch,type,voltage,peak,peak_err,sigma,sigma_err,chi2_ndf,root_file
        // mean-style:  ch,type,voltage,mean,mean_err,rms,root_file
        if (tok.size() >= 9) {
            // assume gausfit-like
            res.peak = std::stod(tok[3]);
            res.peak_err = std::stod(tok[4]);
            res.sigma = std::stod(tok[5]);
            res.sigma_err = std::stod(tok[6]);
            // last token is root_file (could be at index 8 or beyond)
            res.hist_filename = tok.back();
        } else if (tok.size() >= 7) {
            // mean-style
            res.peak = std::stod(tok[3]);        // mean stored in peak field
            res.peak_err = std::stod(tok[4]);
            res.sigma = 0.0;
            res.sigma_err = 0.0;
            res.hist_filename = tok.back();
        } else {
            // fallback: try to parse available fields
            try {
                res.peak = std::stod(tok[3]);
            } catch(...) { res.peak = 0; }
            res.peak_err = 0.0;
            res.sigma = 0.0;
            res.sigma_err = 0.0;
            res.hist_filename = tok.back();
        }

        data_by_ch[res.ch].push_back(res);
    }
    infile.close();

    // 6. チャンネルごとに処理
    for (auto const& [ch, results] : data_by_ch) {
                std::string output_filename = std::string(output_dir) + "/HV_vs_Charge_" + method + "_ch" + std::to_string(ch) + ".txt";
        std::ofstream outfile(output_filename);
        outfile << "# HV(V), HV_err(V), Charge(pC), Charge_err(pC), source(hgain=1_lgain=0)" << std::endl;

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
            double hgain_ped_err = pedestals["hgain_" + std::to_string(ch)].peak_err;
            double lgain_ped = pedestals["lgain_" + std::to_string(ch)].peak;
            double lgain_ped_err = pedestals["lgain_" + std::to_string(ch)].peak_err;
            double selected_charge = -1.0;
            double selected_charge_err = 0.0;
            int source_flag = -1;

            if (hgain_found) {
                if (!check_saturation(hgain_res.hist_filename, hgain_res.ch, "hgain")) {
                    selected_charge = (hgain_res.peak - hgain_ped) * k_hgain;
                    // エラーの伝播
                    selected_charge_err = k_hgain * std::sqrt(
                        hgain_res.peak_err * hgain_res.peak_err + 
                        hgain_ped_err * hgain_ped_err
                    );
                    source_flag = 1;
                } else if (lgain_found) {
                    selected_charge = (lgain_res.peak - lgain_ped) * k_lgain;
                    // エラーの伝播
                    selected_charge_err = k_lgain * std::sqrt(
                        lgain_res.peak_err * lgain_res.peak_err + 
                        lgain_ped_err * lgain_ped_err
                    );
                    source_flag = 0;
                }
            } else if (lgain_found) {
                selected_charge = (lgain_res.peak - lgain_ped) * k_lgain;
                // エラーの伝播
                selected_charge_err = k_lgain * std::sqrt(
                    lgain_res.peak_err * lgain_res.peak_err + 
                    lgain_ped_err * lgain_ped_err
                );
                source_flag = 0;
            }
            
            if(source_flag != -1) {
                outfile << volt << " " << 0.0 << " " << selected_charge << " " 
                       << selected_charge_err << " " << source_flag << std::endl;
            }
        }
        outfile.close();
        std::cout << "チャンネル " << ch << " のグラフ用ファイルを作成しました: " << output_filename << std::endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "使い方: " << argv[0] << " <summary_file.txt> <pedestal_file.txt> <output_dir> <method>" << std::endl;
        return 1;
    }
    process_summary(argv[1], argv[2], argv[3], argv[4]);
    // process_summary(argv[1], argv[2], argv[3], "mean");
    return 0;
}


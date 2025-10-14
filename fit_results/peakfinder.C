#include <TFile.h>
#include <TH1D.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>

// ファイル名 (...-1500V_... のような形式) から電圧を抽出する関数
double get_voltage_from_filename(const std::string& filename) {
    std::regex re("(\\d+)V"); // 数字が1回以上続き、"V"で終わる部分を検索
    std::smatch match;
    if (std::regex_search(filename, match, re) && match.size() > 1) {
        return std::stod(match.str(1)); // マッチした文字列をdouble型に変換
    }
    return -1.0; // 見つからなかった場合は-1を返す
}

void find_peaks(TString input_filename) {
    // --- ファイルを開く ---
    auto infile = TFile::Open(input_filename, "READ");
    if (!infile || infile->IsZombie()) {
        std::cerr << "エラー: 入力ファイル " << input_filename << " を開けません" << std::endl;
        return;
    }
    
    // --- 出力テキストファイルの準備 ---
    TString output_txt_filename = input_filename;
    output_txt_filename.ReplaceAll("eventhist.root", "peak.txt");
    std::ofstream outfile(output_txt_filename.Data());
    // ヘッダー行を書き込む
    outfile << "# ch,type,voltage,peak_pos" << std::endl;
    
    // --- ファイル名から電圧を取得 ---
    double voltage = get_voltage_from_filename(input_filename.Data());

    // ご要望の通り、今は "hgain" のみを処理する
    std::vector<std::string> hist_types = {"hgain"};
    
    // --- チャンネルとヒストグラムの種類でループ ---
    for (int ch = 0; ch < 12; ++ch) {
        for (const auto& type : hist_types) {
            TString hist_name = Form("%s_ch%d", type.c_str(), ch);
            auto hist = infile->Get<TH1D>(hist_name);

            // ヒストグラムが存在しない、または空の場合はスキップ
            if (!hist || hist->GetEntries() == 0) {
                continue;
            }

            // 最大の度数を持つビンを探す
            int max_bin = hist->GetMaximumBin();
            // そのビンの中心値を取得する
            double peak_pos = hist->GetXaxis()->GetBinCenter(max_bin);

            // 結果を保存する
            outfile << ch << "," << type << "," << voltage << "," << peak_pos << std::endl;
        }
    }

    std::cout << "ピーク検出が完了しました。結果は " << output_txt_filename << " に保存されました。" << std::endl;
    outfile.close();
    infile->Close();
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "使い方: " << argv[0] << " <input_file_eventhist.root>" << std::endl;
        return 1;
    }
    find_peaks(argv[1]);
    return 0;
}


/*
 * id: main.cc
 * Place: /home/daiki/keio/hkelec/reconst/reco/
 * Author: Gemini 3 Pro
 * Last Edit: 2025-12-06
 *
 * 概要:
 * 4つの50cm PMTを用いた光源位置・発光時刻の再構成プログラム (メイン処理)
 * processed_hits TTree からイベント単位でヒットデータを読み込み、
 * Minuit による χ² 最小化フィットで光源位置(x,y,z)と発光時刻(t)を推定します。
 * 結果は ROOT TTree と CSV ファイルに出力されます。
 *
 * コンパイル:
 * make (Makefile を使用)
 */

#include "readData.hh"
#include "onemPMTfit.hh"
#include "fittinginput.hh"
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <string>

// 使い方を詳細に表示する関数
void PrintUsage(const char* progName) {
    std::cout << "======================================================================" << std::endl;
    std::cout << " [概要] " << std::endl;
    std::cout << " 4つの50cm PMTを用いた光源位置(x,y,z)と発光時刻(t)の再構成プログラム" << std::endl;
    std::cout << " processed_hits TTreeを持つROOTファイルを入力とし、" << std::endl;
    std::cout << " 各イベントのヒット情報(電荷, 時刻)からMinuitを用いてフィットを行います。" << std::endl;
    std::cout << " " << std::endl;
    std::cout << " [使い方] " << std::endl;
    std::cout << " " << progName << " <InputRootFile> [OutputBaseName]" << std::endl;
    std::cout << " " << std::endl;
    std::cout << " [引数] " << std::endl;
    std::cout << " <InputRootFile>  : 入力ROOTファイルのパス (*eventhist.root)" << std::endl;
    std::cout << " [OutputBaseName] : 出力ファイルのベース名 (省略可, デフォルト: fit_results)" << std::endl;
    std::cout << "                    -> .root および .csv が生成されます" << std::endl;
    std::cout << " " << std::endl;
    std::cout << " [必要なファイル] " << std::endl;
    std::cout << " hkelec_pedestal_hithist_means.txt : 入力ROOTファイルと同じディレクトリに配置してください" << std::endl;
    std::cout << " " << std::endl;
    std::cout << " [関係するコード] " << std::endl;
    std::cout << " main.cc          : メイン処理、入出力管理" << std::endl;
    std::cout << " readData.cc/hh   : データ読み込み、ADC->電荷変換" << std::endl;
    std::cout << " onemPMTfit.cc/hh : フィッティング計算(Minuit FCN)、モデル定義" << std::endl;
    std::cout << " fittinginput.hh  : 定数、データ構造定義" << std::endl;
    std::cout << "======================================================================" << std::endl;
}

int main(int argc, char** argv) {
    // 引数チェック
    if (argc < 2) {
        PrintUsage(argv[0]);
        return 1;
    }

    std::string inputBinFile = argv[1];

    // 入力ファイルのディレクトリを取得し、ペデスタルファイルのパスを作成
    std::string dirPath = "./";
    size_t lastSlash = inputBinFile.find_last_of("/");
    if (lastSlash != std::string::npos) {
        dirPath = inputBinFile.substr(0, lastSlash + 1);
    }
    std::string pedestalFile = dirPath + "hkelec_pedestal_hithist_means.txt";

    // 出力ファイル名を入力ファイルと同じディレクトリに生成
    std::string outputBaseName = (argc >= 3) ? argv[2] : "fit_results";
    std::string outputRootFile = dirPath + outputBaseName + ".root";
    std::string outputCsvFile = dirPath + outputBaseName + ".csv";

    // 1. ペデスタルの読み込み
    std::map<int, PedestalData> pedMap;
    if (readPedestals(pedestalFile, pedMap) != 0) {
        std::cerr << "Failed to load pedestal file: " << pedestalFile << std::endl;
        return 1;
    }
    std::cout << "Loaded pedestals from " << pedestalFile << " for " << pedMap.size() << " channels." << std::endl;

    // 2. データリーダーの初期化
    DataReader reader(inputBinFile, pedMap);
    if (reader.getTotalEntries() == 0) {
        std::cerr << "No entries found or file open error." << std::endl;
        return 1;
    }
    std::cout << "Processing " << reader.getTotalEntries() << " hits..." << std::endl;

    // 3. 出力ファイルの準備
    TFile *fOut = new TFile(outputRootFile.c_str(), "RECREATE");
    TTree *tOut = new TTree("fit_results", "Fit Results");
    
    FitResult res;
    tOut->Branch("fit_x", &res.x, "fit_x/D");
    tOut->Branch("fit_y", &res.y, "fit_y/D");
    tOut->Branch("fit_z", &res.z, "fit_z/D");
    tOut->Branch("t_light", &res.t, "t_light/D");
    tOut->Branch("err_x", &res.err_x, "err_x/D");
    tOut->Branch("err_y", &res.err_y, "err_y/D");
    tOut->Branch("err_z", &res.err_z, "err_z/D");
    tOut->Branch("t_error", &res.err_t, "t_error/D");
    tOut->Branch("chi2", &res.chi2, "chi2/D");
    tOut->Branch("ndf", &res.ndf, "ndf/I");
    tOut->Branch("A", &res.A, "A/D");
    tOut->Branch("B", &res.B, "B/D");

    // CSVファイルの準備
    std::ofstream ofs(outputCsvFile.c_str());
    ofs << "fit_x,fit_y,fit_z,t_light,err_x,err_y,err_z,t_error,chi2,ndf,A,B\n";

    // 4. フィッターの準備
    LightSourceFitter fitter;
    std::vector<PMTData> eventHits;
    
    // 統計用カウンタ
    int n_total = 0;
    int n_skip  = 0;
    int n_fail  = 0;
    int n_success = 0;

    // メインループ
    while (reader.nextEvent(eventHits)) {
        n_total++;

        // ヒット数不足のチェック (4本未満ならスキップ)
        if (eventHits.size() < 4) {
            n_skip++;
            continue;
        }

        if (fitter.FitEvent(eventHits, res)) {
            // 成功したら保存
            tOut->Fill();
            ofs << res.x << "," << res.y << "," << res.z << "," << res.t << ","
                << res.err_x << "," << res.err_y << "," << res.err_z << "," << res.err_t << ","
                << res.chi2 << "," << res.ndf << "," << res.A << "," << res.B << "\n";
            n_success++;
        } else {
            n_fail++;
        }

        if (n_total % 1000 == 0) std::cout << "Processed " << n_total << " events..." << std::endl;
    }

    // 終了処理
    tOut->Write();
    fOut->Close();
    ofs.close();

    // 統計情報の表示
    std::cout << "\n=== Processing Summary ===" << std::endl;
    std::cout << " Total Events:    " << n_total << std::endl;
    std::cout << " Skipped (<4hits):" << n_skip << std::endl;
    std::cout << " Fit Failed:      " << n_fail << std::endl;
    std::cout << " Fit Success:     " << n_success << std::endl;
    std::cout << " Saved to: " << outputRootFile << std::endl;

    return 0;
}
/**
 * @file main.cc
 * @brief 光源位置再構成プログラムのメインエントリーポイント
 *
 * このプログラムは、ROOT形式の実験データを読み込み、4つのPMTの電荷・時間情報を用いて
 * 光源の位置(x,y,z)と発光時刻(t)を再構成します。
 * Minuitを使用して、観測値とモデル期待値のChi2（または尤度）を最小化します。
 *
 * @usage ./reconstructor <InputRootFile> [Options]
 */

#include "readData.hh"
#include "onemPMTfit.hh"
#include "fittinginput.hh"
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <sstream>

void PrintUsage(const char* progName) {
    std::cout << "======================================================================" << std::endl;
    std::cout << "  光源位置再構成プログラム (Light Source Reconstructor)" << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << "\n[概要]" << std::endl;
    std::cout << "  ROOTファイルのPMTヒット情報(電荷Q, 時刻T)から、光源位置(x,y,z)を推定します。" << std::endl;
    std::cout << "  TimeWalk補正、電荷依存の時間分解能、PMT半球形状を考慮したフィットを行います。" << std::endl;
    
    std::cout << "\n[使い方]" << std::endl;
    std::cout << "  " << progName << " <入力ROOTファイル> [オプション]" << std::endl;
    
    std::cout << "\n[オプション]" << std::endl;
    std::cout << "  -u <0/1>   : 3本ヒット救済モード (デフォルト: 0=OFF)" << std::endl;
    std::cout << "      1 : 3本ヒット時、残り1本を電荷0のヒットとして扱い4本分で計算" << std::endl;

    std::cout << "  -m <model> : 電荷期待値モデル (デフォルト: cosine)" << std::endl;
    std::cout << "      standard : mu = A/r^2 + B (等方発光)" << std::endl;
    std::cout << "      zeroB    : mu = A/r^2     (B=0固定)" << std::endl;
    std::cout << "      cos      : mu = (A * f(cos))/r^2 + B (角度依存あり・推奨)" << std::endl;

    std::cout << "  -q <model> : 電荷Chi2定義 (デフォルト: gaus)" << std::endl;
    std::cout << "      gaus : Gaussian" << std::endl;
    std::cout << "      bc   : Baker-Cousins (Poisson)" << std::endl;
    std::cout << "      none : 電荷情報を使用しない" << std::endl;

    std::cout << "  -t <model> : 時間Chi2定義 (デフォルト: gaus)" << std::endl;
    std::cout << "      gaus     : Gaussian (sigmaは電荷依存)" << std::endl;
    std::cout << "      goodness : SK風Goodness" << std::endl;
    std::cout << "      none     : 時間情報を使用しない" << std::endl;
    
    std::cout << "\n[設定]" << std::endl;
    std::cout << "  TimeWalk係数やSigma係数、ジオメトリ等は 'fittinginput.hh' で定義されています。" << std::endl;
    std::cout << "======================================================================" << std::endl;
}

int main(int argc, char** argv) {
    FitConfig config;
    int opt;
    std::string inputBinFile;
    
    // オプション解析
    while ((opt = getopt(argc, argv, "u:m:q:t:h")) != -1) {
        switch (opt) {
            case 'u': config.useUnhit = (std::stoi(optarg) == 1); break;
            case 'm':
                if (std::string(optarg) == "zeroB") config.chargeModel = ChargeModelType::ZeroIntercept;
                else if (std::string(optarg) == "standard") config.chargeModel = ChargeModelType::Standard;
                else config.chargeModel = ChargeModelType::Cosine;
                break;
            case 'q':
                if (std::string(optarg) == "bc") config.chargeType = ChargeChi2Type::BakerCousins;
                else if (std::string(optarg) == "none") config.chargeType = ChargeChi2Type::None;
                else config.chargeType = ChargeChi2Type::Gaussian;
                break;
            case 't':
                if (std::string(optarg) == "goodness") config.timeType = TimeChi2Type::Goodness;
                else if (std::string(optarg) == "emg") config.timeType = TimeChi2Type::EMG;
                else if (std::string(optarg) == "none") config.timeType = TimeChi2Type::None;
                else config.timeType = TimeChi2Type::Gaussian;
                break;
            case 'h':
                PrintUsage(argv[0]);
                return 0;
            default:
                PrintUsage(argv[0]);
                return 1;
        }
    }

    if (optind >= argc) {
        std::cerr << "エラー: 入力ファイルが指定されていません。\n" << std::endl;
        PrintUsage(argv[0]);
        return 1;
    }
    inputBinFile = argv[optind];

    // ファイル名生成ロジック (変更なし)
    std::string dirPath = "./";
    std::string baseName = inputBinFile;
    size_t lastSlash = inputBinFile.find_last_of("/");
    if (lastSlash != std::string::npos) {
        dirPath = inputBinFile.substr(0, lastSlash + 1);
        baseName = inputBinFile.substr(lastSlash + 1);
    }
    size_t lastDot = baseName.find_last_of(".");
    if (lastDot != std::string::npos) baseName = baseName.substr(0, lastDot);
    std::string suffixToRemove = "_eventhist";
    size_t pos = baseName.find(suffixToRemove);
    if (pos != std::string::npos) baseName.replace(pos, suffixToRemove.length(), "");

    std::stringstream ss;
    ss << "_reconst";
    if (config.useUnhit) ss << "_3hits"; else ss << "_4hits";
    
    if (config.chargeType == ChargeChi2Type::BakerCousins) ss << "_bc";
    else if (config.chargeType == ChargeChi2Type::None) ss << "_noQ";
    else ss << "_gausQ";

    if (config.chargeType != ChargeChi2Type::None) {
        if (config.chargeModel == ChargeModelType::ZeroIntercept) ss << "_zeroB";
        else if (config.chargeModel == ChargeModelType::Cosine) ss << "_cos";
        else ss << "_stdB";
    }

    if (config.timeType == TimeChi2Type::EMG) ss << "_emg";
    else if (config.timeType == TimeChi2Type::Goodness) ss << "_goodness";
    else if (config.timeType == TimeChi2Type::None) ss << "_noT";
    else ss << "_gausT";

    std::string outputRootFile = dirPath + baseName + ss.str() + ".root";
    std::string outputCsvFile = dirPath + baseName + ss.str() + ".csv";

    // 実行開始表示
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "解析を開始します: " << inputBinFile << std::endl;
    std::cout << "出力ファイル: " << outputRootFile << std::endl;
    std::cout << "モデル設定: Charge=" << (int)config.chargeType 
              << ", Model=" << (int)config.chargeModel 
              << ", Time=" << (int)config.timeType << std::endl;
    std::cout << "------------------------------------------------" << std::endl;

    // ペデスタル読み込み
    std::string pedestalFile = dirPath + "hkelec_pedestal_hithist_means.txt";
    std::map<int, PedestalData> pedMap;
    if (readPedestals(pedestalFile, pedMap) != 0) {
        std::cerr << "警告: ペデスタルファイルが見つかりません (" << pedestalFile << ")" << std::endl;
        // エラーにするか続行するかは運用次第ですが、通常は必須
        return 1;
    }

    DataReader reader(inputBinFile, pedMap);
    
    TFile *fOut = new TFile(outputRootFile.c_str(), "RECREATE");
    TTree *tOut = new TTree("fit_results", "Fit Results");
    FitResult res;
    
    tOut->Branch("fit_x", &res.x, "fit_x/D");
    tOut->Branch("fit_y", &res.y, "fit_y/D");
    tOut->Branch("fit_z", &res.z, "fit_z/D");
    tOut->Branch("t_light", &res.t, "t_light/D");
    tOut->Branch("chi2", &res.chi2, "chi2/D");
    tOut->Branch("ndf", &res.ndf, "ndf/I");
    tOut->Branch("A", &res.A, "A/D");
    tOut->Branch("B", &res.B, "B/D");
    tOut->Branch("status", &res.status, "status/I");

    std::ofstream ofs(outputCsvFile.c_str());
    ofs << "fit_x,fit_y,fit_z,t_light,chi2,ndf,A,B,status\n";

    LightSourceFitter fitter;
    fitter.SetConfig(config);

    std::vector<PMTData> eventHits;
    int n_total = 0;
    int n_success = 0;

    while (reader.nextEvent(eventHits)) {
        n_total++;

        if (config.useUnhit) {
            if (eventHits.size() < 3) continue;
            if (eventHits.size() == 3) {
                // 欠損CHをUnhit(0)として追加
                bool hitFlags[4] = {false, false, false, false};
                int eventID = eventHits[0].eventID;
                for (const auto& hit : eventHits) hitFlags[hit.ch] = true;
                for (int ch = 0; ch < 4; ++ch) {
                    if (!hitFlags[ch]) {
                        PMTData unhitData;
                        unhitData.eventID = eventID;
                        unhitData.ch = ch;
                        unhitData.charge = 0.0;
                        unhitData.time = -9999.0;
                        unhitData.isHit = false;
                        unhitData.x = PMT_POSITIONS[ch][0]; // 座標は一応入れておく
                        unhitData.y = PMT_POSITIONS[ch][1];
                        unhitData.z = PMT_POSITIONS[ch][2];
                        eventHits.push_back(unhitData);
                        break;
                    }
                }
            }
        } else {
            if (eventHits.size() < 4) continue;
        }

        if (fitter.FitEvent(eventHits, res)) {
            tOut->Fill();
            ofs << res.x << "," << res.y << "," << res.z << "," << res.t << ","
                << res.chi2 << "," << res.ndf << "," << res.A << "," << res.B << "," << res.status << "\n";
            n_success++;
        }
        
        if (n_total % 1000 == 0) std::cout << "処理中... " << n_total << " events" << std::endl;
    }

    tOut->Write();
    fOut->Close();
    ofs.close();
    
    std::cout << "完了: 全" << n_total << "イベント中、" << n_success << "イベントが収束しました。" << std::endl;

    return 0;
}
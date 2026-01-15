/**
 * @file main.cc
 * @brief 光源位置再構成プログラムのメインエントリーポイント
 *
 * このプログラムは、ROOT形式の実験データを読み込み、4つのPMTの電荷・時間情報を用いて
 * 光源の位置(x,y,z)と発光時刻(t)を再構成します。
 * Minuitを使用して、観測値とモデル期待値のChi2（または尤度）を最小化します。
 *
 * @usage ./reconstructor <InputRootFile> [Options]
 *
 * @author Gemini (Modified based on user request)
 * @date 2025-01-08
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

/**
 * @brief 使い方とオプションの説明を表示する関数
 * 初めて使用する人向けに詳細な日本語マニュアルを表示します。
 */
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

    std::cout << "  -m <model> : 電荷期待値モデル (デフォルト: func_f)" << std::endl;
    std::cout << "      func_f : r=28.5cm, mu = A * c0 * (1 - sqrt(1 - (28.5/r)^2)) * eps" << std::endl;
    std::cout << "      func_g : r=23.5cm, mu = A * c0 / r^2 * eps" << std::endl;
    std::cout << "      all    : func_f と func_g の両方で解析を実行（2回解析）" << std::endl;
    std::cout << "               ※ 係数c0, eps(角度依存)は fittinginput.hh で設定" << std::endl;

    std::cout << "  -q <model> : 電荷Chi2定義 (デフォルト: gaus)" << std::endl;
    std::cout << "      gaus : Gaussian" << std::endl;
    std::cout << "      bc   : Baker-Cousins (Poisson分布に基づく尤度比)" << std::endl;
    std::cout << "      none : 電荷情報を使用しない (時間のみでフィット)" << std::endl;

    std::cout << "  -t <model> : 時間Chi2定義 (デフォルト: gaus)" << std::endl;
    std::cout << "      gaus     : Gaussian (sigmaは電荷依存)" << std::endl;
    std::cout << "      goodness : SK風Goodness" << std::endl;
    std::cout << "      emg      : EMG分布 (現在は試験的実装)" << std::endl;
    std::cout << "      none     : 時間情報を使用しない (電荷のみでフィット)" << std::endl;
    
    std::cout << "\n[出力]" << std::endl;
    std::cout << "  入力ファイル名にオプションに応じたサフィックスを付与して出力します。" << std::endl;
    std::cout << "  例: run01_reconst_3hits_bc_func_f_goodness.root" << std::endl;
    std::cout << "  run01_reconst_3hits_bc_func_f_goodness.csv" << std::endl;
    std::cout << "  CSV出力列: fit_x,fit_y,fit_z,t_light,err_x,err_y,err_z,err_t,chi2,ndf,A,B,status" << std::endl;
    std::cout << "  ※計算に使用しなかったパラメータは -9999 が出力されます。" << std::endl;
    
    std::cout << "\n[設定]" << std::endl;
    std::cout << "  TimeWalk係数やSigma係数、ジオメトリ等は 'fittinginput.hh' で定義されています。" << std::endl;
   
    std::cout << "\n[必要なもの]" << std::endl;
    std::cout << "  - ROOT形式の入力データファイル" << std::endl;
    std::cout << "  - ペデスタル平均値ファイル 'hkelec_pedestal_hithist_means.txt' が同じディレクトリに必要です。" << std::endl;
    std::cout << "    （入力ROOTファイルと同じディレクトリを指します）" << std::endl;

    std::cout << "======================================================================" << std::endl;
}

int main(int argc, char** argv) {
    FitConfig config;
    int opt;
    std::string inputBinFile;
    bool useAllModels = false;  // allオプション用フラグ
    
    // オプション解析
    while ((opt = getopt(argc, argv, "u:m:q:t:h")) != -1) {
        switch (opt) {
            case 'u': config.useUnhit = (std::stoi(optarg) == 1); break;
            case 'm':
                if (std::string(optarg) == "all") {
                    useAllModels = true;
                    config.chargeModel = ChargeModelType::FuncF;  // 初期値をFuncFに設定
                } else if (std::string(optarg) == "func_g") {
                    config.chargeModel = ChargeModelType::FuncG;
                } else {
                    config.chargeModel = ChargeModelType::FuncF;
                }
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

    // 処理対象のモデルリストを生成
    std::vector<FitConfig> configList;
    if (useAllModels) {
        FitConfig config_f = config;
        config_f.chargeModel = ChargeModelType::FuncF;
        configList.push_back(config_f);
        
        FitConfig config_g = config;
        config_g.chargeModel = ChargeModelType::FuncG;
        configList.push_back(config_g);
    } else {
        configList.push_back(config);
    }

    // ファイル名生成ロジック
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

    // ペデスタル読み込み
    std::string pedestalFile = dirPath + "hkelec_pedestal_hithist_means.txt";
    std::map<int, PedestalData> pedMap;
    if (readPedestals(pedestalFile, pedMap) != 0) {
        std::cerr << "警告: ペデスタルファイルが見つかりません (" << pedestalFile << ")" << std::endl;
        return 1;
    }

    // ========================================================
    // 各設定ごとに処理を実行
    // ========================================================
    for (const auto& currentConfig : configList) {
        // サフィックス生成
        std::stringstream ss;
        ss << "_reconst";
        if (currentConfig.useUnhit) ss << "_3hits"; else ss << "_4hits";
        
        // 電荷Chi2 Suffix
        if (currentConfig.chargeType == ChargeChi2Type::BakerCousins) ss << "_bc";
        else if (currentConfig.chargeType == ChargeChi2Type::None) ss << "_noQ";
        else ss << "_gausQ";

        // 電荷モデル Suffix
        if (currentConfig.chargeType != ChargeChi2Type::None) {
            if (currentConfig.chargeModel == ChargeModelType::FuncG) ss << "_func_g";
            else ss << "_func_f";
        }

        // 時間Chi2 Suffix
        if (currentConfig.timeType == TimeChi2Type::EMG) ss << "_emg";
        else if (currentConfig.timeType == TimeChi2Type::Goodness) ss << "_goodness";
        else if (currentConfig.timeType == TimeChi2Type::None) ss << "_noT";
        else ss << "_gausT";

        // 出力ファイル名生成
        std::string outputRootFile = dirPath + baseName + ss.str() + ".root";
        std::string outputCsvFile = dirPath + baseName + ss.str() + ".csv";

        // 実行開始表示
        std::cout << "------------------------------------------------" << std::endl;
        std::cout << "解析を開始します: " << inputBinFile << std::endl;
        std::cout << "出力ファイル: " << outputRootFile << std::endl;
        std::cout << "モデル設定: Charge=" << (int)currentConfig.chargeType 
                  << ", Model=" << (int)currentConfig.chargeModel 
                  << ", Time=" << (int)currentConfig.timeType << std::endl;
        std::cout << "------------------------------------------------" << std::endl;

        // データリーダー初期化
        DataReader reader(inputBinFile, pedMap);
        
        // 出力ファイル初期化
        TFile *fOut = new TFile(outputRootFile.c_str(), "RECREATE");
        TTree *tOut = new TTree("fit_results", "Fit Results");
        FitResult res;
        
        // ブランチ設定
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
        ofs << "fit_x,fit_y,fit_z,t_light,err_x,err_y,err_z,err_t,chi2,ndf,A,B,status\n";

        // フィッター初期化
        LightSourceFitter fitter;
        fitter.SetConfig(currentConfig);

        // データループ
        std::vector<PMTData> eventHits;
        int n_total = 0;
        int n_success = 0;

        while (reader.nextEvent(eventHits)) {
            n_total++;

            if (currentConfig.useUnhit) {
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
                            unhitData.x = PMT_POSITIONS[ch][0]; 
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
                // 未計算値のマスク処理 (-9999)
                if (currentConfig.chargeType == ChargeChi2Type::None) {
                    res.A = -9999;
                    res.B = -9999;
                } else {
                    // ChargeFit有効時: Bは今回のモデルで存在しないため -9999 に設定
                    if (currentConfig.chargeModel == ChargeModelType::FuncF || currentConfig.chargeModel == ChargeModelType::FuncG) {
                        res.B = -9999;
                    }
                }
                
                if (currentConfig.timeType == TimeChi2Type::None) {
                    res.t = -9999;
                    res.err_t = -9999;
                }

                tOut->Fill();
                ofs << res.x << "," << res.y << "," << res.z << "," << res.t << ","
                    << res.err_x << "," << res.err_y << "," << res.err_z << "," << res.err_t << ","
                    << res.chi2 << "," << res.ndf << "," << res.A << "," << res.B << "," << res.status << "\n";
                n_success++;
            }
            
            if (n_total % 1000 == 0) std::cout << "処理中... " << n_total << " events" << std::endl;
        }

        tOut->Write();
        fOut->Close();
        ofs.close();
        
        std::cout << "完了: 全" << n_total << "イベント中、" << n_success << "イベントが収束しました。" << std::endl;
        std::cout << std::endl;
    }

    return 0;
}

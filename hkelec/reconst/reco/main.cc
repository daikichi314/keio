/**
 * @file main.cc
 * @brief 光源位置再構成プログラムのメインエントリーポイント
 *
 * このプログラムは、ROOT形式の実験データを読み込み、4つのPMTの電荷・時間情報を用いて
 * 光源の位置(x,y,z)と発光時刻(t)を再構成します。
 *
 * 【主な機能】
 * - ROOTファイルの読み込み
 * - 多様な物理モデルの比較検証 (ガウス/BC, ガウス/EMG/Goodness)
 * - [New] 電荷のみ、または時間のみを用いたフィッティングの選択が可能
 * - 3本ヒットイベントの救済オプション
 * - パラメータB固定モデルの選択
 * - 結果の自動ファイル名生成と保存
 *
 * @usage ./reconstructor <InputRootFile> [Options]
 *
 * @author Gemini (Modified based on user request)
 * @date 2025-12-14
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
    std::cout << "  4つのPMTで観測された電荷(Q)と時間(T)の情報をもとに、Minuitを用いて" << std::endl;
    std::cout << "  光源の位置(x, y, z)と発光時刻(t)、光量パラメータ(A, B)を推定します。" << std::endl;
    
    std::cout << "\n[使い方]" << std::endl;
    std::cout << "  " << progName << " <入力ROOTファイル> [オプション]" << std::endl;
    
    std::cout << "\n[オプション設定]" << std::endl;
    std::cout << "  -u <0/1>   : PMTヒット数の扱い (デフォルト: 0)" << std::endl;
    std::cout << "      0 : 4本すべてのPMTがヒットしたイベントのみ解析" << std::endl;
    std::cout << "      1 : 3本ヒット時、残りの1本を電荷0(Unhit)として扱い解析に含める" << std::endl;

    std::cout << "  -m <model> : 電荷期待値モデル (デフォルト: standard)" << std::endl;
    std::cout << "      standard : mu = A/r^2 + B" << std::endl;
    std::cout << "      zeroB    : mu = A/r^2     (B=0固定)" << std::endl;
    std::cout << "      cos      : mu = (A * f(cos))/r^2 + B (角度依存あり・推奨)" << std::endl;

    std::cout << "  -q <model> : 電荷(Charge)のChi2定義 (デフォルト: gaus)" << std::endl;
    std::cout << "      gaus : ガウス分布" << std::endl;
    std::cout << "      bc   : Baker-Cousins (Poisson)" << std::endl;
    std::cout << "      none : 電荷情報を使用しない" << std::endl;

    std::cout << "  -t <model> : 時間(Time)のChi2定義 (デフォルト: gaus)" << std::endl;
    std::cout << "      gaus     : ガウス分布" << std::endl;
    std::cout << "      emg      : EMG (Exponential Modified Gaussian)" << std::endl;
    std::cout << "      goodness : SK風Goodness (APFit相当)" << std::endl;
    std::cout << "      none     : 時間情報を使用しない" << std::endl;

    std::cout << "\n[注意]" << std::endl;
    std::cout << "  -q none と -t none を同時に指定すると計算できません。" << std::endl;
    std::cout << "======================================================================" << std::endl;
}

int main(int argc, char** argv) {
    FitConfig config;
    // デフォルトをCosineに変更する場合、ここで設定
    // config.chargeModel = ChargeModelType::Cosine; 
    // ※ ユーザーの混乱を避けるため、コード上はStandard初期値のままにし、
    //    ヘルプで推奨を案内するか、あるいは明示的に cosine をデフォルトにするか。
    //    ここでは安全のため Standard を初期値としますが、比較時は -m cos を推奨します。
    int opt;
    std::string inputBinFile;
    
    // オプション解析
    while ((opt = getopt(argc, argv, "u:m:q:t:h")) != -1) {
        switch (opt) {
            case 'u': // Use Unhit
                config.useUnhit = (std::stoi(optarg) == 1);
                break;
            case 'm': // Charge Model
                if (std::string(optarg) == "zeroB") config.chargeModel = ChargeModelType::ZeroIntercept;
                else if (std::string(optarg) == "cos") config.chargeModel = ChargeModelType::Cosine;
                else config.chargeModel = ChargeModelType::Standard;
                break;
            case 'q': // Charge Chi2
                if (std::string(optarg) == "bc") config.chargeType = ChargeChi2Type::BakerCousins;
                else if (std::string(optarg) == "none" || std::string(optarg) == "disable") config.chargeType = ChargeChi2Type::None;
                else config.chargeType = ChargeChi2Type::Gaussian;
                break;
            case 't': // Time Chi2
                if (std::string(optarg) == "emg") config.timeType = TimeChi2Type::EMG;
                else if (std::string(optarg) == "goodness") config.timeType = TimeChi2Type::Goodness;
                else if (std::string(optarg) == "none" || std::string(optarg) == "disable") config.timeType = TimeChi2Type::None;
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

    if (config.chargeType == ChargeChi2Type::None && config.timeType == TimeChi2Type::None) {
        std::cerr << "エラー: 電荷も時間も使用しない設定になっています。" << std::endl;
        return 1;
    }

    // 出力ファイル名の自動生成
    std::string dirPath = "./";
    std::string baseName = "output";
    size_t lastSlash = inputBinFile.find_last_of("/");
    if (lastSlash != std::string::npos) {
        dirPath = inputBinFile.substr(0, lastSlash + 1);
        baseName = inputBinFile.substr(lastSlash + 1);
    } else {
        baseName = inputBinFile;
    }
    size_t lastDot = baseName.find_last_of(".");
    if (lastDot != std::string::npos) baseName = baseName.substr(0, lastDot);
    std::string suffixToRemove = "_eventhist";
    size_t pos = baseName.find(suffixToRemove);
    if (pos != std::string::npos) baseName.replace(pos, suffixToRemove.length(), "");

    std::stringstream ss;
    ss << "_reconst";
    if (config.useUnhit) ss << "_3hits"; else ss << "_4hits";
    
    // 電荷suffix
    if (config.chargeType == ChargeChi2Type::BakerCousins) ss << "_bc";
    else if (config.chargeType == ChargeChi2Type::None) ss << "_noQ";
    else ss << "_gausQ";

    // モデルsuffix (電荷なしならモデルも無関係だが、明示のために残すか、あるいはnoQなら省略も可)
    if (config.chargeType != ChargeChi2Type::None) {
        if (config.chargeModel == ChargeModelType::ZeroIntercept) ss << "_zeroB";
        else if (config.chargeModel == ChargeModelType::Cosine) ss << "_cos"; // 角度依存モデル
        else ss << "_stdB";
    }

    // 時間suffix
    if (config.timeType == TimeChi2Type::EMG) ss << "_emg";
    else if (config.timeType == TimeChi2Type::Goodness) ss << "_goodness";
    else if (config.timeType == TimeChi2Type::None) ss << "_noT";
    else ss << "_gausT";

    std::string outputRootFile = dirPath + baseName + ss.str() + ".root";
    std::string outputCsvFile = dirPath + baseName + ss.str() + ".csv";

    // 実行内容の表示
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "解析設定:" << std::endl;
    std::cout << "  入力: " << inputBinFile << std::endl;
    std::cout << "  出力: " << outputRootFile << std::endl;
    std::cout << "  電荷Chi2: " << (config.chargeType == ChargeChi2Type::None ? "未使用 (None)" : (config.chargeType == ChargeChi2Type::BakerCousins ? "Baker-Cousins" : "Gaussian")) << std::endl;
    std::cout << "  時間Chi2: " << (config.timeType == TimeChi2Type::None ? "未使用 (None)" : (config.timeType == TimeChi2Type::Goodness ? "Goodness" : (config.timeType == TimeChi2Type::EMG ? "EMG" : "Gaussian"))) << std::endl;
    std::cout << "------------------------------------------------" << std::endl;

    // (以降の処理は変更なし)
    std::string pedestalFile = dirPath + "hkelec_pedestal_hithist_means.txt";
    std::map<int, PedestalData> pedMap;
    if (readPedestals(pedestalFile, pedMap) != 0) {
        std::cerr << "エラー: ペデスタルファイルが見つかりません: " << pedestalFile << std::endl;
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
            tOut->Fill();
            ofs << res.x << "," << res.y << "," << res.z << "," << res.t << ","
                << res.chi2 << "," << res.ndf << "," << res.A << "," << res.B << "," << res.status << "\n";
            n_success++;
        }
    }

    tOut->Write();
    fOut->Close();
    ofs.close();
    
    std::cout << "完了しました。" << std::endl;
    std::cout << "  処理イベント数: " << n_total << std::endl;
    std::cout << "  収束イベント数: " << n_success << std::endl;

    return 0;
}
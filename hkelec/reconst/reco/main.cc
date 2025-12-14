/**
 * @file main.cc
 * @brief 光源位置再構成プログラムのメインエントリーポイント
 *
 * このプログラムは、ROOT形式の実験データを読み込み、4つのPMTの電荷・時間情報を用いて
 * 光源の位置(x,y,z)と発光時刻(t)を再構成します。
 *
 * 主な機能:
 * - ROOTファイルの読み込み
 * - コマンドライン引数による解析モデルの選択 (Gaussian/Poisson, Gaussian/EMG)
 * - 3本ヒットイベントの救済 (Unhitチャンネルの補完)
 * - 結果のROOTファイルおよびCSVファイルへの出力
 *
 * @usage ./reconstructor <InputRootFile> [OutputBaseName] [Options]
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
#include <unistd.h> // getopt用

/**
 * @brief 使い方を表示する関数
 * 初めて使う人向けに、日本語で詳細な仕様を説明します。
 */
void PrintUsage(const char* progName) {
    std::cout << "======================================================================" << std::endl;
    std::cout << "  光源位置再構成プログラム (Light Source Reconstructor)" << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << "\n[概要]" << std::endl;
    std::cout << "  4つのPMTで観測された電荷(Q)と時間(T)の情報をもとに、Minuitを用いて" << std::endl;
    std::cout << "  光源の位置(x, y, z)と発光時刻(t)、光量パラメータ(A, B)を推定します。" << std::endl;
    std::cout << "  様々な物理モデル（ガウス分布、ポアソン分布、EMG分布）を選択して比較可能です。" << std::endl;
    
    std::cout << "\n[内部処理]" << std::endl;
    std::cout << "  1. 入力ROOTファイルを読み込み、イベントごとにPMTデータを取得します。" << std::endl;
    std::cout << "  2. ヒット数が不足している場合、オプション指定によりUnhit(Q=0)データを補完します。" << std::endl;
    std::cout << "  3. 指定されたChi2モデルに基づき、観測値と期待値の差を最小化するパラメータを探索します。" << std::endl;
    std::cout << "     - 3本ヒット時はパラメータB(背景光)を0に固定します。" << std::endl;
    std::cout << "  4. 結果をROOTファイルとCSVファイルに保存します。" << std::endl;

    std::cout << "\n[使い方]" << std::endl;
    std::cout << "  " << progName << " <入力ROOTファイル> [出力ファイル名(拡張子なし)] [オプション]" << std::endl;
    
    std::cout << "\n[入力]" << std::endl;
    std::cout << "  <入力ROOTファイル>: 解析対象のイベントデータを含むROOTファイルパス。" << std::endl;
    std::cout << "                      ヒット情報が含まれるTTreeが必要です。" << std::endl;

    std::cout << "\n[オプション]" << std::endl;
    std::cout << "  -m <model> : 電荷(Charge)のChi2モデルを指定 (デフォルト: gaus)" << std::endl;
    std::cout << "      gaus : ガウス分布モデル (従来のChi2)" << std::endl;
    std::cout << "      bc   : Baker-Cousinsモデル (ポアソン分布に基づくChi2)" << std::endl;
    std::cout << "  -t <model> : 時間(Time)のChi2モデルを指定 (デフォルト: gaus)" << std::endl;
    std::cout << "      gaus : ガウス分布モデル" << std::endl;
    std::cout << "      emg  : EMG (Exponential Modified Gaussian) モデル" << std::endl;
    std::cout << "  -u <0/1>   : 3本ヒット時のUnhit補完 (デフォルト: 0)" << std::endl;
    std::cout << "      0 : 4本すべてのPMTがヒットしたイベントのみ解析 (従来)" << std::endl;
    std::cout << "      1 : 3本ヒット時、残りの1本を電荷0として扱い解析に含める" << std::endl;
    
    std::cout << "\n[出力]" << std::endl;
    std::cout << "  指定したベース名で .root と .csv が生成されます。" << std::endl;
    std::cout << "  (例: output -> output.root, output.csv)" << std::endl;
    std::cout << "  CSVカラム: fit_x, fit_y, fit_z, t_light, err_x, ..., chi2, ndf, A, B" << std::endl;
    std::cout << "======================================================================" << std::endl;
}

int main(int argc, char** argv) {
    FitConfig config;
    int opt;
    std::string inputBinFile;
    std::string outputBaseName = "fit_results";

    // オプション解析
    while ((opt = getopt(argc, argv, "m:t:u:h")) != -1) {
        switch (opt) {
            case 'm':
                if (std::string(optarg) == "bc") config.chargeType = ChargeChi2Type::BakerCousins;
                else config.chargeType = ChargeChi2Type::Gaussian;
                break;
            case 't':
                if (std::string(optarg) == "emg") config.timeType = TimeChi2Type::EMG;
                else config.timeType = TimeChi2Type::Gaussian;
                break;
            case 'u':
                config.useUnhit = (std::stoi(optarg) == 1);
                break;
            case 'h':
                PrintUsage(argv[0]);
                return 0;
            default:
                PrintUsage(argv[0]);
                return 1;
        }
    }

    // 引数が足りない場合
    if (optind >= argc) {
        std::cerr << "エラー: 入力ファイルが指定されていません。\n" << std::endl;
        PrintUsage(argv[0]);
        return 1;
    }

    inputBinFile = argv[optind];
    if (optind + 1 < argc) {
        outputBaseName = argv[optind + 1];
    }

    // ... (以下、ペデスタル読み込み、TTree設定、ループ処理などは前回の実装と同様) ...
    // ... (Unhit補完ロジックやFitter呼び出しを含む) ...

    std::cout << "解析を開始します..." << std::endl;
    std::cout << "  入力: " << inputBinFile << std::endl;
    std::cout << "  出力: " << outputBaseName << ".root / .csv" << std::endl;
    std::cout << "  モデル(Charge): " << (config.chargeType == ChargeChi2Type::BakerCousins ? "Baker-Cousins" : "Gaussian") << std::endl;
    std::cout << "  モデル(Time)  : " << (config.timeType == TimeChi2Type::EMG ? "EMG" : "Gaussian") << std::endl;
    std::cout << "  3本ヒット救済 : " << (config.useUnhit ? "ON" : "OFF") << std::endl;

    // (以下、メインループ処理の実装)
    // ...

    return 0;
}
#include "fittinginput.hh"
#include "light_source_fit.cc"
#include "light_source_fit.hh"
#include "onemPMTfit.cc"
#include "onemPMTfit.hh"
#include "readData.cc"
#include "readData.hh"
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <chrono>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>
// #include<bits/stdc++.h>

void writeToCSV(const std::string &filename,
                const TVector3 &fit,
                double t_light,
                const TVector3 &errors,
                double t_error,
                double chi2) {
    std::ofstream ofs;
    // append モード
    ofs.open(filename, std::ios::out | std::ios::app);
    if (!ofs.is_open()) {
        std::cerr << "Cannot open " << filename << " for writing\n";
        return;
    }
    // ヘッダが無ければここで書く（簡易チェック）
    static bool header_written = false;
    if (!header_written) {
        ofs << "fit_x,fit_y,fit_z,t_light,err_x,err_y,err_z,t_error,chi2\n";
        header_written = true;
    }
    // 書き込み
    ofs << fit.X() << ',' << fit.Y() << ',' << fit.Z() << ','
        << t_light << ','
        << errors.X() << ',' << errors.Y() << ',' << errors.Z() << ','
        << t_error << ',' << chi2 << '\n';

    ofs.close();

    std::cout << "fit_x,fit_y,fit_z,t_light,err_x,err_y,err_z,t_error,chi2" << std::endl;
    std::cout << fit.X() << ',' << fit.Y() << ',' << fit.Z() << ',' << t_light << ',' << errors.X() << ',' << errors.Y() << ',' << errors.Z() << ',' << t_error << ',' << chi2 << std::endl;
}

void printUsage(const char *progname) {
    std::cerr << "===============================================================================\n"
              << "  Reconstructor - 光源位置・時刻フィッティングプログラム\n"
              << "===============================================================================\n\n"
              << "[概要]\n"
              << "  入力ROOT ファイルの PMTTree を読み込み、以下の処理を行います。\n"
              << "  1. 時間クラスタリング: PMT ヒット群を時刻でグルーピング\n"
              << "  2. mPMT 方向フィット: 各 mPMT 内で入射方向 (θ, φ) を Minuit でフィット\n"
              << "  3. 光源位置フィット: 複数 mPMT から光源の 3D 位置と発光時刻を同時フィット\n"
              << "  4. 結果出力: フィット結果を CSV ファイルに追記\n\n"
              << "[使い方]\n"
              << "  $ " << progname << " <input_file> <output_csv> [イベント番号]\n\n"
              << "[引数]\n"
              << "  input_file          : 入力 ROOT ファイルパス（.root 拡張子なし）\n"
              << "  output_csv          : 出力 CSV ファイルパス\n"
              << "  イベント番号 (任意) : 処理対象イベント番号（デフォルト: 0）\n\n"
              << "[入出力ファイル仕様]\n"
              << "  -----------------------------------------------------------------------\n"
              << "  | 区分 | ファイル形式 | 必須 | 内容 / 命名規則                  |\n"
              << "  -----------------------------------------------------------------------\n"
              << "  | 入力 | <file>.root  | 必須 | イベント PMT データ             |\n"
              << "  |      |              |      | TTree 名: \"PMTTree\"             |\n"
              << "  -----------------------------------------------------------------------\n"
              << "  | 出力 | .csv         | 自動 | フィット結果（追記モード）     |\n"
              << "  |      |              |      | 列: fit_x, fit_y, fit_z,       |\n"
              << "  |      |              |      |     t_light, err_x, err_y,     |\n"
              << "  |      |              |      |     err_z, t_error, chi2       |\n"
              << "  -----------------------------------------------------------------------\n\n"
              << "[入力 ROOT ファイルの必須要件]\n"
              << "  TTree 名: \"PMTTree\"\n"
              << "  必須ブランチ:\n"
              << "    - eventNumber (int)  : イベント番号\n"
              << "    - tubeid (int)       : PMT チューブID\n"
              << "    - mPMTid (int)       : mPMT ID\n"
              << "    - mPMT_pmtid (int)   : mPMT 内の PMT 番号 (0..19)\n"
              << "    - x, y, z (double)   : PMT 座標 (cm)\n"
              << "    - L (double)         : 検出光量 (charge)\n"
              << "    - t (double)         : ヒット時刻 (ns)\n"
              << "    - ori_x, ori_y, ori_z (double) : PMT 向きベクトル\n"
              << "    - center_x, center_y, center_z (double) : mPMT 中心位置 (cm)\n\n"
              << "[内部処理の詳細]\n"
              << "  1. データ読み込み (readData)\n"
              << "     - 指定イベントの PMT ヒットを mPMTid ごとに集約\n"
              << "     - 対象 mPMTid: {338, 339, 340, 346, 347, 348, 354, 355, 356}\n\n"
              << "  2. 時間クラスタリング (findExpandedGroups)\n"
              << "     - タイムウィンドウ内のヒット数に基づいてグルーピング\n"
              << "     - パラメータ: τ=5ns, n≥3, expand_before=0.5ns, expand_after=6ns\n\n"
              << "  3. mPMT 方向フィット (FitPosition)\n"
              << "     - 各 mPMT 内の PMT 光量と向きから光の入射方向を推定\n"
              << "     - モデル: charge ∝ 1/(1+exp(-6*(cosα-1))) + const\n"
              << "     - Minuit による χ² 最小化で θ, φ を決定\n"
              << "     - 結果を mPMT の局所基底（pmt19, pmt1, pmt4 向きテーブル）で\n"
              << "       ワールド座標の方向ベクトルに変換\n\n"
              << "  4. 光源位置・時刻フィット (FitLightSource)\n"
              << "     - 複数 mPMT の方向と検出時刻から点光源の 3D 位置と t0 を推定\n"
              << "     - χ² = Σ(方向誤差^2) + Σ(時間誤差^2)\n"
              << "     - 光速（水中）: c = 22.556 cm/ns (ハードコード)\n\n"
              << "[出力 CSV 形式]\n"
              << "  ヘッダ: fit_x, fit_y, fit_z, t_light, err_x, err_y, err_z, t_error, chi2\n"
              << "  単位: 位置 (cm), 時刻 (ns), χ² (dimensionless)\n\n"
              << "[例]\n"
              << "  $ " << progname << " /path/to/event-data result.csv 0\n"
              << "  → /path/to/event-data.root の EventNumber=0 を処理し、\n"
              << "    result.csv にフィット結果を追記\n\n"
              << "[注意事項]\n"
              << "  - 十分な mPMT（≥4個）がなければスキップ（戻り値: 2）\n"
              << "  - PMT を含むグループが3未満なら その mPMT はスキップ\n"
              << "  - mPMT_pmtid が 0..19 範囲外の場合、その PMT は無視されます\n"
              << "  - σ_time は全 PMT で 1ns（ハードコード）\n"
              << "===============================================================================" << std::endl;
}

int main(int argc, char *argv[]) {
    // 引数チェック
    if (argc < 3) {
        printUsage(argc > 0 ? argv[0] : "reconstructor");
        return 1;
    }

    std::string inputfilename = argv[1];
    std::string outputcsvfile = argv[2];
    int eventNumber = 0;
    
    if (argc > 3) {
        try {
            eventNumber = std::stoi(argv[3]);
        } catch (const std::exception &e) {
            std::cerr << "Error: イベント番号は整数で指定してください。\n";
            printUsage(argv[0]);
            return 1;
        }
    }

    std::map<int, std::vector<PMTData>> PMTDataMap; // Map to store PMT data keyed by mPMTid
    int status = readData(inputfilename, PMTDataMap);
    if (status == 1) {
        std::cout << "Error reading data for event " << std::endl;
        return 1; // Skip this event if there was an error
    }
    if (status == 2) {
        std::cout << "No data for event " << std::endl;
        return 2;
    }
    if (status == 3) {
        std::cout << "all data was read " << std::endl;
        return 3;
    }
    // std::cout << "Data read successfully." << std::endl;
    //  Perform fitting for each mPMTid
    std::cout << "PMTDataMap size: " << PMTDataMap.size() << std::endl;
    int usePMTnum = 0;

    std::vector<SensorUnit> sensorUnits; // Vector to store fitted sensor units
    std::vector<SensorUnit> PMTsUnits;   // Vector to store 3-inch PMTs
    SensorUnit sensorUnit;
    SensorUnit PMTsUnit;
    for (const auto &pair : PMTDataMap) {
        int mPMTid = pair.first;
        const std::vector<PMTData> &pmtData = pair.second;
        if (pmtData.empty()) {
            std::cout << "\nNo data for mPMTid: " << mPMTid << std::endl;
            continue;
        }

        // filter and grouping 8 cm PMT data in 1 mPMT
        std::vector<PMTData> pmtdata_use = findExpandedGroups(pmtData, 5, 3, 0.5, 6); // Optional: filter or group data
        if (pmtdata_use.size() < 3) {
            continue;
        }

        for (int PMT_i = 0; PMT_i < pmtdata_use.size(); ++PMT_i) {
            PMTsUnit.id = pmtdata_use[PMT_i].tubeid; // PMT id
            PMTsUnit.posx = pmtdata_use[PMT_i].x;    // PMT x position
            PMTsUnit.posy = pmtdata_use[PMT_i].y;    // PMT y position
            PMTsUnit.posz = pmtdata_use[PMT_i].z;    // PMT z position
            PMTsUnit.time = pmtdata_use[PMT_i].t;    // PMT time
            PMTsUnit.L = pmtdata_use[PMT_i].L;
            PMTsUnit.sigma_time = 1;       // PMT time
            PMTsUnits.push_back(PMTsUnit); // Add PMT unit to the vector
        }

        // std::cout << "\nFitting data for mPMTid: " << mPMTid << ", Number of PMTs: " << pmtdata_use.size() << std::endl;

        double fit_theta, fit_phi, err_theta, err_phi;
        FitPosition(pmtdata_use, fit_theta, fit_phi, err_theta, err_phi);
        double errby_phi = err_phi * sin(fit_theta);
        std::cout << "error angle is " << err_theta << " , " << err_phi << std::endl;
        double timesum = 0;
        for (const auto &pmt : pmtdata_use) {
            timesum += pmt.t;
        }

        int target_mPMTid = mPMTid;
        double ori_x, ori_y, ori_z;
        ori_x = pmt19orientations[mPMTid][0];
        ori_y = pmt19orientations[mPMTid][1];
        ori_z = pmt19orientations[mPMTid][2];
        // std::cout << "{" << ori_x << ", " << ori_y << ", " << ori_z << "}," << std::endl;

        double phi0_x, phi0_y, phi0_z;
        phi0_x = pmt1orientations[mPMTid][0];
        phi0_y = pmt1orientations[mPMTid][1];
        phi0_z = pmt1orientations[mPMTid][2];
        // std::cout << " {" << phi0_x << ", " << phi0_y << ", " << phi0_z << "}," << std::endl;

        double phi90_x, phi90_y, phi90_z;
        phi90_x = pmt4orientations[mPMTid][0];
        phi90_y = pmt4orientations[mPMTid][1];
        phi90_z = pmt4orientations[mPMTid][2];
        // std::cout << " {" << phi90_x << ", " << phi90_y << ", " << phi90_z << "}," << std::endl;

        double fit_ori = cos(fit_theta);                  // ori方向成分
        double fit_phi0 = sin(fit_theta) * cos(fit_phi);  // phi0方向成分
        double fit_phi90 = sin(fit_theta) * sin(fit_phi); // phi90方向成分

        double dir_x = fit_ori * ori_x + fit_phi0 * phi0_x + fit_phi90 * phi90_x;
        double dir_y = fit_ori * ori_y + fit_phi0 * phi0_y + fit_phi90 * phi90_y;
        double dir_z = fit_ori * ori_z + fit_phi0 * phi0_z + fit_phi90 * phi90_z;

        sensorUnit.id = pmtdata_use[0].mPMTid;     // mPMTid of the mPMT
        sensorUnit.posx = pmtdata_use[0].center_x; // Assuming center_x is the position of the mPMT
        sensorUnit.posy = pmtdata_use[0].center_y; // Assuming center_y is the position of the mPMT
        sensorUnit.posz = pmtdata_use[0].center_z; // Assuming center_z is the position of the mPMT
        sensorUnit.dirx = dir_x;                   // Direction vector x-component
        sensorUnit.diry = dir_y;                   // Direction vector y-component
        sensorUnit.dirz = dir_z;                   // Direction vector z-component
        sensorUnit.time = timesum / pmtdata_use.size();
        // sensorUnit.sigma_sintheta = sqrt(err_theta * err_theta + errby_phi * errby_phi);
        sensorUnit.sigma_sintheta = err_theta; // Assuming sigma_sintheta is the error in theta
        sensorUnit.sigma_time = 1;
        sensorUnits.push_back(sensorUnit);
        usePMTnum++;
    }
    std::cout << "Fitting " << sensorUnits.size() << " mPMTs." << std::endl;
    std::cout << "Using " << usePMTnum << " PMTs for fitting." << std::endl;

    // Skip event if not enough number of PMTs hit
    if (sensorUnits.size() < 4) {
        // std::cout << "Not enough sensor units for fitting. Skipping event " << ev << std::endl;
        PMTDataMap.clear();  // Clear PMT data map for the next event
        PMTsUnits.clear();   // Clear PMT units for the next event
        sensorUnits.clear(); // Clear sensor units for the next event
        return 2;            // Skip this event if no sensor units were created
    }

    // fitting with light_source_fit
    TVector3 fit_pos, fiterr_pos;
    double fit_time, fiterr_time, chi2;
    FitLightSource(sensorUnits, PMTsUnits, fit_pos, fit_time, fiterr_pos, fiterr_time, chi2);

    // std::cout << "\nTotal time taken: " << duration.count() * 1000 << " ms." << std::endl;
    // std::cout << "Fitted position: (" << fit_pos.X() << ", " << fit_pos.Y() << ", " << fit_pos.Z() << ")" << std::endl;
    // std::cout << "Fitted time: " << fit_time << " ns" << std::endl;
    // std::cout << "Position error: (" << fiterr_pos.X() << ", " << fiterr_pos.Y() << ", " << fiterr_pos.Z() << ")" << std::endl;
    // std::cout << "Time error: " << fiterr_time << " ns" << std::endl;
    // std::cout << "Chi2: " << chi2 << std::endl;
    // std::cout << "Chi2/NDF: " << chi2 / (2 * sensorUnits.size() - 4) << std::endl;

    writeToCSV(
        outputcsvfile,
        fit_pos,
        fit_time,
        fiterr_pos,
        fiterr_time,
        chi2);

    PMTDataMap.clear();  // Clear PMT data map for the next event
    PMTsUnits.clear();   // Clear PMT units for the next event
    sensorUnits.clear(); // Clear sensor units for the next event

    return 0;
}

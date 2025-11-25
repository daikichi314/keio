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

int main() {
    std::string inputfilename = "/home/fukazawa/disk3/workdir_1/MCprod/e-/results/mom350/sub/e-25-all";
    std::string outputcsvfile = "/home/fukazawa/disk3/hoge.csv";

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

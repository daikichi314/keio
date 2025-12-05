// onemPMTfit.cc
// r theta phi 球座標系でフィッティング 現行版
#include "onemPMTfit.hh"
#include "TBox.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TEllipse.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TObject.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TTree.h"
#include "fittinginput.hh"
#include <TMinuit.h>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <stdio.h>
#include <stdlib.h>

// #include "WCSimRootPMT.hh"

// Global constants and variables
static double t_light;
static double error_min, fitting_accuracy;
static double firsthit_time;

static TFile *outputrootfile;
static TTree *outputtree;

// Prepare a vector to store PMT data
static std::vector<PMTData> g_pmtData;

// 2つのベクトルがなす角を計算する関数
double CalculateAngle(double x1, double y1, double z1, double x2, double y2, double z2) {
    // ベクトルの内積
    double dotProduct = x1 * x2 + y1 * y2 + z1 * z2;

    // ベクトルの大きさ
    double magnitude1 = std::sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    double magnitude2 = std::sqrt(x2 * x2 + y2 * y2 + z2 * z2);

    // ベクトルの大きさが0の場合はエラー
    if (magnitude1 == 0 || magnitude2 == 0) {
        std::cerr << "Error: One or both vectors have zero magnitude." << std::endl;
        return -1.0; // エラー値
    }

    // なす角（ラジアン）の計算
    double cosalpha = dotProduct / (magnitude1 * magnitude2);
    // 丸め誤差対策（cosThetaが[-1, 1]の範囲を超えることを防ぐ）
    if (cosalpha > 1.0)
        cosalpha = 1.0;
    if (cosalpha < -1.0)
        cosalpha = -1.0;

    return cosalpha;
}

// 球座標系と直交座標系の変換関数
std::array<double, 3> ConvertToSpherical(
    double x, double y, double z,
    double x0, double y0, double z0 // 球座標の原点
) {
    double dx = x - x0, dy = y - y0, dz = z - z0;
    double r = std::sqrt(dx * dx + dy * dy + dz * dz);
    double theta = std::acos(-dy / r);
    double phi = std::atan2(dz, dx);
    return {r, theta, phi};
}
std::array<double, 3> ConvertToCartesian(
    double r, double theta, double phi,
    double x0, double y0, double z0 // 球座標の原点
) {
    double x = r * std::sin(theta) * std::cos(phi) + x0;
    double y = -r * std::cos(theta) + y0;
    double z = r * std::sin(theta) * std::sin(phi) + z0;
    return {x, y, z};
}

// 角度の正規化関数 
// θを[0, π]に、φを[-π, π]に収める
void NormalizeAngles(double &theta, double &phi) {
    if (theta < 0) {
        theta = -theta;
        phi += TMath::Pi();
    }
    if (theta > TMath::Pi()) {
        theta = 2 * TMath::Pi() - theta;
        phi += TMath::Pi();
    }
    // φ を [-π, π] に戻す
    while (phi <= -TMath::Pi())
        phi += 2 * TMath::Pi();
    while (phi > TMath::Pi())
        phi -= 2 * TMath::Pi();
}

/**
 * @brief PMTヒット群から時間的にまとまったグループを探し、前後に拡張する。
 *
 * @param pmtData        PMTデータのリスト
 * @param tau             グルーピングに使う時間幅（単位：ns）
 * @param n               tauの幅の中でこの数以上ヒットがあればグループとみなす
 * @param expand_before   グループ発見時に前に拡張する時間幅（ns）
 * @param expand_after    グループ発見時に後に拡張する時間幅（ns）
 * @return std::vector<PMTData> グループ化されたPMTデータのリスト
 */


// グルーピングと拡張を行う関数
// 1. ヒットを時刻順にソート
// 2. 各ヒットについて、そのヒットを中心にtau秒以内の他のヒットをカウント
// 3. ヒット数がn以上なら、そのヒットを中心に前後にexpand_before, expand_after秒拡張した範囲のヒットをグループに追加
// 4. グループに追加したヒットは再度グループに使わないようにする        
std::vector<PMTData> findExpandedGroups(const std::vector<PMTData> &pmtData,
                                        double tau,
                                        int n,
                                        double expand_before,
                                        double expand_after) {

    // ヒットを時刻順にpmtdataを並び替え（元の順序は保持されない）
    std::vector<PMTData> sorted_hits = pmtData;
    std::sort(sorted_hits.begin(), sorted_hits.end(),
              [](const PMTData &a, const PMTData &b) { return a.t < b.t; });

    std::set<int> used;             // すでにグループに使われたヒットのインデックス
    std::vector<PMTData> true_hits; // 実際に使うヒット

    // unuseがあるかのチェック変数
    int unuse = 0;
    // ヒットごとにスキャン
    for (size_t i = 0; i < sorted_hits.size(); ++i) {
        //  すでにグループに使われたらスキップ
        if (used.count(i)) {
            // std::cout << "this is used" << std::endl;
            continue;
        }

        double time_i = sorted_hits[i].t;
        // 現在のヒットを中心に、tau秒以内の他のヒットをカウント
        int count = 0;
        for (size_t j = i; j < sorted_hits.size(); ++j) {
            // std::cout << " checking " << j << ", and this time is " << sorted_hits[j].t << std::endl;
            if (sorted_hits[j].t - time_i > tau) {
                // std::cout << "tau is over, tau is " << tau << std::endl;
                continue; // tau秒を超えたら終了
            }
            count++;
            // ヒット数がn未満ならスキップ（グループとみなさない）
        }
        if (count < n) {
            count = 0;
            continue;
        }

        // グループと認められたので、前後に拡張する
        double t_start = time_i - expand_before;
        double t_end = time_i + expand_after;
        // 拡張した時間範囲にあるヒットをグループに追加
        for (size_t j = 0; j < sorted_hits.size(); ++j) {
            if (sorted_hits[j].t >= t_start && sorted_hits[j].t <= t_end) {
                true_hits.push_back(sorted_hits[j]);
                used.insert(j); // グループに使ったヒットを記録
            }
        }
        break;
    }
    return true_hits;
}
// void writeToRoot(
//     TTree *tree,                 // TTree object
//     std::vector<double> &fits,   // fit results
//     double t_light,              // t_light
//     std::vector<double> &errors, // errors of fits
//     double t_error,              // t_error (0 if t is fixed)
//     double minimized,            // -error_min
//     double hitnumber             // how many PMT hits
// ) {
//     double branchFits[3];
//     double branchErrors[3];
//     double branchTLight = t_light;
//     double branchTError = t_error;
//     double branchG = minimized;
//     double branchhitnumber = hitnumber;
//     for (int i = 0; i < 3; ++i) {
//         branchFits[i] = fits[i];
//         branchErrors[i] = errors[i];
//     }

//     double p = TMath::Prob(minimized, hitnumber - 4);

//     // std::cout << "r, theta, phi = " << branchFits[0] << ", " << branchFits[1] << ", " << branchFits[2] << std::endl;
//     // std::cout << "errors = " << branchErrors[0] << ", " << branchErrors[1] << ", " << branchErrors[2] << std::endl;
//     // std::cout << "write :: minimized is " << minimized << std::endl;
//     tree->SetBranchAddress("fits", branchFits);
//     tree->SetBranchAddress("t_light", &branchTLight);
//     tree->SetBranchAddress("errors", branchErrors);
//     tree->SetBranchAddress("t_error", &branchTError);
//     tree->SetBranchAddress("chi2", &branchG);
//     tree->SetBranchAddress("hitnumber", &branchhitnumber);
//     tree->SetBranchAddress("p", &p);
//     tree->Fill();
// }

double CalculateError(double *params) {
    double maxL = 0;
    double chi2 = 0.0;
    for (const auto &pmt : g_pmtData) {
        NormalizeAngles(params[0], params[1]);
        std::array<double, 3> Cartesian = ConvertToCartesian(100, params[0], params[1], pmt.center_x, pmt.center_y, pmt.center_z);
        double dx, dy, dz;
        // 無限遠光源の場合
        dx = std::sin(params[0]) * std::cos(params[1]);
        dy = -std::cos(params[0]);
        dz = std::sin(params[0]) * std::sin(params[1]);
        // dx = -pmt.x + Cartesian[0];
        // dy = -pmt.y + Cartesian[1];
        // dz = -pmt.z + Cartesian[2]; // rを使う場合はこっち
        int id = pmt.mPMT_pmtid;
        // int id =
        if (id < 0 || id > 19) {
            std::cout << "Error: PMT ID out of range." << std::endl;
            continue; // Skip this PMT if the ID is invalid
        }
        double ori_x = pmtDirections[id][0];
        double ori_y = pmtDirections[id][1];
        double ori_z = pmtDirections[id][2];
        double cosalpha = CalculateAngle(ori_x, ori_y, ori_z, dx, dy, dz); // Angle between the PMT orientation and the light direction
        double distance2 = 60;
        double expectedTime;
        double model = 0;
        model = (params[2] / 10000 / (1 + std::exp(-6 * (cosalpha - 1))) + params[3] / 100000);// -6は測定結果によって変更が必要かも N_model
        if (model < 0) {
            model = 0;
        }
        chi2 += (pmt.L - model * distance2) * (pmt.L - model * distance2) / pmt.L;// chi2
        // chi2 += (pmt.L - model) * (pmt.L - model) / pmt.L;
    }
    return chi2;
}

static void FcnForMinuit(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag) {
    fval = CalculateError(par); // Compute the residual sum of squares
}

int fitting(double &fit_theta, double &fit_phi, double &err_theta, double &err_phi,
            double &minimized) {
    TMinuit minuit(4);// rを増やすなら+1
    // Set FCN function (fit function)
    minuit.SetFCN(FcnForMinuit);// Set the function to be minimized
    minuit.SetPrintLevel(-1);// Suppress output

    double maxL = 0;
    for (const auto &pmt : g_pmtData) {
        if (pmt.L > maxL)
            maxL = pmt.L;
    }
    maxL = maxL * 10;
    if (maxL > 800)
        maxL = 800; // Set a maximum limit for maxL to avoid overflow in fitting
    // minuit.DefineParameter(5, "r_light", 60, 1, 3, 400.0);  // rのやつ                                 // 3 to 400
    minuit.DefineParameter(0, "theta_light", 0, 0.1, (-TMath::Pi() / 2), (TMath::Pi() / 2)); // 0 to pi/2
    minuit.DefineParameter(1, "phi_light", 0, 0.3, -TMath::Pi() * 2, TMath::Pi() * 2);       // -pi to pi times 2
    minuit.DefineParameter(2, "A", maxL, 0.1, 0, 1000);
    minuit.DefineParameter(3, "B", 30, 0.1, 0, 100);
    // minuit.DefineParameter(4, "C", 5.5, 0.1, 3, 10);

    // No need for t_light here as it is a loop variable

    // Start the fitting process
    int status = minuit.Migrad();// Perform the minimization
    minuit.GetParameter(0, fit_theta, err_theta);// Retrieve fitted parameters and their errors
    minuit.GetParameter(1, fit_phi, err_phi);// Retrieve fitted parameters and their errors

    double fval, fedm, errdef;
    int nvpar, nparx, istat;
    minuit.mnstat(fval, fedm, errdef, nvpar, nparx, istat);
    minimized = fval;
    // std::cout << "minimized is " << minimized << std::endl;

    minuit.Clear();
    return status;// あんまりいらないかも
}

// Function to perform the fitting process
void FitPosition(std::vector<PMTData> &pmtData,
                 double &fit_theta,
                 double &fit_phi,
                 double &err_theta,
                 double &err_phi) {
    g_pmtData = pmtData; // グローバル変数にデータをセット
    // ori_mPMTを得る
    if (g_pmtData.empty()) {
        std::cerr << "Error: No PMT data available for fitting." << std::endl;
        return;
    }

    double minimized = 0;
    int status = fitting(fit_theta, fit_phi, err_theta, err_phi, minimized);
    if (status != 0) {
         std::cout << "Fitting failed with status: " << status << std::endl;
        return;
    }
    NormalizeAngles(fit_theta, fit_phi);
    // std::cout << "fits is " << fit_theta << " " << fit_phi << std::endl;
    // std::cout << "err is " << err_theta << " " << err_phi << std::endl;

    // outputrootfile->cd();
    // std::vector<double> fits = {0, fit_theta, fit_phi};
    // std::vector<double> errors = {0, err_theta, err_phi};
    // double t_light = 0; // t_light is not used in this context, set to 0
    // writeToRoot(outputtree, fits, t_light, errors, 0, minimized, pmtData.size());
}

int onemPMTfit() {
    outputrootfile = new TFile("fittingResults.root", "RECREATE");
    outputtree = new TTree("fitResults", "Fit results");
    outputtree->Branch("fits", (void *)nullptr, "fits[3]/D");
    outputtree->Branch("t_light", (void *)nullptr, "t_light/D");
    outputtree->Branch("errors", (void *)nullptr, "errors[3]/D");
    outputtree->Branch("t_error", (void *)nullptr, "t_error/D");
    outputtree->Branch("chi2", (void *)nullptr, "G/D");
    outputtree->Branch("hitnumber", (void *)nullptr, "hitnumber/D");
    outputtree->Branch("p", (void *)nullptr, "p/D");

    std::string inputfilename = "/home/fukazawa/disk3/workdir_1/MCprod/e-/results/Spher-30-4-all";
    TFile *file = TFile::Open((inputfilename + ".root").c_str());
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file pmtdata.root" << std::endl;
        return 1;
    }

    TTree *tree = dynamic_cast<TTree *>(file->Get("PMTTree"));
    if (!tree) {
        std::cerr << "Error: Cannot find TTree 'PMTTree' in file." << std::endl;
        return 1;
    }

    // ブランチ変数の定義
    int eventNumber, id, tubeid, mPMTid, mPMT_pmtid;
    double x, y, z, L, t, ori_x, ori_y, ori_z, center_x, center_y, center_z;
    tree->SetBranchAddress("eventNumber", &eventNumber);
    // tree->SetBranchAddress("id", &id);
    //   tree->SetBranchAddress("tubeid", &tubeid);
    tree->SetBranchAddress("mPMTid", &mPMTid);
    tree->SetBranchAddress("mPMT_pmtid", &mPMT_pmtid);
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("z", &z);
    tree->SetBranchAddress("L", &L);
    tree->SetBranchAddress("t", &t);
    tree->SetBranchAddress("ori_x", &ori_x);
    tree->SetBranchAddress("ori_y", &ori_y);
    tree->SetBranchAddress("ori_z", &ori_z);
    tree->SetBranchAddress("center_x", &center_x);
    tree->SetBranchAddress("center_y", &center_y);
    tree->SetBranchAddress("center_z", &center_z);

    Long64_t nEntries = tree->GetEntries();
    if (nEntries == 0) {
        std::cerr << "Error: No entries in TTree." << std::endl;
        return 1;
    }
    tree->GetEntry(nEntries - 1);
    int eventNumber_max = eventNumber + 1;
    PMTData data;
    int data_now = 0;
    int hit_PMT = 0;
    for (int ev = 0; ev < 1000; ev++) {
        hit_PMT = 0;
        for (Long64_t i = data_now; i < nEntries; ++i) {
            tree->GetEntry(i);
            // std::cout << "tree->GetEntry(" << i << ") done" << std::endl;
            if (eventNumber == ev && mPMTid == 347) { // 347 is the mPMT ID for this example
                // if (eventNumber == ev) {
                data_now = i;
                hit_PMT++;
                // data.tubeid = tubeid;
                // data.mPMTid = mPMTid;
                data.mPMT_pmtid = mPMT_pmtid;
                // data.mPMT_pmtid = id - 6574; // Assuming id is the PMT ID in m-PMT
                data.t = t;
                data.x = x;
                data.y = y;
                data.z = z;
                data.L = L;
                data.ori_x = ori_x;
                data.ori_y = ori_y;
                data.ori_z = ori_z;
                data.center_x = center_x;
                data.center_y = center_y;
                data.center_z = center_z;
                g_pmtData.push_back(data);
            }
            if (eventNumber > ev) {
                break;
            }
        }
        if (hit_PMT > 3) {
            // std::cout << "hit_PMT is " << hit_PMT << std::endl;
            g_pmtData = findExpandedGroups(g_pmtData, 5, 3, 0.5, 6);
            double fit_theta, fit_phi;
            double err_theta, err_phi;
            FitPosition(g_pmtData, fit_theta, fit_phi, err_theta, err_phi);
        }
        g_pmtData.clear();
    }
    outputrootfile->cd();
    outputtree->Write();
    //}
    //}
    return 0;
}
#include "readData.hh"
#include "TFile.h"
#include "TTree.h"
#include "fittinginput.hh"
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <set>

int readData(std::string &inputfilename,  std::map<int, std::vector<PMTData>> &PMTDataMap) {
    int ev=0;
    int firstentry=0;
    int lastentry; 
    int eventNumber=0;
    TFile *file = TFile::Open((inputfilename + ".root").c_str());
	std::cout << "open file "<<file<<std::endl;
    if (!file || file->IsZombie()) {
        std::cout << "Error: Cannot open file pmtdata.root" << std::endl;
        return 1;
    }

    TTree *tree = dynamic_cast<TTree *>(file->Get("PMTTree"));
    if (!tree) {
        std::cout << "Error: Cannot find TTree 'PMTTree' in file." << std::endl;
        return 1;
    }

    // ブランチ変数の定義
    int tubeid, mPMTid, mPMT_pmtid;
    double x, y, z, L, t, ori_x, ori_y, ori_z, center_x, center_y, center_z;
    tree->SetBranchAddress("eventNumber", &eventNumber);
    tree->SetBranchAddress("tubeid", &tubeid);
    tree->SetBranchAddress("mPMTid", &mPMTid);
    tree->SetBranchAddress("mPMT_pmtid", &mPMT_pmtid);
    tree->SetBranchAddress("x", &x);// PMTnoiti
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

    // エントリー数の取得
    std::cout << "branch is read" << std::endl;
    Long64_t nEntries = tree->GetEntries();
    if (nEntries == 0) {
        std::cout << "Error: No entries in TTree." << std::endl;
        return 1;
    }
    tree->GetEntry(nEntries - 1);
    int eventNumber_max = eventNumber + 1; // 最大イベント番号を取得
    // もし全イベントを読み終えたら終了
    if (ev >= eventNumber_max) {
        std::cout << "all event is readed, eventNumber_max is " << eventNumber_max << std::endl;
        return 3;
    }
    PMTData data; // PMTData構造体のインスタンスを作成

    // イベントループ データの読み込み 
    for (Long64_t i =0 ; i < nEntries; ++i) {
        tree->GetEntry(i);
        // use mPMT id is {338, 339,340, 346, 347,348, 354,355,356}; the case of the data 30-4-all center is 40,40
        if (eventNumber == ev && (mPMTid == 338 || mPMTid == 339 || mPMTid == 340 || mPMTid == 346 || mPMTid == 347 || mPMTid == 348 || mPMTid == 354 || mPMTid == 355 || mPMTid == 356)) {
            // use uPMT id is {330, 331,332,339, 340,341, 347, 348,349}
            // if (eventNumber == ev && (mPMTid == 330 || mPMTid == 331 || mPMTid == 332 || mPMTid == 339 || mPMTid == 340 || mPMTid == 341 || mPMTid == 347 || mPMTid == 348 || mPMTid == 349)) {
            // if (eventNumber == ev && mPMTid > 317) {
            //  if (eventNumber == ev) {
            data.tubeid = tubeid;
            data.mPMTid = mPMTid;
            data.mPMT_pmtid = mPMT_pmtid;
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
            PMTDataMap[mPMTid].push_back(data);
        }
        if (eventNumber > ev) {
			lastentry = i;
            break;
        }
    }
    // PMTDataMap[i]が空かどうかを確認
    for (const auto &pair : PMTDataMap) {
        if (pair.second.empty()) {
            std::cout << "No data for mPMTid: " << pair.first << std::endl;
        }
    }
    std::cout << "PMTDataMap size: " << PMTDataMap.size() << std::endl;
	file->Close();
	if(PMTDataMap.size() == 0) return 2;
    return 0;
}

int readPMTinfo(std::string &inputfilename, int &target_mPMT_id, int &target_mPMT_pmtid, double &ori_x_return, double &ori_y_return, double &ori_z_return) {
    std::ifstream infile(inputfilename);
    if (!infile) {
        std::cerr << "Error: Cannot open file " << inputfilename << std::endl;
        return 1;
    }

    std::string line;
    int line_num = 0;

    // skip the first 5 lines
    for (int i = 0; i < 5 && std::getline(infile, line); ++i)
        ;

    bool found = false;

    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        int tube_id, mPMT_id, mPMT_pmtid, info;
        double x, y, z, ori_x, ori_y, ori_z;

        if (!(iss >> tube_id >> mPMT_id >> mPMT_pmtid >> x >> y >> z >> ori_x >> ori_y >> ori_z >> info)) {
            std::cerr << "読み取りエラー： " << line << std::endl;
            continue;
        }

        if (mPMT_id == target_mPMT_id && mPMT_pmtid == target_mPMT_pmtid) {
            ori_x_return = ori_x;
            ori_y_return = ori_y;
            ori_z_return = ori_z;
            found = true;
            break; // 最初に一致したものだけを取得する場合
        }
    }

    if (!found) {
        std::cout << "一致するデータが見つかりませんでした。" << std::endl;
    }
    return found ? 0 : 1; // 成功なら0、失敗なら1を返す
}

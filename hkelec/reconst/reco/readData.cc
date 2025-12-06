/*
 * id: readData.cc
 * Place: /home/daiki/keio/hkelec/reconst/reco/
 * Author: Gemini 3 Pro
 * Last Edit: 2025-12-06
 *
 * 概要:
 * データ読み込みとペデスタル補正の実装ファイル
 * ROOT ファイルの processed_hits TTree からヒットデータを読み込み、
 * ペデスタル補正を行って ADC 値を pC 単位の電荷に変換します。
 * イベント単位での逐次読み込み機能を提供します。
 *
 * 主な機能:
 * - ペデスタルファイル読み込み (CSV形式)
 * - ADC → pC 変換 (High Gain 飽和検出、Low Gain 使用切替)
 * - イベント単位でのヒットデータグルーピング
 */

#include "readData.hh"

// ADC -> pC 変換係数
const double K_HGAIN = 0.073; 
const double K_LGAIN = 0.599; 
const double SATURATION_THRESHOLD = 4000.0; // lgainに切り替える閾値

// ペデスタル読み込み
int readPedestals(const std::string &filename, std::map<int, PedestalData> &pedestalMap) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error: Cannot open pedestal file " << filename << std::endl;
        return 1;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        for (char &c : line) if (c == ',') c = ' '; // CSVパース用
        
        std::istringstream iss(line);
        int ch;
        std::string type;
        double mean, err;
        
        if (!(iss >> ch >> type >> mean >> err)) continue;

        if (type == "hgain") pedestalMap[ch].hgain_mean = mean;
        else if (type == "lgain") pedestalMap[ch].lgain_mean = mean;
    }
    return 0;
}

// コンストラクタ
DataReader::DataReader(const std::string &filename, const std::map<int, PedestalData> &pedMap) 
    : pedestalMap(pedMap), currentEntry(0), hasBufferedHit(false) {
    
    file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file " << filename << std::endl;
        tree = nullptr;
        nEntries = 0;
        return;
    }

    tree = dynamic_cast<TTree*>(file->Get("processed_hits"));
    if (!tree) {
        std::cerr << "Error: Cannot find 'processed_hits' tree" << std::endl;
        nEntries = 0;
        return;
    }

    tree->SetBranchAddress("eventID", &b_eventID);
    tree->SetBranchAddress("ch", &b_ch);
    tree->SetBranchAddress("hgain", &b_hgain);
    tree->SetBranchAddress("lgain", &b_lgain);
    tree->SetBranchAddress("tot", &b_tot);
    tree->SetBranchAddress("time_diff", &b_time_diff);

    nEntries = tree->GetEntries();
}

DataReader::~DataReader() {
    if (file) file->Close();
}

// データ変換ロジック
PMTData DataReader::createPMTData() {
    PMTData data;
    data.eventID = b_eventID;
    data.ch = b_ch;
    data.time = b_time_diff;

    // ペデスタル取得
    double ped_h = 0, ped_l = 0;
    auto it = pedestalMap.find(b_ch);
    if (it != pedestalMap.end()) {
        ped_h = it->second.hgain_mean;
        ped_l = it->second.lgain_mean;
    }

    // 電荷計算 (サチュレーション考慮)
    if (b_hgain >= SATURATION_THRESHOLD) {
        data.charge = (b_lgain - ped_l) * K_LGAIN;
    } else {
        data.charge = (b_hgain - ped_h) * K_HGAIN;
    }
    if (data.charge < 0) data.charge = 0;

    // 座標セット
    if (b_ch >= 0 && b_ch < 4) {
        data.x = PMT_POSITIONS[b_ch][0];
        data.y = PMT_POSITIONS[b_ch][1];
        data.z = PMT_POSITIONS[b_ch][2];
        data.dir_x = PMT_DIR[0];
        data.dir_y = PMT_DIR[1];
        data.dir_z = PMT_DIR[2];
    } else {
        // 無効なチャンネルの場合
        data.x = 0; data.y = 0; data.z = 0;
    }

    return data;
}

// 次のイベント読み出し
bool DataReader::nextEvent(std::vector<PMTData> &eventHits) {
    eventHits.clear();
    if (!tree) return false;

    // バッファ分を処理
    if (hasBufferedHit) {
        b_eventID = buf_eventID;
        b_ch = buf_ch;
        b_hgain = buf_hgain;
        b_lgain = buf_lgain;
        b_tot = buf_tot;
        b_time_diff = buf_time_diff;
        
        // 【修正】有効なチャンネル(0-3)のみ追加
        if (b_ch >= 0 && b_ch < 4) {
            eventHits.push_back(createPMTData());
        }
        hasBufferedHit = false;
    }

    int currentEventID = -1;
    if (!eventHits.empty()) currentEventID = eventHits[0].eventID;

    while (currentEntry < nEntries) {
        tree->GetEntry(currentEntry);
        currentEntry++;

        // まだヒットリストが空なら、現在のIDをイベントIDとする
        if (eventHits.empty()) currentEventID = b_eventID;

        // IDが変わったらバッファしてリターン
        if (b_eventID != currentEventID) {
            hasBufferedHit = true;
            buf_eventID = b_eventID;
            buf_ch = b_ch;
            buf_hgain = b_hgain;
            buf_lgain = b_lgain;
            buf_tot = b_tot;
            buf_time_diff = b_time_diff;
            return true;
        }

        // 【修正】有効なチャンネル(0-3)のみ追加
        if (b_ch >= 0 && b_ch < 4) {
            eventHits.push_back(createPMTData());
        }
    }

    return !eventHits.empty();
}
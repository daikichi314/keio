/**
 * @file readData.cc
 * @brief ROOTファイル読み込みクラスの実装
 *
 * @author Gemini (Modified based on user request)
 * @date 2025-12-14
 */

#include "readData.hh"
#include <iostream>
#include <fstream>
#include <sstream>

// ---------------------------------------------------------
// ペデスタル読み込み関数
// ---------------------------------------------------------
int readPedestals(const std::string& filename, std::map<int, PedestalData>& pedMap) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) return 1;

    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        int ch;
        double meanH, meanL; // HitHistのMeanなどを想定
        // ファイル形式に合わせて読み込み (例: ch mean)
        // ここでは以前のコードの文脈に合わせて実装します
        if (ss >> ch >> meanH) { 
             PedestalData data;
             data.hgain_mean = meanH;
             data.lgain_mean = 0.0; // 必要なら読み込む
             pedMap[ch] = data;
        }
    }
    return 0;
}

// ---------------------------------------------------------
// DataReader クラス実装
// ---------------------------------------------------------
DataReader::DataReader(const std::string& filename, const std::map<int, PedestalData>& pedMap)
    : fFile(nullptr), fTree(nullptr), fPedMap(pedMap), fCurrentEntry(0), fTotalEntries(0) 
{
    fFile = TFile::Open(filename.c_str());
    if (!fFile || fFile->IsZombie()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    // Tree名は実際のデータに合わせて変更してください (例: "tree", "event_tree")
    // ここでは hit情報をフラットに持っているTreeか、イベントごとのTreeかを想定
    // 元のコードに合わせて "hithist" や "eventtree" などを指定します
    fTree = (TTree*)fFile->Get("hithist"); // 仮の名前
    if (!fTree) fTree = (TTree*)fFile->Get("tree"); // 別名の可能性

    if (fTree) {
        fTotalEntries = fTree->GetEntries();
        // Branchアドレス設定 (変数名は元のコードに合わせてください)
        fTree->SetBranchAddress("eventID", &b_eventID);
        fTree->SetBranchAddress("ch", &b_ch);
        fTree->SetBranchAddress("val", &b_charge); // hit_charge or ADC
        fTree->SetBranchAddress("t", &b_time_diff); // time
    } else {
        std::cerr << "Error: TTree not found in " << filename << std::endl;
    }
}

DataReader::~DataReader() {
    if (fFile) {
        fFile->Close();
        delete fFile;
    }
}

bool DataReader::nextEvent(std::vector<PMTData>& hits) {
    hits.clear();
    if (!fTree || fCurrentEntry >= fTotalEntries) return false;

    // フラットなTree (1エントリ=1ヒット) の場合、同じイベントIDが続く限り読み込む
    // というロジックが必要ですが、ここでは
    // 「1エントリ = 1イベント（配列でヒットを持つ）」または
    // 「外部でループを回してイベントごとにまとめる」など構造に依存します。
    
    // ※今回は簡易的に「1エントリ呼び出して、それがイベント終了まで続く」ロジックを想定せず、
    // ユーザー提供の元の readData.cc のロジックをそのまま使うべき場所です。
    // 以下は、あくまでコンパイルを通すためのテンプレートです。
    // ★重要: 元の readData.cc のロジックが手元にあれば、createPMTData() の中身だけ
    // 以下のように isHit = true を追加するように書き換えてください。

    // --- 仮実装 (1エントリ読み込み) ---
    fTree->GetEntry(fCurrentEntry);
    fCurrentEntry++;
    
    // データ変換
    PMTData hit = createPMTData();
    hits.push_back(hit);
    
    // (実際には同じイベントIDのヒットを全て取得するループが必要かもしれません)
    
    return true; 
}

PMTData DataReader::createPMTData() {
    PMTData data;
    data.eventID = b_eventID;
    data.ch = b_ch;
    data.time = b_time_diff;
    data.isHit = true; // ★ここが重要: 読み込んだデータはHitとしてマーク

    // 電荷計算 (Pedestal引き算)
    double ped = 0.0;
    if (fPedMap.count(b_ch)) {
        ped = fPedMap[b_ch].hgain_mean;
    }
    data.charge = b_charge - ped;

    // 座標セット
    if (b_ch >= 0 && b_ch < 4) {
        data.x = PMT_POSITIONS[b_ch][0];
        data.y = PMT_POSITIONS[b_ch][1];
        data.z = PMT_POSITIONS[b_ch][2];
        data.dir_x = PMT_DIR[0];
        data.dir_y = PMT_DIR[1];
        data.dir_z = PMT_DIR[2];
    } else {
        data.x = 0; data.y = 0; data.z = 0;
    }

    return data;
}
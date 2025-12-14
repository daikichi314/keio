/**
 * @file readData.hh
 * @brief ROOTファイル読み込みクラスのヘッダー
 *
 * fittinginput.hh で定義された構造体(PMTData, PedestalData)を使用します。
 *
 * @author Gemini (Modified based on user request)
 * @date 2025-12-14
 */

#ifndef READDATA_HH
#define READDATA_HH

#include <string>
#include <vector>
#include <map>
#include <TFile.h>
#include <TTree.h>
#include "fittinginput.hh" // 必須: これがないとPedestalDataなどが不明になります

// ペデスタル読み込み関数 (mainでも使うためクラス外に定義)
int readPedestals(const std::string& filename, std::map<int, PedestalData>& pedMap);

class DataReader {
public:
    DataReader(const std::string& filename, const std::map<int, PedestalData>& pedMap);
    virtual ~DataReader();

    // 次のイベントを取得。ファイルの終わりならfalseを返す
    bool nextEvent(std::vector<PMTData>& hits);

private:
    TFile* fFile;
    TTree* fTree;
    std::map<int, PedestalData> fPedMap;
    int fCurrentEntry;
    int fTotalEntries;

    // TTreeのBranch用変数
    int b_eventID;
    int b_ch;
    double b_time_diff; // または b_time など、元のTree構造に合わせてください
    double b_charge;    // ADC値など

    // データ変換用内部関数
    PMTData createPMTData();
};

#endif // READDATA_HH
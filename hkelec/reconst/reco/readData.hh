/*
 * id: readData.hh
 * Place: /home/daiki/keio/hkelec/reconst/reco/
 * Author: Gemini 3 Pro
 * Last Edit: 2025-12-06
 *
 * 概要:
 * データ読み込み機能のヘッダーファイル
 * DataReader クラスとペデスタル読み込み関数の宣言を含みます。
 * イベント単位での逐次的なデータアクセスを提供します。
 */

#ifndef READ_DATA_HH
#define READ_DATA_HH // インクルードガード

#include "fittinginput.hh" // 共通のデータ構造定義を読み込む
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <TFile.h> // ROOTファイルの操作用
#include <TTree.h> // TTreeの操作用

// ペデスタル情報をテキストファイルから読み込み、マップに格納する関数
int readPedestals(const std::string &filename, std::map<int, PedestalData> &pedestalMap);

// 逐次処理を行うためのデータリーダークラスの定義
class DataReader {
public:
    // コンストラクタ: ファイル名とペデスタルマップを受け取って初期化
    DataReader(const std::string &filename, const std::map<int, PedestalData> &pedMap);
    // デストラクタ: ファイルを閉じるなどの後処理
    ~DataReader();

    // 次のイベントデータを取得する関数
    // eventHitsベクトルにデータを詰め込み、成功ならtrueを返す
    bool nextEvent(std::vector<PMTData> &eventHits);

    // 全エントリー数を返す関数
    long getTotalEntries() const { return nEntries; }
    // 現在の読み込み位置を返す関数
    long getCurrentEntry() const { return currentEntry; }

private:
    TFile *file; // ROOTファイルポインタ
    TTree *tree; // TTreeポインタ
    long nEntries; // TTreeの総エントリー数
    long currentEntry; // 現在読んでいるエントリー番号
    
    // TTreeから読み込むためのブランチ変数
    int b_eventID;
    int b_ch;
    double b_hgain;
    double b_lgain;
    double b_tot;
    double b_time_diff; // ns単位

    // ペデスタルマップへの参照 (コピーせず参照を持つ)
    const std::map<int, PedestalData> &pedestalMap;

    // 次のイベント判定のために先読みしたデータを一時保存する変数群
    bool hasBufferedHit; // バッファにデータがあるかフラグ
    int buf_eventID;
    int buf_ch;
    double buf_hgain;
    double buf_lgain;
    double buf_tot;
    double buf_time_diff;

    // 現在のブランチ変数の値からPMTData構造体を作成する内部関数
    PMTData createPMTData();
};

#endif // READ_DATA_HH
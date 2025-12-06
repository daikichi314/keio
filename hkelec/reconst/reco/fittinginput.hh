/*
 * id: fittinginput.hh
 * Place: /home/daiki/keio/hkelec/reconst/reco/
 * Author: Gemini 3 Pro
 * Last Edit: 2025-12-06
 *
 * 概要:
 * 共通定数とデータ構造の定義ファイル
 * 光速、PMT 位置・向き、時刻補正値などの物理定数を定義し、
 * PMTData、PedestalData、FitResult などのデータ構造体を提供します。
 */

#ifndef FITTINGINPUT_HH
#define FITTINGINPUT_HH

#include <vector>
#include <cmath>

// 定数定義
const double C_LIGHT = 29.970255; // cm/ns (n=1.0003 in air)

// PMTごとの時間補正値 (ns)
// 将来的にはPMT個別の値と誤差を入れる
const double TIME_CORRECTION_VAL[4] = {200.0, 200.0, 200.0, 200.0};
const double TIME_CORRECTION_ERR[4] = {0.0, 0.0, 0.0, 0.0}; // 将来的な誤差考慮用

// 時刻分解能を取得する関数 (電荷依存性を考慮予定)
// 実装は onemPMTfit.cc に記述
double GetSigmaTime(int ch, double charge);

// PMT座標定義 (4 PMTs)
// 向きはすべて真上 (0, 0, 1) と仮定
const double PMT_POSITIONS[4][3] = {
    {-35.0,  35.0, 80.5}, // PMT#1 (CH0)
    { 35.0,  35.0, 80.5}, // PMT#2 (CH1)
    {-35.0, -35.0, 80.5}, // PMT#3 (CH2)
    { 35.0, -35.0, 80.5}  // PMT#4 (CH3)
};

const double PMT_DIR[3] = {0.0, 0.0, 1.0}; // すべて真上を向いている

struct PMTData {
    int eventID;
    int ch;            // チャンネル (0-3)
    double time;       // 補正済み時刻 (ns)
    double charge;     // 変換済み電荷 (pC)
    double x, y, z;    // PMT位置
    double dir_x, dir_y, dir_z; // PMT向き
};

struct PedestalData {
    double hgain_mean;
    double lgain_mean;
};

// フィッティング結果を格納する構造体
struct FitResult {
    double x, y, z, t;
    double err_x, err_y, err_z, err_t;
    double A, B; // 光量モデルパラメータ
    double chi2;
    int ndf;
    int status; // Minuit status
};

#endif // FITTINGINPUT_HH
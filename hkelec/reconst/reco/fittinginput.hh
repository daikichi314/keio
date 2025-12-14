/**
 * @file fittinginput.hh
 * @brief フィッティングに必要な定数、データ構造、設定構造体を定義するヘッダーファイル
 *
 * @author Gemini (Modified based on user request)
 * @date 2025-12-14
 */

#ifndef FITTINGINPUT_HH
#define FITTINGINPUT_HH

#include <vector>
#include <cmath>
#include <string>

// =========================================================
// 定数定義
// =========================================================
const double C_LIGHT = 29.970255; // cm/ns

// PMTごとの時間補正値
const double TIME_CORRECTION_VAL[4] = {200.0, 200.0, 200.0, 200.0};
const double TIME_CORRECTION_ERR[4] = {0.0, 0.0, 0.0, 0.0};

// PMT座標定義 (4 PMTs)
const double PMT_POSITIONS[4][3] = {
    {-35.0,  35.0, 80.5}, // PMT#1 (CH0)
    { 35.0,  35.0, 80.5}, // PMT#2 (CH1)
    {-35.0, -35.0, 80.5}, // PMT#3 (CH2)
    { 35.0, -35.0, 80.5}  // PMT#4 (CH3)
};

// PMTの向き (全PMT共通でZ軸方向を向いていると仮定)
// ※ 個別に設定が必要な場合はPMTData構造体内の値を使用してください
const double PMT_DIR[3] = {0.0, 0.0, 1.0}; 

// =========================================================
// 関数プロトタイプ
// =========================================================
double GetEMG_Sigma(int ch, double charge);
double GetEMG_Tau(int ch, double charge);
double GetTimeResolution(int ch, double charge);

// =========================================================
// データ構造体
// =========================================================
struct PMTData {
    int eventID;
    int ch;
    double time;
    double charge;
    double x, y, z;
    double dir_x, dir_y, dir_z; // PMTの向き
    bool isHit;
};

struct PedestalData {
    double hgain_mean;
    double lgain_mean;
};

// =========================================================
// フィッティング設定
// =========================================================

enum class ChargeChi2Type {
    Gaussian,
    BakerCousins,
    None
};

/**
 * @enum ChargeModelType
 * @brief 電荷の期待値モデル
 */
enum class ChargeModelType {
    Standard,       // 従来: mu = A/r^2 + B (等方発光仮定)
    ZeroIntercept,  // B固定: mu = A/r^2 (B=0)
    Cosine          // 新機能: mu = (A * f(cos_alpha)) / r^2 + B (角度依存あり)
};

enum class TimeChi2Type {
    Gaussian,
    EMG,
    Goodness,
    None
};

struct FitConfig {
    ChargeChi2Type chargeType = ChargeChi2Type::Gaussian;
    ChargeModelType chargeModel = ChargeModelType::Cosine; // デフォルトを角度依存モデルに変更推奨
    TimeChi2Type timeType = TimeChi2Type::Gaussian;
    bool useUnhit = false;
};

struct FitResult {
    double x, y, z, t;
    double err_x, err_y, err_z, err_t;
    double A, B;
    double chi2;
    int ndf;
    int status;
};

#endif // FITTINGINPUT_HH
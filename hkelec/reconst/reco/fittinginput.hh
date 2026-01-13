/**
 * @file fittinginput.hh
 * @brief フィッティングに必要な定数、データ構造、設定構造体を定義するヘッダーファイル
 *
 * このファイルには以下の要素が含まれます：
 * 1. 物理定数 (光速など)
 * 2. 検出器ジオメトリ定義 (PMT座標、向き)
 * 3. 補正パラメータ (TimeWalk, TimeOffset, 電荷モデル係数)
 * 4. 入力データの構造体 (PMTData)
 * 5. 解析設定用の列挙型と構造体 (FitConfig)
 *
 * 新しい電荷モデル (SolidAngle) 用の係数はこのファイルで設定してください。
 *
 * @author Gemini (Modified based on user request)
 * @date 2025-01-08
 */

#ifndef FITTINGINPUT_HH
#define FITTINGINPUT_HH

#include <vector>
#include <cmath>
#include <string>

// =========================================================
// 物理定数・PMT補正値
// =========================================================
const double C_LIGHT = 29.970255; // 光速 [cm/ns] (空気中の屈折率 n=1.0003 を考慮)

// PMTごとの時間補正値 (ns) : ケーブル遅延やT0オフセット
// t_expected の計算に使用されます: t_exp = t0 + tof + TW + CORRECTION
const double TIME_CORRECTION_VAL[4] = {200.0, 200.0, 200.0, 200.0};
const double TIME_CORRECTION_ERR[4] = {0.0, 0.0, 0.0, 0.0};

// =========================================================
// TimeWalk (TW) および 時間分解能 (Sigma_t) のパラメータ
// =========================================================
// 関数形: f(q) = c0 * q^{-1/2} + c1 + c2 * q + c3 * q^2
// 配列の並び: {c0, c1, c2, c3}

// TimeWalk補正係数 [ch][param]
const double TW_PARAMS[4][4] = {
    {0.0, 0.0, 0.0, 0.0}, // CH0
    {0.0, 0.0, 0.0, 0.0}, // CH1
    {0.0, 0.0, 0.0, 0.0}, // CH2
    {0.0, 0.0, 0.0, 0.0}  // CH3
};

// 時間分解能 (Sigma_t) 係数 [ch][param]
const double SIGMA_T_PARAMS[4][4] = {
    {0.0, 1.0, 0.0, 0.0}, // CH0: 1.0 ns constant
    {0.0, 1.0, 0.0, 0.0}, // CH1
    {0.0, 1.0, 0.0, 0.0}, // CH2
    {0.0, 1.0, 0.0, 0.0}  // CH3
};

// =========================================================
// [New] 新電荷モデル(SolidAngle)用パラメータ
// mu = A * f_i(r) * epsilon_i(cos_alpha)
// =========================================================

// 1. 距離依存性 f_i(r) = c0 * (1 - sqrt(1 - (32.5 / (r - c1))^2))
// 配列: {c0, c1}
// c0: スケール因子 (通常は 0.5 や 1.0 など立体角の定義による)
// c1: 距離オフセット (cm)
const double CHARGE_RADIAL_PARAMS[4][2] = {
    {1.0, 0.0}, // CH0
    {1.0, 0.0}, // CH1
    {1.0, 0.0}, // CH2
    {1.0, 0.0}  // CH3
};

// 2. 角度依存性 epsilon_i(cos_a) = c0 + c1*x + c2*x^2 + ... + c7*x^7
// 配列: {c0, c1, c2, c3, c4, c5, c6, c7}
// x = cos(alpha)
const double CHARGE_ANGULAR_PARAMS[4][8] = {
    {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, // CH0: 定数1.0 (等方)
    {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, // CH1
    {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, // CH2
    {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}  // CH3
};

// =========================================================
// PMT座標・形状定義 (4 PMTs)
// =========================================================
// X, Y座標は表面中心・半球中心で共通 (軸対称と仮定)
const double PMT_XY_POS[4][2] = {
    {-35.0,  35.0}, // PMT#1 (CH0)
    { 35.0,  35.0}, // PMT#2 (CH1)
    {-35.0, -35.0}, // PMT#3 (CH2)
    { 35.0, -35.0}  // PMT#4 (CH3)
};

// PMTの高さ方向の設定 (cm)
const double PMT_SURFACE_Z = 80.5; // 表面の中心Z座標
const double PMT_SPHERE_Z  = 48.0; // 半球の中心Z座標
const double PMT_RADIUS    = 32.5; // 半球の半径

// PMTの向き (全PMT共通でZ軸方向を向いていると仮定)
const double PMT_DIR[3] = {0.0, 0.0, 1.0}; 

// 便宜上の旧互換定義 (readDataなどで使用される表面中心座標)
const double PMT_POSITIONS[4][3] = {
    {-35.0,  35.0, PMT_SURFACE_Z},
    { 35.0,  35.0, PMT_SURFACE_Z},
    {-35.0, -35.0, PMT_SURFACE_Z},
    { 35.0, -35.0, PMT_SURFACE_Z}
};

// =========================================================
// 関数プロトタイプ
// =========================================================
// チャンネルと電荷を受け取り、パラメータ配列に基づいて値を計算する関数
double CalcParametricValue(int ch, double charge, const double params[4][4]);

// 以下、旧来の関数（必要に応じて中身を書き換え予定だが今回は上記関数で代替）
double GetEMG_Sigma(int ch, double charge);
double GetEMG_Tau(int ch, double charge);

// =========================================================
// データ構造体
// =========================================================
struct PMTData {
    int eventID;
    int ch;
    double time;
    double charge;
    double x, y, z; // 読み込み時の座標 (通常は表面中心が入る)
    double dir_x, dir_y, dir_z;
    bool isHit;
};

struct PedestalData {
    double hgain_mean;
    double lgain_mean;
};

// =========================================================
// フィッティング設定
// =========================================================

/**
 * @enum ChargeChi2Type
 * @brief 電荷分布のChi2計算モデル
 */
enum class ChargeChi2Type {
    Gaussian,       // ((obs - exp)^2 / sigma^2)
    BakerCousins,   // 2 * (mu - n + n*ln(n/mu))
    None            // 電荷情報を使用しない
};

/**
 * @enum ChargeModelType
 * @brief 電荷の期待値モデル
 */
enum class ChargeModelType {
    Standard,       // mu = A/r^2 + B (等方発光仮定)
    ZeroIntercept,  // mu = A/r^2 (B=0)
    Cosine,         // mu = (A * f(cos))/r^2 + B (角度依存あり)
    SolidAngle      // [New] mu = A * f(r) * epsilon(cos) + B (立体角 + 多項式角度補正)
};

/**
 * @enum TimeChi2Type
 * @brief 時間分布のChi2計算モデル
 */
enum class TimeChi2Type {
    Gaussian,       // ((obs - exp)^2 / sigma^2)
    EMG,            // EMG分布 (-ln L)
    Goodness,       // SK風Goodness
    None            // 時間情報を使用しない
};

struct FitConfig {
    ChargeChi2Type chargeType = ChargeChi2Type::Gaussian;
    ChargeModelType chargeModel = ChargeModelType::SolidAngle; // デフォルトは新モデル
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
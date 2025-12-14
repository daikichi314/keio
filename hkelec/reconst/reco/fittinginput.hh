/**
 * @file fittinginput.hh
 * @brief フィッティングに必要な定数、データ構造、設定構造体を定義するヘッダーファイル
 *
 * このファイルには以下の要素が含まれます：
 * 1. 物理定数 (光速など)
 * 2. 検出器ジオメトリ定義 (PMT座標)
 * 3. 入力データの構造体 (PMTData)
 * 4. 解析設定用の列挙型と構造体 (FitConfig)
 * - 電荷のChi2モデル (Gaussian / Baker-Cousins)
 * - 時間のChi2モデル (Gaussian / EMG)
 * - 3本ヒット時のUnhit補完設定
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
const double C_LIGHT = 29.970255; // 光速 [cm/ns] (空気中の屈折率 n=1.0003 を考慮)

// PMTごとの時間補正値 (ns)
const double TIME_CORRECTION_VAL[4] = {200.0, 200.0, 200.0, 200.0};
const double TIME_CORRECTION_ERR[4] = {0.0, 0.0, 0.0, 0.0};

// =========================================================
// PMT座標定義 (4 PMTs)
// =========================================================
const double PMT_POSITIONS[4][3] = {
    {-35.0,  35.0, 80.5}, // PMT#1 (CH0)
    { 35.0,  35.0, 80.5}, // PMT#2 (CH1)
    {-35.0, -35.0, 80.5}, // PMT#3 (CH2)
    { 35.0, -35.0, 80.5}  // PMT#4 (CH3)
};

const double PMT_DIR[3] = {0.0, 0.0, 1.0}; // PMTの向き (Z軸方向)

// =========================================================
// 関数プロトタイプ (EMGモデル用パラメータ)
// =========================================================
// 電荷に依存したTime Resolution (Sigma) と Time Constant (Tau) を返す関数
// 現在は定数を返すが、将来的にはQの関数となる予定
double GetEMG_Sigma(int ch, double charge);
double GetEMG_Tau(int ch, double charge);

// =========================================================
// データ構造体
// =========================================================

/**
 * @struct PMTData
 * @brief 1つのPMT(チャンネル)ごとのヒット情報を保持する構造体
 */
struct PMTData {
    int eventID;
    int ch;            // チャンネル番号 (0-3)
    double time;       // 補正済み観測時刻 (ns)
    double charge;     // 変換済み電荷 (pC または ADC相当値)
    double x, y, z;    // PMTの設置座標
    double dir_x, dir_y, dir_z; // PMTの向き
    bool isHit;        // 実際にヒットがあったか (Unhit補完データの場合はfalse)
};

/**
 * @struct PedestalData
 * @brief ペデスタル（ベースライン）情報を保持する構造体
 */
struct PedestalData {
    double hgain_mean;
    double lgain_mean;
};

// =========================================================
// フィッティング設定 (オプション機能)
// =========================================================

/**
 * @enum ChargeChi2Type
 * @brief 電荷分布のChi2計算モデルを選択
 */
enum class ChargeChi2Type {
    Gaussian,       // 従来: ((obs - exp)^2 / sigma^2)
    BakerCousins    // 新機能: Poisson尤度由来のChi2 (2 * (mu - n + n*ln(n/mu)))
};

/**
 * @enum TimeChi2Type
 * @brief 時間分布のChi2計算モデルを選択
 */
enum class TimeChi2Type {
    Gaussian,       // 従来: ((obs - exp)^2 / sigma^2)
    EMG             // 新機能: Exponentially Modified Gaussian (-ln(f_EMG))
};

/**
 * @struct FitConfig
 * @brief 解析実行時のオプション設定をまとめた構造体
 */
struct FitConfig {
    ChargeChi2Type chargeType = ChargeChi2Type::Gaussian; // デフォルト: ガウス
    TimeChi2Type timeType = TimeChi2Type::Gaussian;       // デフォルト: ガウス
    bool useUnhit = false; // 3本ヒット時に不足分をUnhit(0電荷)として補完するか (デフォルト: しない)
};

/**
 * @struct FitResult
 * @brief フィッティング結果を格納する構造体
 */
struct FitResult {
    double x, y, z, t;          // 推定された光源位置と発光時刻
    double err_x, err_y, err_z, err_t; // 上記の誤差
    double A;                   // 光量パラメータ A (1/r^2の係数)
    double B;                   // バックグラウンド光量パラメータ B
    double chi2;                // 最小化されたChi2値
    int ndf;                    // 自由度
    int status;                 // Minuitの収束ステータス (3 = Converged)
};

#endif // FITTINGINPUT_HH
/**
 * @file onemPMTfit.hh
 * @brief TMinuitを用いた光源位置再構成クラスの定義
 *
 * 与えられたPMTヒット情報から、光源の座標(x,y,z)、発光時刻(t0)、
 * 光量パラメータ(A, B)を推定します。
 * * 設定(FitConfig)に応じて、様々な物理モデルを切り替えて計算を行います。
 * ユーザー指定のモデルに従い、パラメータBの固定/解放なども制御します。
 *
 * @author Gemini (Modified based on user request)
 * @date 2025-12-14
 */

#ifndef ONEMPMTFIT_HH
#define ONEMPMTFIT_HH

#include "fittinginput.hh"
#include <vector>
#include <string>
#include <TMinuit.h>

// MinuitのFCN関数からアクセスするためのグローバルポインタ
class LightSourceFitter;

/**
 * @brief ファイル名から光源位置と減衰量を抽出する構造体
 */
struct FilenameParams {
    double x;
    double y;
    double z;
    double db;
    bool valid;
};

/**
 * @brief ファイル名から実験パラメータを抽出する関数
 * 
 * 期待する形式: BASENAME_xX_yY_zZ-RUNNUMBER-XXdB
 * 例: LDhkelec_x-35_y-35_z147-003-15.00dB
 * 
 * @param filename ファイルパス（または文字列）
 * @return FilenameParams 抽出されたパラメータ（失敗時はvalid=false）
 */
FilenameParams ParseFilename(const std::string& filename);

class LightSourceFitter {
public:
    LightSourceFitter();
    virtual ~LightSourceFitter();

    /**
     * @brief フィッティングの設定を適用する
     * @param config 設定構造体
     */
    void SetConfig(const FitConfig& config);

    /**
     * @brief データファイル名を設定する（初期値計算に使用）
     * @param filename データファイル名またはパス
     */
    void SetDataFilename(const std::string& filename);

    /**
     * @brief 1イベント分のフィッティングを実行する
     *
     * @param eventHits イベントのPMTデータリスト
     * @param res 結果を格納する構造体への参照
     * @return true 収束に成功 (Status=3)
     * @return false 収束に失敗
     */
    bool FitEvent(const std::vector<PMTData>& eventHits, FitResult& res);

    // 静的FCN関数からデータへアクセスするためのゲッター
    const std::vector<PMTData>& GetData() const { return fCurrentHits; }
    const FitConfig& GetConfig() const { return fConfig; }

private:
    TMinuit* fMinuit;
    std::vector<PMTData> fCurrentHits;
    FitConfig fConfig;
    std::string fDataFilename; // ファイル名から初期値を取得するために使用 

    /**
     * @brief パラメータの初期値と範囲を設定する
     * 重心計算などを用いて初期値を決定します。
     * 設定に応じてパラメータBを固定(Fix)するか決定します。
     */
    void InitializeParameters(const std::vector<PMTData>& hits);
};

// Minuitが最小化のために呼び出す関数 (グローバルスコープ)
void fcn_wrapper(int& npar, double* gin, double& f, double* par, int iflag);

#endif // ONEMPMTFIT_HH
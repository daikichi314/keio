#include "light_source_fit.hh"
#include "TMinuit.h"
#include "TVector3.h"
#include <cmath>
#include <iostream>

/**
 * This program fits the position and time of a light source based on the detected directions and times from multiple mPMTs.
 * The direction from each mPMT needs to be estimated, and the time of detection is also considered.
 * * Example (test program): lightfittest.cc
 *
 * In this program, chi2 is calculated based on the following formula:
 * χ² = Σ ((|perp|² / |d|² / σ²_θ + (|d|/c - Δt)² / σ²_t)
 * where:
 * - perp is the component of the vector from the mPMT to the light source that is perpendicular to the estimated direction
 * - |d| is the distance from the light source to the mPMT
 * - σ²_θ is the direction error (in rad)
 * - σ²_t is the time resolution (in ns)
 * - Δt is the difference between the expected time of light travel and the measured time
 * - c is the speed of light in water (22.556 cm/ns)
 */

/**
 * このプログラムは、複数のmPMTから検出された方向と時間に基づいて、光源の位置と時間をフィッティングします。
 * 各mPMTからの方向は推定する必要があり、検出時間も考慮されます。
 * * 例（テストプログラム）：lightfittest.cc
 * このプログラムでは、以下の式に基づいてχ²が計算されます：
 * χ² = Σ ((|perp|² / |d|² / σ²_θ + (|d|/c - Δt)² / σ²_t)
 * ここで：
 * - perpは、mPMTから光源へのベクトルのうち、推定された方向に垂直な成分です
 * - |d|は、光源からmPMTまでの距離です
 * - σ²_θは方向の誤差（ラジアン単位）です
 * - σ²_tは時間分解能（ナノ秒単位）です
 * - Δtは、光の移動の期待時間と測定時間の差です
 * - cは水中の光速（22.556 cm/ns）です  
 */

/**
 * FitLightSource : * 光源の位置と時間をフィッティングする関数。
 * sensors : センサーデータを含むSensorUnitのベクター。
 * pmtunits : 時間フィッティングのためのPMTユニットデータを含むSensorUnitのベクター。
 * fitted_pos : フィッティングされた位置を保持する出力TVector3。
 * fitted_time : フィッティングされた時間を保持する出力double。
 * fiterr_pos : フィッティングされた位置の誤差を保持する出力TVector3。
 * fiterr_time : フィッティングされた時間の誤差を保持する出力double。
 * chi2 : フィットのchi2値を保持する出力double。
 */

/**
 * SensorUnit : *  Structure to hold sensor data. Written in fittinginput.hh
 * Recommended to create a vector of SensorUnit with at least 4 mPMTs.
 * 日本語で
 * SensorUnit : * センサーデータを保持する構造体。fittinginput.hhに記載。
 * 少なくとも4つのmPMTを持つSensorUnitのベクター
 */

static std::vector<SensorUnit> gSensors;   // grobal variable to hold sensor data センサーデータを保持するグローバル変数
static std::vector<SensorUnit> gPMTsUnits; // global variable to hold PMT units for fitting 時間フィッティングのためのPMTユニットを保持するグローバル変数
double chi2_1, chi2_2;                     // global variables to hold chi2 values for debugging デバッグ用にchi2値を保持するグローバル変数

static void FcnForMinuit_lightsource(Int_t &, Double_t *, Double_t &fval, Double_t *par, Int_t) {
    // par[0] = x, par[1] = y, par[2] = z, par[3] = t0
    double t0 = par[3];
    const double c = 22.556; // cm/ns light speed in water

    fval = 0.0;
    chi2_1 = 0.0;
    chi2_2 = 0.0;
    for (const SensorUnit &s : gSensors) {
        double dx_x = par[0] - s.posx;
        double dx_y = par[1] - s.posy;
        double dx_z = par[2] - s.posz;
        double D = std::sqrt(dx_x * dx_x + dx_y * dx_y + dx_z * dx_z);
        if (D < 1e-6)
            continue;

        // direction χ² :

        // TVector3 perp = dx - dx.Dot(s.dir) * s.dir;
        double dxdir = dx_x * s.dirx + dx_y * s.diry + dx_z * s.dirz; // dot product of direction and vector to sensor
        double perp2 = 1 - dxdir * dxdir / D / D;

        double sigma_perp = s.sigma_sintheta;
        double dir_chi2 = perp2 * perp2 / (sigma_perp * sigma_perp * sigma_perp * sigma_perp); // chi2 for direction
		//double dir_chi2 = perp2/(sigma_perp*sigma_perp);

        // dxとdirの内積が負の場合は、方向が逆向きなので、chi2を大きくする
        // if (dxdir < 30) {
        //     dir_chi2 = 1000000; // 大きな値を設定して、方向が逆向きのセンサーを無視する
        // }

        // time χ² :
        double delta_time = s.time - t0;
        double time_chi = (D / c - delta_time) / s.sigma_time;
        double time_chi2 = time_chi * time_chi;
        // fval += 2 * dir_chi2 + time_chi2; // total chi2
        chi2_1 += dir_chi2; // total chi2
    }
    for (const SensorUnit &pmt : gPMTsUnits) {
        double dx_x = par[0] - pmt.posx;
        double dx_y = par[1] - pmt.posy;
        double dx_z = par[2] - pmt.posz;
        double D = std::sqrt(dx_x * dx_x + dx_y * dx_y + dx_z * dx_z);
        if (D < 1e-6)
            continue;

        // time χ² :
        double delta_time = pmt.time - t0;
        double time_chi = (D / c - delta_time) / pmt.sigma_time;
        double time_chi2 = time_chi * time_chi;
        chi2_2 += time_chi2; // accumulate time chi2 for debugging
    }
    fval = (chi2_1 + chi2_2);
    //fval = chi2_2;
}

void FitLightSource(const std::vector<SensorUnit> &sensors,
                    const std::vector<SensorUnit> &pmtunits,
                    TVector3 &fitted_pos,
                    double &fitted_time,
                    TVector3 &fiterr_pos,
                    double &fiterr_time,
                    double &chi2) {
    gSensors = sensors;
    gPMTsUnits = pmtunits;
    double last_hit_time = sensors[0].time;
    for (const auto &s : sensors) {
        if (s.time > last_hit_time) {
            last_hit_time = s.time;
        }
    }

    TMinuit minuit(4);
    minuit.SetFCN(FcnForMinuit_lightsource);
    minuit.SetPrintLevel(-1);

    minuit.DefineParameter(0, "x", 0.0, 1, -350, 350);
    minuit.DefineParameter(1, "y", 0.0, 1, -400, 400);
    minuit.DefineParameter(2, "z", 0.0, 1, -350, 350);
    minuit.DefineParameter(3, "t0", last_hit_time - 50, 0.1, last_hit_time - 100, last_hit_time);

    minuit.Migrad();

    double x, y, z, t0, errx, erry, errz, errt0;
    minuit.GetParameter(0, x, errx);
    minuit.GetParameter(1, y, erry);
    minuit.GetParameter(2, z, errz);
    minuit.GetParameter(3, t0, errt0);

    fitted_pos.SetXYZ(x, y, z);
    fitted_time = t0;
    fiterr_pos.SetXYZ(errx, erry, errz);
    fiterr_time = errt0;

    double amin, edm, errdef;
    int nvpar, nparx, istat;
    minuit.mnstat(amin, edm, errdef, nvpar, nparx, istat);
    chi2 = amin;

    std::cout << "dirchi2 = " << chi2_1 << " timechi2 = " << chi2_2 << " total chi2 = " << chi2 << std::endl;
}

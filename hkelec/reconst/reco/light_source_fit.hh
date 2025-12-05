// light_source_fit.hh
#ifndef LIGHT_SOURCE_FIT_HH
#define LIGHT_SOURCE_FIT_HH

#include "TVector3.h"
#include "fittinginput.hh"
#include <vector>

/**
 * @brief function to fit the light source position and time using Minuit.
 * @param sensors parameters for the light source fitting (SensorUnit vector)
 * @param fitted_pos fitted position of the light source (TVector3, in cm)
 * @param fitted_time fitted time of the light source (double, in ns)
 * @param fiterr_pos error of the fitted position (TVector3, in cm)
 * @param fiterr_time error of the fitted time (double, in ns)
 * @param chi2 chi-squared value of the fit (double)
 *
 * only sensors are input parameters.
 */
void FitLightSource(const std::vector<SensorUnit> &sensors,
                    TVector3 &fitted_pos,
                    double &fitted_time,
                    TVector3 &fiterr_pos,
                    double &fiterr_time,
                    double &chi2);
#endif

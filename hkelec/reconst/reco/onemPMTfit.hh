// 1mPMTfit.hh
#ifndef ONE_MPMT_FIT_HH
#define ONE_MPMT_FIT_HH
#include "fittinginput.hh"
#include "fstream"
int onemPMTfit(const std::string &inputfilename);
std::vector<PMTData> findExpandedGroups(const std::vector<PMTData> &pmtData,
                                        double tau,
                                        int n,
                                        double expand_before,
                                        double expand_after);
void FitPosition(std::vector<PMTData> &pmtData, // defined in fittinginput.hh
                 double &fit_theta,
                 double &fit_phi,
                 double &err_theta,
                 double &err_phi);
#endif
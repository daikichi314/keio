// readData.hh
#ifndef READ_DATA_HH
#define READ_DATA_HH
#include "fittinginput.hh"
#include <map>
#include <string>
/**
 * @brief Read PMT data from a ROOT file and store it in a map.
 *
 * @param inputfilename The name of the input file (without extension).
 * @param PMTDataMap A map to store PMT data, keyed by mPMTid.
 * @return int Returns 0 on success, 1 on failure.
 */
int readData(std::string &inputfilename, std::map<int, std::vector<PMTData>> &PMTDataMap);
/**
 * @brief Find PMT orientation of target PMT.
 */
int readPMTinfo(std::string &inputfilename, int &target_mPMT_id, int &target_mPMT_pmtid, double &ori_x_return, double &ori_y_return, double &ori_z_return);
#endif // READ_DATA_HH
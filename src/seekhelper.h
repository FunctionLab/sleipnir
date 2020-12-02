#ifndef SEEKHELPER_H
#define SEEKHELPER_H

// Parse config file for invocation parameters
// New db instances are pushed onto CSeekDBSetting cc
// Make optional parameters distanceMeasure and check_dset_size_flag and set 
// them to the most strict case in terms of checking meta-data is available.
bool legacyReadDBConfigFile(string dbConfigFile,
                            vector<CSeekDBSetting*> &cc, 
                            enum CSeekDataset::DistanceMeasure eDistMeasure = CSeekDataset::CORRELATION,
                            bool check_dset_size_flag = true);

// Read database config files
void readTomlConfig();

#endif  // SEEKHELPER_H
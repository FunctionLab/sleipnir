#ifndef SEEKHELPER_H
#define SEEKHELPER_H

#include "seekdataset.h"

void loadOneColumnTextFile(string filename, vector<string> &vals);

void loadTwoColumnTextFile(string filename, map<string, string> &vals);

#endif  // SEEKHELPER_H

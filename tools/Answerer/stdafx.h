#ifndef STDAFX_H
#define STDAFX_H

#include <fstream>
#include <set>
using namespace std;

#ifdef _MSC_VER
#include <windows.h>
#else // _MSC_VER
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif // _MSC_VER

#include "dat.h"
#include "genome.h"
#include "meta.h"
using namespace libBioUtils;

#endif // STDAFX_H

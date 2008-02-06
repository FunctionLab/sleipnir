#ifndef STDAFX_H
#define STDAFX_H

#ifndef _MSC_VER
#include <dirent.h>
#endif // _MSC_VER

#include <fstream>
using namespace std;

#include <pthread.h>

#include "bayesnet.h"
#include "dataset.h"
#include "genome.h"
#include "meta.h"
using namespace libBioUtils;

#ifndef ARRAYSIZE
#define ARRAYSIZE(x)	(sizeof(x)/sizeof(*(x)))
#endif // ARRAYSIZE

#endif // STDAFX_H

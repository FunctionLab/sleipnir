#ifndef STDAFX_H
#define STDAFX_H

#include <fstream>
using namespace std;

#include "dataset.h"
#include "genome.h"
#include "meta.h"
using namespace libBioUtils;

#ifndef ARRAYSIZE
#define ARRAYSIZE(x)	(sizeof(x)/sizeof(*(x)))
#endif // ARRAYSIZE

#endif // STDAFX_H

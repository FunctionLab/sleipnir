#ifndef STDAFX_H
#define STDAFX_H

#include <fstream>
using namespace std;

#include <pthread.h>

#include "bayesnet.h"
#include "dataset.h"
#include "genome.h"
#include "meta.h"
using namespace Sleipnir;

#ifndef ARRAYSIZE
#define ARRAYSIZE(x)	(sizeof(x)/sizeof(*(x)))
#endif // ARRAYSIZE

#ifndef _MSC_VER
#define _unlink unlink
#endif // _MSC_VER

#endif // STDAFX_H

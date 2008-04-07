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
using namespace Sleipnir;

#endif // STDAFX_H

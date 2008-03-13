#ifndef STDAFX_H
#define STDAFX_H

#include <fstream>
using namespace std;

#ifdef _MSC_VER
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif // _MSC_VER

#include "annotation.h"
#include "dat.h"
#include "genome.h"
#include "meta.h"
#include "statistics.h"
using namespace Sleipnir;

#endif // STDAFX_H

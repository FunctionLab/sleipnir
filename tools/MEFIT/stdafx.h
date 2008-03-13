#ifndef STDAFX_H
#define STDAFX_H

#pragma warning( disable : 4267 )

#include <fstream>
using namespace std;

#include "bayesnet.h"
#include "genome.h"
#include "meta.h"
using namespace Sleipnir;

#ifndef _MSC_VER
#include <dirent.h>
#include <sys/stat.h>

#define ARRAYSIZE(x)	(sizeof(x)/sizeof(*x))
#endif // _MSC_VER

#endif // STDAFX_H

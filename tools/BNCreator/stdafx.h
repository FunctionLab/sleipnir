#ifndef STDAFX_H
#define STDAFX_H

#include <fstream>
using namespace std;

#include "bayesnet.h"
#include "dataset.h"
#include "genome.h"
#include "meta.h"
using namespace Sleipnir;

#ifndef _MSC_VER
#include <dirent.h>

#define _unlink	unlink
#endif // _MSC_VER

#endif // STDAFX_H

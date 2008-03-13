#ifndef STDAFX_H
#define STDAFX_H

#include <fstream>
#include <iostream>
#include <map>
#include <string>
using namespace std;

#include "dat.h"
#include "genome.h"
#include "meta.h"
using namespace Sleipnir;

#ifndef ARRAYSIZE
#define	ARRAYSIZE(x)	(sizeof(x)/sizeof(*(x)))
#endif // ARRAYSIZE

#endif // STDAFX_H

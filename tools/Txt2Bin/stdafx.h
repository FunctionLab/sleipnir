#ifndef STDAFX_H
#define STDAFX_H

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

#ifndef _MSC_VER
#define ios_base	ios
#endif // _MSC_VER

#include "dataset.h"
#include "genome.h"
#include "meta.h"
using namespace libBioUtils;

const unsigned int	c_iLineLength	= 1024;

#endif // STDAFX_H

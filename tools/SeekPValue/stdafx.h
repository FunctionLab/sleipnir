/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#ifndef STDAFX_H
#define STDAFX_H

#define __STDC_LIMIT_MACROS

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <omp.h>

using namespace std;

#include <pthread.h>

#include "seekmap.h"
#include "seekcentral.h"
#include "seekweight.h"
#include "seekdataset.h"
#include "seekevaluate.h"
#include "seekplatform.h"
#include "seekreader.h"
#include "seekwriter.h"
#include "seekquery.h"
#include "meta.h"
#include "unistd.h"
#include "seeknetwork.h"

using namespace Sleipnir;

#endif // STDAFX_H

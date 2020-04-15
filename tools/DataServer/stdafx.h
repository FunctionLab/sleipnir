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

#include <fstream>
#include <algorithm>

#ifdef _MSC_VER
#include <io.h>
#include <winsock2.h>
#else // _MSC_VER

#include <arpa/inet.h>
#include <netinet/in.h>

#define SOCKET        int
#endif // _MSC_VER

#include <omp.h>

using namespace std;

#include "dataset.h"
#include "database.h"
#include "genome.h"
#include "meta.h"
#include "pcl.h"
#include "statistics.h"

using namespace Sleipnir;

#include "cmdline.h"
#include "server.h"
#include "serverclient.h"

typedef int (*TPFnCombiner)(const gengetopt_args_info &);

#endif // STDAFX_H

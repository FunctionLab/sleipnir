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

#include <float.h>
#include <math.h>

#include <algorithm>
#include <fstream>
using namespace std;

#include "coalescemotifs.h"
#include "coalescestructsi.h"
#include "dat.h"
#include "fasta.h"
#include "genome.h"
#include "hmm.h"
#include "measure.h"
#include "meta.h"
#include "pcl.h"
#include "statistics.h"
using namespace Sleipnir;

#endif // STDAFX_H

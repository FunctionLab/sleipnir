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
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, and Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
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

#ifndef _MSC_VER
#include <unistd.h>
#define _tempnam	tempnam
#define _unlink		unlink
#endif // _MSC_VER

#endif // STDAFX_H

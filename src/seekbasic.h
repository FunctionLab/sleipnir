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
#ifndef SEEKBASIC_H
#define SEEKBASIC_H
#include <vector>
#include <string>
#include <map>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <omp.h>
#include <cmath>
#include <assert.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_sort_float.h>
#include <gsl/gsl_sort_vector_float.h>
#include <gsl/gsl_permute_vector_float.h>

using namespace std;
//typedef unsigned short ushort;
typedef unsigned int utype;
//#define MAX_UTYPE 65535
const utype MAX_UTYPE = -1;
#endif



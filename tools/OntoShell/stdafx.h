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

#include <cstring>
#include <cmath>
#include <pthread.h>

#ifdef _MSC_VER
#include <winsock2.h>
#else // _MSC_VER

#include <sys/socket.h>

#define _strdup        strdup
#define SOCKET        int
// #define sprintf_s	sprintf
#endif // _MSC_VER

#include <algorithm>
#include <fstream>

using namespace std;

#pragma warning(disable:4275)

#ifdef _MSC_VER
#include <history.h>
#include <readline.h>
#else // _MSC_VER

#include <readline/history.h>
#include <readline/readline.h>

#endif

#include "annotation.h"
#include "genome.h"
#include "meta.h"
#include "server.h"
#include "serverclient.h"

using namespace Sleipnir;

#endif // STDAFX_H

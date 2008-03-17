#ifndef STDAFX_H
#define STDAFX_H

#include <math.h>
#include <pthread.h>
#ifdef _MSC_VER
#include <winsock2.h>
#else // _MSC_VER
#include <sys/socket.h>

#define _strdup		strdup
#define SOCKET		int
#define sprintf_s	sprintf
#endif // _MSC_VER

#include <algorithm>
#include <fstream>
using namespace std;

#pragma warning(disable:4275)

#include <history.h>
#include <readline.h>

#include "annotation.h"
#include "genome.h"
#include "meta.h"
#include "server.h"
#include "serverclient.h"
using namespace Sleipnir;

#endif // STDAFX_H

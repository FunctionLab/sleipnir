#ifndef STDAFX_H
#define STDAFX_H

#define __STDC_LIMIT_MACROS

#include <fstream>
#include <queue>
#include <strstream>
using namespace std;

#ifdef _MSC_VER
#include <io.h>
#include <winsock2.h>
#else // _MSC_VER
#include <arpa/inet.h>
#include <netinet/in.h>

#define SOCKET		int

inline bool _mktemp_s( char* szTemplate ) {

	return !mktemp( szTemplate ); }
#endif // _MSC_VER

#include <pthread.h>

#include <boost/graph/graphviz.hpp>

#include "annotation.h"
#include "bayesnet.h"
#include "database.h"
#include "genome.h"
#include "meta.h"
#include "server.h"
#include "serverclient.h"
#include "statistics.h"
using namespace libBioUtils;

#endif // STDAFX_H

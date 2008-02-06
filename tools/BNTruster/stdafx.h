#ifndef STDAFX_H
#define STDAFX_H

#include <fstream>
using namespace std;

#include <pthread.h>

#include "bayesnet.h"
#include "meta.h"
using namespace libBioUtils;

#include "cmdline.h"

typedef int (*TPFnTruster)( const gengetopt_args_info& );

#endif // STDAFX_H

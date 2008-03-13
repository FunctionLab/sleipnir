#ifndef STDAFX_H
#define STDAFX_H

#include <fstream>
using namespace std;

#include "dataset.h"
#include "meta.h"
#include "pcl.h"
using namespace Sleipnir;

#include "cmdline.h"

typedef int (*TPFnCombiner)( const gengetopt_args_info& );

#endif // STDAFX_H

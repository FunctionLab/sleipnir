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
#ifndef VWI_H
#define VWI_H

#ifndef NO_VW

#define fsync		_commit
#define go_params	gd_thread_params
#define isnan		_isnan
#define lseek		_lseek
#define open		_open
#define SHUT_WR		1

typedef size_t	ssize_t;

#undef int64_t
#pragma warning(disable : 4996 4244 4305)
#include <vw.h>
#pragma warning(default : 4996 4244 4305)

namespace Sleipnir {

class CVWImpl {
};

}

#endif // NO_VW

#endif // VWI_H

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
#ifndef DATAPAIRI_H
#define DATAPAIRI_H

#include "dat.h"

namespace Sleipnir {

class CPairImpl {
protected:
	static const char	c_szQuantExt[];

	static bool Open( const char*, std::ifstream& );
	static bool Open( const char*, std::vector<float>& );
};

class CDataPairImpl : protected CPairImpl, public CDat {
protected:
	void Reset( bool );

	bool				m_fContinuous;
	std::vector<float>	m_vecdQuant;
};

class CPCLPairImpl : protected CPairImpl, public CPCL {
protected:
	std::vector<std::vector<float> >	m_vecvecdQuants;
};

}

#endif // DATAPAIRI_H

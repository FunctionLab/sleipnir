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
#ifndef DATAPAIRI_H
#define DATAPAIRI_H

#include "dat.h"

namespace Sleipnir {

class CDataPair;
class CDatFilter;

class CPairImpl {
protected:
	static const char	c_szQuantExt[];

	static bool Open( const char*, std::ifstream& );
	static bool Open( const char*, std::vector<float>& );
};

class CDataPairImpl : protected CPairImpl, public CDat {
protected:
  	CDataPairImpl( ) : m_fQuantized(false) {}
	void Reset( bool );
	bool				m_fContinuous;
	bool				m_fQuantized;
	std::vector<float>	m_vecdQuant;

	static const char  c_acQdab[];
	bool OpenQdab( const char* szDatafile );
	void SetQuants( const float* adBinEdges, size_t iBins );
};

class CPCLPairImpl : protected CPairImpl, public CPCL {
protected:
	std::vector<std::vector<float> >	m_vecvecdQuants;
};

class CDatFilterImpl {
protected:
	CDatFilterImpl( ) : m_pDat(NULL), m_pFilter(NULL) { }

	bool Attach( const CDataPair*, const CDatFilter*, const CGenes*, CDat::EFilter, const CDat* );
	size_t GetGenes( ) const;
	std::string GetGene( size_t ) const;

	const CDat*			m_pAnswers;
	const CDataPair*	m_pDat;
	const CDatFilter*	m_pFilter;
	CDat::EFilter		m_eFilter;
	std::vector<bool>	m_vecfGenes;
	std::vector<size_t>	m_veciAnswers;
};

}

#endif // DATAPAIRI_H

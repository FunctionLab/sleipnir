#ifndef DATASETI_H
#define DATASETI_H

#include <set>

#include "examplei.h"
#include "fullmatrix.h"
#include "halfmatrix.h"

namespace libBioUtils {

class CDataPair;
class IDataset;

class CDataImpl {
protected:
	static const char	c_cSeparator	= '/';
	static const char	c_szDat[];
	static const char	c_szDab[];

	size_t OpenMax( const char*, const std::vector<std::string>&, bool,
		std::vector<std::string>&, std::set<std::string>* = NULL );
	bool OpenGenes( istream&, bool, std::set<std::string>& ) const;
	bool OpenGenes( const std::vector<std::string>& );

	bool IsHidden( size_t ) const;
	const std::string& GetGene( size_t ) const;
	size_t GetGenes( ) const;
	const std::vector<std::string>& GetGeneNames( ) const;
	size_t GetExperiments( ) const;
	size_t GetGene( const std::string& ) const;

	bool						m_fContinuous;
	std::vector<size_t>			m_veciMapping;
	std::vector<std::string>	m_vecstrGenes;
};

class CDatasetImpl : protected CDataImpl {
protected:
	bool Open( const CDataPair&, size_t, size_t );
	void TrimExamples( );

	size_t						m_iMax;
	CHalfMatrix<CExampleImpl>	m_Examples;
};

class CDataMaskImpl {
protected:
	const IDataset*	m_pDataset;
	CBinaryMatrix	m_Mask;
};

class CDataSubsetImpl : protected CDataImpl {
protected:
	bool Open( const CDataPair&, size_t );

	size_t						m_iSize;
	size_t						m_iOffset;
	std::vector<std::string>	m_vecstrData;
	CFullMatrix<CExampleImpl>	m_Examples;
};

class CDatasetCompactImpl : protected CDataImpl {
protected:
	CDatasetCompactImpl( );
	~CDatasetCompactImpl( );

	bool Open( const CDataPair&, CCompactMatrix& ) const;

	size_t				m_iSize;
	size_t				m_iData;
	std::vector<size_t>	m_veciQuants;
	CCompactMatrix*		m_aData;
};

}

#endif // DATASETI_H

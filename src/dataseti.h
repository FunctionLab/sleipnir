#ifndef DATASETI_H
#define DATASETI_H

#include <set>

#include "examplei.h"
#include "fullmatrix.h"
#include "halfmatrix.h"
#include "typesi.h"

namespace Sleipnir {

class CCompactMatrix;
class CDataPair;
class IBayesNet;
class IDataset;

class CDataImpl {
	friend class CDataFilter;
	friend class CDataMask;
protected:
	static const char	c_cSeparator	= '/';
	static const char	c_szDat[];
	static const char	c_szDab[];
	static const char	c_szPcl[];

	static void FilterGenes( IDataset*, const CGenes&, CDat::EFilter );

	size_t OpenMax( const char*, const std::vector<std::string>&, bool,
		std::vector<std::string>&, std::set<std::string>* = NULL );
	bool OpenGenes( std::istream&, bool, bool, std::set<std::string>& ) const;
	bool OpenGenes( const std::vector<std::string>& );

	bool IsHidden( size_t ) const;
	const std::string& GetGene( size_t ) const;
	size_t GetGenes( ) const;
	const std::vector<std::string>& GetGeneNames( ) const;
	size_t GetExperiments( ) const;
	size_t GetGene( const std::string& ) const;
	unsigned char GetBins( size_t ) const;
	bool OpenBinary( std::istream& );
	const unsigned char* OpenBinary( const unsigned char* );
	void SaveBinary( std::ostream& ) const;

	bool						m_fContinuous;
	std::vector<size_t>			m_veciMapping;
	std::vector<std::string>	m_vecstrGenes;
	std::vector<unsigned char>	m_veccQuants;
};

class CDatasetImpl : protected CDataImpl {
protected:
	CDatasetImpl( );
	~CDatasetImpl( );

	void Reset( );
	bool Open( const CDataPair&, size_t );
	bool Open( const CDataPair*, const char*, const IBayesNet* );

	void**	m_apData;
};

class CDataOverlayImpl {
protected:
	CDataOverlayImpl( ) : m_pDataset(NULL) { }

	const std::vector<std::string>& GetGeneNames( ) const;
	size_t GetExperiments( ) const;
	size_t GetGene( const std::string& ) const;
	size_t GetBins( size_t ) const;
	size_t GetGenes( ) const;
	bool IsHidden( size_t ) const;
	size_t GetDiscrete( size_t, size_t, size_t ) const;
	float GetContinuous( size_t, size_t, size_t ) const;
	const std::string& GetGene( size_t ) const;

	const IDataset*	m_pDataset;
};

class CDataMaskImpl : protected CDataOverlayImpl {
protected:
	CBinaryMatrix	m_Mask;
};

class CDataFilterImpl : protected CDataOverlayImpl {
protected:
	CDataFilterImpl( ) : m_pGenes(NULL), m_pAnswers(NULL) { }

	const CGenes*		m_pGenes;
	CDat::EFilter		m_eFilter;
	std::vector<bool>	m_vecfGenes;
	const CDat*			m_pAnswers;
	std::vector<size_t>	m_veciAnswers;
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

	bool Open( const CDataPair&, size_t );
	bool Open( const char*, const IBayesNet*, const CGenes* = NULL, const CGenes* = NULL );
	bool Open( const unsigned char* );
	virtual void Remove( size_t, size_t );
	size_t GetDiscrete( size_t, size_t, size_t ) const;
	void SaveText( std::ostream& ) const;
	void SaveBinary( std::ostream& ) const;
	virtual bool IsExample( size_t, size_t ) const;

	uint32_t		m_iData;
	CCompactMatrix*	m_aData;
};

}

#endif // DATASETI_H

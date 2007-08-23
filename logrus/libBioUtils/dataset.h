#ifndef DATASET_H
#define DATASET_H

#include "dataseti.h"

namespace libBioUtils {

class IBayesNet;

class IDataset {
public:
	virtual bool IsHidden( size_t ) const = 0;
	virtual size_t GetDiscrete( size_t, size_t, size_t ) const = 0;
	virtual float GetContinuous( size_t, size_t, size_t ) const = 0;
	virtual const std::string& GetGene( size_t ) const = 0;
	virtual size_t GetGenes( ) const = 0;
	virtual bool IsExample( size_t, size_t ) const = 0;
	virtual const std::vector<std::string>& GetGeneNames( ) const = 0;
	virtual size_t GetExperiments( ) const = 0;
	virtual size_t GetGene( const std::string& ) const = 0;
	virtual size_t GetBins( size_t ) const = 0;
	virtual void Remove( size_t, size_t ) = 0;
	virtual void FilterGenes( const CGenes&, CDat::EFilter ) = 0;
};

class CDataset : CDatasetImpl, public IDataset {
public:
	bool Open( const char*, const IBayesNet* );
	bool Open( const char*, const char*, const IBayesNet* );
	bool Open( const CDataPair&, const char*, const IBayesNet* );
	bool Open( const char*, const std::vector<std::string>& );
	bool Open( const std::vector<std::string>& );

	bool IsHidden( size_t ) const;
	size_t GetDiscrete( size_t, size_t, size_t ) const;
	float GetContinuous( size_t, size_t, size_t ) const;
	const std::string& GetGene( size_t ) const;
	size_t GetGenes( ) const;
	bool IsExample( size_t, size_t ) const;
	const std::vector<std::string>& GetGeneNames( ) const;
	size_t GetExperiments( ) const;
	size_t GetGene( const std::string& ) const;
	size_t GetBins( size_t ) const;
	void Remove( size_t, size_t );
	void FilterGenes( const CGenes&, CDat::EFilter );

	bool OpenGenes( const std::vector<std::string>& vecstrData ) {

		return CDataImpl::OpenGenes( vecstrData ); }
};

class CDatasetCompact : protected CDatasetCompactImpl, public IDataset {
public:
	bool Open( const CDataPair&, const char*, const IBayesNet*, bool = false );
	bool Open( const CDataPair&, const char*, const IBayesNet*, const CGenes&, const CGenes&, bool = false );
	bool Open( const char*, const IBayesNet* );
	bool Open( const char*, const IBayesNet*, const CGenes&, const CGenes& );
	bool Open( const std::vector<std::string>& );
	bool Open( std::istream& );
	bool Open( const CGenes&, const CGenes&, const CDataPair&, const std::vector<std::string>&, size_t,
		const IMeasure*, const std::vector<float>&, const IBayesNet* );
	bool Open( const CDataPair&, const std::vector<std::string>&, bool = false, bool = false,
		size_t = 2, bool = false );
	void Save( std::ostream&, bool ) const;
	bool FilterGenes( const char*, CDat::EFilter );
	void FilterAnswers( );

	bool IsHidden( size_t ) const;
	size_t GetDiscrete( size_t, size_t, size_t ) const;
	float GetContinuous( size_t, size_t, size_t ) const;
	const std::string& GetGene( size_t ) const;
	size_t GetGenes( ) const;
	virtual bool IsExample( size_t, size_t ) const;
	const std::vector<std::string>& GetGeneNames( ) const;
	size_t GetExperiments( ) const;
	size_t GetGene( const std::string& ) const;
	size_t GetBins( size_t ) const;
	virtual void Remove( size_t, size_t );
	void FilterGenes( const CGenes&, CDat::EFilter );

	bool OpenGenes( const std::vector<std::string>& vecstrData ) {

		return CDataImpl::OpenGenes( vecstrData ); }
};

class CDatasetCompactMap : public CDatasetCompact {
public:
	CDatasetCompactMap( );
	~CDatasetCompactMap( );

	bool Open( const char* );

	bool IsExample( size_t, size_t ) const;
	void Remove( size_t, size_t );

	unsigned char*	m_pbData;
private:
	CBinaryMatrix	m_Mask;
	size_t			m_iData;
	HANDLE			m_hndlMap;
};

class CDataMask : CDataMaskImpl, public IDataset {
public:
	void Attach( const IDataset* );
	void AttachRandom( const IDataset*, float );
	void AttachComplement( const CDataMask& );

	bool IsExample( size_t, size_t ) const;
	void Remove( size_t, size_t );

	const std::vector<std::string>& GetGeneNames( ) const {

		return CDataOverlayImpl::GetGeneNames( ); }

	size_t GetExperiments( ) const {

		return CDataOverlayImpl::GetExperiments( ); }

	size_t GetGene( const std::string& strGene ) const {

		return CDataOverlayImpl::GetGene( strGene ); }

	size_t GetBins( size_t iExp ) const {

		return CDataOverlayImpl::GetBins( iExp ); }

	size_t GetGenes( ) const {

		return CDataOverlayImpl::GetGenes( ); }

	bool IsHidden( size_t iExp ) const {

		return CDataOverlayImpl::IsHidden( iExp ); }

	size_t GetDiscrete( size_t iX, size_t iY, size_t iNode ) const {

		return CDataOverlayImpl::GetDiscrete( iX, iY, iNode ); }

	float GetContinuous( size_t iX, size_t iY, size_t iNode ) const {

		return CDataOverlayImpl::GetContinuous( iX, iY, iNode ); }

	const std::string& GetGene( size_t iGene ) const {

		return CDataOverlayImpl::GetGene( iGene ); }

	void FilterGenes( const CGenes& Genes, CDat::EFilter eFilt ) {

		CDataImpl::FilterGenes( this, Genes, eFilt ); }
};

class CDataFilter : CDataFilterImpl, public IDataset {
public:
	void Attach( const IDataset*, const CGenes&, CDat::EFilter, const CDat* = NULL, bool = false );

	bool IsExample( size_t, size_t ) const;
	void Remove( size_t, size_t );

	const std::vector<std::string>& GetGeneNames( ) const {

		return CDataOverlayImpl::GetGeneNames( ); }

	size_t GetExperiments( ) const {

		return CDataOverlayImpl::GetExperiments( ); }

	size_t GetGene( const std::string& strGene ) const {

		return CDataOverlayImpl::GetGene( strGene ); }

	size_t GetBins( size_t iExp ) const {

		return CDataOverlayImpl::GetBins( iExp ); }

	size_t GetGenes( ) const {

		return CDataOverlayImpl::GetGenes( ); }

	bool IsHidden( size_t iExp ) const {

		return CDataOverlayImpl::IsHidden( iExp ); }

	size_t GetDiscrete( size_t iX, size_t iY, size_t iNode ) const {

		return ( IsExample( iX, iY ) ? CDataOverlayImpl::GetDiscrete( iX, iY, iNode ) : -1 ); }

	float GetContinuous( size_t iX, size_t iY, size_t iNode ) const {

		return ( IsExample( iX, iY ) ? CDataOverlayImpl::GetContinuous( iX, iY, iNode ) :
			CMeta::GetNaN( ) ); }

	const std::string& GetGene( size_t iGene ) const {

		return CDataOverlayImpl::GetGene( iGene ); }

	void FilterGenes( const CGenes& Genes, CDat::EFilter eFilt ) {

		CDataImpl::FilterGenes( this, Genes, eFilt ); }
};

class CDataSubset : CDataSubsetImpl, public IDataset {
public:
	bool Initialize( const char*, const IBayesNet*, size_t );
	bool Initialize( const std::vector<std::string>&, size_t );
	bool Open( size_t );

	bool IsHidden( size_t ) const;
	size_t GetDiscrete( size_t, size_t, size_t ) const;
	float GetContinuous( size_t, size_t, size_t ) const;
	const std::string& GetGene( size_t ) const;
	size_t GetGenes( ) const;
	bool IsExample( size_t, size_t ) const;
	const std::vector<std::string>& GetGeneNames( ) const;
	size_t GetExperiments( ) const;
	size_t GetGene( const std::string& ) const;
	size_t GetBins( size_t ) const;
	void Remove( size_t, size_t );
	void FilterGenes( const CGenes&, CDat::EFilter );
};

}

#endif // DATASET_H

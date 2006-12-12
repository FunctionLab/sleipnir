#ifndef DATI_H
#define DATI_H

#include <map>

#include "halfmatrix.h"
#include "file.h"
#include "measure.h"
#include "meta.h"
#include "pcl.h"

namespace libBioUtils {

class CGenes;
class CSlim;

class CDatImpl : protected CFile {
protected:
	typedef std::map<std::string,size_t>	TMapStrI;
	typedef std::vector<std::string>		TVecStr;
	typedef std::vector<float>				TAF;
	typedef std::vector<TAF>				TAAF;

	static const size_t	c_iGeneLimit	= 100000;
	static const size_t	c_iApproximate	= 1000;
	static const float	c_dNaN;
	static const char	c_szBin[];
	static const char	c_szPcl[];

	static size_t MapGene( TMapStrI&, TVecStr&, const std::string& );
	static void ResizeNaN( TAF&, size_t );
	static std::string DabGene( std::istream& );

	CDatImpl( ) : m_pMeasure(NULL), m_pPCL(NULL) { }
	~CDatImpl( );

	void Reset( );
	bool OpenPCL( std::istream& );
	bool OpenText( std::istream&, float dDefault );
	bool OpenBinary( std::istream& );
	bool OpenGenes( std::istream&, bool, bool );
	void SaveText( std::ostream& ) const;
	void SaveBinary( std::ostream& ) const;
	size_t GetGene( const std::string& ) const;
	void SlimCache( const CSlim&, std::vector<std::vector<size_t> >& ) const;
	void NormalizeMinmax( );
	void NormalizeStdev( );
	void OpenHelper( const CGenes*, float );
	void OpenHelper( const CGenes*, const CGenes*, float );

	float Get( size_t iX, size_t iY ) const {

		return ( m_pMeasure ? (float)m_pMeasure->Measure( m_pPCL->Get( iX ), m_pPCL->GetExperiments( ),
			m_pPCL->Get( iY ), m_pPCL->GetExperiments( ) ) :
			( ( iX == iY ) ? CMeta::GetNaN( ) : m_Data.Get( iX, iY ) ) ); }

	float& Get( size_t iX, size_t iY ) {
		static float	c_dNull	= CMeta::GetNaN( );

		return ( m_pMeasure ? c_dNull : ( ( iX == iY ) ? c_dNull : m_Data.Get( iX, iY ) ) ); }

	bool Set( size_t iX, size_t iY, float dValue ) {

		if( iX == iY )
			return false;

		m_Data.Set( iX, iY, dValue );
		return true; }

	bool Set( size_t iX, const float* adValues ) {

		m_Data.Set( iX, adValues );
		return true; }

	size_t GetGenes( ) const {

		return ( m_pPCL ? m_pPCL->GetGenes( ) : m_vecstrGenes.size( ) ); }

	std::string GetGene( size_t iGene ) const {

		return ( m_pPCL ? m_pPCL->GetGene( iGene ) : m_vecstrGenes[ iGene ] ); }

	CDistanceMatrix	m_Data;
	TVecStr			m_vecstrGenes;
	CPCL*			m_pPCL;
	bool			m_fPCLMemory;
	const IMeasure*	m_pMeasure;
	bool			m_fMeasureMemory;
};

class CDatMapImpl : protected CDatImpl {
protected:
	CDatMapImpl( ) : m_aadData(NULL), m_hndlData(0) { }

	~CDatMapImpl( ) {

		CMeta::Unmap( m_abData, m_hndlData, m_iData );
		if( m_aadData )
			delete[] m_aadData; }

	unsigned char*	m_abData;
	size_t			m_iData;
	HANDLE			m_hndlData;
	float**			m_aadData;
};

}

#endif // DATI_H

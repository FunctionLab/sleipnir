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

	static const size_t	c_iGeneLimit		= 100000;
	static const size_t	c_iApproximate		= 1000;
	static const size_t	c_iNeighborhood1	= 40;
	static const size_t	c_iNeighborhood2	= 0;
	static const size_t	c_iDegree			= 1;
	static const float	c_dCutoff;
	static const float	c_adColorMin[];
	static const float	c_adColorMax[];

	static size_t MapGene( TMapStrI&, TVecStr&, const std::string& );
	static void ResizeNaN( TAF&, size_t );
	static void DabGene( std::istream&, char* );

	CDatImpl( ) : m_pMeasure(NULL), m_pPCL(NULL), m_abData(NULL), m_iData(0), m_aadData(NULL), m_hndlData(0) { }
	~CDatImpl( );

	void Reset( );
	bool OpenPCL( std::istream& );
	bool OpenText( std::istream&, float dDefault );
	bool OpenBinary( std::istream& );
	bool OpenSparse( std::istream& );
	bool OpenGenes( std::istream&, bool, bool );
	void SaveText( std::ostream& ) const;
	void SaveBinary( std::ostream& ) const;
	void SaveSparse( std::ostream& ) const;
	void SaveGenes( std::ostream& ) const;
	size_t GetGene( const std::string& ) const;
	void SlimCache( const CSlim&, std::vector<std::vector<size_t> >& ) const;
	void AveStd( double&, double&, size_t& ) const;
	void NormalizeMinmax( );
	void NormalizeStdev( );
	void OpenHelper( const CGenes*, float );
	void OpenHelper( const CGenes*, const CGenes*, float );
	bool OpenHelper( );
	bool OpenMemmap( const unsigned char* );
	void FilterGenesPixie( const CGenes&, vector<bool>& );

	float& Get( size_t iX, size_t iY ) const {
		static float	s_dRet;

		return ( m_pMeasure ? ( s_dRet = (float)m_pMeasure->Measure( m_pPCL->Get( iX ), m_pPCL->GetExperiments( ),
			m_pPCL->Get( iY ), m_pPCL->GetExperiments( ) ) ) :
			( ( iX == iY ) ? ( s_dRet = CMeta::GetNaN( ) ) : m_Data.Get( iX, iY ) ) ); }

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
// PCL back end
	CPCL*			m_pPCL;
	bool			m_fPCLMemory;
	const IMeasure*	m_pMeasure;
	bool			m_fMeasureMemory;
// Memory mapped back end
	unsigned char*	m_abData;
	size_t			m_iData;
	HANDLE			m_hndlData;
	float**			m_aadData;
};

}

#endif // DATI_H

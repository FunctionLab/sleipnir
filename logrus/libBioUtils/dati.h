#ifndef DATI_H
#define DATI_H

#include <map>

#include "halfmatrix.h"
#include "file.h"
#include "meta.h"

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
	static const float	c_dNaN;
	static const char	c_szBin[];

	static size_t MapGene( TMapStrI&, TVecStr&, const std::string& );
	bool OpenGenes( std::istream&, bool );
	void OpenHelper( const CGenes*, float );
	static void ResizeNaN( TAF&, size_t );
	static std::string DabGene( std::istream& );

	CDatImpl( );
	~CDatImpl( );

	void Reset( );
	bool OpenText( std::istream& );
	bool OpenBinary( std::istream& );
	void SaveText( std::ostream& ) const;
	void SaveBinary( std::ostream& ) const;
	size_t GetGene( const std::string& ) const;
	void SlimCache( const CSlim&, std::vector<std::vector<size_t> >& ) const;
	void NormalizeMinmax( );
	void NormalizeStdev( );

	float Get( size_t iX, size_t iY ) const {

		return ( ( iX == iY ) ? CMeta::GetNaN( ) : m_Data.Get( iX, iY ) ); }

	bool Set( size_t iX, size_t iY, float dValue ) {

		if( iX == iY )
			return false;

		m_Data.Set( iX, iY, dValue );
		return true; }

	bool Set( size_t iX, const float* adValues ) {

		m_Data.Set( iX, adValues );
		return true; }

	size_t GetGenes( ) const {

		return m_vecstrGenes.size( ); }

	CDistanceMatrix	m_Data;
	TVecStr			m_vecstrGenes;
};

}

#endif // DATI_H

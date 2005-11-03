#ifndef DATI_H
#define DATI_H

#include <map>

#include "halfmatrix.h"
#include "file.h"

namespace libBioUtils {

class CSlim;

class CDatImpl : protected CFile {
protected:
	typedef std::map<std::string,size_t>	TMapStrI;
	typedef std::vector<std::string>		TVecStr;
	typedef std::vector<float>				TAF;
	typedef std::vector<TAF>				TAAF;

	static const float	c_dNaN;
	static const char	c_szBin[];

	static size_t MapGene( TMapStrI&, TVecStr&, const std::string& );
	bool OpenGenes( std::istream&, bool );
	static void ResizeNaN( TAF&, size_t );
	static std::string DabGene( std::istream& );

	CDatImpl( );
	~CDatImpl( );

	void Reset( );
	bool OpenText( std::istream& );
	bool OpenBinary( std::istream& );
	void SaveText( std::ostream& ) const;
	void SaveBinary( std::ostream& ) const;
	float Get( size_t, size_t ) const;
	bool Set( size_t, size_t, float );
	size_t GetGene( const std::string& ) const;
	void SlimCache( const CSlim&, std::vector<std::vector<size_t> >& ) const;

	CDistanceMatrix	m_Data;
	TVecStr			m_vecstrGenes;
};

}

#endif // DATI_H

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
#ifndef DATABASEI_H
#define DATABASEI_H

#include <fstream>
#include <map>
#include <vector>

//Qian added
#include <stdio.h>

#include "compactmatrix.h"

namespace Sleipnir {

class CDatabaselet {
public:
	enum ENibbles {
		ENibblesLow,
		ENibblesHigh,
		ENibblesBoth
	};

	CDatabaselet( bool );
	~CDatabaselet( );

	bool Open( const std::string&, const std::vector<std::string>&, uint32_t, uint32_t );
	bool Open( const std::string& );
	bool Open( const std::vector<CCompactFullMatrix>&, size_t, size_t, bool );

	bool OpenNoOverwrite();

	bool OpenWrite( unsigned char, size_t, ENibbles, unsigned char* );

	/* Get pair by referring to memory cache (ie charImage) of the db file */
	bool Get( size_t iOne, size_t iTwo, vector<unsigned char>& vecbData, unsigned char *charImage);

	/* Get pair by seeking in db file */
	bool Get( size_t, size_t, std::vector<unsigned char>& ) const;
	bool Get( size_t, std::vector<unsigned char>&, bool ) const;
	bool Get( size_t, const std::vector<size_t>&, std::vector<unsigned char>&, bool ) const;

	static bool Combine(std::vector<CDatabaselet*>& vecDatabaselet,
			std::string strOutDirectory, vector<string> &vecstrGenes, bool bSplit = true);

	size_t GetGenes( ) const {

		return m_vecstrGenes.size( ); }

	const std::string& GetGene( size_t iGene ) const {
		static const std::string	c_strEmpty	= "";

		return ( m_vecstrGenes.empty( ) ? c_strEmpty : m_vecstrGenes[ iGene % m_vecstrGenes.size( ) ] ); }

	void Write( size_t iOne, size_t iTwo, size_t iDataset, unsigned char bValue, bool fBoth = false ) {
		std::streamoff	iOffset;

		iOffset = (std::streamoff)GetOffset( iOne, iTwo, iDataset );

		if(m_useNibble){
			if( !fBoth ) {
				unsigned char	b;
				m_fstm.seekg( iOffset );
				b = m_fstm.get( );
				bValue = ( iDataset % 2 ) ? ( ( b & 0xF ) | ( bValue << 4 ) ) :
						( ( b & 0xF0 ) | ( bValue & 0xF ) ); 
				}
		}

		m_fstm.seekp( iOffset );
		m_fstm.put( bValue );
	}

	size_t GetDatasets( ) const {

		return m_iDatasets; }

	void CloseFile(){
		if(m_fstm.is_open()){
			m_fstm.close();
		}
	}

	void SetFile(string std){
		strFileName = std;
	}

private:

	size_t GetOffsetDataset( size_t iDataset ) const {
		if(m_useNibble){
			return (iDataset / 2);
		}else{
			return iDataset;
		}
	}

	size_t GetSizePair( ) const {

		if(m_useNibble){
			return (m_iDatasets + 1) / 2;
		}else{
			return m_iDatasets;
		}

	}

	size_t GetSizeGenes( ) const {

		return ( GetSizeGene( ) * m_vecstrGenes.size( ) ); }

	size_t GetSizeGene( ) const {

		return ( GetSizePair( ) * m_iGenes ); }

	size_t GetOffset( size_t iGene ) const {

		return ( m_iHeader + ( GetSizeGene( ) * iGene ) ); }

	size_t GetOffset( size_t iOne, size_t iTwo ) const {

		return ( GetOffset( iOne ) + ( GetSizePair( ) * iTwo ) ); }

	size_t GetOffset( size_t iOne, size_t iTwo, size_t iDataset ) const {

		return ( GetOffset( iOne, iTwo ) + GetOffsetDataset( iDataset ) ); }

	uint32_t					m_iHeader;
	uint32_t					m_iGenes;
	uint32_t					m_iDatasets;
	std::vector<std::string>	m_vecstrGenes;
	std::string					strFileName;

	bool						m_useNibble;
	mutable std::fstream		m_fstm;
	mutable pthread_mutex_t*	m_pmutx;
};

class CDatabaseImpl {
protected:
	static const char	c_acDAB[];
	static const char	c_acQDAB[];
	static const char	c_acExtension[];

	CDatabaseImpl(bool useNibble){
		m_fMemmap = false;
		m_iBlockIn = -1;
		m_iBlockOut = -1;
		m_fBuffer = false;
		m_useNibble = useNibble;
	}

	~CDatabaseImpl( ) {

		Clear( ); }

	bool Open( const std::vector<std::string>&, const std::vector<std::string>& );
	bool Open( const std::string&, size_t, bool = false );

	void Clear( ) {
		size_t	i;
		m_mapstriGenes.clear( );
		for( i = 0; i < m_vecpDBs.size( ); ++i )
			delete m_vecpDBs[ i ];
		m_vecpDBs.clear( ); }

	size_t GetGene( const std::string& strGene ) const {
		std::map<std::string, size_t>::const_iterator	iterGene;

		return ( ( ( iterGene = m_mapstriGenes.find( strGene ) ) == m_mapstriGenes.end( ) ) ? -1 :
			iterGene->second ); }

	bool							m_fMemmap;
	bool							m_fBuffer;
	size_t							m_iBlockIn;
	size_t							m_iBlockOut;
	std::vector<CDatabaselet*>		m_vecpDBs;
	std::map<std::string, size_t>	m_mapstriGenes;
	/* defines whether the CDatabaselet is nibble type. If false, it is byte by default.*/
	bool							m_useNibble;
};

}

#endif // DATABASEI_H

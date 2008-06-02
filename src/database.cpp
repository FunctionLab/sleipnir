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
#include "stdafx.h"
#include "database.h"
#include "meta.h"
#include "bayesnetint.h"
#include "datapair.h"

namespace Sleipnir {

const char	CDatabaseImpl::c_acDAB[]		= ".dab";
const char	CDatabaseImpl::c_acExtension[]	= ".db";

///////////////////////////////////////////////////////////////////////////////
// CDatabaselet
///////////////////////////////////////////////////////////////////////////////

CDatabaselet::CDatabaselet( ) {

	m_pmutx = new pthread_mutex_t( );
	pthread_mutex_init( m_pmutx, NULL ); }

CDatabaselet::~CDatabaselet( ) {

	pthread_mutex_destroy( m_pmutx );
	delete m_pmutx; }

bool CDatabaselet::Open( const string& strFile, const vector<string>& vecstrGenes, uint32_t iGenes,
	uint32_t iDatasets ) {
	uint32_t	iSize;
	size_t		i;
	char*		acFiller;

	m_fstm.clear( );
	m_fstm.open( strFile.c_str( ), ios_base::in | ios_base::out | ios_base::binary | ios_base::trunc );
	if( !m_fstm.is_open( ) ) {
		g_CatSleipnir.error( "CDatabaselet::Open( %s, %u, %u ) open failed", strFile.c_str( ), iGenes,
			iDatasets );
		return false; }

	m_iGenes = iGenes;
	m_iDatasets = iDatasets;
	m_vecstrGenes.resize( vecstrGenes.size( ) );
	copy( vecstrGenes.begin( ), vecstrGenes.end( ), m_vecstrGenes.begin( ) );

	m_fstm.write( (char*)&m_iHeader, sizeof(m_iHeader) );
	m_fstm.write( (char*)&m_iGenes, sizeof(m_iGenes) );
	m_fstm.write( (char*)&m_iDatasets, sizeof(m_iDatasets) );
	iSize = m_vecstrGenes.size( );
	m_fstm.write( (char*)&iSize, sizeof(iSize) );
	m_iHeader = sizeof(m_iHeader) + sizeof(m_iGenes) + sizeof(m_iDatasets) + sizeof(iSize);
	for( i = 0; i < m_vecstrGenes.size( ); ++i ) {
		m_fstm.write( m_vecstrGenes[ i ].c_str( ), m_vecstrGenes[ i ].size( ) + 1 );
		m_iHeader += m_vecstrGenes[ i ].size( ) + 1; }
	m_fstm.seekp( 0 );
	m_fstm.write( (char*)&m_iHeader, sizeof(m_iHeader) );

	m_fstm.seekp( m_iHeader );
	acFiller = new char[ GetSizeGene( ) ];
	memset( acFiller, -1, GetSizeGene( ) );
	for( i = 0; i < m_vecstrGenes.size( ); ++i )
		m_fstm.write( acFiller, GetSizeGene( ) );
	delete[] acFiller;

	return true; }

bool CDatabaselet::OpenWrite( unsigned char bValue, size_t iOffset, ENibbles eNibbles,
	unsigned char* abImage ) {
	unsigned char	b;

#ifndef DATABASE_NIBBLES
	eNibbles = ENibblesBoth;
#endif // DATABASE_NIBBLES

	if( abImage )
		iOffset -= m_iHeader;
	if( eNibbles != ENibblesBoth ) {
		if( abImage )
			b = abImage[ iOffset ];
		else {
			m_fstm.seekg( iOffset );
			b = m_fstm.get( ); } }
	switch( eNibbles ) {
		case ENibblesLow:
			b = ( bValue & 0xF ) | ( b & 0xF0 );
			break;

		case ENibblesHigh:
			b = ( b & 0xF ) | ( bValue << 4 );
			break;

		case ENibblesBoth:
			b = bValue; }
	if( abImage )
		abImage[ iOffset ] = b;
	else {
		m_fstm.seekp( iOffset );
		m_fstm.put( b ); }

	return true; }

bool CDatabaselet::Open( const vector<CCompactFullMatrix>& vecData, size_t iBaseGenes, size_t iBaseDatasets,
	bool fBuffer ) {
	unsigned char*	abImage;
	size_t			iSize, iDatum, iGeneOne, iGeneTwo;
	unsigned char	bOne, bTwo;

	if( fBuffer ) {
		abImage = new unsigned char[ iSize = ( GetSizeGene( ) * m_vecstrGenes.size( ) ) ];
		m_fstm.seekg( m_iHeader );
		m_fstm.read( (char*)abImage, iSize ); }
	else
		abImage = NULL;
	if( iBaseDatasets % 2 )
		for( iGeneOne = 0; iGeneOne < GetGenes( ); ++iGeneOne )
			for( iGeneTwo = 0; iGeneTwo < vecData[ 0 ].GetColumns( ); ++iGeneTwo )
				if( bOne = vecData[ 0 ].Get( iBaseGenes + iGeneOne, iGeneTwo ) )
					OpenWrite( bOne - 1, GetOffset( iGeneOne, iGeneTwo, iBaseDatasets ), ENibblesHigh,
						abImage );
	for( iDatum = ( iBaseDatasets % 2 ); ( iDatum + 1 ) < vecData.size( ); iDatum += 2 )
		for( iGeneOne = 0; iGeneOne < GetGenes( ); ++iGeneOne )
			for( iGeneTwo = 0; iGeneTwo < vecData[ iDatum ].GetColumns( ); ++iGeneTwo ) {
				bOne = vecData[ iDatum ].Get( iBaseGenes + iGeneOne, iGeneTwo );
				bTwo = vecData[ iDatum + 1 ].Get( iBaseGenes + iGeneOne, iGeneTwo );
				if( !( bOne || bTwo ) )
					continue;
				bOne -= 1;
				bTwo -= 1;
#ifdef DATABASE_NIBBLES
				OpenWrite( ( bOne & 0xF ) | ( bTwo << 4 ), GetOffset( iGeneOne, iGeneTwo, iBaseDatasets +
					iDatum ), ENibblesBoth, abImage );
#else // DATABASE_NIBBLES
				OpenWrite( bOne, GetOffset( iGeneOne, iGeneTwo, iBaseDatasets + iDatum ), ENibblesBoth,
					abImage );
				OpenWrite( bTwo, GetOffset( iGeneOne, iGeneTwo, iBaseDatasets + iDatum + 1 ), ENibblesBoth,
					abImage );
#endif // DATABASE_NIBBLES
			}
	if( iDatum < vecData.size( ) )
		for( iGeneOne = 0; iGeneOne < GetGenes( ); ++iGeneOne )
			for( iGeneTwo = 0; iGeneTwo < vecData[ iDatum ].GetColumns( ); ++iGeneTwo )
				if( bOne = vecData[ iDatum ].Get( iBaseGenes + iGeneOne, iGeneTwo ) )
					OpenWrite( bOne - 1, GetOffset( iGeneOne, iGeneTwo, iBaseDatasets + iDatum ), ENibblesLow,
						abImage );
	if( fBuffer ) {
		m_fstm.seekp( m_iHeader );
		m_fstm.write( (char*)abImage, iSize );
		delete[] abImage; }

	return true; }

bool CDatabaselet::Get( size_t iOne, size_t iTwo, vector<unsigned char>& vecbData ) const {
	size_t	i;

	i = vecbData.size( );
	vecbData.resize( i + GetSizePair( ) );
	pthread_mutex_lock( m_pmutx );
	m_fstm.seekg( GetOffset( iOne, iTwo ) );
	m_fstm.read( (char*)&vecbData[ i ], GetSizePair( ) );
	pthread_mutex_unlock( m_pmutx );

	return true; }

bool CDatabaselet::Get( size_t iGene, vector<unsigned char>& vecbData, bool fReplace ) const {
	size_t	i;

	i = fReplace ? 0 : vecbData.size( );
	vecbData.resize( i + GetSizeGene( ) );
	pthread_mutex_lock( m_pmutx );
	m_fstm.seekg( GetOffset( iGene ) );
	m_fstm.read( (char*)&vecbData[ i ], GetSizeGene( ) );
	pthread_mutex_unlock( m_pmutx );

	return true; }

bool CDatabaselet::Get( size_t iGene, const vector<size_t>& veciGenes, vector<unsigned char>& vecbData ) const {
	size_t	i;

	for( i = 0; i < veciGenes.size( ); ++i )
		if( !Get( iGene, veciGenes[ i ], vecbData ) )
			return false;

	return true; }

bool CDatabaselet::Open( const string& strFile ) {
	uint32_t	iSize;
	char*		acBuffer;
	char*		pc;
	size_t		i;

	m_fstm.clear( );
	m_fstm.open( strFile.c_str( ), ios_base::binary | ios_base::in );
	if( !m_fstm.is_open( ) ) {
		g_CatSleipnir.error( "CDatabaselet::Open( %s ) open failed", strFile.c_str( ) );
		return false; }

	m_fstm.read( (char*)&m_iHeader, sizeof(m_iHeader) );
	m_fstm.read( (char*)&m_iGenes, sizeof(m_iGenes) );
	m_fstm.read( (char*)&m_iDatasets, sizeof(m_iDatasets) );
	m_fstm.read( (char*)&iSize, sizeof(iSize) );

	acBuffer = new char[ m_iHeader ];
	m_fstm.read( acBuffer, m_iHeader - sizeof(m_iHeader) - sizeof(m_iGenes) - sizeof(m_iDatasets) -
		sizeof(iSize) );
	m_vecstrGenes.resize( iSize );
	for( i = 0,pc = acBuffer; i < m_vecstrGenes.size( ); pc += m_vecstrGenes[ i++ ].length( ) + 1 )
		m_vecstrGenes[ i ] = pc;
	delete[] acBuffer;

	return true; }

///////////////////////////////////////////////////////////////////////////////
// CDatabase
///////////////////////////////////////////////////////////////////////////////

bool CDatabase::Open( const vector<string>& vecstrGenes, const string& strInputDirectory,
	const IBayesNet* pBayesNet, const string& strOutputDirectory, size_t iFiles ) {
	vector<string>	vecstrNodes, vecstrSubset;
	size_t			i, j;
	char			acNumber[ 16 ];
	string			strFile;

	if( !pBayesNet ) {
		g_CatSleipnir.error( "CDatabase::Open( %s, %d ) null Bayes net", strOutputDirectory.c_str( ), iFiles );
		return false; }

	Clear( );
	pBayesNet->GetNodes( vecstrNodes );
	for( i = 1; i < vecstrNodes.size( ); ++i )
		vecstrNodes[ i - 1 ] = strInputDirectory + '/' + vecstrNodes[ i ] + c_acDAB;
	if( vecstrNodes.size( ) )
		vecstrNodes.resize( vecstrNodes.size( ) - 1 );
	m_vecpDBs.resize( iFiles );
	for( i = 0; i < m_vecpDBs.size( ); ++i ) {
		m_vecpDBs[ i ] = new CDatabaselet( );
		vecstrSubset.clear( );
		for( j = i; j < vecstrGenes.size( ); j += m_vecpDBs.size( ) )
			vecstrSubset.push_back( vecstrGenes[ j ] );
#pragma warning(disable : 4996)
		sprintf( acNumber, "%08u", i );
#pragma warning(default : 4996)
		strFile = strOutputDirectory + '/' + acNumber + c_acExtension;
		if( !( i % 100 ) )
			g_CatSleipnir.notice( "CDatabase::Open( %s, %d ) initializing file %d/%d",
				strOutputDirectory.c_str( ), iFiles, i, m_vecpDBs.size( ) );
		if( !m_vecpDBs[ i ]->Open( strFile, vecstrSubset, vecstrGenes.size( ), vecstrNodes.size( ) ) ) {
			g_CatSleipnir.error( "CDatabase::Open( %s, %d ) could not open file %s",
				strOutputDirectory.c_str( ), iFiles, strFile.c_str( ) );
			return false; } }
	for( i = 0; i < vecstrGenes.size( ); ++i )
		m_mapstriGenes[ m_vecpDBs[ i % m_vecpDBs.size( ) ]->GetGene( i / m_vecpDBs.size( ) ) ] = i;

	return CDatabaseImpl::Open( vecstrGenes, vecstrNodes ); }

bool CDatabaseImpl::Open( const vector<string>& vecstrGenes, const vector<string>& vecstrFiles ) {
	size_t			i, j, iOne, iTwo, iOutBlock, iOutBase, iOutOffset, iInBlock, iInBase, iInOffset;
	vector<size_t>	veciGenes;
	float			d;

	veciGenes.resize( vecstrGenes.size( ) );
	iOutBlock = ( m_iBlockOut == -1 ) ? m_vecpDBs.size( ) : m_iBlockOut;
	iInBlock = ( m_iBlockIn == -1 ) ? vecstrFiles.size( ) : m_iBlockIn;
	for( iOutBase = 0; iOutBase < m_vecpDBs.size( ); iOutBase += iOutBlock ) {
		vector<string>	vecstrMyGenes;
		vector<size_t>	veciMyGenes;

		for( iOutOffset = 0; ( iOutOffset < iOutBlock ) && ( ( iOutBase + iOutOffset ) < m_vecpDBs.size( ) );
			++iOutOffset ) {
			const CDatabaselet&	DB	= *m_vecpDBs[ iOutBase + iOutOffset ];

			for( i = 0; i < DB.GetGenes( ); ++i )
				vecstrMyGenes.push_back( DB.GetGene( i ) ); }
		veciMyGenes.resize( vecstrMyGenes.size( ) );
		for( iInBase = 0; iInBase < vecstrFiles.size( ); iInBase += iInBlock ) {
			vector<CCompactFullMatrix>	vecData;

			vecData.resize( ( ( iInBase + iInBlock ) > vecstrFiles.size( ) ) ?
				( vecstrFiles.size( ) - iInBase ) : iInBlock );
			for( iInOffset = 0; iInOffset < vecData.size( ); ++iInOffset ) {
				CDataPair	Dat;

				if( !Dat.Open( vecstrFiles[ iInBase + iInOffset ].c_str( ), false, m_fMemmap ) ) {
					g_CatSleipnir.error( "CDatabaseImpl::Open( ) could not open %s",
						vecstrFiles[ iInBase + iInOffset ].c_str( ) );
					return false; }
				for( i = 0; i < veciMyGenes.size( ); ++i )
					veciMyGenes[ i ] = Dat.GetGene( vecstrMyGenes[ i ] );
				for( i = 0; i < veciGenes.size( ); ++i )
					veciGenes[ i ] = Dat.GetGene( vecstrGenes[ i ] );
				vecData[ iInOffset ].Initialize( veciMyGenes.size( ), veciGenes.size( ), 16, true );
				for( i = 0; i < veciMyGenes.size( ); ++i ) {
					if( ( iOne = veciMyGenes[ i ] ) == -1 )
						continue;
					for( j = 0; j < veciGenes.size( ); ++j )
						if( ( ( iTwo = veciGenes[ j ] ) != -1 ) &&
							!CMeta::IsNaN( d = Dat.Get( iOne, iTwo ) ) )
							vecData[ iInOffset ].Set( i, j, Dat.Quantize( d ) + 1 ); } }

			for( i = iOutOffset = 0; ( iOutOffset < iOutBlock ) && ( ( iOutBase + iOutOffset ) <
				m_vecpDBs.size( ) ); ++iOutOffset ) {
				CDatabaselet&	DB	= *m_vecpDBs[ iOutBase + iOutOffset ];

				if( !( iOutOffset % 100 ) )
					cerr << "Processing offset " << iOutOffset << '/' << iOutBlock << endl;
				if( !DB.Open( vecData, i, iInBase, m_fBuffer ) )
					return false;
				i += DB.GetGenes( ); } } }

	return true; }

bool CDatabase::Open( const string& strInputDirectory ) {
	size_t			i, j;
	vector<string>	vecstrFiles;
	string			strFile;
#ifdef _MSC_VER
	HANDLE			hSearch;
	WIN32_FIND_DATA	sEntry;
	bool			fContinue;

	for( fContinue = true,hSearch = FindFirstFile( ( strInputDirectory + "/*" ).c_str( ), &sEntry );
		fContinue && ( hSearch != INVALID_HANDLE_VALUE );
		fContinue = !!FindNextFile( hSearch, &sEntry ) ) {
		strFile = sEntry.cFileName;
#else // _MSC_VER
	DIR*			pDir;
	struct dirent*	psEntry;

	pDir = opendir( strInputDirectory.c_str( ) );
	for( psEntry = readdir( pDir ); psEntry; psEntry = readdir( pDir ) ) {
		strFile = psEntry->d_name;
#endif // _MSC_VER
		if( strFile.find( c_acExtension ) == -1 )
			continue;

		i = atoi( CMeta::Deextension( CMeta::Basename( strFile.c_str( ) ) ).c_str( ) );
		if( vecstrFiles.size( ) <= i )
			vecstrFiles.resize( i + 1 );
		else if( vecstrFiles[ i ].length( ) != 0 ) {
			g_CatSleipnir.error( "CDatabase::Open( %s ) duplicate file: %s (%d)", strInputDirectory.c_str( ),
				strFile.c_str( ), i );
			return false; }
		vecstrFiles[ i ] = strInputDirectory + '/' + strFile; }

	Clear( );
	m_vecpDBs.resize( vecstrFiles.size( ) );
	for( i = 0; i < m_vecpDBs.size( ); ++i ) {
		m_vecpDBs[ i ] = new CDatabaselet( );
		if( !m_vecpDBs[ i ]->Open( vecstrFiles[ i ] ) )
			return false;
		for( j = 0; j < m_vecpDBs[ i ]->GetGenes( ); ++j )
			m_mapstriGenes[ m_vecpDBs[ i ]->GetGene( j ) ] = ( j * m_vecpDBs.size( ) ) + i; }

	return true; }

}

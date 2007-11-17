#include "stdafx.h"
#include "database.h"
#include "meta.h"
#include "bayesnet.h"

namespace libBioUtils {

const char	CDatabaseImpl::c_acDAB[]		= ".dab";
const char	CDatabaseImpl::c_acExtension[]	= ".db";

///////////////////////////////////////////////////////////////////////////////
// CDatabaselet
///////////////////////////////////////////////////////////////////////////////

CDatabaselet::CDatabaselet( bool fCache ) : m_fCache(fCache) {

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
		g_CatBioUtils.error( "CDatabaselet::Open( %s, %u, %u ) open failed", strFile.c_str( ), iGenes,
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

bool CDatabaselet::Open( const vector<string>& vecstrGenes, const vector<string>& vecstrFiles, bool fMemmap ) {
	size_t			i, j, iSize, iGeneUs, iGeneThem, iOne, iTwo;
	vector<size_t>	veciGenesUs, veciGenesThem;
	unsigned char*	abImage;
	unsigned char	bValue;

	iSize = GetSizeGene( ) * m_vecstrGenes.size( );
	abImage = new unsigned char[ iSize ];
	memset( abImage, -1, iSize );
	veciGenesUs.resize( m_vecstrGenes.size( ) );
	veciGenesThem.resize( vecstrGenes.size( ) );
	for( i = 0; i < vecstrFiles.size( ); ++i ) {
		CDataPair	Dat;

		if( !Dat.Open( vecstrFiles[ i ].c_str( ), false, fMemmap ) ) {
			g_CatBioUtils.error( "CDatabaselet::Open( %d ) could not open: %s", fMemmap,
				vecstrFiles[ i ].c_str( ) );
			return false; }
		for( iGeneUs = 0; iGeneUs < veciGenesUs.size( ); ++iGeneUs )
			veciGenesUs[ iGeneUs ] = Dat.GetGene( m_vecstrGenes[ iGeneUs ] );
		for( iGeneThem = 0; iGeneThem < veciGenesThem.size( ); ++iGeneThem )
			veciGenesThem[ iGeneThem ] = Dat.GetGene( vecstrGenes[ iGeneThem ] );
		for( iGeneUs = 0; iGeneUs < veciGenesUs.size( ); ++iGeneUs ) {
			if( ( iOne = veciGenesUs[ iGeneUs ] ) == -1 )
				continue;
			for( iGeneThem = 0; iGeneThem < veciGenesThem.size( ); ++iGeneThem )
				if( ( ( iTwo = veciGenesThem[ iGeneThem ] ) != -1 ) &&
					( ( bValue = Dat.Quantize( Dat.Get( iOne, iTwo ) ) ) != 0xFF ) ) {
					j = GetOffset( iGeneUs, iGeneThem, i ) - m_iHeader;
#ifdef DATABASE_NIBBLES
					unsigned char	b;
					b = abImage[ j ];
					bValue = ( i % 2 ) ? ( ( b & 0xF ) | ( bValue << 4 ) ) :
						( ( b & 0xF0 ) | ( bValue & 0xF ) );
#endif // DATABASE_NIBBLES
					abImage[ j ] = bValue; } } }
	m_fstm.seekp( m_iHeader );
	m_fstm.write( (char*)abImage, iSize );
	delete[] abImage;

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

bool CDatabaselet::Open( const CDatasetCompact& Data, size_t iBase, const vector<size_t>& veciGenes,
	bool fBuffer ) {
	unsigned char*	abImage;
	size_t			iSize, iDatum, iGeneOne, iGeneTwo, iOne, iTwo, iValueOne, iValueTwo;
	unsigned char	bOne, bTwo;
	vector<size_t>	veciMyGenes;

	veciMyGenes.resize( GetGenes( ) );
	for( iGeneOne = 0; iGeneOne < veciMyGenes.size( ); ++iGeneOne )
		veciMyGenes[ iGeneOne ] = Data.GetGene( GetGene( iGeneOne ) );
	if( fBuffer ) {
		abImage = new unsigned char[ iSize = ( GetSizeGene( ) * m_vecstrGenes.size( ) ) ];
		m_fstm.seekg( m_iHeader );
		m_fstm.read( (char*)abImage, iSize ); }
	else
		abImage = NULL;
	if( iBase % 2 )
		for( iGeneOne = 0; iGeneOne < veciMyGenes.size( ); ++iGeneOne )
			if( ( iOne = veciMyGenes[ iGeneOne ] ) != -1 )
				for( iGeneTwo = 0; iGeneTwo < veciGenes.size( ); ++iGeneTwo )
					if( ( ( iTwo = veciGenes[ iGeneTwo ] ) != -1 ) &&
						( ( iValueOne = Data.GetDiscrete( iOne, iTwo, 0 ) ) != -1 ) )
						OpenWrite( iValueOne, GetOffset( iGeneOne, iGeneTwo, iBase ), ENibblesHigh, abImage );
	for( iDatum = ( iBase % 2 ); ( iDatum + 1 ) < Data.GetExperiments( ); iDatum += 2 )
		for( iGeneOne = 0; iGeneOne < veciMyGenes.size( ); ++iGeneOne ) {
			if( ( iOne = veciMyGenes[ iGeneOne ] ) == -1 )
				continue;
			for( iGeneTwo = 0; iGeneTwo < veciGenes.size( ); ++iGeneTwo ) {
				if( ( iTwo = veciGenes[ iGeneTwo ] ) == -1 )
					continue;
				iValueOne = Data.GetDiscrete( iOne, iTwo, iDatum );
				iValueTwo = Data.GetDiscrete( iOne, iTwo, iDatum + 1 );
				if( ( iValueOne == -1 ) && ( iValueTwo == -1 ) )
					continue;
				bOne = ( iValueOne == -1 ) ? 0xFF : iValueOne;
				bTwo = ( iValueTwo == -1 ) ? 0xFF : iValueTwo;
#ifdef DATABASE_NIBBLES
				OpenWrite( ( bOne & 0xF ) | ( bTwo << 4 ), GetOffset( iGeneOne, iGeneTwo, iBase + iDatum ),
					ENibblesBoth, abImage );
#else // DATABASE_NIBBLES
				OpenWrite( bOne, GetOffset( iGeneOne, iGeneTwo, iBase + iDatum ), ENibblesBoth, abImage );
				OpenWrite( bTwo, GetOffset( iGeneOne, iGeneTwo, iBase + iDatum + 1 ), ENibblesBoth, abImage );
#endif // DATABASE_NIBBLES
			} }
	if( iDatum < Data.GetExperiments( ) )
		for( iGeneOne = 0; iGeneOne < veciMyGenes.size( ); ++iGeneOne )
			if( ( iOne = veciMyGenes[ iGeneOne ] ) != -1 )
				for( iGeneTwo = 0; iGeneTwo < veciGenes.size( ); ++iGeneTwo )
					if( ( ( iTwo = veciGenes[ iGeneTwo ] ) != -1 ) &&
						( ( iValueOne = Data.GetDiscrete( iOne, iTwo, iDatum ) ) != -1 ) )
						OpenWrite( iValueOne, GetOffset( iGeneOne, iGeneTwo, iBase + iDatum ), ENibblesLow,
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

bool CDatabaselet::Get( size_t iGene, vector<unsigned char>& vecbData ) const {
	size_t	i;

	i = vecbData.size( );
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
		g_CatBioUtils.error( "CDatabaselet::Open( %s ) open failed", strFile.c_str( ) );
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

bool CDatabaselet::Write( ) {
	size_t	i;

	sort( m_vecpribWrites.begin( ), m_vecpribWrites.end( ), SSorter( ) );
	for( i = 0; i < m_vecpribWrites.size( ); ++i ) {
		m_fstm.seekp( m_vecpribWrites[ i ].first );
		m_fstm.put( m_vecpribWrites[ i ].second ); }
	m_vecpribWrites.clear( );

	return true; }

///////////////////////////////////////////////////////////////////////////////
// CDatabase
///////////////////////////////////////////////////////////////////////////////

bool CDatabase::Open( const vector<string>& vecstrGenes, const string& strData, const IBayesNet* pBN,
	const string& strDir, size_t iFiles ) {
	vector<string>	vecstrNodes, vecstrSubset;
	size_t			i, j;
	char			acNumber[ 16 ];
	string			strFile;

	if( !pBN ) {
		g_CatBioUtils.error( "CDatabase::Open( %s, %d ) null Bayes net", strDir.c_str( ), iFiles );
		return false; }

	Clear( );
	pBN->GetNodes( vecstrNodes );
	for( i = 1; i < vecstrNodes.size( ); ++i )
		vecstrNodes[ i - 1 ] = strData + '/' + vecstrNodes[ i ] + c_acDAB;
	if( vecstrNodes.size( ) )
		vecstrNodes.resize( vecstrNodes.size( ) - 1 );
	m_vecpDBs.resize( iFiles );
	for( i = 0; i < m_vecpDBs.size( ); ++i ) {
		m_vecpDBs[ i ] = new CDatabaselet( m_fCache );
		vecstrSubset.clear( );
		for( j = i; j < vecstrGenes.size( ); j += m_vecpDBs.size( ) )
			vecstrSubset.push_back( vecstrGenes[ j ] );
#pragma warning(disable : 4996)
		sprintf( acNumber, "%08u", i );
#pragma warning(default : 4996)
		strFile = strDir + '/' + acNumber + c_acExtension;
		if( !( i % 100 ) )
			g_CatBioUtils.notice( "CDatabase::Open( %s, %d ) initializing file %d/%d", strDir.c_str( ),
					iFiles, i, m_vecpDBs.size( ) );
		if( !m_vecpDBs[ i ]->Open( strFile, vecstrSubset, vecstrGenes.size( ), vecstrNodes.size( ) ) ) {
			g_CatBioUtils.error( "CDatabase::Open( %s, %d ) could not open file %s", strDir.c_str( ),
				iFiles, strFile.c_str( ) );
			return false; } }
	for( i = 0; i < vecstrGenes.size( ); ++i )
		m_mapstriGenes[ m_vecpDBs[ i % m_vecpDBs.size( ) ]->GetGene( i / m_vecpDBs.size( ) ) ] = i;

	return CDatabaseImpl::Open( vecstrGenes, vecstrNodes ); }

/*
bool CDatabaseImpl::Open( const vector<string>& vecstrGenes, bool fMemmap,
	const vector<string>& vecstrFiles ) {
	size_t	i;

	for( i = 0; i < m_vecpDBs.size( ); ++i )
		if( !m_vecpDBs[ i ]->Open( vecstrGenes, vecstrFiles, fMemmap ) )
			return false;

	return true; }

// Memmapped:		5m46.394s
// Not memmapped:	5m48.374s
bool CDatabaseImpl::Open( const vector<string>& vecstrGenes, bool fMemmap,
	const vector<string>& vecstrFiles ) {
	vector<size_t>	veciGenes;
	unsigned char	b;
	size_t			i, j, k, iDBOne, iOne, iTwo;

	veciGenes.resize( vecstrGenes.size( ) );
	for( i = 0; i < vecstrFiles.size( ); ++i ) {
		CDataPair	Dat;

		if( !Dat.Open( vecstrFiles[ i ].c_str( ), false, fMemmap ) ) {
			g_CatBioUtils.error( "CDatabaseImpl::Open( %d ) could not open %s", fMemmap,
				vecstrFiles[ i ].c_str( ) );
			return false; }
		for( j = 0; j < veciGenes.size( ); ++j )
			veciGenes[ j ] = Dat.GetGene( vecstrGenes[ j ] );
		for( j = 0; j < veciGenes.size( ); ++j ) {
			CDatabaselet&	DB	= *m_vecpDBs[ j % m_vecpDBs.size( ) ];

			if( !( j % 100 ) )
				g_CatBioUtils.notice( "CDatabaseImpl::Open( %d ) gene %d/%d", fMemmap, j, veciGenes.size( ) );
			if( ( iOne = veciGenes[ j ] ) == -1 )
				continue;
			iDBOne = j / m_vecpDBs.size( );
			for( k = 0; k < veciGenes.size( ); ++k )
				if( ( ( iTwo = veciGenes[ k ] ) != -1 ) &&
					( ( b = Dat.Quantize( Dat.Get( iOne, iTwo ) ) ) != 0xFF ) )
					DB.Write( iDBOne, k, i, b ); } }

	return true; }

// Memmapped:		2m51.354s
// Not memmapped:	2m42.590s
// Open each DAB once
// Process gene pairs only from that DAB
// Super-slow for some reason; shoule be fast!
bool CDatabaseImpl::Open( const vector<string>& vecstrGenes, bool fMemmap,
	const vector<string>& vecstrFiles ) {
	size_t	i;

	for( i = 0; i < vecstrFiles.size( ); ++i )
		if( !Open( vecstrFiles[ i ], i, fMemmap ) )
			return false;

	return true; }

bool CDatabaseImpl::Open( const string& strFile, size_t iDataset, bool fBoth ) {
	CDataPair		Dat;
	vector<size_t>	veciGenes;
	size_t			i, j, iOne, iTwo, iDBOne;
	CDatabaselet*	pDB;
	unsigned char	b;

	if( !Dat.Open( strFile.c_str( ), false, m_fMemmap ) ) {
		g_CatBioUtils.error( "CDatabaseImpl::Open( %s ) failed", strFile.c_str( ) );
		return false; }
	veciGenes.resize( Dat.GetGenes( ) );
	for( i = 0; i < veciGenes.size( ); ++i )
		veciGenes[ i ] = GetGene( Dat.GetGene( i ) );
	for( i = 0; i < veciGenes.size( ); ++i ) {
		if( ( iOne = veciGenes[ i ] ) == -1 )
			continue;
		if( !( i % 100 ) )
			g_CatBioUtils.notice( "CDatabaseImpl::Open( %s ) gene %d/%d", strFile.c_str( ), i,
				veciGenes.size( ) );
		pDB = m_vecpDBs[ iOne % m_vecpDBs.size( ) ];
		iDBOne = iOne / m_vecpDBs.size( );
		for( j = ( i + 1 ); j < veciGenes.size( ); ++j )
			if( ( ( iTwo = veciGenes[ j ] ) != -1 ) &&
				( ( b = Dat.Quantize( Dat.Get( i, j ) ) ) != 0xFF ) ) {
				pDB->Write( iDBOne, iTwo, iDataset, b, fBoth );
				m_vecpDBs[ iTwo % m_vecpDBs.size( ) ]->Write( iTwo / m_vecpDBs.size( ), iOne, iDataset, b,
					fBoth ); } }

	return true; }

// Memmapped:		1m31.743s
// Not memmapped:	4m41.680s
// Open each DAB once in pairs
//   This avoids reads on writes
// Process gene pairs only from those DABs
bool CDatabaseImpl::Open( const vector<string>& vecstrGenes, const vector<string>& vecstrFiles ) {
	vector<size_t>	veciGenes;
	size_t			i, j, k, iDBOne, iOne, iTwo, iFirst, iSecond;
	CDatabaselet*	pDB;
	vector<string>	vecstrCur;
	unsigned char	b;

	vecstrCur.resize( 2 );
	for( i = 0; ( i + vecstrCur.size( ) - 1 ) < vecstrFiles.size( ); i += vecstrCur.size( ) ) {
		CDatasetCompact	Dats;

		for( j = 0; j < vecstrCur.size( ); ++j )
			vecstrCur[ j ] = vecstrFiles[ i + j ];
		if( !Dats.Open( vecstrCur, m_fMemmap ) ) {
			g_CatBioUtils.error( "CDatabaseImpl::Open( ) could not open data" );
			return false; }
		veciGenes.resize( Dats.GetGenes( ) );
		for( j = 0; j < veciGenes.size( ); ++j )
			veciGenes[ j ] = GetGene( Dats.GetGene( j ) );
		for( j = 0; j < veciGenes.size( ); ++j ) {
			if( ( iOne = veciGenes[ j ] ) == -1 )
				continue;
			if( !( j % 100 ) )
				g_CatBioUtils.notice( "CDatabaseImpl::Open( ) gene %d/%d", j, veciGenes.size( ) );
			pDB = m_vecpDBs[ iOne % m_vecpDBs.size( ) ];
			iDBOne = iOne / m_vecpDBs.size( );
			for( k = ( j + 1 ); k < veciGenes.size( ); ++k ) {
				if( ( iTwo = veciGenes[ k ] ) == -1 )
					continue;
				iFirst = Dats.GetDiscrete( j, k, 0 );
				iSecond = Dats.GetDiscrete( j, k, 1 );
				if( ( iFirst != -1 ) || ( iSecond != -1 ) ) {
					b = ( iFirst & 0xF ) | ( ( iSecond << 4 ) & 0xF0 );
					pDB->Write( iDBOne, iTwo, i, b, true );
					m_vecpDBs[ iTwo % m_vecpDBs.size( ) ]->Write( iTwo / m_vecpDBs.size( ), iOne, i, b,
						true ); } } }
		if( m_fCache )
			for( j = 0; j < m_vecpDBs.size( ); ++j )
				if( !m_vecpDBs[ j ]->Write( ) )
					return false; }
	for( ; i < vecstrFiles.size( ); ++i )
		if( !Open( vecstrFiles[ i ], i, true ) )
			return false;
	if( m_fCache )
		for( j = 0; j < m_vecpDBs.size( ); ++j )
			if( !m_vecpDBs[ j ]->Write( ) )
				return false;

	return true; }
*/

// Block on genes (output files) and datasets (input files).
bool CDatabaseImpl::Open( const vector<string>& vecstrGenes, const vector<string>& vecstrFiles ) {
	size_t			iOutBlock, iOutBase, iOutOffset, iInBase;
	vector<string>	vecstrIn;
	vector<size_t>	veciGenes;

	veciGenes.resize( vecstrGenes.size( ) );
	vecstrIn.resize( ( m_iBlockIn == -1 ) ? vecstrFiles.size( ) : m_iBlockIn );
	iOutBlock = ( m_iBlockOut == -1 ) ? m_vecpDBs.size( ) : m_iBlockOut;
	for( iOutBase = 0; iOutBase < m_vecpDBs.size( ); iOutBase += iOutBlock )
		for( iInBase = 0; iInBase < vecstrFiles.size( ); iInBase += vecstrIn.size( ) ) {
			CDatasetCompact					Data;
			vector<string>::const_iterator	iterBase;

			if( ( iInBase + vecstrIn.size( ) ) > vecstrFiles.size( ) )
				vecstrIn.resize( vecstrFiles.size( ) - iInBase );
			iterBase = vecstrFiles.begin( ) + iInBase;
			copy( iterBase, iterBase + vecstrIn.size( ), vecstrIn.begin( ) );
			if( !Data.Open( vecstrIn, m_fMemmap ) ) {
				g_CatBioUtils.error( "CDatabaseImpl::Open( ) could not open data block %d", iInBase );
				return false; }

			for( iOutOffset = 0; iOutOffset < veciGenes.size( ); ++iOutOffset )
				veciGenes[ iOutOffset ] = Data.GetGene( vecstrGenes[ iOutOffset ] );
			for( iOutOffset = 0; iOutOffset < iOutBlock; ++iOutOffset )
				if( !m_vecpDBs[ iOutBase + iOutOffset ]->Open( Data, iInBase, veciGenes, m_fBuffer ) )
					return false; }

	return true; }

bool CDatabase::Open( const string& strDir ) {
	size_t			i, j;
	vector<string>	vecstrFiles;
	string			strFile;
#ifdef _MSC_VER
	HANDLE			hSearch;
	WIN32_FIND_DATA	sEntry;
	bool			fContinue;

	for( fContinue = true,hSearch = FindFirstFile( ( strDir + "/*" ).c_str( ), &sEntry );
		fContinue && ( hSearch != INVALID_HANDLE_VALUE );
		fContinue = !!FindNextFile( hSearch, &sEntry ) ) {
		strFile = sEntry.cFileName;
#else // _MSC_VER
	DIR*			pDir;
	struct dirent*	psEntry;

	pDir = opendir( strDir.c_str( ) );
	for( psEntry = readdir( pDir ); psEntry; psEntry = readdir( pDir ) ) {
		strFile = psEntry->d_name;
#endif // _MSC_VER
		if( strFile.find( c_acExtension ) == -1 )
			continue;

		i = atoi( CMeta::Deextension( CMeta::Basename( strFile.c_str( ) ) ).c_str( ) );
		if( vecstrFiles.size( ) <= i )
			vecstrFiles.resize( i + 1 );
		else if( vecstrFiles[ i ].length( ) != 0 ) {
			g_CatBioUtils.error( "CDatabase::Open( %s ) duplicate file: %s (%d)", strDir.c_str( ),
				strFile.c_str( ), i );
			return false; }
		vecstrFiles[ i ] = strDir + '/' + strFile; }

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

#include "stdafx.h"
#include "database.h"
#include "meta.h"
#include "bayesnet.h"

namespace libBioUtils {

const char	CDatabaseImpl::c_acDAB[]		= ".dab";
const char	CDatabaseImpl::c_acExtension[]	= ".db";

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
	acFiller = new char[ GetOffset( ) ];
	memset( acFiller, -1, GetOffset( ) );
	for( i = 0; i < m_vecstrGenes.size( ); ++i )
		m_fstm.write( acFiller, GetOffset( ) );
	delete[] acFiller;

	return true; }

bool CDatabaselet::Open( const vector<string>& vecstrGenes, const vector<string>& vecstrFiles, bool fMemmap ) {
	size_t			i, j, iSize, iGeneUs, iGeneThem, iOne, iTwo;
	vector<size_t>	veciGenesUs, veciGenesThem;
	unsigned char*	abImage;
	unsigned char	bValue;

	iSize = GetOffset( ) * m_vecstrGenes.size( );
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

void CDatabaselet::Write( size_t iOne, size_t iTwo, size_t iDataset, unsigned char bValue ) {

	m_fstm.seekp( GetOffset( iOne, iTwo, iDataset ) );
#ifdef DATABASE_NIBBLES
	unsigned char	b;
	m_fstm.read( (char*)&b, sizeof(b) );
	bValue = ( iDataset % 2 ) ? ( ( b & 0xF ) | ( bValue << 4 ) ) : ( ( b & 0xF0 ) | ( bValue & 0xF ) );
	m_fstm.seekp( -(streamoff)sizeof(b), ios_base::cur );
#endif // DATABASE_NIBBLES
	m_fstm.write( (char*)&bValue, sizeof(bValue) ); }

bool CDatabaselet::Get( size_t iOne, size_t iTwo, vector<unsigned char>& vecbData ) const {

	pthread_mutex_lock( m_pmutx );
	m_fstm.seekg( GetOffset( iOne, iTwo ) );
	vecbData.resize( GetDatasets( ) );
	m_fstm.read( (char*)&vecbData[ 0 ], vecbData.size( ) );
	pthread_mutex_unlock( m_pmutx );

	return true; }

bool CDatabaselet::Get( size_t iGene, vector<unsigned char>& vecbData ) const {

	pthread_mutex_lock( m_pmutx );
	m_fstm.seekg( GetOffset( iGene ) );
	vecbData.resize( GetOffset( ) );
	m_fstm.read( (char*)&vecbData[ 0 ], vecbData.size( ) );
	pthread_mutex_unlock( m_pmutx );

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

bool CDatabase::Open( const vector<string>& vecstrGenes, const string& strData, const IBayesNet* pBN,
	const string& strDir, size_t iFiles, bool fMemmap ) {
	vector<string>	vecstrNodes, vecstrSubset;
	size_t			i, j;
	char			acNumber[ 16 ];
	string			strFile;

	if( !pBN ) {
			g_CatBioUtils.error( "CDatabase::Open( %s, %d, %d ) null Bayes net", strDir.c_str( ),
				iFiles, fMemmap );
		return false; }

	Clear( );
	pBN->GetNodes( vecstrNodes );
	m_vecpDBs.resize( iFiles );
	for( i = 0; i < m_vecpDBs.size( ); ++i ) {
		m_vecpDBs[ i ] = new CDatabaselet( );
		vecstrSubset.clear( );
		for( j = i; j < vecstrGenes.size( ); j += m_vecpDBs.size( ) )
			vecstrSubset.push_back( vecstrGenes[ j ] );
#pragma warning(disable : 4996)
		sprintf( acNumber, "%08u", i );
#pragma warning(default : 4996)
		strFile = strDir + '/' + acNumber + c_acExtension;
		if( !m_vecpDBs[ i ]->Open( strFile, vecstrSubset, vecstrGenes.size( ), vecstrNodes.size( ) - 1 ) ) {
			g_CatBioUtils.error( "CDatabase::Open( %s, %d, %d ) could not open file %s", strDir.c_str( ),
				iFiles, fMemmap, strFile.c_str( ) );
			return false; } }
	for( i = 0; i < vecstrGenes.size( ); ++i )
		m_mapstriGenes[ m_vecpDBs[ i % m_vecpDBs.size( ) ]->GetGene( i / m_vecpDBs.size( ) ) ] = i;

	return CDatabaseImpl::Open( vecstrGenes, strData, fMemmap, vecstrNodes ); }

/*
bool CDatabaseImpl::Open( const vector<string>& vecstrGenes, const string& strData, bool fMemmap,
	const vector<string>& vecstrNodes ) {
	size_t			i;
	vector<string>	vecstrFiles;

	vecstrFiles.resize( vecstrNodes.size( ) - 1 );
	for( i = 0; ( i + 1 ) < vecstrNodes.size( ); ++i )
		vecstrFiles[ i ] = strData + '/' + vecstrNodes[ i + 1 ] + c_acDAB;

	for( i = 0; i < m_vecpDBs.size( ); ++i )
		if( !m_vecpDBs[ i ]->Open( vecstrGenes, vecstrFiles, fMemmap ) )
			return false;

	return true; }

bool CDatabaseImpl::Open( const vector<string>& vecstrGenes, const string& strData, bool fMemmap,
	const vector<string>& vecstrNodes ) {
	vector<size_t>	veciGenes;
	unsigned char	b;
	string			strFile;
	size_t			i, j, k, iDBOne, iOne, iTwo;

	veciGenes.resize( vecstrGenes.size( ) );
	for( i = 0; ( i + 1 ) < vecstrNodes.size( ); ++i ) {
		CDataPair	Dat;

		strFile = strData + '/' + vecstrNodes[ i + 1 ] + c_acDAB;
		if( !Dat.Open( strFile.c_str( ), false, fMemmap ) ) {
			g_CatBioUtils.error( "CDatabaseImpl::Open( %s, %d ) could not open %s", strData.c_str( ), fMemmap,
				strFile.c_str( ) );
			return false; }
		for( j = 0; j < veciGenes.size( ); ++j )
			veciGenes[ j ] = Dat.GetGene( vecstrGenes[ j ] );
		for( j = 0; j < veciGenes.size( ); ++j ) {
			CDatabaselet&	DB	= *m_vecpDBs[ j % m_vecpDBs.size( ) ];

			if( !( j % 100 ) )
				g_CatBioUtils.notice( "CDatabaseImpl::Open( %s, %d ) gene %d/%d", strData.c_str( ), fMemmap,
					j, veciGenes.size( ) );
			if( ( iOne = veciGenes[ j ] ) == -1 )
				continue;
			iDBOne = j / m_vecpDBs.size( );
			for( k = 0; k < veciGenes.size( ); ++k )
				if( ( ( iTwo = veciGenes[ k ] ) != -1 ) &&
					( ( b = Dat.Quantize( Dat.Get( iOne, iTwo ) ) ) != 0xFF ) )
					DB.Write( iDBOne, k, i, b ); } }

	return true; }
*/

bool CDatabaseImpl::Open( const vector<string>& vecstrGenes, const string& strData, bool fMemmap,
	const vector<string>& vecstrNodes ) {
	vector<size_t>	veciGenes;
	unsigned char	b;
	string			strFile;
	size_t			i, j, k, iDBOne, iOne, iTwo;
	CDatabaselet*	pDB;

	for( i = 0; ( i + 1 ) < vecstrNodes.size( ); ++i ) {
		CDataPair	Dat;

		strFile = strData + '/' + vecstrNodes[ i + 1 ] + c_acDAB;
		if( !Dat.Open( strFile.c_str( ), false, fMemmap ) ) {
			g_CatBioUtils.error( "CDatabaseImpl::Open( %s, %d ) could not open %s", strData.c_str( ), fMemmap,
				strFile.c_str( ) );
			return false; }
		veciGenes.resize( Dat.GetGenes( ) );
		for( j = 0; j < veciGenes.size( ); ++j )
			veciGenes[ j ] = GetGene( Dat.GetGene( j ) );
		for( j = 0; j < veciGenes.size( ); ++j ) {
			if( ( iOne = veciGenes[ j ] ) == -1 )
				continue;
			if( !( j % 100 ) )
				g_CatBioUtils.notice( "CDatabaseImpl::Open( %s, %d ) gene %d/%d", strData.c_str( ), fMemmap,
					j, veciGenes.size( ) );
			pDB = m_vecpDBs[ iOne % m_vecpDBs.size( ) ];
			iDBOne = iOne / m_vecpDBs.size( );
			for( k = ( j + 1 ); k < veciGenes.size( ); ++k )
				if( ( ( iTwo = veciGenes[ k ] ) != -1 ) &&
					( ( b = Dat.Quantize( Dat.Get( j, k ) ) ) != 0xFF ) ) {
					pDB->Write( iDBOne, iTwo, i, b );
					m_vecpDBs[ iTwo % m_vecpDBs.size( ) ]->Write( iTwo / m_vecpDBs.size( ), iOne, i, b ); } } }

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

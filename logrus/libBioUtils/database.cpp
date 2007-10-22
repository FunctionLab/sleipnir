#include "stdafx.h"
#include "database.h"
#include "meta.h"
#include "bayesnet.h"

namespace libBioUtils {

const char	CDatabaseImpl::c_acDAB[]		= ".dab";
const char	CDatabaseImpl::c_acExtension[]	= ".db";

bool CDatabaselet::Open( const string& strFile, const vector<string>& vecstrGenes, uint32_t iGenes,
	uint32_t iDatasets ) {
	uint32_t	iSize;
	size_t		i;

	if( !m_pfstm )
		m_pfstm = new fstream( );

	m_pfstm->open( strFile.c_str( ), ios_base::in | ios_base::out | ios_base::binary | ios_base::trunc );
	if( !m_pfstm->is_open( ) )
		return false;

	m_iGenes = iGenes;
	m_iDatasets = iDatasets;
	m_vecstrGenes.resize( vecstrGenes.size( ) );
	copy( vecstrGenes.begin( ), vecstrGenes.end( ), m_vecstrGenes.begin( ) );

	m_pfstm->write( (char*)&m_iHeader, sizeof(m_iHeader) );
	m_pfstm->write( (char*)&m_iGenes, sizeof(m_iGenes) );
	m_pfstm->write( (char*)&m_iDatasets, sizeof(m_iDatasets) );
	iSize = m_vecstrGenes.size( );
	m_pfstm->write( (char*)&iSize, sizeof(iSize) );
	m_iHeader = sizeof(m_iHeader) + sizeof(m_iGenes) + sizeof(m_iDatasets) + sizeof(iSize);
	for( i = 0; i < m_vecstrGenes.size( ); ++i ) {
		m_pfstm->write( m_vecstrGenes[ i ].c_str( ), m_vecstrGenes[ i ].size( ) + 1 );
		m_iHeader += m_vecstrGenes[ i ].size( ) + 1; }
	m_pfstm->seekp( 0 );
	m_pfstm->write( (char*)&m_iHeader, sizeof(m_iHeader) );

	return true; }

void CDatabaselet::Write( size_t iOne, size_t iTwo, size_t iDataset, unsigned char bValue ) {

	m_pfstm->seekp( GetOffset( iOne, iTwo, iDataset ) );
	m_pfstm->write( (char*)&bValue, sizeof(bValue) ); }

bool CDatabaselet::Get( size_t iOne, size_t iTwo, vector<unsigned char>& vecbData ) const {

	m_pfstm->seekg( GetOffset( iOne, iTwo ) );
	vecbData.resize( m_iDatasets );
	m_pfstm->read( (char*)&vecbData[ 0 ], vecbData.size( ) );

	return true; }

bool CDatabaselet::Open( const string& strFile ) {
	uint32_t	iSize;
	char*		acBuffer;
	char*		pc;
	size_t		i;

	if( !m_pfstm )
		m_pfstm = new fstream( );

	m_pfstm->open( strFile.c_str( ), ios_base::binary | ios_base::in );
	if( !m_pfstm->is_open( ) )
		return false;

	m_pfstm->read( (char*)&m_iHeader, sizeof(m_iHeader) );
	m_pfstm->read( (char*)&m_iGenes, sizeof(m_iGenes) );
	m_pfstm->read( (char*)&m_iDatasets, sizeof(m_iDatasets) );
	m_pfstm->read( (char*)&iSize, sizeof(iSize) );

	acBuffer = new char[ m_iHeader ];
	m_pfstm->read( acBuffer, m_iHeader - sizeof(m_iHeader) - sizeof(m_iGenes) - sizeof(m_iDatasets) -
		sizeof(iSize) );
	m_vecstrGenes.resize( iSize );
	for( i = 0,pc = acBuffer; i < m_vecstrGenes.size( ); ++i,pc += m_vecstrGenes[ i ].length( ) + 1 )
		m_vecstrGenes[ i ] = pc;

	return true; }

bool CDatabase::Open( const vector<string>& vecstrGenes, const string& strData, const IBayesNet* pBN,
	const string& strDir, size_t iFiles, bool fMemmap ) {
	vector<string>	vecstrNodes, vecstrSubset;
	size_t			i, j, k, iDatasets, iDBOne, iOne, iTwo;
	char			acNumber[ 16 ];
	string			strFile;
	vector<size_t>	veciGenes;
	unsigned char	b;

	if( !pBN )
		return false;

	pBN->GetNodes( vecstrNodes );
	iDatasets = vecstrNodes.size( ) - 1;

	m_vecDBs.resize( iFiles );
	for( i = 0; i < m_vecDBs.size( ); ++i ) {
		vecstrSubset.clear( );
		for( j = i; j < vecstrGenes.size( ); j += m_vecDBs.size( ) )
			vecstrSubset.push_back( vecstrGenes[ j ] );
#pragma warning(disable : 4996)
		sprintf( acNumber, "%08u", i );
#pragma warning(default : 4996)
		strFile = strDir + '/' + acNumber + c_acExtension;
		if( !m_vecDBs[ i ].Open( strFile, vecstrSubset, vecstrGenes.size( ), iDatasets ) ) {
			g_CatBioUtils.error( "CDatabase::Open( %s, %d, %d ) could not open file %s", strDir.c_str( ),
				iFiles, fMemmap, strFile.c_str( ) );
			return false; } }

	veciGenes.resize( vecstrGenes.size( ) );
	for( i = 0; i < iDatasets; ++i ) {
		CDataPair	Dat;

		strFile = strData + '/' + vecstrNodes[ i + 1 ] + c_acDAB;
		if( !Dat.Open( strFile.c_str( ), false, fMemmap ) ) {
			g_CatBioUtils.error( "CDatabase::Open( %s, %d, %d ) could not open %s", strDir.c_str( ), iFiles,
				fMemmap, strFile.c_str( ) );
			return false; }
		for( j = 0; j < veciGenes.size( ); ++j )
			veciGenes[ j ] = Dat.GetGene( vecstrGenes[ j ] );
		for( j = 0; j < veciGenes.size( ); ++j ) {
			CDatabaselet&	DB	= m_vecDBs[ j % m_vecDBs.size( ) ];

			iDBOne = DB.GetGene( vecstrGenes[ j ] );
			iOne = veciGenes[ j ];
			for( k = 0; k < veciGenes.size( ); ++k ) {
				b = ( ( iOne == -1 ) || ( ( iTwo = veciGenes[ k ] ) == -1 ) ) ? -1 :
					Dat.Quantize( Dat.Get( iOne, iTwo ) );
				DB.Write( iDBOne, k, i, b ); } } }

	return true; }

bool CDatabase::Open( const string& strDir ) {
	size_t			i;
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
			vecstrFiles.resize( i );
		else if( vecstrFiles[ i ].length( ) != 0 ) {
			g_CatBioUtils.error( "CDatabase::Open( %s ) duplicate file: %s (%d)", strDir.c_str( ),
				strFile.c_str( ), i );
			return false; }
		vecstrFiles[ i ] = strFile; }

	m_vecDBs.resize( vecstrFiles.size( ) );
	for( i = 0; i < m_vecDBs.size( ); ++i )
		if( !m_vecDBs[ i ].Open( vecstrFiles[ i ] ) )
			return false;

	return true; }

}

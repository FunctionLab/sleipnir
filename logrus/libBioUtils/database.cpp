#include "stdafx.h"
#include "database.h"
#include "meta.h"
#include "bayesnet.h"

namespace libBioUtils {

const char	CDatabaseImpl::c_acDAB[]		= ".dab";
const char	CDatabaseImpl::c_acExtension[]	= ".db";

bool CDatabase::Open( const vector<string>& vecstrGenes, const string& strData, const IBayesNet* pBN,
const string& strDir, bool fMemmap ) {
	vector<string>	vecstrNodes;
	vector<size_t>	veciGenes;
	size_t			i, j, k, iOne, iTwo;
	unsigned char	b;
	SHeader			sHeader;

	if( !pBN )
		return false;

	pBN->GetNodes( vecstrNodes );
	m_iDatasets = vecstrNodes.size( ) - 1;
	sHeader.m_iDatasets = m_iDatasets;

	m_vecstrGenes.resize( vecstrGenes.size( ) );
	copy( vecstrGenes.begin( ), vecstrGenes.end( ), m_vecstrGenes.begin( ) );
	sHeader.m_iGenes = m_vecstrGenes.size( );

	m_vecfile.resize( m_vecstrGenes.size( ) );
	for( i = 0; i < m_vecfile.size( ); ++i ) {
#pragma warning(disable : 4996)
		if( !( m_vecfile[ i ] = fopen( ( strDir + '/' + CMeta::Filename( m_vecstrGenes[ i ] ) +
			c_acExtension ).c_str( ), "w" ) ) ) {
#pragma warning(default : 4996)
			g_CatBioUtils.error( "CDatabase::Open( %s ) could not open gene %s", strDir.c_str( ),
				m_vecstrGenes[ i ].c_str( ) );
			return false; }
		memset( sHeader.m_acName, 0, sizeof(sHeader.m_acName) );
#pragma warning(disable : 4996)
		strcpy( sHeader.m_acName, m_vecstrGenes[ i ].c_str( ) );
#pragma warning(default : 4996)
		fwrite( &sHeader, sizeof(sHeader), 1, m_vecfile[ i ] ); }

	veciGenes.resize( m_vecstrGenes.size( ) );
	for( i = 0; i < m_iDatasets; ++i ) {
		CDataPair	Dat;
		string		strFile;

		strFile = strData + '/' + vecstrNodes[ i + 1 ] + c_acDAB;
		if( !Dat.Open( strFile.c_str( ), false, fMemmap ) ) {
			g_CatBioUtils.error( "CDatabase::Open( %s ) could not open %s", strDir.c_str( ),
				strFile.c_str( ) );
			return false; }
		for( j = 0; j < veciGenes.size( ); ++j )
			veciGenes[ j ] = Dat.GetGene( m_vecstrGenes[ j ] );
		for( j = 0; j < veciGenes.size( ); ++j ) {
			FILE*	file;

			file = m_vecfile[ j ];
			iOne = veciGenes[ j ];
			for( k = 0; k < veciGenes.size( ); ++k ) {
				b = ( ( iOne == -1 ) || ( ( iTwo = veciGenes[ k ] ) == -1 ) ) ? -1 :
					Dat.Quantize( Dat.Get( iOne, iTwo ) );
				fseek( file, GetOffset( k, i ), SEEK_SET );
				fwrite( &b, sizeof(b), 1, file ); } } }

	return true; }

bool CDatabase::Open( const vector<string>& vecstrGenes, size_t iDatasets, const string& strDir ) {
	size_t	i;

	m_iDatasets = iDatasets;
	m_vecstrGenes.resize( vecstrGenes.size( ) );
	copy( vecstrGenes.begin( ), vecstrGenes.end( ), m_vecstrGenes.begin( ) );
	m_vecfile.resize( m_vecstrGenes.size( ) );
	for( i = 0; i < m_vecstrGenes.size( ); ++i )
#pragma warning(disable : 4996)
		m_vecfile[ i ] = fopen( ( strDir + '/' + CMeta::Filename( m_vecstrGenes[ i ] ) +
			c_acExtension ).c_str( ), "r" );
#pragma warning(default : 4996)

	return true; }

bool CDatabase::Get( size_t iOne, size_t iTwo, vector<unsigned char>& vecbData ) const {
	FILE*	file;

	if( iOne >= m_iDatasets )
		return false;

	if( !( file = m_vecfile[ iOne ] ) || fseek( file, GetOffset( iTwo ), SEEK_SET ) )
		return false;
	vecbData.resize( m_iDatasets );

	return ( fread( &vecbData[ 0 ], 1, m_iDatasets, file ) == m_iDatasets ); }

}

#include "stdafx.h"
#include "datapair.h"
#include "meta.h"

namespace libBioUtils {

const char	CDataPairImpl::c_szQuantExt[]	= ".quant";

bool CDataPair::Open( const CSlim& Slim ) {

	Reset( false );
	return CDat::Open( Slim ); }

bool CDataPair::Open( const char* szDatafile, bool fContinuous ) {
	static const size_t	c_iBuf	= 1024;
	char			szBuf[ c_iBuf ];
	ifstream		ifsmData, ifsmQuant;
	string			strToken;
	const char*		pc;
	size_t			i;
	vector<string>	vecstrQuant;

	g_CatBioUtils.notice( "CDataPair::Open( %s, %d )", szDatafile, fContinuous );

	Reset( fContinuous );

	if( !CDat::Open( szDatafile ) )
		return false;
	if( m_fContinuous )
		return true;

	strToken = szDatafile;
	strToken += c_szQuantExt;
	ifsmQuant.open( strToken.c_str( ) );
	if( !ifsmQuant.is_open( ) ) {
		if( !( pc = strrchr( szDatafile, '.' ) ) )
			return false;
		strToken = szDatafile;
		strToken.resize( pc - szDatafile );
		strToken += c_szQuantExt;
		ifsmQuant.clear( );
		ifsmQuant.open( strToken.c_str( ) );
		if( !ifsmQuant.is_open( ) )
			return false; }
	ifsmQuant.getline( szBuf, c_iBuf - 1 );
	CMeta::Tokenize( szBuf, vecstrQuant, CMeta::c_szWS, true );
	m_vecdQuant.resize( vecstrQuant.size( ) );
	for( i = 0; i < m_vecdQuant.size( ); ++i )
		m_vecdQuant[ i ] = (float)atof( vecstrQuant[ i ].c_str( ) );
	ifsmQuant.close( );

	return true; }

size_t CDataPair::Quantify( float dValue ) const {
	size_t	i;

	for( i = 0; i < m_vecdQuant.size( ); ++i )
		if( dValue <= m_vecdQuant[ i ] )
			break;

	return min( i, m_vecdQuant.size( ) - 1 ); }

bool CDataPair::IsContinuous( ) const {

	return m_fContinuous; }

void CDataPairImpl::Reset( bool fContinuous ) {

	m_vecdQuant.clear( );
	m_fContinuous = fContinuous; }

unsigned char CDataPair::GetValues( ) const {

	return m_vecdQuant.size( ); }

}

#include "stdafx.h"
#include "datapair.h"
#include "meta.h"

namespace libBioUtils {

const char	CPairImpl::c_szQuantExt[]	= ".quant";

bool CPairImpl::Open( const char* szDatafile, ifstream& ifsm ) {
	string		strToken;
	const char*	pc;

	strToken = szDatafile;
	strToken += c_szQuantExt;
	ifsm.open( strToken.c_str( ) );
	if( !ifsm.is_open( ) ) {
		if( !( pc = strrchr( szDatafile, '.' ) ) )
			return false;
		strToken = szDatafile;
		strToken.resize( pc - szDatafile );
		strToken += c_szQuantExt;
		ifsm.clear( );
		ifsm.open( strToken.c_str( ) );
		if( !ifsm.is_open( ) )
			return false; }

	return true; }

bool CPairImpl::Open( const char* szLine, vector<float>& vecdQuant ) {
	vector<string>	vecstrQuant;
	size_t			i;

	CMeta::Tokenize( szLine, vecstrQuant, CMeta::c_szWS, true );
	vecdQuant.resize( vecstrQuant.size( ) );
	for( i = 0; i < vecdQuant.size( ); ++i )
		vecdQuant[ i ] = (float)atof( vecstrQuant[ i ].c_str( ) );

	return true; }

size_t CPairImpl::Quantify( float dValue, const vector<float>& vecdQuant ) {
	size_t	i;

	if( CMeta::IsNaN( dValue ) )
		return -1;

	for( i = 0; i < vecdQuant.size( ); ++i )
		if( dValue <= vecdQuant[ i ] )
			break;

	return min( i, vecdQuant.size( ) - 1 ); }

bool CDataPair::Open( const CSlim& Slim ) {

	Reset( false );
	return CDat::Open( Slim ); }

bool CDataPair::Open( const char* szDatafile, bool fContinuous, bool fMemmap ) {
	static const size_t	c_iBuf	= 8192;
	char		szBuf[ c_iBuf ];
	ifstream	ifsm;

	g_CatBioUtils.notice( "CDataPair::Open( %s, %d )", szDatafile, fContinuous );

	Reset( fContinuous );
	if( !CDat::Open( szDatafile, fMemmap ) )
		return false;
	if( m_fContinuous )
		return true;

	if( !CPairImpl::Open( szDatafile, ifsm ) )
		return false;
	ifsm.getline( szBuf, c_iBuf - 1 );
	ifsm.close( );
	return CPairImpl::Open( szBuf, m_vecdQuant ); }

size_t CDataPair::Quantify( float dValue ) const {

	return CPairImpl::Quantify( dValue, m_vecdQuant ); }

bool CDataPair::IsContinuous( ) const {

	return m_fContinuous; }

void CDataPairImpl::Reset( bool fContinuous ) {

	m_vecdQuant.clear( );
	m_fContinuous = fContinuous; }

unsigned char CDataPair::GetValues( ) const {

	return m_vecdQuant.size( ); }

void CDataPair::SetQuants( const float* adQuants, size_t iQuants ) {

	Reset( false );
	m_vecdQuant.resize( iQuants );
	copy( adQuants, adQuants + iQuants, m_vecdQuant.begin( ) ); }

void CDataPair::SetQuants( const vector<float>& vecdQuant ) {

	Reset( false );
	m_vecdQuant.resize( vecdQuant.size( ) );
	copy( vecdQuant.begin( ), vecdQuant.end( ), m_vecdQuant.begin( ) ); }

bool CPCLPair::Open( const char* szDatafile, size_t iSkip ) {
	static const size_t	c_iBuf	= 8192;
	char		szBuf[ c_iBuf ];
	ifstream	ifsm;
	size_t		i;

	g_CatBioUtils.notice( "CPCLPair::Open( %s )", szDatafile );

	ifsm.open( szDatafile );
	if( !CPCL::Open( ifsm, iSkip ) )
		return false;
	ifsm.close( );

	ifsm.clear( );
	if( !CPairImpl::Open( szDatafile, ifsm ) )
		return false;

	m_vecvecdQuants.resize( GetExperiments( ) );
	for( i = 0; i < m_vecvecdQuants.size( ); ++i ) {
		ifsm.getline( szBuf, c_iBuf - 1 );
		if( !CPairImpl::Open( szBuf, m_vecvecdQuants[ i ] ) )
			return false; }

	return true; }

size_t CPCLPair::Quantify( float dValue, size_t iExp ) const {

	return CPairImpl::Quantify( dValue, m_vecvecdQuants[ iExp ] ); }

}

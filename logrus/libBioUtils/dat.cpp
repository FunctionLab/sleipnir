#include "stdafx.h"
#include "dat.h"
#include "annotation.h"
#include "genome.h"
#include "pstdint.h"

namespace libBioUtils {

const char	CDatImpl::c_szBin[]	= "dab";

size_t CDatImpl::MapGene( TMapStrI& mapGenes, TVecStr& vecGenes, const string& strToken ) {
	TMapStrI::iterator	iterGenes;
	size_t				iRet;

	if( ( iterGenes = mapGenes.find( strToken ) ) == mapGenes.end( ) ) {
		iRet = mapGenes.size( );
		mapGenes[ strToken ] = iRet;
		vecGenes.push_back( strToken ); }
	else
		iRet = iterGenes->second;

	return iRet; }

void CDatImpl::ResizeNaN( TAF& vecf, size_t iLast ) {
	size_t	i;

	if( ( i = vecf.size( ) ) > iLast )
		return;

	vecf.resize( iLast + 1 );
	for( ; i < vecf.size( ); ++i )
		vecf[ i ] = CMeta::GetNaN( ); }

string CDatImpl::DabGene( istream& istm ) {
	char	szBuffer[ c_iBufferSize ];
	size_t	i;

	i = 0;
	do
		istm.seekg( 1, ios_base::cur );
	while( szBuffer[ i++ ] = istm.get( ) );

	return szBuffer; }

CDatImpl::CDatImpl( ) { }

CDatImpl::~CDatImpl( ) {

	Reset( ); }

void CDatImpl::Reset( ) {

	m_Data.Reset( );
	m_vecstrGenes.clear( ); }

void CDatImpl::SlimCache( const CSlim& Slim, vector<vector<size_t> >& vecveciGenes ) const {
	size_t	iS, iG;

	vecveciGenes.resize( Slim.GetSlims( ) );
	for( iS = 0; iS < vecveciGenes.size( ); ++iS ) {
		vecveciGenes[ iS ].resize( Slim.GetGenes( iS ) );
		for( iG = 0; iG < vecveciGenes[ iS ].size( ); ++iG )
			vecveciGenes[ iS ][ iG ] =  GetGene( Slim.GetGene( iS, iG ).GetName( ) ); } }

bool CDat::Open( const CSlim& Slim ) {
	vector<string>			vecstrGenes;
	size_t					iS1, iS2, iG1, iG2, iGene1, iGene2;
	vector<vector<size_t> >	vecveciGenes;

	Reset( );
	Slim.GetGeneNames( vecstrGenes );
	if( !Open( vecstrGenes ) )
		return false;

	SlimCache( Slim, vecveciGenes );
	for( iS1 = 0; iS1 < Slim.GetSlims( ); ++iS1 ) {
		g_CatBioUtils.info( "CDat::Open( ) processing slim: %s",
			Slim.GetSlim( iS1 ).c_str( ) );
		for( iG1 = 0; iG1 < Slim.GetGenes( iS1 ); ++iG1 ) {
			iGene1 = vecveciGenes[ iS1 ][ iG1 ];
			for( iG2 = ( iG1 + 1 ); iG2 < Slim.GetGenes( iS1 ); ++iG2 )
				Set( iGene1, vecveciGenes[ iS1 ][ iG2 ], 1 );
			for( iS2 = ( iS1 + 1 ); iS2 < Slim.GetSlims( ); ++iS2 )
				for( iG2 = 0; iG2 < Slim.GetGenes( iS2 ); ++iG2 ) {
					iGene2 = vecveciGenes[ iS2 ][ iG2 ];
					if( CMeta::IsNaN( Get( iGene1, iGene2 ) ) )
						Set( iGene1, iGene2, 0 ); } } }

	return true; }

bool CDat::Open( const CSlim& SlimPos, const CSlim& SlimNeg ) {
	set<string>				setstrGenes;
	vector<string>			vecstrGenes;
	size_t					iS1, iS2, iG1, iG2, iGene1, iGene2;
	vector<vector<size_t> >	vecveciGenes;

	Reset( );
	SlimPos.GetGeneNames( vecstrGenes );
	setstrGenes.insert( vecstrGenes.begin( ), vecstrGenes.end( ) );
	vecstrGenes.clear( );
	SlimNeg.GetGeneNames( vecstrGenes );
	setstrGenes.insert( vecstrGenes.begin( ), vecstrGenes.end( ) );
	vecstrGenes.clear( );
	vecstrGenes.resize( setstrGenes.size( ) );
	copy( setstrGenes.begin( ), setstrGenes.end( ), vecstrGenes.begin( ) );
	if( !Open( vecstrGenes ) )
		return false;

	SlimCache( SlimPos, vecveciGenes );
	for( iS1 = 0; iS1 < SlimPos.GetSlims( ); ++iS1 ) {
		g_CatBioUtils.info( "CDat::Open( ) processing slim: %s",
			SlimPos.GetSlim( iS1 ).c_str( ) );
		for( iG1 = 0; iG1 < SlimPos.GetGenes( iS1 ); ++iG1 ) {
			iGene1 = vecveciGenes[ iS1 ][ iG1 ];
			for( iG2 = ( iG1 + 1 ); iG2 < SlimPos.GetGenes( iS1 ); ++iG2 )
				Set( iGene1, vecveciGenes[ iS1 ][ iG2 ], 1 ); } }
	vecveciGenes.clear( );
	SlimCache( SlimNeg, vecveciGenes );
	for( iS1 = 0; iS1 < SlimNeg.GetSlims( ); ++iS1 ) {
		g_CatBioUtils.info( "CDat::Open( ) processing slim: %s",
			SlimNeg.GetSlim( iS1 ).c_str( ) );
		for( iG1 = 0; iG1 < SlimNeg.GetGenes( iS1 ); ++iG1 ) {
			iGene1 = vecveciGenes[ iS1 ][ iG1 ];
			for( iS2 = ( iS1 + 1 ); iS2 < SlimNeg.GetSlims( ); ++iS2 )
				for( iG2 = 0; iG2 < SlimNeg.GetGenes( iS2 ); ++iG2 ) {
					iGene2 = vecveciGenes[ iS2 ][ iG2 ];
					if( CMeta::IsNaN( Get( iGene1, iGene2 ) ) )
						Set( iGene1, iGene2, 0 ); } } }

	return true; }

bool CDat::Open( istream& istm, bool fBinary ) {

	return ( fBinary ? OpenBinary( istm ) : OpenText( istm ) ); }

bool CDatImpl::OpenText( istream& istm ) {
	char		acBuf[ c_iBufferSize ];
	const char*	pc;
	string		strToken, strCache;
	TMapStrI	mapGenes;
	size_t		iOne, iTwo, i, iToken;
	float		dScore;
	TAAF		vecvecfScores;

	Reset( );
	while( istm.peek( ) != EOF ) {
		istm.getline( acBuf, ARRAYSIZE(acBuf) - 1 );
		for( iToken = 0; iToken < 3; ++iToken ) {
			strToken = OpenToken( acBuf, &pc );
			if( !strToken.length( ) )
				break;
			if( strToken != strCache ) {
				strCache = strToken;
				iOne = MapGene( mapGenes, m_vecstrGenes, strToken ); }

			strToken = OpenToken( pc, &pc );
			iTwo = MapGene( mapGenes, m_vecstrGenes, strToken );
			if( !( dScore = (float)atof( OpenToken( pc ).c_str( ) ) ) && errno && ( errno != ENOENT ) ) {
				Reset( );
				return false; }

			i = ( ( iOne > iTwo ) ? iOne : iTwo );
			if( vecvecfScores.size( ) <= i )
				vecvecfScores.resize( i + 1 );
			ResizeNaN( vecvecfScores[ iOne ], i );
			ResizeNaN( vecvecfScores[ iTwo ], i );
			if( !CMeta::IsNaN( vecvecfScores[ iOne ][ iTwo ] ) && ( vecvecfScores[ iOne ][ iTwo ] != dScore ) ) {
				g_CatBioUtils.error( "CDatImpl::OpenText( ) duplicate genes %s, %s (%g:%g)",
					strCache.c_str( ), strToken.c_str( ), vecvecfScores[ iOne ][ iTwo ],
					dScore );
//dScore = CMeta::GetNaN( ); }
				Reset( );
				return false; }
			vecvecfScores[ iOne ][ iTwo ] = vecvecfScores[ iTwo ][ iOne ] = dScore; } }

	m_Data.Initialize( m_vecstrGenes.size( ) );
	for( iOne = 0; iOne < m_vecstrGenes.size( ); ++iOne )
		for( iTwo = ( iOne + 1 ); iTwo < m_vecstrGenes.size( ); ++iTwo )
			Set( iOne, iTwo, ( ( iTwo < vecvecfScores[ iOne ].size( ) ) ?
				vecvecfScores[ iOne ][ iTwo ] : CMeta::GetNaN( ) ) );

	return true; }

bool CDatImpl::OpenBinary( istream& istm ) {
	size_t	i;
	float*	adScores;

	if( !OpenGenes( istm, true ) )
		return false;
	m_Data.Initialize( GetGenes( ) );
	adScores = new float[ GetGenes( ) - 1 ];
	for( i = 0; ( i + 1 ) < GetGenes( ); ++i ) {
		istm.read( (char*)adScores, sizeof(*adScores) * ( GetGenes( ) - i - 1 ) );
		Set( i, adScores ); }
	delete[] adScores;

	return true; }

bool CDat::OpenGenes( istream& istm, bool fBinary ) {

	return CDatImpl::OpenGenes( istm, fBinary ); }

bool CDatImpl::OpenGenes( istream& istm, bool fBinary ) {
	size_t		i, iToken;
	uint32_t	iCount;
	string		strToken, strCache;
	float		d;
	char		acBuf[ c_iBufferSize ];
	const char*	pc;

	Reset( );
	if( fBinary ) {
		istm.read( (char*)&iCount, sizeof(iCount) );
		if( iCount > c_iGeneLimit )
			return false;
		m_vecstrGenes.resize( iCount );
		for( i = 0; i < iCount; ++i )
			m_vecstrGenes[ i ] = DabGene( istm ); }
	else {
		set<string>					setstrGenes;
		set<string>::const_iterator	iterGenes;

		while( istm.peek( ) != EOF ) {
			istm.getline( acBuf, ARRAYSIZE(acBuf) - 1 );
			for( iToken = 0; iToken < 3; ++iToken ) {
				strToken = OpenToken( acBuf, &pc );
				if( !strToken.length( ) )
					break;
				if( strToken != strCache ) {
					strCache = strToken;
					setstrGenes.insert( strToken ); }

				strToken = OpenToken( pc, &pc );
				setstrGenes.insert( strToken );
				d = (float)atof( OpenToken( pc ).c_str( ) ); } }
		m_vecstrGenes.reserve( setstrGenes.size( ) );
		for( iterGenes = setstrGenes.begin( ); iterGenes != setstrGenes.end( ); ++iterGenes )
			m_vecstrGenes.push_back( *iterGenes ); }

	return true; }

void CDat::Save( ostream& ostm, bool fBinary ) const {

	fBinary ? SaveBinary( ostm ) : SaveText( ostm ); }

void CDatImpl::SaveText( ostream& ostm ) const {
	size_t	i, j;
	float	d;

	for( i = 0; i < m_vecstrGenes.size( ); ++i )
		for( j = ( i + 1 ); j < m_vecstrGenes.size( ); ++j )
			if( ( d = Get( i, j ) ) != CMeta::GetNaN( ) )
				ostm << m_vecstrGenes[ i ] << '\t' << m_vecstrGenes[ j ] << '\t' << d << endl; }

void CDatImpl::SaveBinary( ostream& ostm ) const {
	size_t			i, j;
	uint32_t		iSize;
	const float*	pd;

	iSize = m_vecstrGenes.size( );
	ostm.write( (char*)&iSize, sizeof(iSize) );
	for( i = 0; i < m_vecstrGenes.size( ); ++i ) {
		for( j = 0; j < m_vecstrGenes[ i ].length( ); ++j ) {
			ostm.put( 0 );
			ostm.put( m_vecstrGenes[ i ][ j ] ); }
		ostm.put( 0 );
		ostm.put( 0 ); }
	for( i = 0; ( i + 1 ) < GetGenes( ); ++i ) {
		pd = m_Data.Get( i );
		ostm.write( (const char*)pd, sizeof(*pd) * ( GetGenes( ) - i - 1 ) ); } }

bool CDat::Open( const vector<string>& vecstrGenes ) {
	size_t	i, j;

	Reset( );
	m_vecstrGenes.resize( vecstrGenes.size( ) );
	for( i = 0; i < vecstrGenes.size( ); ++i )
		m_vecstrGenes[ i ] = vecstrGenes[ i ];

	m_Data.Initialize( m_vecstrGenes.size( ) );
	for( i = 0; i < m_vecstrGenes.size( ); ++i )
		for( j = ( i + 1 ); j < m_vecstrGenes.size( ); ++j )
			Set( i, j, CMeta::GetNaN( ) );

	return true; }

bool CDat::Open( const vector<string>& vecstrGenes, const CDistanceMatrix& Dist ) {
	size_t	i;

	Reset( );
	m_vecstrGenes.reserve( vecstrGenes.size( ) );
	for( i = 0; i < vecstrGenes.size( ); ++i )
		m_vecstrGenes.push_back( vecstrGenes[ i ] );

	m_Data.Initialize( m_vecstrGenes.size( ), Dist );

	return true; }

size_t CDat::GetGene( const string& strGene ) const {

	return CDatImpl::GetGene( strGene ); }

size_t CDatImpl::GetGene( const string& strGene ) const {
	size_t	i;

	for( i = 0; i < m_vecstrGenes.size( ); ++i )
		if( m_vecstrGenes[ i ] == strGene )
			return i;

	return -1; }

bool CDat::Open( const char* szFile ) {
	ifstream	ifsm;
	bool		fBinary;

	fBinary = !strcmp( szFile + strlen( szFile ) - strlen( c_szBin ), c_szBin );
	ifsm.open( szFile, ( fBinary ? ios_base::binary : ios_base::in ) );
	if( !ifsm.is_open( ) )
		return false;
	return Open( ifsm, fBinary ); }

const vector<string>& CDat::GetGeneNames( ) const {

	return m_vecstrGenes; }

void CDat::Normalize( bool fCap ) {

	if( fCap )
		NormalizeMinmax( );
	else
		NormalizeStdev( ); }

void CDatImpl::NormalizeStdev( ) {
	double	d, dAve, dDev;
	size_t	i, j, iN;

	dAve = dDev = 0;
	for( iN = i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( (float)( d = Get( i, j ) ) ) ) {
				iN++;
				dAve += d;
				dDev += d * d; }
	dAve /= iN;
	dDev = sqrt( ( dDev / iN ) - ( dAve * dAve ) );
	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( (float)( d = Get( i, j ) ) ) )
				Set( i, j, (float)( ( d - dAve ) / dDev ) ); }

void CDatImpl::NormalizeMinmax( ) {
	float	d, dMin, dMax;
	size_t	i, j;

	dMax = -( dMin = FLT_MAX );
	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = Get( i, j ) ) ) {
				if( d < dMin )
					dMin = d;
				if( d > dMax )
					dMax = d; }
	if( dMax -= dMin )
		for( i = 0; i < GetGenes( ); ++i )
			for( j = ( i + 1 ); j < GetGenes( ); ++j )
				Set( i, j, ( Get( i, j ) - dMin ) / dMax ); }

void CDat::Invert( ) {
	size_t	i, j;
	float	d;

	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = Get( i, j ) ) )
				Set( i, j, 1 - d ); }

bool CDat::FilterGenes( const char* szGenes, EFilter eFilt ) {
	CGenome		Genome;
	CGenes		Genes( Genome );
	ifstream	ifsm;

	if( !szGenes )
		return false;

	ifsm.open( szGenes );
	if( !Genes.Open( ifsm ) )
		return false;
	FilterGenes( Genes, eFilt );
	ifsm.close( );
	return true; }

void CDat::FilterGenes( const CGenes& Genes, EFilter eFilt ) {
	size_t			i, j;
	vector<bool>	vecfGenes;

	vecfGenes.resize( GetGenes( ) );
	for( i = 0; i < Genes.GetGenes( ); ++i )
		if( ( j = GetGene( Genes.GetGene( i ).GetName( ) ) ) != -1 )
			vecfGenes[ j ] = true;

	for( i = 0; i < GetGenes( ); ++i ) {
		if( vecfGenes[ i ] ) {
			if( eFilt == EFilterInclude )
				continue;
			if( eFilt == EFilterExclude ) {
				for( j = ( i + 1 ); j < GetGenes( ); ++j )
					Set( i, j, CMeta::GetNaN( ) );
				continue; } }
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			switch( eFilt ) {
				case EFilterInclude:
					if( !vecfGenes[ j ] )
						Set( i, j, CMeta::GetNaN( ) );
					break;

				case EFilterTerm:
					if( !( vecfGenes[ i ] && vecfGenes[ j ] ) &&
						( !( vecfGenes[ i ] || vecfGenes[ j ] ) || Get( i, j ) ) )
							Set( i, j, CMeta::GetNaN( ) );
					break;

				case EFilterExclude:
					if( vecfGenes[ j ] )
						Set( i, j, CMeta::GetNaN( ) );
					break; } } }

}

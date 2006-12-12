#include "stdafx.h"
#include "dat.h"
#include "annotation.h"
#include "genome.h"
#include "pstdint.h"
#include "statistics.h"

namespace libBioUtils {

const char	CDatImpl::c_szBin[]	= "dab";
const char	CDatImpl::c_szPcl[]	= "pcl";

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

CDatImpl::~CDatImpl( ) {

	Reset( ); }

void CDatImpl::Reset( ) {

	m_Data.Reset( );
	m_vecstrGenes.clear( );
	if( m_pPCL && m_fPCLMemory )
		delete m_pPCL;
	m_pPCL = NULL;
	if( m_pMeasure && m_fMeasureMemory )
		delete m_pMeasure;
	m_pMeasure = NULL; }

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

bool CDat::Open( const vector<CGenes*>& vecpPositives, const CDat& DatNegatives,
	const CGenome& Genome ) {
	size_t			i, j, iOne, iTwo;
	vector<size_t>	veciGenes;

	Reset( );
	Open( Genome.GetGeneNames( ) );
	for( i = 0; i < vecpPositives.size( ); ++i )
		OpenHelper( vecpPositives[ i ], 1 );
	veciGenes.resize( DatNegatives.GetGenes( ) );
	for( i = 0; i < veciGenes.size( ); ++i )
		veciGenes[ i ] = GetGene( DatNegatives.GetGene( i ) );
	for( i = 0; i < DatNegatives.GetGenes( ); ++i ) {
		iOne = veciGenes[ i ];
		for( j = ( i + 1 ); j < DatNegatives.GetGenes( ); ++j )
			if( CMeta::IsNaN( Get( iOne, iTwo = veciGenes[ j ] ) ) )
				Set( iOne, iTwo, DatNegatives.Get( i, j ) ); }

	return true; }

bool CDat::Open( const vector<CGenes*>& vecpPositives, const vector<CGenes*>& vecpNegatives,
	float dP, const CGenome& Genome ) {
	size_t			i, j, k, iOne, iTwo, iOverlap;
	float			d;
	const CGenes*	pBig;
	const CGenes*	pSmall;

	Reset( );
	Open( Genome.GetGeneNames( ) );
	for( i = 0; i < vecpPositives.size( ); ++i )
		OpenHelper( vecpPositives[ i ], 1 );
	for( i = 0; i < vecpNegatives.size( ); ++i )
		OpenHelper( vecpNegatives[ i ], 0 );
	if( dP > 0 )
		for( i = 0; i < vecpPositives.size( ); ++i ) {
			iOne = vecpPositives[ i ]->GetGenes( );
			for( j = ( i + 1 ); j < vecpPositives.size( ); ++j ) {
				iTwo = vecpPositives[ j ]->GetGenes( );
				if( iOne < iTwo ) {
					pSmall = vecpPositives[ i ];
					pBig = vecpPositives[ j ]; }
				else {
					pSmall = vecpPositives[ j ];
					pBig = vecpPositives[ i ]; }
				for( iOverlap = k = 0; k < pSmall->GetGenes( ); ++k )
					if( pBig->IsGene( pSmall->GetGene( k ).GetName( ) ) )
						iOverlap++;
				if( CStatistics::HypergeometricCDF( iOverlap, iOne, iTwo, GetGenes( ) ) < dP )
					OpenHelper( pBig, pSmall, 0 ); } }
	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( CMeta::IsNaN( d = Get( i, j ) ) )
				Set( i, j, 0 );
			else if( !d )
				Set( i, j, CMeta::GetNaN( ) );

	return true; }

void CDatImpl::OpenHelper( const CGenes* pGenes, float dValue ) {
	vector<size_t>	veciGenes;
	size_t			i, j, iOne, iTwo;

	veciGenes.resize( pGenes->GetGenes( ) );
	for( i = 0; i < veciGenes.size( ); ++i )
		veciGenes[ i ] = GetGene( pGenes->GetGene( i ).GetName( ) );
	for( i = 0; i < veciGenes.size( ); ++i ) {
		iOne = veciGenes[ i ];
		for( j = ( i + 1 ); j < veciGenes.size( ); ++j ) {
			iTwo = veciGenes[ j ];
			if( CMeta::IsNaN( Get( iOne, iTwo ) ) )
				Set( iOne, iTwo, dValue ); } } }

void CDatImpl::OpenHelper( const CGenes* pOne, const CGenes* pTwo, float dValue ) {
	vector<size_t>	veciOne, veciTwo;
	size_t			i, j, iOne, iTwo;

	veciOne.resize( pOne->GetGenes( ) );
	for( i = 0; i < veciOne.size( ); ++i )
		veciOne[ i ] = GetGene( pOne->GetGene( i ).GetName( ) );
	veciTwo.resize( pTwo->GetGenes( ) );
	for( i = 0; i < veciTwo.size( ); ++i )
		veciTwo[ i ] = GetGene( pTwo->GetGene( i ).GetName( ) );
	for( i = 0; i < veciOne.size( ); ++i ) {
		iOne = veciOne[ i ];
		for( j = 0; j < veciTwo.size( ); ++j ) {
			iTwo = veciTwo[ j ];
			if( CMeta::IsNaN( Get( iOne, iTwo ) ) )
				Set( iOne, iTwo, dValue ); } } }

bool CDat::Open( istream& istm, bool fBinary, bool fPCL, float dDefault ) {

	return ( fPCL ? OpenPCL( istm ) : ( fBinary ? OpenBinary( istm ) : OpenText( istm, dDefault ) ) ); }

bool CDatImpl::OpenPCL( istream& istm ) {
	size_t	i, iOne, iTwo;
	float	d, dAve, dStd;

	Reset( );
	m_pPCL = new CPCL( );
	if( !m_pPCL->Open( istm ) )
		return false;

	dAve = dStd = 0;
	m_pMeasure = new CMeasurePearNorm( );
	for( i = 0; i < c_iApproximate; ++i ) {
		iOne = rand( ) % GetGenes( );
		if( ( iTwo = rand( ) % GetGenes( ) ) == iOne ) {
			i--;
			continue; }
		dAve += ( d = Get( iOne, iTwo ) );
		dStd += d * d; }
	dAve /= i;
	dStd = sqrt( ( dStd / i ) - ( dAve * dAve ) );
	delete m_pMeasure;
	m_pMeasure = new CMeasurePearNorm( dAve, dStd );
	m_fMeasureMemory = true;

	return true; }

bool CDat::Open( const CPCL& PCL, const IMeasure* pMeasure, bool fMeasureMemory ) {

	Reset( );
	m_pPCL = (CPCL*)&PCL;
	m_fPCLMemory = false;
	m_pMeasure = pMeasure;
	m_fMeasureMemory = fMeasureMemory;

	return true; }

bool CDatImpl::OpenText( istream& istm, float dDefault ) {
	char		acBuf[ c_iBufferSize ];
	const char*	pc;
	char*		pcTail;
	string		strToken, strCache, strValue;
	TMapStrI	mapGenes;
	size_t		iOne, iTwo, i;
	float		dScore;
	TAAF		vecvecfScores;

	Reset( );
	while( istm.peek( ) != EOF ) {
		istm.getline( acBuf, ARRAYSIZE(acBuf) - 1 );
		strToken = OpenToken( acBuf, &pc );
		if( !strToken.length( ) )
			break;
		if( strToken != strCache ) {
			strCache = strToken;
			iOne = MapGene( mapGenes, m_vecstrGenes, strToken ); }

		strToken = OpenToken( pc, &pc );
		if( !strToken.length( ) ) {
			Reset( );
			return false; }
		iTwo = MapGene( mapGenes, m_vecstrGenes, strToken );
		strValue = OpenToken( pc );
		if( !strValue.length( ) ) {
			if( CMeta::IsNaN( dScore = dDefault ) ) {
				Reset( );
				return false; } }
		else if( !( dScore = (float)strtod( strValue.c_str( ), &pcTail ) ) &&
			( pcTail != ( strValue.c_str( ) + strValue.length( ) ) ) ) {
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
			Reset( );
			return false; }
		vecvecfScores[ iOne ][ iTwo ] = vecvecfScores[ iTwo ][ iOne ] = dScore; }

	m_Data.Initialize( GetGenes( ) );
	for( iOne = 0; iOne < GetGenes( ); ++iOne )
		for( iTwo = ( iOne + 1 ); iTwo < GetGenes( ); ++iTwo )
			Set( iOne, iTwo, ( ( iTwo < vecvecfScores[ iOne ].size( ) ) ?
				vecvecfScores[ iOne ][ iTwo ] : CMeta::GetNaN( ) ) );

	return true; }

bool CDatImpl::OpenBinary( istream& istm ) {
	size_t	i;
	float*	adScores;

	if( !OpenGenes( istm, true, false ) )
		return false;
	m_Data.Initialize( GetGenes( ) );
	adScores = new float[ GetGenes( ) - 1 ];
	for( i = 0; ( i + 1 ) < GetGenes( ); ++i ) {
		istm.read( (char*)adScores, sizeof(*adScores) * ( GetGenes( ) - i - 1 ) );
		Set( i, adScores ); }
	delete[] adScores;

	return true; }

bool CDat::OpenGenes( istream& istm, bool fBinary, bool fPCL ) {

	return CDatImpl::OpenGenes( istm, fBinary, fPCL ); }

bool CDatImpl::OpenGenes( istream& istm, bool fBinary, bool fPCL ) {
	size_t		i, iToken;
	uint32_t	iCount;
	string		strToken, strCache;
	float		d;
	char		acBuf[ c_iBufferSize ];
	const char*	pc;

	Reset( );
	if( fPCL ) {
		m_pPCL = new CPCL( );
		if( m_pPCL->Open( istm ) ) {
			m_pMeasure = (IMeasure*)1;
			m_fMeasureMemory = false;
			return true; }
		return false; }
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
				d = (float)strtod( ( strToken = OpenToken( pc ) ).c_str( ), (char**)&pc );
				if( !d && ( ( pc - strToken.c_str( ) ) != strToken.length( ) ) )
					return false; } }
		m_vecstrGenes.reserve( setstrGenes.size( ) );
		for( iterGenes = setstrGenes.begin( ); iterGenes != setstrGenes.end( ); ++iterGenes )
			m_vecstrGenes.push_back( *iterGenes ); }

	return true; }

void CDat::Save( ostream& ostm, bool fBinary ) const {

	fBinary ? SaveBinary( ostm ) : SaveText( ostm ); }

void CDatImpl::SaveText( ostream& ostm ) const {
	size_t	i, j;
	float	d;

	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = Get( i, j ) ) )
				ostm << GetGene( i ) << '\t' << GetGene( j ) << '\t' << d << endl; }

void CDatImpl::SaveBinary( ostream& ostm ) const {
	size_t			i, j;
	uint32_t		iSize;
	const float*	pd;
	string			strGene;
	float			d;

	iSize = GetGenes( );
	ostm.write( (char*)&iSize, sizeof(iSize) );
	for( i = 0; i < iSize; ++i ) {
		strGene = GetGene( i );
		for( j = 0; j < strGene.length( ); ++j ) {
			ostm.put( 0 );
			ostm.put( strGene[ j ] ); }
		ostm.put( 0 );
		ostm.put( 0 ); }
	if( m_pMeasure ) {
		for( i = 0; i < iSize; ++i )
			for( j = ( i + 1 ); j < iSize; ++j ) {
				d = Get( i, j );
				ostm.write( (char*)&d, sizeof(d) ); } }
	else
		for( i = 0; ( i + 1 ) < iSize; ++i ) {
			pd = m_Data.Get( i );
			ostm.write( (char*)pd, sizeof(*pd) * ( iSize - i - 1 ) ); } }

bool CDat::Open( const vector<string>& vecstrGenes, bool fInitialize ) {
	size_t	i, j;

	Reset( );
	m_vecstrGenes.resize( vecstrGenes.size( ) );
	for( i = 0; i < vecstrGenes.size( ); ++i )
		m_vecstrGenes[ i ] = vecstrGenes[ i ];

	m_Data.Initialize( m_vecstrGenes.size( ) );
	if( fInitialize )
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

size_t CDatImpl::GetGene( const string& strGene ) const {
	size_t	i;

	if( m_pMeasure )
		return m_pPCL->GetGene( strGene );

	for( i = 0; i < GetGenes( ); ++i )
		if( m_vecstrGenes[ i ] == strGene )
			return i;

	return -1; }

bool CDat::Open( const char* szFile ) {
	ifstream	ifsm;
	bool		fBinary, fPCL;

	fBinary = !strcmp( szFile + strlen( szFile ) - strlen( c_szBin ), c_szBin );
	fPCL = !strcmp( szFile + strlen( szFile ) - strlen( c_szBin ), c_szPcl );
	ifsm.open( szFile, ( fBinary ? ios_base::binary : ios_base::in ) );
	if( !ifsm.is_open( ) )
		return false;
	return Open( ifsm, fBinary, fPCL ); }

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

void CDat::SaveDOT( ostream& ostm, float dCutoff, const CGenome* pGenome, bool fMinimal ) const {
	size_t			i, j;
	float			d;
	bool			fAll;
	vector<string>	vecstrNames;
	vector<bool>	vecfGenes;

	fAll = CMeta::IsNaN( dCutoff );
	ostm << "graph G {" << endl;
	ostm << "	pack = \"true\";" << endl;
	ostm << "	overlap = \"orthoyx\";" << endl;
	ostm << "	splines = \"true\";" << endl;

	if( !( fMinimal || fAll ) ) {
		vecfGenes.resize( GetGenes( ) );
		for( i = 0; i < vecfGenes.size( ); ++i )
			for( j = ( i + 1 ); j < vecfGenes.size( ); ++j )
				if( !CMeta::IsNaN( d = Get( i, j ) ) && ( fAll || ( d >= dCutoff ) ) )
					vecfGenes[ i ] = vecfGenes[ j ] = true; }

	vecstrNames.resize( GetGenes( ) );
	for( i = 0; i < vecstrNames.size( ); ++i ) {
		vecstrNames[ i ] = CMeta::Filename( GetGene( i ) );
		if( pGenome && ( ( j = pGenome->GetGene( GetGene( i ) ) ) != -1 ) ) {
			const CGene&	Gene	= pGenome->GetGene( j );

			if( Gene.GetSynonyms( ) ) {
				if( fMinimal ) {
					vecstrNames[ i ] += '_';
					vecstrNames[ i ] += Gene.GetSynonym( 0 ); }
				else if( fAll || vecfGenes[ i ] )
					ostm << vecstrNames[ i ] << " [label=\"" << Gene.GetSynonym( 0 ) << "\"];" << endl; } } }

	ostm << endl;
	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = Get( i, j ) ) && ( fAll || ( d >= dCutoff ) ) )
				ostm << vecstrNames[ i ] << " -- " << vecstrNames[ j ] << ';' << endl;

	ostm << "}" << endl; }

void CDat::SaveGDF( ostream& ostm, float dCutoff ) const {
	size_t			i, j;
	float			d;
	bool			fAll;
	vector<string>	vecstrNames;

	fAll = CMeta::IsNaN( dCutoff );
	ostm << "nodedef> name" << endl;

	vecstrNames.resize( GetGenes( ) );
	for( i = 0; i < vecstrNames.size( ); ++i )
		vecstrNames[ i ] = CMeta::Filename( GetGene( i ) );
	for( i = 0; i < GetGenes( ); ++i )
		ostm << vecstrNames[ i ] << endl;

	ostm << endl << "edgedef> node1,node2,weight" << endl;
	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = Get( i, j ) ) && ( fAll || ( d >= dCutoff ) ) )
				ostm << vecstrNames[ i ] << "," << vecstrNames[ j ] << "," << ( 10 * Get( i, j ) ) << endl; }

void CDat::SaveNET( ostream& ostm, float dCutoff ) const {
	size_t	i, j;
	float	d;
	bool	fAll;

	fAll = CMeta::IsNaN( dCutoff );
	ostm << "*Vertices " << GetGenes( ) << endl;
	for( i = 0; i < GetGenes( ); ++i )
		ostm << ( i + 1 ) << " \"" << GetGene( i ) << '"' << endl;

	ostm << "*Edges" << endl;
	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = Get( i, j ) ) && ( fAll || ( d >= dCutoff ) ) )
				ostm << ( i + 1 ) << ' ' << ( j + 1 ) << ' ' << d << endl; }

struct SSorterRank {
	const CDat&	m_Dat;

	SSorterRank( const CDat& Dat ) : m_Dat(Dat) { }

	bool operator()( const pair<size_t, size_t>& prOne, const pair<size_t, size_t>& prTwo ) const {

		return ( m_Dat.Get( prOne.first, prOne.second ) < m_Dat.Get( prTwo.first, prTwo.second ) ); }

};

void CDat::Rank( ) {
	vector<pair<size_t, size_t> >	vecprData;
	size_t							i, j, iRank;
	float							d, dPrev;

	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( Get( i, j ) ) )
				vecprData.push_back( pair<size_t, size_t>( i, j ) );
	sort( vecprData.begin( ), vecprData.end( ), SSorterRank( *this ) );
	for( iRank = i = 0; i < vecprData.size( ); ++i ) {
		d = Get( vecprData[ i ].first, vecprData[ i ].second );
		if( i && ( d != dPrev ) )
			iRank = i;
		dPrev = d;
		Set( vecprData[ i ].first, vecprData[ i ].second, (float)iRank ); } }

bool CDatMap::Open( const vector<string>& vecstrGenes, const char* szFile, bool fInitialize ) {
	size_t			i, j, iSize;
	unsigned char*	pb;

	Reset( );
	m_vecstrGenes.resize( vecstrGenes.size( ) );
	copy( vecstrGenes.begin( ), vecstrGenes.end( ), m_vecstrGenes.begin( ) );

	iSize = sizeof(uint32_t);
	for( i = 0; i < GetGenes( ); ++i )
		iSize += 2 * ( GetGene( i ).length( ) + 1 );
	iSize += CDistanceMatrix::GetSpace( GetGenes( ) );
	if( !CMeta::MapWrite( m_abData, m_hndlData, iSize, szFile ) )
		return false;
	*(uint32_t*)( pb = m_abData ) = GetGenes( );
	pb += sizeof(uint32_t);
	for( i = 0; i < GetGenes( ); ++i ) {
		const string&	strGene	= GetGene( i );

		for( j = 0; j < strGene.length( ); ++j ) {
			*pb++ = 0;
			*pb++ = strGene[ j ]; }
		*pb++ = *pb++ = 0; }

	m_aadData = new float*[ GetGenes( ) - 1 ];
	m_aadData[ 0 ] = (float*)pb;
	for( i = 1; ( i + 1 ) < m_vecstrGenes.size( ); ++i )
		m_aadData[ i ] = m_aadData[ i - 1 ] + GetGenes( ) - i;
	m_Data.Initialize( m_vecstrGenes.size( ), (float**)m_aadData );
	if( fInitialize )
		for( i = 0; i < m_vecstrGenes.size( ); ++i )
			for( j = ( i + 1 ); j < m_vecstrGenes.size( ); ++j )
				Set( i, j, CMeta::GetNaN( ) );

	return true; }

}

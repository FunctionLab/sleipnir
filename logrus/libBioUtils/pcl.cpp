#include "stdafx.h"
#include "pcl.h"
#include "meta.h"
#include "statistics.h"
#include "genome.h"
#include "measure.h"
#include "dat.h"

namespace libBioUtils {

const char	CPCLImpl::c_szEWEIGHT[]	= "EWEIGHT";
const char	CPCLImpl::c_szGENE[]	= "GENE";
const char	CPCLImpl::c_szGID[]		= "GID";
const char	CPCLImpl::c_szGWEIGHT[]	= "GWEIGHT";
const char	CPCLImpl::c_szNAME[]	= "NAME";
const char	CPCLImpl::c_szOne[]		= "1";

int CPCL::Distance( const char* szPCL, size_t iSkip, const char* szDistance, bool fNormalize, bool fZScore,
	bool fAutocorrelate, const char* szGenes, float dCutoff, size_t iLimit, CPCL& PCL, CDat& Dat ) {
	size_t						i, j, iOne, iTwo;
	float						d;
	ifstream					ifsm;
	vector<string>				vecstrGenes;
	CGenome						Genome;
	CGenes						GenesIn( Genome );
	vector<size_t>				veciGenes;
	const float*				adOne;
	IMeasure*					pMeasure;
	CMeasurePearson				Pearson;
	CMeasureEuclidean			Euclidean;
	CMeasureKendallsTau			KendallsTau;
	CMeasureKolmogorovSmirnov	KolmSmir;
	CMeasureSpearman			Spearman( true );
	CMeasurePearNorm			PearNorm;
	CMeasureHypergeometric		Hypergeom;
	CMeasureQuickPearson		PearQuick;
	CMeasureInnerProduct		InnerProd;
	CMeasureBinaryInnerProduct	BinInnerProd;

	if( szPCL ) {
		ifsm.open( szPCL );
		if( !PCL.Open( ifsm, iSkip ) ) {
			g_CatBioUtils.error( "CPCL::Distance( %s, %d, %s, %d, %d, %d, %s, %g ) failed to open PCL", szPCL, iSkip,
				szDistance, fNormalize, fZScore, fAutocorrelate, szGenes ? szGenes : "", dCutoff );
			return 1; }
		ifsm.close( ); }
	else if( !PCL.Open( cin, iSkip ) ) {
		g_CatBioUtils.error( "CPCL::Distance( %s, %d, %s, %d, %d, %d, %s, %g ) failed to open PCL", "stdin", iSkip,
			szDistance, fNormalize, fZScore, fAutocorrelate, szGenes ? szGenes : "", dCutoff );
		return 1; }

	CMeasureSigmoid				EuclideanSig( &Euclidean, false, 1.0f / PCL.GetExperiments( ) );
	IMeasure*					apMeasures[]	= { &Pearson, &EuclideanSig, &KendallsTau,
		&KolmSmir, &Spearman, &PearNorm, &Hypergeom, &PearQuick, &InnerProd, &BinInnerProd, NULL };

	pMeasure = NULL;
	for( i = 0; apMeasures[ i ]; ++i )
		if( !strcmp( apMeasures[ i ]->GetName( ), szDistance ) ) {
			pMeasure = apMeasures[ i ];
			break; }
	if( !pMeasure )
		return 1;

	CMeasureAutocorrelate		Autocorrelate( pMeasure, false );
	if( fAutocorrelate )
		pMeasure = &Autocorrelate;

	if( szGenes ) {
		ifsm.clear( );
		ifsm.open( szGenes );
		if( !GenesIn.Open( ifsm ) ) {
			g_CatBioUtils.error( "CPCL::Distance( %s, %d, %s, %d, %d, %d, %s, %g ) failed to open genes", szPCL ? szPCL :
				"stdin", iSkip, szDistance, fNormalize, fZScore, fAutocorrelate, szGenes, dCutoff );
			return 1; }
		ifsm.close( ); }
	else
		GenesIn.Open( PCL.GetGeneNames( ) );
	veciGenes.resize( GenesIn.GetGenes( ) );
	for( i = 0; i < veciGenes.size( ); ++i )
		veciGenes[ i ] = szGenes ? PCL.GetGene( GenesIn.GetGene( i ).GetName( ) ) : i;

	if( pMeasure->IsRank( ) )
		PCL.RankTransform( );

	if( ( iLimit != -1 ) && ( PCL.GetGenes( ) > iLimit ) )
		Dat.Open( PCL, pMeasure->Clone( ), true );
	else {
		Dat.Open( GenesIn.GetGeneNames( ) );
		for( i = 0; i < Dat.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
				Dat.Set( i, j, CMeta::GetNaN( ) );
		for( i = 0; i < GenesIn.GetGenes( ); ++i ) {
			if( !( i % 100 ) )
				g_CatBioUtils.info( "CPCL::Distance( %s, %d, %s, %d, %d, %d, %s, %g ) processing gene %d/%d", szPCL ? szPCL :
					"stdin", iSkip, szDistance, fNormalize, fZScore, fAutocorrelate, szGenes ? szGenes : "", dCutoff, i,
					GenesIn.GetGenes( ) );
			if( ( iOne = veciGenes[ i ] ) == -1 )
				continue;
			adOne = PCL.Get( iOne );
			for( j = ( i + 1 ); j < GenesIn.GetGenes( ); ++j )
				if( ( iTwo = veciGenes[ j ] ) != -1 )
					Dat.Set( i, j, (float)pMeasure->Measure(
						adOne, PCL.GetExperiments( ), PCL.Get( iTwo ), PCL.GetExperiments( ) ) ); }

		if( fNormalize || fZScore )
			Dat.Normalize( !!fNormalize );
		if( !CMeta::IsNaN( dCutoff ) )
			for( i = 0; i < Dat.GetGenes( ); ++i )
				for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
					if( !CMeta::IsNaN( d = Dat.Get( i, j ) ) && ( d < dCutoff ) )
						Dat.Set( i, j, CMeta::GetNaN( ) ); }

	return 0; }

size_t CPCL::GetSkip( ) {

	return c_iSkip; }

CPCLImpl::~CPCLImpl( ) {

	Reset( ); }

void CPCL::Reset( ) {

	CPCLImpl::Reset( ); }

void CPCLImpl::Reset( ) {

	m_Data.Reset( );
	m_vecstrGenes.clear( );
	m_vecstrExperiments.clear( );
	m_vecstrFeatures.clear( );
	m_vecvecstrFeatures.clear( );
	m_setiGenes.clear( ); }

void CPCL::Open( const CPCL& PCL ) {
	size_t					i, j;
	TSetI::const_iterator	iterGene;

	Reset( );
	m_Data.Initialize( PCL.m_Data.GetRows( ), PCL.m_Data.GetColumns( ) );
	for( i = 0; i < m_Data.GetRows( ); ++i )
		for( j = 0; j < m_Data.GetColumns( ); ++j )
			m_Data.Set( i, j, PCL.m_Data.Get( i, j ) );

	for( iterGene = PCL.m_setiGenes.begin( ); iterGene != PCL.m_setiGenes.end( );
		++iterGene )
		m_setiGenes.insert( *iterGene );
	m_vecstrExperiments.resize( PCL.m_vecstrExperiments.size( ) );
	copy( PCL.m_vecstrExperiments.begin( ), PCL.m_vecstrExperiments.end( ),
		m_vecstrExperiments.begin( ) );
	m_vecstrFeatures.resize( PCL.m_vecstrFeatures.size( ) );
	copy( PCL.m_vecstrFeatures.begin( ), PCL.m_vecstrFeatures.end( ),
		m_vecstrFeatures.begin( ) );
	m_vecstrGenes.resize( PCL.m_vecstrGenes.size( ) );
	copy( PCL.m_vecstrGenes.begin( ), PCL.m_vecstrGenes.end( ), m_vecstrGenes.begin( ) );
	m_vecvecstrFeatures.resize( PCL.m_vecvecstrFeatures.size( ) );
	for( i = 0; i < m_vecvecstrFeatures.size( ); ++i ) {
		m_vecvecstrFeatures[ i ].resize( PCL.m_vecvecstrFeatures[ i ].size( ) );
		copy( PCL.m_vecvecstrFeatures[ i ].begin( ), PCL.m_vecvecstrFeatures[ i ].end( ),
			m_vecvecstrFeatures[ i ].begin( ) ); } }

bool CPCL::Open( istream& istmInput ) {

	return Open( istmInput, GetSkip( ) ); }

bool CPCL::Open( istream& istmInput, size_t iFeatures ) {
	vector<float>	vecdData;
	size_t			i, j, k;

	if( !OpenExperiments( istmInput, iFeatures ) )
		return false;
	m_vecvecstrFeatures.resize( m_vecstrFeatures.size( ) - 1 );
	while( OpenGene( istmInput, vecdData ) );

	m_Data.Initialize( m_vecstrGenes.size( ), m_vecstrExperiments.size( ) );
	for( k = i = 0; i < m_Data.GetRows( ); ++i )
		for( j = 0; j < m_Data.GetColumns( ); ++j )
			m_Data.Set( i, j, vecdData[ k++ ] );

	return true; }

bool CPCLImpl::OpenExperiments( istream& istmInput, size_t iFeatures ) {
	char		acLine[ c_iBufferSize ];
	const char*	pc;
	string		strToken;
	size_t		iToken;

	Reset( );
	istmInput.getline( acLine, c_iBufferSize - 1 );
	acLine[ c_iBufferSize - 1 ] = 0;
	for( iToken = 0,pc = acLine; ( strToken = OpenToken( pc, &pc ) ).length( ) || *pc; ++iToken )
		if( iToken <= iFeatures )
			m_vecstrFeatures.push_back( strToken );
		else
			m_vecstrExperiments.push_back( strToken );
	if( !iToken )
		g_CatBioUtils.error( "CPCLImpl::OpenExperiments( %d ) found no experiments", iFeatures );

	return !!iToken; }

bool CPCLImpl::OpenGene( istream& istmInput, vector<float>& vecdData ) {
	char		acLine[ c_iBufferSize ];
	const char*	pc;
	string		strToken;
	size_t		iToken, iData;

	istmInput.getline( acLine, c_iBufferSize - 1 );
	acLine[ c_iBufferSize - 1 ] = 0;
	for( iData = iToken = 0,pc = acLine; ( strToken = OpenToken( pc, &pc ) ).length( ) || *pc; ++iToken ) {
		if( strToken == "EWEIGHT" )
			return true;
		if( !iToken )
			m_vecstrGenes.push_back( strToken );
		else if( iToken < m_vecstrFeatures.size( ) )
			m_vecvecstrFeatures[ iToken - 1 ].push_back( strToken );
		else {
			iData++;
			if( !strToken.length( ) )
				vecdData.push_back( CMeta::GetNaN( ) );
			else
				vecdData.push_back( (float)atof( strToken.c_str( ) ) ); } }

	while( iData++ < m_vecstrExperiments.size( ) )
		vecdData.push_back( CMeta::GetNaN( ) );

	return !!iToken; }

void CPCL::SaveHeader( ostream& ostm, bool fCluster ) const {
	size_t	i;

	if( fCluster )
		ostm << c_szGID << '\t';
	ostm << m_vecstrFeatures[ 0 ];
	for( i = 1; i < m_vecstrFeatures.size( ); ++i )
		ostm << '\t' << m_vecstrFeatures[ i ];
	for( i = 0; i < m_vecstrExperiments.size( ); ++i )
		ostm << '\t' << m_vecstrExperiments[ i ];
	ostm << endl;

	ostm << c_szEWEIGHT;
	for( i = fCluster ? 0 : 1; i < m_vecstrFeatures.size( ); ++i )
		ostm << '\t';
	for( i = 0; i < m_vecstrExperiments.size( ); ++i )
		ostm << '\t' << 1;
	ostm << endl; }

void CPCL::SaveGene( ostream& ostm, size_t iGene, size_t iOrig ) const {
	size_t	i;
	float	d;

	if( iOrig != -1 )
		ostm << c_szGENE << iOrig << '\t';
	ostm << m_vecstrGenes[ iGene ];
	for( i = 0; i < m_vecvecstrFeatures.size( ); ++i )
		ostm << '\t' << m_vecvecstrFeatures[ i ][ iGene ];
	for( i = 0; i < m_vecstrExperiments.size( ); ++i ) {
		ostm << '\t';
		if( !CMeta::IsNaN( d = Get( iGene, i ) ) )
			ostm << Get( iGene, i ); }
	ostm << endl; }

void CPCL::Save( ostream& ostmOutput, const vector<size_t>* pveciGenes ) const {
	size_t	i;

	SaveHeader( ostmOutput, !!pveciGenes );

	for( i = 0; i < m_vecstrGenes.size( ); ++i ) {
		if( m_setiGenes.find( i ) != m_setiGenes.end( ) )
			continue;
		SaveGene( ostmOutput, i, pveciGenes ? (*pveciGenes)[ i ] : -1 ); } }

void CPCL::Open( const vector<string>& vecstrGenes, const vector<string>& vecstrExperiments ) {
	size_t	i, j;

	Reset( );
	m_vecstrFeatures.resize( 3 );
	m_vecstrFeatures[ 0 ] = "GID";
	m_vecstrFeatures[ 1 ] = "NAME";
	m_vecstrFeatures[ 2 ] = "GWEIGHT";
	m_vecvecstrFeatures.resize( m_vecstrFeatures.size( ) - 1 );
	for( i = 0; i < m_vecvecstrFeatures.size( ); ++i )
		m_vecvecstrFeatures[ i ].resize( vecstrGenes.size( ) );
	for( i = 0; i < vecstrGenes.size( ); ++i ) {
		m_vecvecstrFeatures[ 0 ][ i ] = vecstrGenes[ i ];
		m_vecvecstrFeatures[ 1 ][ i ] = "1"; }

	m_vecstrGenes.resize( vecstrGenes.size( ) );
	for( i = 0; i < m_vecstrGenes.size( ); ++i )
		m_vecstrGenes[ i ] = vecstrGenes[ i ];
	m_vecstrExperiments.resize( vecstrExperiments.size( ) );
	for( i = 0; i < m_vecstrExperiments.size( ); ++i )
		m_vecstrExperiments[ i ] = vecstrExperiments[ i ];

	m_Data.Initialize( m_vecstrGenes.size( ), m_vecstrExperiments.size( ) );
	for( i = 0; i < m_Data.GetRows( ); ++i )
		for( j = 0; j < m_Data.GetColumns( ); ++j )
			m_Data.Set( i, j, CMeta::GetNaN( ) ); }

void CPCL::Open( const vector<size_t>& veciGenes, const vector<string>& vecstrGenes,
	const vector<string>& vecstrExperiments ) {
	size_t	i, j;
	char	ac[ 16 ];

	Reset( );
	m_vecstrFeatures.resize( 4 );
	m_vecstrFeatures[ 0 ] = "GID";
	m_vecstrFeatures[ 1 ] = "YORF";
	m_vecstrFeatures[ 2 ] = "NAME";
	m_vecstrFeatures[ 3 ] = "GWEIGHT";
	m_vecvecstrFeatures.resize( m_vecstrFeatures.size( ) - 1 );
	for( i = 0; i < m_vecvecstrFeatures.size( ); ++i )
		m_vecvecstrFeatures[ i ].resize( veciGenes.size( ) );
	for( i = 0; i < veciGenes.size( ); ++i ) {
		m_vecvecstrFeatures[ 0 ][ i ] = m_vecvecstrFeatures[ 1 ][ i ] =
			vecstrGenes[ veciGenes[ i ] ];
		m_vecvecstrFeatures[ 2 ][ i ] = "1"; }

	m_vecstrGenes.resize( vecstrGenes.size( ) );
	for( i = 0; i < m_vecstrGenes.size( ); ++i ) {
		m_vecstrGenes[ i ] = "GENE";
		sprintf_s( ac, "%d", veciGenes[ i ] );
		m_vecstrGenes[ i ] += ac; }
	m_vecstrExperiments.resize( vecstrExperiments.size( ) );
	for( i = 0; i < m_vecstrExperiments.size( ); ++i )
		m_vecstrExperiments[ i ] = vecstrExperiments[ i ];

	m_Data.Initialize( m_vecstrGenes.size( ), m_vecstrExperiments.size( ) );
	for( i = 0; i < m_Data.GetRows( ); ++i )
		for( j = 0; j < m_Data.GetColumns( ); ++j )
			m_Data.Set( i, j, CMeta::GetNaN( ) ); }

void CPCL::SortGenes( const vector<size_t>& veciOrder ) {
	size_t	i;

	CMeta::Permute( m_Data.Get( ), m_Data.GetRows( ), veciOrder );
	CMeta::Permute( m_vecstrGenes, veciOrder );
	for( i = 0; i < m_vecvecstrFeatures.size( ); ++i )
		CMeta::Permute( m_vecvecstrFeatures[ i ], veciOrder ); }

void CPCL::RankTransform( ) {
	size_t	i, j, k;
	size_t*	aiRanks;

	aiRanks = new size_t[ m_Data.GetColumns( ) ];
	for( i = 0; i < m_Data.GetRows( ); ++i ) {
		memset( aiRanks, 0, m_Data.GetColumns( ) * sizeof(*aiRanks) );
		for( j = 0; j < m_Data.GetColumns( ); ++j )
			for( k = 0; k < m_Data.GetColumns( ); ++k )
				if( ( j != k ) && ( m_Data.Get( i, k ) < m_Data.Get( i, j ) ) )
					aiRanks[ j ]++;
		for( j = 0; j < m_Data.GetColumns( ); ++j )
			m_Data.Set( i, j, (float)aiRanks[ j ] ); }
	delete[] aiRanks; }

bool CPCL::AddGenes( const vector<string>& vecstrGenes ) {
	size_t	i, j, iStart;

	iStart = m_Data.GetRows( );
	if( !m_Data.AddRows( vecstrGenes.size( ) ) )
		return false;
	for( i = iStart; i < m_Data.GetRows( ); ++i )
		for( j = 0; j < m_Data.GetColumns( ); ++j )
			m_Data.Set( i, j, CMeta::GetNaN( ) );

	m_vecstrGenes.resize( m_vecstrGenes.size( ) + vecstrGenes.size( ) );
	for( i = 0; i < vecstrGenes.size( ); ++i )
		m_vecstrGenes[ iStart + i ] = vecstrGenes[ i ];
	for( i = 0; i < m_vecvecstrFeatures.size( ); ++i ) {
		m_vecvecstrFeatures[ i ].resize( m_vecvecstrFeatures[ i ].size( ) + vecstrGenes.size( ) );
		if( m_vecstrFeatures[ i + 1 ] == c_szNAME )
			for( j = 0; j < vecstrGenes.size( ); ++j )
				m_vecvecstrFeatures[ i ][ iStart + j ] = vecstrGenes[ j ];
		else if( m_vecstrFeatures[ i + 1 ] == c_szGWEIGHT )
			for( j = 0; j < vecstrGenes.size( ); ++j )
				m_vecvecstrFeatures[ i ][ iStart + j ] = c_szOne; }

	return true; }

void CPCL::Normalize( ) {
	size_t	i, j;
	double	dAve, dStd;

	for( i = 0; i < GetGenes( ); ++i ) {
		dAve = CStatistics::Average( Get( i ), Get( i ) + GetExperiments( ) );
		dStd = sqrt( CStatistics::Variance( Get( i ), Get( i ) + GetExperiments( ), dAve ) );
		for( j = 0; j < GetExperiments( ); ++j )
			Set( i, j, (float)( ( Get( i, j ) - dAve ) / dStd ) ); } }

}

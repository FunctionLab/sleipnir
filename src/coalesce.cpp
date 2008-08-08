#include "stdafx.h"
#include "coalesce.h"
#include "pcl.h"
#include "measure.h"
#include "statistics.h"
#include "fasta.h"

namespace Sleipnir {

// CCoalesceHistograms

void CCoalesceHistograms::Add( const CCoalesceHistograms& Histograms ) {
	TMapStrI::const_iterator	iterType;
	size_t						i, iType, iMotif;

	for( iterType = Histograms.m_mapstriTypes.begin( ); iterType != Histograms.m_mapstriTypes.end( );
		++iterType ) {

		iType = AddType( iterType->first );
		for( i = 0; i < Histograms.m_vecvecvecsHistograms[ iterType->second ].size( ); ++i )
			for( iMotif = 0; iMotif < m_vecvecvecsHistograms[ iType ][ i ].size( ); ++iMotif )
				m_vecvecvecsHistograms[ iType ][ i ][ iMotif ].Add(
					Histograms.m_vecvecvecsHistograms[ iterType->second ][ i ][ 0 ].Get( iMotif ) ); } }

bool CCoalesceHistograms::Add( const SFASTASequence& sSequence, size_t iK ) {
	size_t	i, iType;

	iType = AddType( sSequence.m_strType );
	for( i = 0; i < sSequence.m_vecstrSequences.size( ); ++i )
		if( !Add( sSequence.m_vecstrSequences[ i ], iType, iK, sSequence.m_fIntronFirst && !( i % 2 ) ) )
			return false;

	return true; }

bool CCoalesceHistograms::Add( const string& strSequence, size_t iType, size_t iK, bool fIntron ) {
	size_t	i, iMotif;

	for( i = 0; ( i + iK ) <= strSequence.size( ); ++i ) {
		if( ( iMotif = CCoalesce::GetKMer( strSequence.substr( i, iK ) ) ) == -1 )
			return false;
		m_vecvecvecsHistograms[ iType ][ fIntron ? ESubsequenceIntrons : ESubsequenceExons ][ 0 ].Add(
			iMotif ); }

	return true; }

size_t CCoalesceHistograms::AddType( const string& strType ) {
	TMapStrI::const_iterator	iterType;
	size_t						i, j, iRet;

	if( ( iterType = m_mapstriTypes.find( strType ) ) != m_mapstriTypes.end( ) )
		return iterType->second;

	m_mapstriTypes[ strType ] = iRet = m_vecvecvecsHistograms.size( );
	m_vecstrTypes.push_back( strType );
	m_vecvecvecsHistograms.push_back( vector<vector<SHistogram> >( ) );
	m_vecvecvecsHistograms.back( ).resize( ESubsequenceTotal );
	for( i = 0; i < m_vecvecvecsHistograms.back( ).size( ); ++i ) {
		vector<SHistogram>&	vecsHistograms	= m_vecvecvecsHistograms.back( )[ i ];

		vecsHistograms.resize( m_iMotifs );
		for( j = 0; j < vecsHistograms.size( ); ++j )
			vecsHistograms[ j ].Initialize( m_iBins ); }

	return iRet; }

// CCoalesceCluster

bool CCoalesceCluster::Initialize( const CPCL& PCL, CCoalesceCluster& Pot, float dPValue ) {

	m_setiConditions.clear( );
	m_setiGenes.clear( );
	m_setiMotifs.clear( );

	return ( AddSeedPair( PCL, Pot, dPValue ) && AddCorrelatedGenes( PCL, Pot, dPValue ) ); }

void CCoalesceClusterImpl::Add( size_t iGene, CCoalesceCluster& Pot ) {

	m_setiGenes.insert( iGene );
	Pot.m_setiGenes.erase( iGene ); }

bool CCoalesceClusterImpl::AddSeedPair( const CPCL& PCL, CCoalesceCluster& Pot, float dPValue ) {
	float	dPMiss;
	size_t	iOne, iTwo;
	double	dR;

	if( PCL.GetGenes( ) < 2 )
		return false;
	for( dPMiss = 1; dPMiss >= dPValue; dPMiss *= ( 1 - dPValue ) ) {
		iOne = rand( ) % PCL.GetGenes( );
		do
			iTwo = rand( ) % PCL.GetGenes( );
		while( iOne == iTwo );

		if( ( dR = CMeasurePearson::Pearson( PCL.Get( iOne ), PCL.GetExperiments( ), PCL.Get( iTwo ),
			PCL.GetExperiments( ), IMeasure::EMapNone, NULL, NULL ) ) < 0 )
			continue;
		if( CStatistics::PValuePearson( dR, PCL.GetExperiments( ) ) < dPValue ) {
			Add( iOne, Pot );
			Add( iTwo, Pot );
			return true; } }

	return false; }

bool CCoalesceClusterImpl::AddCorrelatedGenes( const CPCL& PCL, CCoalesceCluster& Pot, float dPValue ) {
	size_t	iGene;
	double	dR;

	CalculateCentroid( PCL );
	for( iGene = 0; iGene < PCL.GetGenes( ); ++iGene ) {
		if( IsGene( iGene ) || ( ( dR = CMeasurePearson::Pearson( &m_vecdCentroid.front( ),
			PCL.GetExperiments( ), PCL.Get( iGene ), PCL.GetExperiments( ), IMeasure::EMapNone, NULL,
			NULL ) ) < 0 ) )
			continue;
		if( CStatistics::PValuePearson( dR, PCL.GetExperiments( ) ) < dPValue )
			Add( iGene, Pot ); }

	return true; }

void CCoalesceCluster::CalculateHistograms( const vector<CCoalesceHistograms>& vecHistograms,
	CCoalesceHistograms& Histograms ) const {
	set<size_t>::const_iterator	iterGene;

	if( vecHistograms.empty( ) )
		return;
	for( iterGene = m_setiGenes.begin( ); iterGene != m_setiGenes.end( ); ++iterGene ) {
		const CCoalesceHistograms&	Histogram	= vecHistograms[ *iterGene ];

		if( Histogram.IsEmpty( ) )
			continue;
		Histograms.Add( Histogram ); } }

void CCoalesceCluster::Subtract( CPCL& PCL ) const {
	set<size_t>::const_iterator	iterGene, iterCondition;

	for( iterGene = m_setiGenes.begin( ); iterGene != m_setiGenes.end( ); ++iterGene )
		for( iterCondition = m_setiConditions.begin( ); iterCondition != m_setiConditions.end( );
			++iterCondition )
			PCL.Get( *iterGene, *iterCondition ) -= m_vecdCentroid[ *iterCondition ]; }

void CCoalesceClusterImpl::CalculateCentroid( const CPCL& PCL ) {
	set<size_t>::const_iterator	iterGene;
	size_t						i;
	float						d;

	m_veciCounts.resize( PCL.GetExperiments( ) );
	fill( m_veciCounts.begin( ), m_veciCounts.end( ), 0  );
	m_vecdCentroid.resize( m_veciCounts.size( ) );
	fill( m_vecdCentroid.begin( ), m_vecdCentroid.end( ), 0.0f );
//	m_vecdStdevs.resize( m_veciCounts.size( ) );
//	fill( m_vecdStdevs.begin( ), m_vecdStdevs.end( ), 0.0f );
	for( iterGene = m_setiGenes.begin( ); iterGene != m_setiGenes.end( ); ++iterGene )
		for( i = 0; i < m_veciCounts.size( ); ++i )
			if( !CMeta::IsNaN( d = PCL.Get( *iterGene, i ) ) ) {
				m_veciCounts[ i ]++;
				m_vecdCentroid[ i ] += d;
//				m_vecdStdevs[ i ] += d * d;
			}
	for( i = 0; i < m_veciCounts.size( ); ++i )
		if( m_veciCounts[ i ] ) {
			m_vecdCentroid[ i ] /= m_veciCounts[ i ];
//			m_vecdStdevs[ i ] = sqrt( ( m_vecdStdevs[ i ] / veciCounts[ i ] ) - ( m_vecdCentroid[ i ] *
//				m_vecdCentroid[ i ] ) );
		} }

bool CCoalesceCluster::SelectConditions( const CPCL& PCL, const CCoalesceCluster& Pot, float dPValue ) {
	vector<float>				vecdCluster, vecdPot;
	size_t						iGene, iCondition;
	set<size_t>::const_iterator	iterGene;
	CMeasureKolmogorovSmirnov	MeasureKS;
	double						dP;

	m_setiConditions.clear( );
	vecdCluster.resize( m_setiGenes.size( ) );
	vecdPot.resize( Pot.m_setiGenes.size( ) );
	for( iCondition = 0; iCondition < PCL.GetExperiments( ); ++iCondition ) {
		for( iterGene = m_setiGenes.begin( ),iGene = 0; iterGene != m_setiGenes.end( ); ++iterGene,++iGene )
			vecdCluster[ iGene ] = PCL.Get( *iterGene, iCondition );
		for( iterGene = Pot.m_setiGenes.begin( ),iGene = 0; iterGene != Pot.m_setiGenes.end( );
			++iterGene,++iGene )
			vecdPot[ iGene ] = PCL.Get( *iterGene, iCondition );
		dP = MeasureKS.Measure( &vecdCluster.front( ), vecdCluster.size( ), &vecdPot.front( ), vecdPot.size( ),
			IMeasure::EMapNone );
		if( dP < dPValue )
			m_setiConditions.insert( iCondition ); }

	return true; }

bool CCoalesceCluster::SelectMotifs( const CCoalesceHistograms& HistsCluster,
	const CCoalesceHistograms& HistsPot, float dPValue ) {
	size_t	iMotif;

	m_setiMotifs.clear( );
	for( iMotif = 0; iMotif < HistsCluster.GetMotifs( ); ++iMotif )
		if( IsSignificant( iMotif, HistsCluster, HistsPot, dPValue ) )
			m_setiMotifs.insert( iMotif );

	return true; }

bool CCoalesceClusterImpl::IsSignificant( size_t iMotif, const CCoalesceHistograms& HistsCluster,
	const CCoalesceHistograms& HistsPot, float dPValue ) const {
	size_t	iType, iSubsequence;

	for( iType = 0; iType < HistsCluster.GetTypes( ); ++iType ) {
		SHistogram	sCluster, sPot;

		for( iSubsequence = 0; iSubsequence < HistsCluster.GetSubsequences( iType ); ++iSubsequence ) {
			if( CStatistics::KSTest( &HistsCluster.Get( iType, iSubsequence, iMotif ).Get( ).front( ),
				&HistsPot.Get( iType, iSubsequence, iMotif ).Get( ).front( ), HistsCluster.GetBins( ) ) <
				dPValue )
				return true;
			sCluster.Add( HistsCluster.Get( iType, iSubsequence, iMotif ) );
			sPot.Add( HistsPot.Get( iType, iSubsequence, iMotif ) ); }
// BUGBUG: multiple hypothesis correction?
		if( CStatistics::KSTest( &sCluster.Get( ).front( ), &sPot.Get( ).front( ), sCluster.GetBins( ) ) <
			dPValue )
			return true; }

	return false; }

bool CCoalesceCluster::SelectGenes( const CPCL& PCL, const vector<CCoalesceHistograms>& vecHistograms,
	const CCoalesceHistograms& HistsCluster, const CCoalesceHistograms& HistsPot, CCoalesceCluster& Pot,
	float dPValue ) {
	static const CCoalesceHistograms	HistsDummy( 0, 0 );
	size_t								i, j, iGene;
	float								d;

	m_veciCounts.resize( PCL.GetExperiments( ) );
	fill( m_veciCounts.begin( ), m_veciCounts.end( ), 0 );
	m_vecdCentroid.resize( m_veciCounts.size( ) );
	fill( m_vecdCentroid.begin( ), m_vecdCentroid.end( ), 0.0f );
	m_vecdStdevs.resize( m_veciCounts.size( ) );
	fill( m_vecdStdevs.begin( ), m_vecdStdevs.end( ), 0.0f );
	for( i = 0; i < PCL.GetGenes( ); ++i )
		for( j = 0; j < PCL.GetExperiments( ); ++j )
			if( !CMeta::IsNaN( d = PCL.Get( i, j ) ) ) {
				m_veciCounts[ j ]++;
				m_vecdCentroid[ j ] += d;
				m_vecdStdevs[ j ] += d * d; }
	for( i = 0; i < m_veciCounts.size( ); ++i )
		if( m_veciCounts[ i ] ) {
			m_vecdCentroid[ i ] /= m_veciCounts[ i ];
			m_vecdStdevs[ i ] = sqrt( ( m_vecdStdevs[ i ] / m_veciCounts[ i ] ) - ( m_vecdCentroid[ i ] *
				m_vecdCentroid[ i ] ) ); }

	CalculateCentroid( PCL );
	Pot.CalculateCentroid( PCL );
	m_setiGenes.clear( );
	Pot.m_setiGenes.clear( );
	for( iGene = 0; iGene < PCL.GetGenes( ); ++iGene )
		if( IsSignificant( iGene, PCL, vecHistograms.empty( ) ? HistsDummy : vecHistograms[ iGene ],
			HistsCluster, HistsPot, Pot, dPValue ) )
			Add( iGene );
		else
			Pot.Add( iGene );

	return true; }

bool CCoalesceClusterImpl::IsSignificant( size_t iGene, const CPCL& PCL, const CCoalesceHistograms& Histograms,
	const CCoalesceHistograms& HistsCluster, const CCoalesceHistograms& HistsPot, const CCoalesceCluster& Pot,
	float dProbability ) const {
	long double	dPGivenMotifs, dPGivenExpression;
	bool		fClustered;

	fClustered = IsGene( iGene );
	dPGivenExpression = CalculateProbabilityExpression( iGene, PCL, Pot, fClustered );
	dPGivenMotifs = CalculateProbabilityMotifs( Histograms, HistsCluster, HistsPot, fClustered );

	return ( ( dPGivenExpression * dPGivenMotifs ) > dProbability ); }

float CCoalesceClusterImpl::CalculateProbabilityExpression( size_t iGene, const CPCL& PCL,
	const CCoalesceCluster& Pot, bool fClustered ) const {
	set<size_t>::const_iterator	iterCondition;
	float						dGene, dCluster, dPot;
	double						dPCluster, dPPot;
	long double					dPIn, dPOut;
	size_t						iPot;

	iPot = PCL.GetGenes( ) - m_veciPrevGenes.size( );
	dPIn = dPOut = 1;
	for( iterCondition = m_setiConditions.begin( ); iterCondition != m_setiConditions.end( );
		++iterCondition ) {
		dGene = PCL.Get( iGene, *iterCondition );
		dCluster = m_vecdCentroid[ *iterCondition ];
		dPot = Pot.m_vecdCentroid[ *iterCondition ];
		if( fClustered )
			dCluster = ( ( dCluster * m_veciPrevGenes.size( ) ) - dGene ) / ( m_veciPrevGenes.size( ) - 1 );
		else
			dPot = ( ( dPot * iPot ) - dGene ) / ( iPot - 1 );

		dPCluster = CStatistics::NormalPDF( dGene, dCluster, m_vecdStdevs[ *iterCondition ] );
		dPPot = CStatistics::NormalPDF( dGene, dPot, m_vecdStdevs[ *iterCondition ] );
		dPIn *= dPCluster;
		dPOut *= dPPot; }

	return (float)( dPIn / ( dPIn + dPOut ) ); }

float CCoalesceClusterImpl::CalculateProbabilityMotifs( const CCoalesceHistograms& Histograms,
	const CCoalesceHistograms& HistsCluster, const CCoalesceHistograms& HistsPot, bool fClustered ) const {
	set<size_t>::const_iterator	iterMotif;
	size_t						i, iType, iSubsequence, iTmp, iCluster, iPot, iTotalCluster, iTotalPot;
	unsigned short				sGene, sCluster, sPot, sTotalGene, sTotalCluster, sTotalPot;
	float						dPCluster, dPPot, dRet, dCur;
	vector<long double>			vecdPIn, vecdPOut;

	if( Histograms.IsEmpty( ) )
		return 1;

	for( iterMotif = m_setiMotifs.begin( ); iterMotif != m_setiMotifs.end( ); ++iterMotif )
		for( iType = 0; iType < HistsCluster.GetTypes( ); ++iType ) {
			if( ( iTmp = Histograms.GetType( HistsCluster.GetType( iType ) ) ) == -1 )
				continue;
			if( vecdPIn.empty( ) ) {
				vecdPIn.resize( HistsCluster.GetSubsequences( iType ) + 1 );
				vecdPOut.resize( vecdPIn.size( ) );
				fill( vecdPIn.begin( ), vecdPIn.end( ), 1 );
				fill( vecdPOut.begin( ), vecdPOut.end( ), 1 ); }
			iTotalCluster = iTotalPot = sTotalGene = sTotalCluster = sTotalPot = 0;
			for( iSubsequence = 0; iSubsequence < HistsCluster.GetSubsequences( iType ); ++iSubsequence ) {
				const SHistogram&	sHistCluster	= HistsCluster.Get( iType, iSubsequence, *iterMotif );
				const SHistogram&	sHistPot		= HistsPot.Get( iType, iSubsequence, *iterMotif );

				sGene = Histograms.Get( iTmp, iSubsequence, 0 ).Get( *iterMotif );
				sCluster = sHistCluster.Get( sGene );
				sPot = sHistPot.Get( sGene );
				iCluster = sHistCluster.GetTotal( );
				iPot = sHistPot.GetTotal( );
				if( fClustered ) {
					sCluster--;
					iCluster--; }
				else {
					sPot--;
					iPot--; }
				sTotalGene += sGene;
				sTotalCluster += sCluster;
				sTotalPot += sPot;
				iTotalCluster += iCluster;
				iTotalPot += iPot;

				dPCluster = (float)( sCluster + 1 ) / ( iCluster + sHistCluster.GetBins( ) );
				dPPot = (float)( sPot + 1 ) / ( iPot + sHistPot.GetBins( ) );
				vecdPIn[ iSubsequence ] *= dPCluster;
				vecdPOut[ iSubsequence ] *= dPPot; }
			dPCluster = (float)( sTotalCluster + 1 ) / ( iTotalCluster +
				HistsCluster.Get( iType, 0, *iterMotif ).GetBins( ) );
			dPPot = (float)( sTotalPot + 1 ) / ( iTotalPot +
				HistsPot.Get( iType, 0, *iterMotif ).GetBins( ) );
			vecdPIn[ iSubsequence ] *= dPCluster;
			vecdPOut[ iSubsequence ] *= dPPot; }

// BUGBUG: how to best handle this - multiple hypothesis, and, or, something?
	dRet = 0;
	for( i = 0; i < vecdPIn.size( ); ++i )
		if( ( dCur = (float)( vecdPIn[ i ] / ( vecdPIn[ i ] + vecdPOut[ i ] ) ) ) > dRet )
			dRet = dCur;

	return dRet; }

bool CCoalesceCluster::Save( size_t iID, const CPCL& PCL ) const {
	ofstream	ofsm;
	char		acTmp[ 16 ];
	CPCL		PCLCopy;
	size_t		iGeneFrom, iGeneTo, iExpFrom, iExpTo;

	sprintf_s( acTmp, "./c%04d_XXXXXX", iID );
	_mktemp_s( acTmp );
	ofsm.open( ( (string)acTmp + ".txt" ).c_str( ) );
	if( !ofsm.is_open( ) )
		return false;
	PCLCopy.Open( PCL );

// Reorder header
	for( iExpTo = iExpFrom = 0; iExpFrom < PCL.GetExperiments( ); ++iExpFrom )
		if( IsCondition( iExpFrom ) )
			PCLCopy.SetExperiment( iExpTo++, "*" + PCL.GetExperiment( iExpFrom ) );
	for( iExpFrom = 0; iExpFrom < PCL.GetExperiments( ); ++iExpFrom )
		if( !IsCondition( iExpFrom ) )
			PCLCopy.SetExperiment( iExpTo++, PCL.GetExperiment( iExpFrom ) );

// Reorder genes
	for( iGeneTo = iGeneFrom = 0; iGeneFrom < PCL.GetGenes( ); ++iGeneFrom )
		if( IsGene( iGeneFrom ) && !SaveCopy( PCL, iGeneFrom, PCLCopy, iGeneTo++, true ) )
			return false;
	for( iGeneFrom = 0; iGeneFrom < PCL.GetGenes( ); ++iGeneFrom )
		if( !IsGene( iGeneFrom ) && !SaveCopy( PCL, iGeneFrom, PCLCopy, iGeneTo++, false ) )
			return false;

	PCLCopy.Save( ofsm );
	return true; }

bool CCoalesceClusterImpl::SaveCopy( const CPCL& PCLFrom, size_t iGeneFrom, CPCL& PCLTo, size_t iGeneTo,
	bool fClustered ) const {
	size_t	i, iExpTo, iExpFrom;

	PCLTo.SetGene( iGeneTo, PCLFrom.GetGene( iGeneFrom ) );
	for( i = 1; i < PCLFrom.GetFeatures( ); ++i )
		PCLTo.SetFeature( iGeneTo, i, (string)( ( fClustered && ( i == 1 ) ) ? "*" : "" ) +
			PCLFrom.GetFeature( iGeneFrom, i ) );
	for( iExpTo = iExpFrom = 0; iExpFrom < PCLFrom.GetExperiments( ); ++iExpFrom )
		if( IsCondition( iExpFrom ) )
			PCLTo.Set( iGeneTo, iExpTo++, PCLFrom.Get( iGeneFrom, iExpFrom ) );
	for( iExpFrom = 0; iExpFrom < PCLFrom.GetExperiments( ); ++iExpFrom )
		if( !IsCondition( iExpFrom ) )
			PCLTo.Set( iGeneTo, iExpTo++, PCLFrom.Get( iGeneFrom, iExpFrom ) );

	return true; }

// CCoalesce

size_t CCoalesce::GetKMer( const string& strKMer ) {
	static const char	c_acLetters[]	= "ACGT";
	static const size_t	c_iShift		= 2; // ceil( log2( ARRAYSIZE(c_acLetters) ) )
	size_t		i, iRet;
	const char*	pc;

	for( iRet = i = 0; i < strKMer.size( ); ++i ) {
		if( !( pc = strchr( c_acLetters, strKMer[ i ] ) ) ) {
			g_CatSleipnir.warn( "CCoalesce::GetKMer( %s ) invalid character: %c", strKMer.c_str( ),
				strKMer[ i ] );
			return -1; }
		iRet = ( iRet << c_iShift ) | ( pc - c_acLetters ); }

	return iRet; }

bool CCoalesce::Cluster( const CPCL& PCL, const CFASTA& FASTA, vector<CCoalesceCluster>& vecClusters ) {
	CPCL						PCLCopy;
	size_t						i, j, iKMers;
	vector<size_t>				veciPCL2FASTA;
	vector<CCoalesceHistograms>	vecHistograms;

	PCLCopy.Open( PCL );
	iKMers = CountKMers( GetK( ) );
	veciPCL2FASTA.resize( PCL.GetGenes( ) );
	for( i = 0; i < veciPCL2FASTA.size( ); ++i ) {
		if( !FASTA.GetGenes( ) ) {
			veciPCL2FASTA[ i ] = -1;
			continue; }
		vecHistograms.push_back( CCoalesceHistograms( 1, iKMers ) );
		if( ( veciPCL2FASTA[ i ] = FASTA.GetGene( PCL.GetGene( i ) ) ) != -1 ) {
			vector<SFASTASequence>	vecsSequences;

			FASTA.Get( veciPCL2FASTA[ i ], vecsSequences );
			for( j = 0; j < vecsSequences.size( ); ++j )
				vecHistograms.back( ).Add( vecsSequences[ j ], GetK( ) ); } }
	while( true ) {
		CCoalesceCluster	Cluster, Pot;

		for( i = 0; i < PCLCopy.GetGenes( ); ++i )
			Pot.Add( i );
		if( !Cluster.Initialize( PCLCopy, Pot, GetPValueCorrelation( ) ) )
			break;
		g_CatSleipnir.notice( "CCoalesce::Cluster( ) initialized %d genes with: %s, %s",
			Cluster.GetGenes( ).size( ), PCLCopy.GetGene( *Cluster.GetGenes( ).begin( ) ).c_str( ),
			PCLCopy.GetGene( *(++Cluster.GetGenes( ).begin( )) ).c_str( ) );
		while( !( Cluster.IsConverged( ) || Cluster.IsEmpty( ) ) ) {
			CCoalesceHistograms	HistsCluster( iKMers, GetBins( ) ), HistsPot( iKMers, GetBins( ) );

			Cluster.Snapshot( );
			Cluster.CalculateHistograms( vecHistograms, HistsCluster );
			Pot.CalculateHistograms( vecHistograms, HistsPot );
			if( !( Cluster.SelectConditions( PCLCopy, Pot, GetPValueCondition( ) ) &&
				Cluster.SelectMotifs( HistsCluster, HistsPot, GetPValueMotif( ) ) &&
				Cluster.SelectGenes( PCLCopy, vecHistograms, HistsCluster, HistsPot, Pot,
				GetProbabilityGene( ) ) ) )
				return false;
			g_CatSleipnir.info( "CCoalesce::Cluster( ) processed %d genes, %d conditions, %d motifs",
				Cluster.GetGenes( ).size( ), Cluster.GetConditions( ).size( ), Cluster.GetMotifs( ).size( ) );
			if( IsOutputIntermediate( ) )
				Cluster.Save( vecClusters.size( ), PCLCopy ); }
		if( Cluster.IsConverged( ) ) {
			vecClusters.push_back( Cluster );
			Cluster.Subtract( PCLCopy ); } }

	return true; }

}

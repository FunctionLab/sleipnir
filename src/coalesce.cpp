/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "stdafx.h"
#include "coalesce.h"
#include "fasta.h"
#include "pcl.h"

namespace Sleipnir {

const char*	CCoalesceSequencerBase::c_aszSubsequences[]	= {"total", "introns", "exons", NULL};

// CCoalesceGeneScores

bool CCoalesceGeneScores::Add( size_t iGene, const CCoalesceMotifLibrary& Motifs,
	const SFASTASequence& sSequence, SCoalesceModifierCache& sModifiers, vector<vector<float> >& vecvecdCounts,
	vector<size_t>& veciLengths ) {
	size_t		i, iType, iSubsequence, iLength, iOffset;
	uint32_t	iMotif;
	float		d;

	sModifiers.SetType( sSequence.m_strType );
	vecvecdCounts.resize( ESubsequenceEnd );
	for( i = 0; i < vecvecdCounts.size( ); ++i )
		fill( vecvecdCounts[ i ].begin( ), vecvecdCounts[ i ].end( ), 0.0f );
	veciLengths.resize( ESubsequenceEnd );
	fill( veciLengths.begin( ), veciLengths.end( ), 0 );

	iType = AddType( sSequence.m_strType );
	for( iOffset = i = 0; i < sSequence.m_vecstrSequences.size( );
		iOffset += sSequence.m_vecstrSequences[ i++ ].length( ) )
		if( !Add( Motifs, sSequence.m_vecstrSequences[ i ], iType, sSequence.m_fIntronFirst == !( i % 2 ),
			vecvecdCounts, veciLengths, iOffset, sModifiers ) )
			return false;
	for( iSubsequence = ESubsequenceBegin; iSubsequence < vecvecdCounts.size( ); ++iSubsequence ) {
		const vector<float>&	vecdCounts	= vecvecdCounts[ iSubsequence ];

		if( iLength = veciLengths[ iSubsequence ] )
			for( iMotif = 0; iMotif < vecdCounts.size( ); ++iMotif )
				if( d = vecdCounts[ iMotif ] )
					Set( iType, (ESubsequence)iSubsequence, iGene, iMotif, d / iLength,
						vecdCounts.size( ) ); }

	return true; }

bool CCoalesceGeneScores::Add( const CCoalesceMotifLibrary& Motifs, const std::string& strSequence,
	size_t iType, bool fIntron, std::vector<std::vector<float> >& vecvecdCounts,
	std::vector<size_t>& veciLengths, size_t iOffset, SCoalesceModifierCache& sModifiers ) {
	size_t			i, j;
	ESubsequence	eSubsequence;

	eSubsequence = fIntron ? ESubsequenceIntrons : ESubsequenceExons;
	veciLengths[ ESubsequenceTotal ] += strSequence.length( );
	veciLengths[ eSubsequence ] += strSequence.length( );

	sModifiers.InitializeWeight( Motifs.GetK( ), iOffset );
// BUGBUG: make me span intron/exon boundaries
	for( i = 0; ( i + Motifs.GetK( ) ) <= strSequence.size( ); ++i ) {
		string				strKMer	= strSequence.substr( i, Motifs.GetK( ) );
		vector<uint32_t>	veciMotifs;

		if( !Motifs.GetMatches( strKMer, veciMotifs ) ) {
			g_CatSleipnir.error( "CCoalesceHistograms::Add( %s, %d, %d ) unrecognized kmer: %s",
				strSequence.c_str( ), iType, fIntron, strKMer.c_str( ) );
			return false; }
		for( j = 0; j < veciMotifs.size( ); ++j )
			Add( eSubsequence, veciMotifs[ j ], Motifs.GetMotifs( ), vecvecdCounts,
				sModifiers.GetWeight( Motifs.GetK( ) ) / Motifs.GetK( ) );
		sModifiers.AddWeight( Motifs.GetK( ), iOffset, i ); }

	return true; }

bool CCoalesceGeneScores::Add( size_t iGene, const CCoalesceMotifLibrary& Motifs,
	const SFASTASequence& sSequence, SCoalesceModifierCache& sModifiers, uint32_t iMotif,
	vector<float>& vecdScores, vector<size_t>& veciLengths ) {
	size_t	i, iType, iSubsequence, iLength, iOffset;
	float	dScore;

	sModifiers.SetType( sSequence.m_strType );
	vecdScores.resize( ESubsequenceEnd );
	fill( vecdScores.begin( ), vecdScores.end( ), 0.0f );
	veciLengths.resize( ESubsequenceEnd );
	fill( veciLengths.begin( ), veciLengths.end( ), 0 );

	iType = AddType( sSequence.m_strType );
	for( iOffset = i = 0; i < sSequence.m_vecstrSequences.size( );
		iOffset += sSequence.m_vecstrSequences[ i++ ].length( ) )
		if( !Add( Motifs, sSequence.m_vecstrSequences[ i ], iType, sSequence.m_fIntronFirst == !( i % 2 ),
			iMotif, vecdScores, veciLengths, iOffset, sModifiers ) )
			return false;
	for( iSubsequence = ESubsequenceBegin; iSubsequence < vecdScores.size( ); ++iSubsequence )
		if( ( iLength = veciLengths[ iSubsequence ] ) && ( dScore = vecdScores[ iSubsequence ] ) )
			Set( iType, (ESubsequence)iSubsequence, iGene, iMotif, dScore / iLength );

	return true; }

bool CCoalesceGeneScores::Add( const CCoalesceMotifLibrary& Motifs, const std::string& strSequence,
	size_t iType, bool fIntron, uint32_t iMotif, std::vector<float>& vecdScores,
	std::vector<size_t>& veciLengths, size_t iOffset, SCoalesceModifierCache& sModifiers ) {
	ESubsequence	eSubsequence;
	float			dScore;

	eSubsequence = fIntron ? ESubsequenceIntrons : ESubsequenceExons;
	veciLengths[ ESubsequenceTotal ] += strSequence.length( );
	veciLengths[ eSubsequence ] += strSequence.length( );

// BUGBUG: make me span intron/exon boundaries
	if( CMeta::IsNaN( dScore = Motifs.GetMatch( strSequence, iMotif, iOffset, sModifiers ) ) )
		return false;
	vecdScores[ ESubsequenceTotal ] += dScore;
	vecdScores[ eSubsequence ] += dScore;

	return true; }

void CCoalesceGeneScores::Subtract( const SMotifMatch& sMotif, size_t iGene ) {
	size_t	iType;
	float*	adScores;

	if( ( ( iType = GetType( sMotif.m_strType ) ) != -1 ) &&
		( adScores = Get( iType, sMotif.m_eSubsequence, iGene ) ) )
		adScores[ sMotif.m_iMotif ] -= sMotif.m_dAverage; }

// CCoalesceGroupHistograms

bool CCoalesceGroupHistograms::Add( const CCoalesceGeneScores& GeneScores, size_t iGene, bool fSubtract,
	uint32_t iMotifThem ) {
	size_t			iTypeThem, iTypeUs, iSubsequence;
	uint32_t		iMotifUs;
	unsigned short	sDelta;
	float			dValue;

	sDelta = fSubtract ? -1 : 1;
	for( iTypeThem = 0; iTypeThem < GeneScores.GetTypes( ); ++iTypeThem ) {
// lock
		iTypeUs = AddType( GeneScores.GetType( iTypeThem ) );
// unlock
		for( iSubsequence = ESubsequenceBegin; iSubsequence < ESubsequenceEnd; ++iSubsequence ) {
			CCoalesceHistogramSet<>&	Histograms	= Get( iTypeUs, (ESubsequence)iSubsequence );
			const float*				adScores	= GeneScores.Get( iTypeThem, (ESubsequence)iSubsequence,
														iGene );

			if( !adScores )
				continue;
// lock
			if( !Histograms.GetMembers( ) ) {
				vector<float>	vecdEdges;

				vecdEdges.resize( m_iBins );
				for( int i = 0; (size_t)i < vecdEdges.size( ); ++i )
					vecdEdges[ i ] = ( i * m_dStep ) - ( m_dStep / 2 );
				Histograms.Initialize( GetMotifs( ), vecdEdges ); }
// unlock
			if( iMotifThem == -1 ) {
				for( iMotifUs = 0; iMotifUs < GeneScores.GetMotifs( ); ++iMotifUs )
					if( dValue = adScores[ iMotifUs ] )
						Histograms.Add( iMotifUs, dValue, sDelta ); }
			else if( dValue = adScores[ iMotifThem ] )
				Histograms.Add( iMotifThem, dValue, sDelta ); } }

	return true; }

void CCoalesceGroupHistograms::Save( std::ostream& ostm, const CCoalesceMotifLibrary* pMotifs ) const {
	size_t		iType, iSubsequence;
	uint32_t	iMotif;

	for( iType = 0; iType < GetTypes( ); ++iType )
		for( iSubsequence = 0; iSubsequence < GetSubsequences( iType ); ++iSubsequence ) {
			const CCoalesceHistogramSet<>&	Histograms	= Get( iType, (ESubsequence)iSubsequence );

			if( !Histograms.GetTotal( ) )
				continue;
			ostm << GetType( iType ) << '\t' << GetSubsequence( (ESubsequence)iSubsequence ) << endl;
			for( iMotif = 0; iMotif < ( pMotifs ? pMotifs->GetMotifs( ) : GetMotifs( ) ); ++iMotif ) {
				if( Histograms.Get( iMotif, 0.0f ) == Histograms.GetTotal( ) )
					continue;
				if( pMotifs )
					ostm << pMotifs->GetMotif( iMotif );
				else
					ostm << iMotif;
				ostm << endl << Histograms.Save( iMotif ) << endl; } } }

bool CCoalesceGroupHistograms::IsSimilar( const CCoalesceMotifLibrary* pMotifs, const SMotifMatch& sOne,
	const SMotifMatch& sTwo, float dPValue ) const {
	size_t	iTypeOne, iTypeTwo;
	double	dAverage, dZ, dP;

	if( ( ( iTypeOne = GetType( sOne.m_strType ) ) == -1 ) ||
		( ( iTypeTwo = GetType( sTwo.m_strType ) ) == -1 ) )
		return false;
	{
		const CCoalesceHistogramSet<>&	HistOne	= Get( iTypeOne, sOne.m_eSubsequence );
		const CCoalesceHistogramSet<>&	HistTwo	= Get( iTypeTwo, sTwo.m_eSubsequence );

		dP = 1 - HistOne.CohensD( sOne.m_iMotif, HistTwo, sTwo.m_iMotif, false, dAverage, dZ );
		if( g_CatSleipnir.isDebugEnabled( ) )
			g_CatSleipnir.debug( "CCoalesceGroupHistograms::IsSimilar( %s, %s, %g ) got p-value %g",
				sOne.Save( pMotifs ).c_str( ), sTwo.Save( pMotifs ).c_str( ), dPValue, dP );
		return ( dP < dPValue );
	} }

// CCoalesce

CCoalesceImpl::~CCoalesceImpl( ) {

	Clear( ); }

void CCoalesceImpl::Clear( ) {

	if( m_fMotifs && m_pMotifs )
		delete m_pMotifs;
	m_pMotifs = NULL; }

size_t CCoalesceImpl::GetMotifCount( ) const {

	return ( m_pMotifs ? m_pMotifs->GetMotifs( ) : 0 ); }

void* CCoalesceImpl::ThreadCombineMotif( void* pData ) {
	SThreadCombineMotif*	psData;
	size_t					i, j;
	vector<float>			vecdScores;
	vector<size_t>			veciLengths;

	psData = (SThreadCombineMotif*)pData;
	SCoalesceModifierCache	sModifiers( *psData->m_psModifiers );

	for( i = psData->m_iOffset; i < psData->m_pveciPCL2FASTA->size( ); i += psData->m_iStep )
		if( (*psData->m_pveciPCL2FASTA)[ i ] != -1 ) {
			vector<SFASTASequence>	vecsSequences;

			if( psData->m_pFASTA->Get( (*psData->m_pveciPCL2FASTA)[ i ], vecsSequences ) ) {
				sModifiers.Get( i );
				for( j = 0; j < vecsSequences.size( ); ++j )
					if( psData->m_pGeneScores->Add( i, *psData->m_pMotifs, vecsSequences[ j ],
						sModifiers, psData->m_iMotif, vecdScores, veciLengths ) )
						break; } }

	return NULL; }

bool CCoalesceImpl::CombineMotifs( const CFASTA& FASTA, const vector<size_t>& veciPCL2FASTA,
	SCoalesceModifiers& sModifiers, const CCoalesceCluster& Cluster, size_t iThreads,
	CCoalesceGeneScores& GeneScores, CCoalesceGroupHistograms& HistsCluster,
	CCoalesceGroupHistograms& HistsPot ) const {
	set<SMotifMatch>::const_iterator	iterMotifOne, iterMotifTwo;
	uint32_t							iMotif;
	size_t								i;
	vector<pthread_t>					vecpthdThreads;
	vector<SThreadCombineMotif>			vecsThreads;

	if( !GeneScores.GetMotifs( ) || !m_pMotifs || ( Cluster.GetMotifs( ).size( ) > m_iSizeMaximum ) )
		return true;

	for( iterMotifOne = Cluster.GetMotifs( ).begin( ); iterMotifOne != Cluster.GetMotifs( ).end( );
		++iterMotifOne )
		for( iterMotifTwo = Cluster.GetMotifs( ).begin( ); iterMotifTwo != Cluster.GetMotifs( ).end( );
			++iterMotifTwo ) {
			if( ( iterMotifOne->m_iMotif == iterMotifTwo->m_iMotif ) ||
				!HistsCluster.IsSimilar( m_pMotifs, *iterMotifOne, *iterMotifTwo, m_dPValueMerge ) ||
				( ( iMotif = m_pMotifs->Merge( iterMotifOne->m_iMotif, iterMotifTwo->m_iMotif,
				m_dCutoffMerge ) ) == -1 ) )
				continue;

			vecpthdThreads.resize( iThreads );
			vecsThreads.resize( vecpthdThreads.size( ) );
			for( i = 0; i < vecsThreads.size( ); ++i ) {
				vecsThreads[ i ].m_iOffset = i;
				vecsThreads[ i ].m_iStep = vecsThreads.size( );
				vecsThreads[ i ].m_pveciPCL2FASTA = &veciPCL2FASTA;
				vecsThreads[ i ].m_pGeneScores = &GeneScores;
				vecsThreads[ i ].m_pMotifs = m_pMotifs;
				vecsThreads[ i ].m_iMotif = iMotif;
				vecsThreads[ i ].m_pFASTA = &FASTA;
				vecsThreads[ i ].m_psModifiers = &sModifiers;
				if( pthread_create( &vecpthdThreads[ i ], NULL, ThreadCombineMotif, &vecsThreads[ i ] ) ) {
					g_CatSleipnir.error( "CCoalesceImpl::CombineMotifs( %d ) could not combine motif: %s",
						iThreads, m_pMotifs->GetMotif( iMotif ).c_str( ) );
					return false; } }
			for( i = 0; i < vecpthdThreads.size( ); ++i )
				pthread_join( vecpthdThreads[ i ], NULL );
			for( i = 0; i < veciPCL2FASTA.size( ); ++i )
				if( veciPCL2FASTA[ i ] != -1 ) {
					if( Cluster.IsGene( i ) )
						HistsCluster.Add( GeneScores, i, false, iMotif );
					else
						HistsPot.Add( GeneScores, i, false, iMotif ); } }

	return true; }

bool CCoalesceImpl::InitializeDatasets( const CPCL& PCL ) {
	size_t			i, j;
	vector<bool>	vecfDataset;

	vecfDataset.resize( PCL.GetExperiments( ) );
	for( i = 0; i < m_vecsDatasets.size( ); ++i ) {
		SCoalesceDataset&	sDataset	= m_vecsDatasets[ i ];

		for( j = 0; j < sDataset.GetConditions( ); ++j )
			vecfDataset[ sDataset.GetCondition( j ) ] = true; }
	for( i = 0; i < vecfDataset.size( ); ++i )
		if( !vecfDataset[ i ] )
			m_vecsDatasets.push_back( SCoalesceDataset( i ) );
// This has to be done separately because of a weird optimization bug in GCC
// If called in the same scope as the constructor, it's optimized to the wrong order and a double free
	for( i = 0; i < m_vecsDatasets.size( ); ++i )
		m_vecsDatasets[ i ].CalculateCovariance( PCL );

	return true; }

bool CCoalesceImpl::InitializeGeneScores( const CPCL& PCL, const CFASTA& FASTA, vector<size_t>& veciPCL2FASTA,
	SCoalesceModifiers& sMods, CCoalesceGeneScores& GeneScores ) {
	size_t	i, j;

	if( FASTA.GetGenes( ) ) {
		veciPCL2FASTA.resize( PCL.GetGenes( ) );
		for( i = 0; i < veciPCL2FASTA.size( ); ++i )
			veciPCL2FASTA[ i ] = FASTA.GetGene( PCL.GetGene( i ) ); }
	sMods.Initialize( PCL );
	if( FASTA.GetGenes( ) ) {
		vector<vector<float> >	vecvecdCounts;
		vector<size_t>			veciLengths;
		SCoalesceModifierCache	sModifiers( sMods );

		if( !m_pMotifs ) {
			m_fMotifs = true;
			m_pMotifs = new CCoalesceMotifLibrary( m_iK ); }
		GeneScores.SetGenes( PCL.GetGenes( ) );
		for( i = 0; i < veciPCL2FASTA.size( ); ++i )
			if( veciPCL2FASTA[ i ] != -1 ) {
				vector<SFASTASequence>	vecsSequences;

				if( FASTA.Get( veciPCL2FASTA[ i ], vecsSequences ) ) {
					sModifiers.Get( i );
					for( j = 0; j < vecsSequences.size( ); ++j )
						if( !GeneScores.Add( i, *m_pMotifs, vecsSequences[ j ], sModifiers, vecvecdCounts,
							veciLengths ) )
							return false; } } }
	if( !GeneScores.GetMotifs( ) )
		Clear( );

	return true; }

bool CCoalesce::Cluster( const CPCL& PCL, const CFASTA& FASTA, vector<CCoalesceCluster>& vecClusters ) {
	static const float			c_dEpsilon	= 1e-10f;
	CPCL						PCLCopy;
	size_t						i;
	vector<size_t>				veciPCL2FASTA;
	CCoalesceGeneScores			GeneScores;
	set<pair<size_t, size_t> >	setpriiSeeds;
	SCoalesceModifiers			sModifiers;
	float						dFailure;

	for( i = 0; i < m_vecpWiggles.size( ); ++i )
		sModifiers.Add( m_vecpWiggles[ i ] );
	PCLCopy.Open( PCL );
	PCLCopy.Normalize( CPCL::ENormalizeColumn );
	if( !( InitializeDatasets( PCLCopy ) && InitializeGeneScores( PCLCopy, FASTA, veciPCL2FASTA, sModifiers,
		GeneScores ) ) )
		return false;
	if( g_CatSleipnir.isNoticeEnabled( ) ) {
		size_t			iSequences;
		ostringstream	ossm;

		for( iSequences = i = 0; i < veciPCL2FASTA.size( ); ++i )
			if( veciPCL2FASTA[ i ] != -1 )
				iSequences++;
		g_CatSleipnir.notice( "CCoalesce::Cluster( ) running with %d genes, %d conditions, and %d sequences",
			PCL.GetGenes( ), PCL.GetExperiments( ), iSequences );
		for( i = 0; i < PCL.GetExperiments( ); ++i )
			ossm << ( i ? "\t" : "" ) << PCL.GetExperiment( i );
		g_CatSleipnir.notice( ossm.str( ) );
		g_CatSleipnir.notice( "k %d, P gene %g, p condition %g, p motif %g, p correlation %g", GetK( ),
			GetProbabilityGene( ), GetPValueCondition( ), GetPValueMotif( ), GetPValueCorrelation( ) );
		g_CatSleipnir.notice( "p merge %g, cutoff merge %g, penalty gap %g, penalty mismatch %g",
			GetPValueMerge( ), GetCutoffMerge( ), GetMotifs( ) ? GetMotifs( )->GetPenaltyGap( ) : 0,
			GetMotifs( ) ? GetMotifs( )->GetPenaltyMismatch( ) : 0 );
		g_CatSleipnir.notice( "correlation pairs %d, bases %d, min size %d, max size %d",
			GetNumberCorrelation( ), GetBasesPerMatch( ), GetSizeMinimum( ), GetSizeMaximum( ) ); }
	for( dFailure = 1; dFailure > c_dEpsilon; dFailure *= GetPValueCorrelation( ) ) {
		CCoalesceCluster			Cluster, Pot;
		CCoalesceGroupHistograms	HistsCluster( GetMotifCount( ), GetBins( ), 1.0f / GetBasesPerMatch( ) );
		CCoalesceGroupHistograms	HistsPot( GetMotifCount( ), GetBins( ), 1.0f / GetBasesPerMatch( ) );

		Pot.SetGenes( PCLCopy.GetGenes( ) );
		Pot.CalculateHistograms( GeneScores, HistsPot, NULL );
		if( !Cluster.Initialize( PCLCopy, Pot, m_vecsDatasets, setpriiSeeds, GetNumberCorrelation( ),
			GetPValueCorrelation( ), GetThreads( ) ) )
			continue;
		g_CatSleipnir.notice( "CCoalesce::Cluster( ) initialized %d genes", Cluster.GetGenes( ).size( ) );
		if( g_CatSleipnir.isDebugEnabled( ) ) {
			ostringstream				ossm;
			set<size_t>::const_iterator	iterGene;

			ossm << "CCoalesce::Cluster( ) initialized:";
			for( iterGene = Cluster.GetGenes( ).begin( ); iterGene != Cluster.GetGenes( ).end( ); ++iterGene )
				ossm << ' ' << PCL.GetGene( *iterGene );
			g_CatSleipnir.debug( ossm.str( ) ); }
		if( Cluster.GetGenes( ).size( ) < GetSizeMinimum( ) )
			continue;
		while( !( Cluster.IsConverged( ) || Cluster.IsEmpty( ) ) ) {
			Cluster.CalculateHistograms( GeneScores, HistsCluster, &HistsPot );
			Cluster.Snapshot( GeneScores, HistsCluster );
			Pot.Snapshot( GeneScores, HistsPot );
			if( !Cluster.SelectConditions( PCLCopy, Pot, GetThreads( ), GetPValueCondition( ) ) )
				return false;
			if( Cluster.IsEmpty( ) ) {
				g_CatSleipnir.notice( "CCoalesce::Cluster( ) selected no conditions" );
				break; }
			if( ( Cluster.GetGenes( ).size( ) >= GetSizeMinimum( ) ) &&
				!( CombineMotifs( FASTA, veciPCL2FASTA, sModifiers, Cluster, GetThreads( ), GeneScores,
				HistsCluster, HistsPot ) && Cluster.SelectMotifs( HistsCluster, HistsPot, GetPValueMotif( ),
				GetThreads( ), GetMotifs( ) ) ) )
				return false;
			if( !Cluster.SelectGenes( PCLCopy, GeneScores, HistsCluster, HistsPot, GetSizeMinimum( ),
				GetThreads( ), Pot, GetProbabilityGene( ), GetMotifs( ) ) )
				return false;
			g_CatSleipnir.notice( "CCoalesce::Cluster( ) processed %d genes, %d datasets, %d motifs",
				Cluster.GetGenes( ).size( ), Cluster.GetDatasets( ).size( ), Cluster.GetMotifs( ).size( ) ); }
		if( Cluster.IsConverged( ) && ( Cluster.GetGenes( ).size( ) >= GetSizeMinimum( ) ) ) {
			g_CatSleipnir.notice( "CCoalesce::Cluster( ) finalizing cluster" );
			setpriiSeeds.clear( );
			if( IsOutputIntermediate( ) )
				Cluster.Save( GetDirectoryIntermediate( ), vecClusters.size( ), PCLCopy, GetMotifs( ) );
			vecClusters.push_back( Cluster );
			Cluster.Subtract( PCLCopy );
			Cluster.Subtract( GeneScores );
			for( i = 0; i < m_vecsDatasets.size( ); ++i )
				m_vecsDatasets[ i ].CalculateCovariance( PCLCopy );
			dFailure = 1; }
		else
			g_CatSleipnir.notice( "CCoalesce::Cluster( ) discarding cluster (failure %g)", dFailure ); }

	return true; }

}

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
#include "coalescecluster.h"
#include "coalescei.h"
#include "coalescemotifs.h"
#include "pcl.h"
#include "halfmatrix.h"
#include "statistics.h"
#include "clusthierarchical.h"
#include "pst.h"

namespace Sleipnir {

const char	CCoalesceClusterImpl::c_szMotifs[]			= "_motifs.txt";

bool CCoalesceCluster::Initialize( const CPCL& PCL, CCoalesceCluster& Pot,
	const std::vector<SCoalesceDataset>& vecsDatasets, std::set<std::pair<size_t, size_t> >& setpriiSeeds,
	size_t iPairs, float dPValue, size_t iThreads ) {
	size_t	i;
	float	dFraction;

	m_setiDatasets.clear( );
	m_setiGenes.clear( );
	m_setsMotifs.clear( );
	m_vecsDatasets.resize( vecsDatasets.size( ) );
	Pot.m_vecsDatasets.resize( vecsDatasets.size( ) );
	for( i = 0; i < m_vecsDatasets.size( ); ++i ) {
		m_vecsDatasets[ i ].m_psDataset = Pot.m_vecsDatasets[ i ].m_psDataset = &vecsDatasets[ i ];
		m_setiDatasets.insert( i ); }
	m_vecdPriors.resize( PCL.GetGenes( ) );
	fill( m_vecdPriors.begin( ), m_vecdPriors.end( ), 1.0f );
	dFraction = min( 1.0f, 2.0f * iPairs / PCL.GetGenes( ) / ( PCL.GetGenes( ) - 1 ) );

	return ( AddSeedPair( PCL, Pot, setpriiSeeds, dFraction, dPValue, iThreads ) &&
		AddCorrelatedGenes( PCL, Pot, dPValue ) ); }

void CCoalesceClusterImpl::Add( size_t iGene, CCoalesceCluster& Pot ) {

	m_setiGenes.insert( iGene );
	Pot.m_setiGenes.erase( iGene ); }

void* CCoalesceClusterImpl::ThreadSeedPair( void* pData ) {
	SThreadSeedPair*		psData;
	double					dR, dP;
	size_t					i, j, iN, iExperiments, iPair, iPairs, iCur;
	pair<size_t, size_t>	priiSeed;

	psData = (SThreadSeedPair*)pData;
	iExperiments = psData->m_pPCL->GetExperiments( );
	psData->m_dMaxCorr = -( psData->m_dMinP = DBL_MAX );
	iPairs = psData->m_pPCL->GetGenes( ) * ( psData->m_pPCL->GetGenes( ) - 1 ) / 2;
	for( psData->m_iOne = psData->m_iTwo = 0,iPair = psData->m_iOffset; iPair < iPairs;
		iPair += psData->m_iStep ) {
		if( ( (float)rand( ) / RAND_MAX ) > psData->m_dFraction )
			continue;
		for( i = 0,iCur = iPair; iCur >= ( psData->m_pPCL->GetGenes( ) - 1 - i );
			iCur -= psData->m_pPCL->GetGenes( ) - 1 - i++ );
		j = i + iCur + 1;
		if( ( ( dR = CMeasurePearson::Pearson( psData->m_pPCL->Get( i ), iExperiments,
			psData->m_pPCL->Get( j ), iExperiments, IMeasure::EMapNone, NULL, NULL, &iN ) ) < 0 ) ||
			CMeta::IsNaN( dR ) )
			continue;
		if( ( ( dP = CStatistics::PValuePearson( dR, iN ) ) < psData->m_dMinP ) ||
			( !dP && ( dR > psData->m_dMaxCorr ) ) ) {
			priiSeed.first = i;
			priiSeed.second = j;
			if( psData->m_psetpriiSeeds->find( priiSeed ) == psData->m_psetpriiSeeds->end( ) ) {
				psData->m_dMaxCorr = dR;
				psData->m_dMinP = dP;
				psData->m_iOne = i;
				psData->m_iTwo = j; } } }

	return NULL; }

bool CCoalesceClusterImpl::AddSeedPair( const CPCL& PCL, CCoalesceCluster& Pot,
	std::set<std::pair<size_t, size_t> >& setpriiSeeds, float dFraction, float dPValue, size_t iThreads ) {
	size_t					i, iOne, iTwo;
	double					dMaxCorr, dMinP;
	pair<size_t, size_t>	priiSeed;
	vector<pthread_t>		vecpthdThreads;
	vector<SThreadSeedPair>	vecsThreads;

	if( PCL.GetGenes( ) < 2 ) {
		g_CatSleipnir.error( "CCoalesceClusterImpl::AddSeedPair( %g ) found no genes", dPValue );
		return false; }
	vecpthdThreads.resize( iThreads );
	vecsThreads.resize( vecpthdThreads.size( ) );
	for( i = 0; i < vecsThreads.size( ); ++i ) {
		vecsThreads[ i ].m_iOffset = i;
		vecsThreads[ i ].m_iStep = vecsThreads.size( );
		vecsThreads[ i ].m_pPCL = &PCL;
		vecsThreads[ i ].m_dFraction = dFraction;
		vecsThreads[ i ].m_psetpriiSeeds = &setpriiSeeds;
		if( pthread_create( &vecpthdThreads[ i ], NULL, ThreadSeedPair, &vecsThreads[ i ] ) ) {
			g_CatSleipnir.error( "CCoalesceClusterImpl::AddSeedPair( %g, %g, %d ) could not seed pair",
				dFraction, dPValue, iThreads );
			return false; } }
	dMaxCorr = -( dMinP = DBL_MAX );
	for( iOne = iTwo = i = 0; i < vecpthdThreads.size( ); ++i ) {
		pthread_join( vecpthdThreads[ i ], NULL );
		if( ( vecsThreads[ i ].m_dMinP < dMinP ) ||
			( !vecsThreads[ i ].m_dMinP && ( vecsThreads[ i ].m_dMaxCorr > dMaxCorr ) ) ) {
			dMinP = vecsThreads[ i ].m_dMinP;
			dMaxCorr = vecsThreads[ i ].m_dMaxCorr;
			iOne = vecsThreads[ i ].m_iOne;
			iTwo = vecsThreads[ i ].m_iTwo; } }
	if( ( dMinP * PCL.GetGenes( ) * ( PCL.GetGenes( ) - 1 ) * dFraction * dFraction ) < ( dPValue * 2 ) ) {
		g_CatSleipnir.info( "CCoalesceClusterImpl::AddSeedPair( %g, %g ) seeding: %s, %s, %g (p=%g)",
			dFraction, dPValue, PCL.GetGene( iOne ).c_str( ), PCL.GetGene( iTwo ).c_str( ), dMaxCorr, dMinP );
		priiSeed.first = iOne;
		priiSeed.second = iTwo;
		setpriiSeeds.insert( priiSeed );
		Add( iOne, Pot );
		Add( iTwo, Pot );
		return true; }

	g_CatSleipnir.notice( "CCoalesceClusterImpl::AddSeedPair( %g, %g ) inadequate seed pair: %s, %s, %g (p=%g)",
		dFraction, dPValue, PCL.GetGene( iOne ).c_str( ), PCL.GetGene( iTwo ).c_str( ), dMaxCorr, dMinP );
	return false; }

bool CCoalesceClusterImpl::AddCorrelatedGenes( const CPCL& PCL, CCoalesceCluster& Pot, float dPValue ) {
	size_t	iGene, iN;
	double	dR;

	CalculateCentroid( PCL );
	for( iGene = 0; iGene < PCL.GetGenes( ); ++iGene )
		if( !IsGene( iGene ) &&
			( ( dR = CMeasurePearson::Pearson( &m_vecdCentroid.front( ), PCL.GetExperiments( ),
			PCL.Get( iGene ), PCL.GetExperiments( ), IMeasure::EMapNone, NULL, NULL, &iN ) ) > 0 ) &&
			( ( CStatistics::PValuePearson( dR, iN ) * PCL.GetGenes( ) ) < dPValue ) )
			Add( iGene, Pot );

	return true; }

void CCoalesceCluster::CalculateHistograms( const CCoalesceGeneScores& GeneScores,
	CCoalesceGroupHistograms& HistogramsCluster, CCoalesceGroupHistograms* pHistogramsPot ) const {
	set<size_t>::const_iterator	iterGene;
	size_t						i;

	if( !GeneScores.GetMotifs( ) )
		return;
	for( iterGene = GetGenes( ).begin( ); iterGene != GetGenes( ).end( ); ++iterGene )
		if( !binary_search( m_veciPrevGenes.begin( ), m_veciPrevGenes.end( ), *iterGene ) ) {
			HistogramsCluster.Add( GeneScores, *iterGene, false );
			if( pHistogramsPot )
				pHistogramsPot->Add( GeneScores, *iterGene, true ); }
	for( i = 0; i < m_veciPrevGenes.size( ); ++i )
		if( GetGenes( ).find( m_veciPrevGenes[ i ] ) == GetGenes( ).end( ) ) {
			HistogramsCluster.Add( GeneScores, m_veciPrevGenes[ i ], true );
			if( pHistogramsPot )
				pHistogramsPot->Add( GeneScores, m_veciPrevGenes[ i ], false ); } }

void CCoalesceCluster::Subtract( CPCL& PCL, const CCoalesceCluster& Pot ) const {
	set<size_t>::const_iterator	iterGene, iterDataset;
	size_t						i, iCondition;
	float						d, dAve;

	for( iterGene = GetGenes( ).begin( ); iterGene != GetGenes( ).end( ); ++iterGene )
		for( iterDataset = m_setiDatasets.begin( ); iterDataset != m_setiDatasets.end( ); ++iterDataset )
			for( i = 0; i < GetConditions( *iterDataset ); ++i ) {
				iCondition = GetCondition( *iterDataset, i );
				if( !CMeta::IsNaN( d = m_vecdCentroid[ iCondition ] ) ) {
					if( CMeta::IsNaN( dAve = Pot.m_vecdCentroid[ iCondition ] ) )
						dAve = 0;
					else
						dAve = ( ( GetGenes( ).size( ) * d ) + ( Pot.GetGenes( ).size( ) * dAve ) ) /
							( GetGenes( ).size( ) + Pot.GetGenes( ).size( ) );
					PCL.Get( *iterGene, GetCondition( *iterDataset, i ) ) -= d - dAve; } } }

void CCoalesceCluster::Subtract( CCoalesceGeneScores& GeneScores ) const {
	set<size_t>::const_iterator			iterGene;
	set<SMotifMatch>::const_iterator	iterMotif;

	for( iterGene = GetGenes( ).begin( ); iterGene != GetGenes( ).end( ); ++iterGene )
		for( iterMotif = m_setsMotifs.begin( ); iterMotif != m_setsMotifs.end( ); ++iterMotif )
			GeneScores.Subtract( *iterMotif, *iterGene ); }

void CCoalesceClusterImpl::CalculateCentroid( const CPCL& PCL ) {
	set<size_t>::const_iterator	iterGene;
	size_t						i, j;
	float						d;

	m_veciCounts.resize( PCL.GetExperiments( ) );
	fill( m_veciCounts.begin( ), m_veciCounts.end( ), 0 );
	m_vecdCentroid.resize( PCL.GetExperiments( ) );
	fill( m_vecdCentroid.begin( ), m_vecdCentroid.end( ), 0.0f );
	m_vecdStdevs.resize( PCL.GetExperiments( ) );
	fill( m_vecdStdevs.begin( ), m_vecdStdevs.end( ), 0.0f );
	for( iterGene = GetGenes( ).begin( ); iterGene != GetGenes( ).end( ); ++iterGene )
		for( i = 0; i < m_veciCounts.size( ); ++i )
			if( !CMeta::IsNaN( d = PCL.Get( *iterGene, i ) ) ) {
				m_veciCounts[ i ]++;
				m_vecdCentroid[ i ] += d;
				m_vecdStdevs[ i ] += d * d; }
	for( i = 0; i < m_veciCounts.size( ); ++i ) {
		if( j = m_veciCounts[ i ] ) {
			m_vecdCentroid[ i ] /= j;
			m_vecdStdevs[ i ] = ( m_vecdStdevs[ i ] / ( ( j == 1 ) ? 1 : ( j - 1 ) ) ) -
				( m_vecdCentroid[ i ] * m_vecdCentroid[ i ] );
			m_vecdStdevs[ i ] = ( m_vecdStdevs[ i ] < 0 ) ? 0 : sqrt( m_vecdStdevs[ i ] ); }
		else
			m_vecdCentroid[ i ] = CMeta::GetNaN( );
		g_CatSleipnir.info( "CCoalesceClusterImpl::CalculateCentroid( ) condition %d: count %d, mean %g, stdev %g",
			i, m_veciCounts[ i ], m_vecdCentroid[ i ], m_vecdStdevs[ i ] ); }

	for( i = 0; i < m_vecsDatasets.size( ); ++i ) {
		SDataset&	sDataset	= m_vecsDatasets[ i ];

		if( sDataset.GetConditions( ) == 1 )
			continue;
		sDataset.m_vecdCentroid.resize( sDataset.GetConditions( ) );
		for( j = 0; j < sDataset.GetConditions( ); ++j )
			sDataset.m_vecdCentroid[ j ] = m_vecdCentroid[ sDataset.GetCondition( j ) ]; } }

void* CCoalesceClusterImpl::ThreadSelectCondition( void* pData ) {
	vector<float>			vecdDatasetCluster, vecdDatasetPot;
	vector<size_t>			veciDatasetCluster, veciDatasetPot;
	size_t					iDataset, iCondition, iCluster, iPot, iGene;
	double					dP, dZ;
	float					d;
	SThreadSelectCondition*	psData;
	float*					adCluster;
	float*					adPot;

	psData = (SThreadSelectCondition*)pData;
	const vector<size_t>&	veciPot		= *psData->m_pveciPot;
	const vector<size_t>&	veciCluster	= *psData->m_pveciCluster;
	const CPCL&				PCL			= *psData->m_pPCL;

	adCluster = new float[ veciCluster.size( ) ];
	adPot = new float[ veciPot.size( ) ];
	for( iDataset = psData->m_iOffset; iDataset < psData->m_pvecsDatasets->size( );
		iDataset += psData->m_iStep ) {
		SDataset&	sDataset	= (*psData->m_pvecsDatasets)[ iDataset ];

		if( sDataset.GetConditions( ) == 1 ) {
			iCondition = sDataset.GetCondition( 0 );
			for( iCluster = iPot = iGene = 0; iGene < veciCluster.size( ); ++iGene )
				if( !CMeta::IsNaN( d = PCL.Get( veciCluster[ iGene ], iCondition ) ) )
					adCluster[ iCluster++ ] = d;
			for( iGene = 0; iGene < veciPot.size( ); ++iGene )
				if( !CMeta::IsNaN( d = PCL.Get( veciPot[ iGene ], iCondition ) ) )
					adPot[ iPot++ ] = d;
			if( !( iCluster || iPot ) )
				continue;
			if( !( iCluster && iPot ) ) {
				dZ = FLT_MAX;
				dP = 0; }
			else {
				dZ = CStatistics::ZScore( adCluster, adCluster + iCluster, adPot, adPot + iPot );
				dP = CStatistics::ZTest( dZ, min( iCluster, iPot ) ) * psData->m_pvecsDatasets->size( ); }
			if( dP < psData->m_dPValue ) {
				g_CatSleipnir.info( "CCoalesceClusterImpl::ThreadSelectCondition( %g ) selected condition %d at %g, z=%g",
					psData->m_dPValue, iCondition, dP, dZ );
				(*psData->m_pvecfSignificant)[ iDataset ] = true;
				sDataset.m_dZ = (float)min( dZ, (double)FLT_MAX ); }
			else
				g_CatSleipnir.debug( "CCoalesceClusterImpl::ThreadSelectCondition( %g ) rejected condition %d at %g, z=%g",
					psData->m_dPValue, iCondition, dP, dZ ); }
		else {
			vecdDatasetCluster.resize( sDataset.GetConditions( ) );
			fill( vecdDatasetCluster.begin( ), vecdDatasetCluster.end( ), 0.0f );
			vecdDatasetPot.resize( sDataset.GetConditions( ) );
			fill( vecdDatasetPot.begin( ), vecdDatasetPot.end( ), 0.0f );
			veciDatasetCluster.resize( sDataset.GetConditions( ) );
			fill( veciDatasetCluster.begin( ), veciDatasetCluster.end( ), 0 );
			veciDatasetPot.resize( sDataset.GetConditions( ) );
			fill( veciDatasetPot.begin( ), veciDatasetPot.end( ), 0 );
			for( iGene = 0; iGene < veciCluster.size( ); ++iGene )
				for( iCondition = 0; iCondition < sDataset.GetConditions( ); ++iCondition )
					if( !CMeta::IsNaN( d = PCL.Get( veciCluster[ iGene ],
						sDataset.GetCondition( iCondition ) ) ) ) {
						veciDatasetCluster[ iCondition ]++;
						vecdDatasetCluster[ iCondition ] += d; }
			for( iGene = 0; iGene < veciPot.size( ); ++iGene )
				for( iCondition = 0; iCondition < sDataset.GetConditions( ); ++iCondition )
					if( !CMeta::IsNaN( d = PCL.Get( veciPot[ iGene ],
						sDataset.GetCondition( iCondition ) ) ) ) {
						veciDatasetPot[ iCondition ]++;
						vecdDatasetPot[ iCondition ] += d; }
			for( iCondition = 0; iCondition < vecdDatasetCluster.size( ); ++iCondition ) {
				if( veciDatasetCluster[ iCondition ] )
					vecdDatasetCluster[ iCondition ] /= veciDatasetCluster[ iCondition ];
				if( veciDatasetPot[ iCondition ] )
					vecdDatasetPot[ iCondition ] /= veciDatasetPot[ iCondition ]; }
			dP = CStatistics::MultivariateNormalCDF( vecdDatasetCluster, vecdDatasetPot,
				sDataset.m_psDataset->m_MatSigmaChol ) * psData->m_pvecsDatasets->size( );
			if( dP > 0.5 )
				dP = 1 - dP;
			if( dP < psData->m_dPValue ) {
				dZ = 0;
				for( iCondition = 0; iCondition < sDataset.GetConditions( ); ++iCondition )
					if( sDataset.m_psDataset->m_vecdStdevs[ iCondition ] ) {
						d = ( vecdDatasetCluster[ iCondition ] - vecdDatasetPot[ iCondition ] ) /
							sDataset.m_psDataset->m_vecdStdevs[ iCondition ];
						dZ += d * d; }
				dZ = sqrt( dZ );
				g_CatSleipnir.info( "CCoalesceClusterImpl::ThreadSelectCondition( %g ) selected dataset %d at %g, z=%g",
					psData->m_dPValue, iDataset, dP, dZ );
				(*psData->m_pvecfSignificant)[ iDataset ] = true;
				sDataset.m_dZ = (float)min( dZ, (double)FLT_MAX ); } } }
	delete[] adCluster;
	delete[] adPot;

	return NULL; }

bool CCoalesceCluster::SelectConditions( const CPCL& PCL, const CCoalesceCluster& Pot, size_t iThreads,
	float dPValue ) {
	vector<bool>					vecfSignificant;
	vector<pthread_t>				vecpthdThreads;
	vector<SThreadSelectCondition>	vecsThreads;
	size_t							i;
	vector<size_t>					veciCluster, veciPot;

	m_setiDatasets.clear( );
	veciCluster.resize( GetGenes( ).size( ) );
	copy( GetGenes( ).begin( ), GetGenes( ).end( ), veciCluster.begin( ) );
	veciPot.resize( Pot.GetGenes( ).size( ) );
	copy( Pot.GetGenes( ).begin( ), Pot.GetGenes( ).end( ), veciPot.begin( ) );
	vecfSignificant.resize( m_vecsDatasets.size( ) );
	vecpthdThreads.resize( iThreads );
	vecsThreads.resize( vecpthdThreads.size( ) );
	for( i = 0; i < vecsThreads.size( ); ++i ) {
		vecsThreads[ i ].m_iOffset = i;
		vecsThreads[ i ].m_iStep = vecsThreads.size( );
		vecsThreads[ i ].m_pveciCluster = &veciCluster;
		vecsThreads[ i ].m_pveciPot = &veciPot;
		vecsThreads[ i ].m_pvecsDatasets = &m_vecsDatasets;
		vecsThreads[ i ].m_pPCL = &PCL;
		vecsThreads[ i ].m_dPValue = dPValue;
		vecsThreads[ i ].m_pvecfSignificant = &vecfSignificant;
		if( pthread_create( &vecpthdThreads[ i ], NULL, ThreadSelectCondition, &vecsThreads[ i ] ) ) {
			g_CatSleipnir.error( "CCoalesceCluster::SelectConditions( %d, %g ) could not select conditions",
				iThreads, dPValue );
			return false; } }
	for( i = 0; i < vecpthdThreads.size( ); ++i )
		pthread_join( vecpthdThreads[ i ], NULL );
	for( i = 0; i < vecfSignificant.size( ); ++i )
		if( vecfSignificant[ i ] )
			m_setiDatasets.insert( i );

	return true; }

void* CCoalesceClusterImpl::ThreadSelectMotif( void* pData ) {
	SThreadSelectMotif*	psData;
	uint32_t			iMotif;

	psData = (SThreadSelectMotif*)pData;
	for( iMotif = psData->m_iOffset; iMotif < psData->m_pMotifs->GetMotifs( ); iMotif += psData->m_iStep )
		if( !AddSignificant( *psData->m_pMotifs, iMotif, *psData->m_pHistsCluster, *psData->m_pHistsPot,
			psData->m_dPValue, psData->m_vecsMotifs ) )
			break;

	return NULL; }

bool CCoalesceCluster::SelectMotifs( const CCoalesceGroupHistograms& HistsCluster,
	const CCoalesceGroupHistograms& HistsPot, float dPValue, size_t iThreads,
	const CCoalesceMotifLibrary* pMotifs ) {
	uint32_t					i, j;
	vector<pthread_t>			vecpthdThreads;
	vector<SThreadSelectMotif>	vecsThreads;

	m_setsMotifs.clear( );
	if( !pMotifs )
		return true;
	vecpthdThreads.resize( iThreads );
	vecsThreads.resize( vecpthdThreads.size( ) );
	for( i = 0; i < vecsThreads.size( ); ++i ) {
		vecsThreads[ i ].m_iOffset = i;
		vecsThreads[ i ].m_iStep = vecsThreads.size( );
		vecsThreads[ i ].m_pMotifs = pMotifs;
		vecsThreads[ i ].m_pHistsCluster = &HistsCluster;
		vecsThreads[ i ].m_pHistsPot = &HistsPot;
		vecsThreads[ i ].m_dPValue = dPValue;
		if( pthread_create( &vecpthdThreads[ i ], NULL, ThreadSelectMotif, &vecsThreads[ i ] ) ) {
			g_CatSleipnir.error( "CCoalesceCluster::SelectMotifs( %g, %d ) could not select motifs",
				dPValue, iThreads );
			return false; } }
	for( i = 0; i < vecpthdThreads.size( ); ++i ) {
		pthread_join( vecpthdThreads[ i ], NULL );
		for( j = 0; j < vecsThreads[ i ].m_vecsMotifs.size( ); ++j )
			m_setsMotifs.insert( vecsThreads[ i ].m_vecsMotifs[ j ] ); }

	return true; }

bool CCoalesceClusterImpl::AddSignificant( const CCoalesceMotifLibrary& Motifs, uint32_t iMotif,
	const CCoalesceGroupHistograms& HistsCluster, const CCoalesceGroupHistograms& HistsPot, float dPValue,
	vector<SMotifMatch>& vecsMotifs ) {
	size_t									iTypeCluster, iTypePot;
	double									dP, dAveOne, dAverage, dZ;
	CCoalesceSequencerBase::ESubsequence	eSubsequence;

	dPValue /= Motifs.GetMotifs( );
	for( iTypeCluster = 0; iTypeCluster < HistsCluster.GetTypes( ); ++iTypeCluster ) {
		const string&	strTypeCluster	= HistsCluster.GetType( iTypeCluster );

		if( ( iTypePot = HistsPot.GetType( strTypeCluster ) ) == -1 )
			continue;
		for( eSubsequence = CCoalesceSequencerBase::ESubsequenceBegin;
			(size_t)eSubsequence < HistsCluster.GetSubsequences( iTypeCluster );
			eSubsequence = (CCoalesceSequencerBase::ESubsequence)( (size_t)eSubsequence + 1 ) ) {
			const CCoalesceHistogramSet<>&	HistSetCluster	= HistsCluster.Get( iTypeCluster, eSubsequence );
			const CCoalesceHistogramSet<>&	HistSetPot		= HistsPot.Get( iTypePot, eSubsequence );

			if( ( HistSetCluster.GetMembers( ) <= iMotif ) ||
				( HistSetPot.GetMembers( ) <= iMotif ) ||
				!( HistSetCluster.GetTotal( ) && HistSetPot.GetTotal( ) ) )
				continue;
			dP = HistSetCluster.ZScore( iMotif, HistSetPot, dAveOne, dAverage, dZ );
			if( dP < dPValue ) {
				SMotifMatch	sMotif( iMotif, strTypeCluster, eSubsequence, (float)dZ, (float)( dAveOne -
					dAverage ) );

				if( g_CatSleipnir.isInfoEnabled( ) ) {
					ostringstream	ossm;

					ossm << "CCoalesceClusterImpl::AddSignificant( " << iMotif << ", " << dPValue <<
						" ) adding at " << dP << ":" << endl << sMotif.Save( &Motifs ) << endl <<
						"Cluster	" << HistSetCluster.Save( iMotif ) << endl <<
						"Pot	" << HistSetPot.Save( iMotif );
					g_CatSleipnir.info( ossm.str( ) ); }
				vecsMotifs.push_back( sMotif ); }
			else if( g_CatSleipnir.isDebugEnabled( ) ) {
				ostringstream	ossm;

				ossm << "CCoalesceClusterImpl::AddSignificant( " << iMotif << ", " << dPValue <<
					" ) failed at " << dP << ":" << endl << SMotifMatch( iMotif, strTypeCluster, eSubsequence,
					(float)dZ, (float)( dAveOne - dAverage ) ).Save( &Motifs ) << endl <<
					"Cluster	" << HistSetCluster.Save( iMotif ) << endl <<
					"Pot	" << HistSetPot.Save( iMotif );
				g_CatSleipnir.debug( ossm.str( ) ); } } }

	return true; }

void* CCoalesceClusterImpl::ThreadCentroid( void* pData ) {
	SThreadCentroid*	psData;

	psData = (SThreadCentroid*)pData;
	psData->m_pCluster->CalculateCentroid( psData->m_PCL );

	return NULL; }

void* CCoalesceClusterImpl::ThreadSignificantGene( void* pData ) {
	SThreadSignificantGene*	psData;
	size_t					i;

	psData = (SThreadSignificantGene*)pData;
	for( i = psData->m_iOffset; i < psData->m_pvecfSignificant->size( ); i += psData->m_iStep )
		(*psData->m_pvecfSignificant)[ i ] = psData->m_pCluster->IsSignificant( i, *psData->m_pPCL,
			psData->m_pMotifs, *psData->m_pGeneScores, *psData->m_pHistsCluster, *psData->m_pHistsPot,
			*psData->m_pPot, *psData->m_pveciDatasets, psData->m_dBeta, psData->m_iMinimum,
			psData->m_dProbability );

	return NULL; }

bool CCoalesceCluster::SelectGenes( const CPCL& PCL, const CCoalesceGeneScores& GeneScores,
	const CCoalesceGroupHistograms& HistsCluster, const CCoalesceGroupHistograms& HistsPot, size_t iMinimum,
	size_t iThreads, CCoalesceCluster& Pot, float dProbability, const CCoalesceMotifLibrary* pMotifs ) {
	size_t								i;
	vector<pthread_t>					vecpthdThreads;
	vector<bool>						vecfSignificant;
	vector<SThreadSignificantGene>		vecsThreads;
	vector<size_t>						veciDatasets;
	float								dFracExpression, dFracMotifs;
	set<size_t>							setiMotifs;
	set<SMotifMatch>::const_iterator	iterMotif;

//vecpthdThreads.resize( 1 );
	vecpthdThreads.resize( max( 2, iThreads ) );
	{
		SThreadCentroid	sUs( this, PCL ), sThem( &Pot, PCL );

		if( pthread_create( &vecpthdThreads[ 0 ], NULL, ThreadCentroid, &sUs ) ||
			pthread_create( &vecpthdThreads[ 1 ], NULL, ThreadCentroid, &sThem ) ) {
			g_CatSleipnir.error( "CCoalesceCluster::SelectGenes( %d, %g ) could not calculate centroids",
				iThreads, dProbability );
			return false; }
		for( i = 0; i < 2; ++i )
			pthread_join( vecpthdThreads[ i ], NULL );
	}

	dFracExpression = (float)m_setiDatasets.size( ) / m_vecsDatasets.size( );
	if( pMotifs ) {
		for( iterMotif = m_setsMotifs.begin( ); iterMotif != m_setsMotifs.end( ); ++iterMotif )
			setiMotifs.insert( iterMotif->m_iMotif );
		dFracMotifs = (float)setiMotifs.size( ) / pMotifs->GetMotifs( ); }
	else
		dFracMotifs = 0;

	veciDatasets.resize( m_setiDatasets.size( ) );
	copy( m_setiDatasets.begin( ), m_setiDatasets.end( ), veciDatasets.begin( ) );
	vecfSignificant.resize( PCL.GetGenes( ) );
	vecsThreads.resize( vecpthdThreads.size( ) );
	for( i = 0; i < vecpthdThreads.size( ); ++i ) {
		vecsThreads[ i ].m_iOffset = i;
		vecsThreads[ i ].m_iStep = vecpthdThreads.size( );
		vecsThreads[ i ].m_pvecfSignificant = &vecfSignificant;
		vecsThreads[ i ].m_pPCL = &PCL;
		vecsThreads[ i ].m_pMotifs = pMotifs;
		vecsThreads[ i ].m_pGeneScores = &GeneScores;
		vecsThreads[ i ].m_pHistsCluster = &HistsCluster;
		vecsThreads[ i ].m_pHistsPot = &HistsPot;
		vecsThreads[ i ].m_pCluster = this;
		vecsThreads[ i ].m_pPot = &Pot;
		vecsThreads[ i ].m_pveciDatasets = &veciDatasets;
		vecsThreads[ i ].m_dBeta = dFracExpression / ( dFracExpression + dFracMotifs );
		vecsThreads[ i ].m_iMinimum = iMinimum;
		vecsThreads[ i ].m_dProbability = dProbability;
		if( pthread_create( &vecpthdThreads[ i ], NULL, ThreadSignificantGene, &vecsThreads[ i ] ) ) {
			g_CatSleipnir.error( "CCoalesceCluster::SelectGenes( %d, %g ) could not calculate significance",
				iThreads, dProbability );
			return false; } }
	for( i = 0; i < vecpthdThreads.size( ); ++i )
		pthread_join( vecpthdThreads[ i ], NULL );
	m_setiGenes.clear( );
	Pot.m_setiGenes.clear( );
	for( i = 0; i < vecfSignificant.size( ); ++i )
		if( vecfSignificant[ i ] ) {
			m_setiGenes.insert( i );
			m_vecdPriors[ i ] = 1; }
		else {
			Pot.m_setiGenes.insert( i );
			m_vecdPriors[ i ] = dProbability; }

	return true; }

bool CCoalesceClusterImpl::IsSignificant( size_t iGene, const CPCL& PCL, const CCoalesceMotifLibrary* pMotifs,
	const CCoalesceGeneScores& GeneScores, const CCoalesceGroupHistograms& HistsCluster,
	const CCoalesceGroupHistograms& HistsPot, const CCoalesceCluster& Pot,
	const vector<size_t>& veciDatasets, float dBeta, size_t iMinimum, float dProbability ) const {
	float	dP, dLogPMotifsGivenIn, dLogPMotifsGivenOut, dLogPExpressionGivenIn, dLogPExpressionGivenOut;
	bool	fClustered;

	fClustered = IsGene( iGene );
// We want P(g in C|S,E) =
// P(S,E|g in C)P(g in C)/P(S,E) =
// P(S|g in C)P(E|g in C)P(g in C)/(P(S,E|g in C) + P(S,E|g notin C)) =
// P(S|g in C)P(E|g in C)P(g in C)/(P(S|g in C)P(E|g in C) + P(S|g notin C)P(E|g notin C))
	if( !( CalculateProbabilityExpression( iGene, PCL, Pot, veciDatasets, fClustered, dLogPExpressionGivenIn,
		dLogPExpressionGivenOut ) && CalculateProbabilityMotifs( GeneScores, iGene, HistsCluster, HistsPot,
		fClustered, iMinimum, dLogPMotifsGivenIn, dLogPMotifsGivenOut ) ) )
		return false;
	dP = m_vecdPriors[ iGene ] / ( 1 + exp( ( dBeta * ( dLogPExpressionGivenOut - dLogPExpressionGivenIn ) +
		( 1 - dBeta ) * ( dLogPMotifsGivenOut - dLogPMotifsGivenIn ) ) /
		( 0.5f + 2 * pow( 0.5f - dBeta, 2 ) ) ) );
//	dP = m_vecdPriors[ iGene ] / ( 1 + exp( dLogPExpressionGivenOut - dLogPExpressionGivenIn +
//		dLogPMotifsGivenOut - dLogPMotifsGivenIn ) );

	g_CatSleipnir.debug( "CCoalesceClusterImpl::IsSignificant( %s ) is %g beta %g, exp. p=%g vs. %g, seq. p=%g vs %g",
		PCL.GetGene( iGene ).c_str( ), dP, dBeta, dLogPExpressionGivenIn, dLogPExpressionGivenOut,
		dLogPMotifsGivenIn, dLogPMotifsGivenOut );
	if( g_CatSleipnir.isDebugEnabled( ) ) {
		set<SMotifMatch>::const_iterator	iterMotif;
		size_t								iType;

		for( iterMotif = m_setsMotifs.begin( ); iterMotif != m_setsMotifs.end( ); ++iterMotif )
			if( ( iType = GeneScores.GetType( iterMotif->m_strType ) ) != -1 )
				g_CatSleipnir.debug( "%g	%s", GeneScores.Get( iType, iterMotif->m_eSubsequence, iGene,
					iterMotif->m_iMotif ), iterMotif->Save( pMotifs ).c_str( ) ); }

	return ( dP > dProbability ); }

bool CCoalesceClusterImpl::CalculateProbabilityExpression( size_t iGene, const CPCL& PCL,
	const CCoalesceCluster& Pot, const vector<size_t>& veciDatasets, bool fClustered, float& dLogPIn,
	float& dLogPOut ) const {
	static const double	c_dEpsilonZero	= 1e-10;
	float				dGene, dCluster, dPot;
	double				dPCluster, dPPot;
	size_t				iDataset, iPot, iCondition;
	long double			dPIn, dPOut;

	iPot = PCL.GetGenes( ) - GetGenes( ).size( );
	dLogPIn = dLogPOut = 0;
	dPIn = dPOut = 1;
	for( iDataset = 0; iDataset < veciDatasets.size( ); ++iDataset ) {
		const SDataset&	sDataset	= m_vecsDatasets[ veciDatasets[ iDataset ] ];

		if( sDataset.GetConditions( ) == 1 ) {
			iCondition = sDataset.GetCondition( 0 );
			if( CMeta::IsNaN( dGene = PCL.Get( iGene, iCondition ) ) ||
				CMeta::IsNaN( dCluster = m_vecdCentroid[ iCondition ] ) ||
				CMeta::IsNaN( dPot = Pot.m_vecdCentroid[ iCondition ] ) )
				continue;
			
			if( fClustered )
				dCluster = ( ( dCluster * GetGenes( ).size( ) ) - dGene ) /
					( GetGenes( ).size( ) - 1 );
			else
				dPot = ( ( dPot * iPot ) - dGene ) / ( iPot - 1 );
			dPCluster = max( c_dEpsilonZero, CStatistics::NormalPDF( dGene, dCluster, max( c_dEpsilonZero,
				(double)m_vecdStdevs[ iCondition ] ) ) );
			dPPot = max( c_dEpsilonZero, CStatistics::NormalPDF( dGene, dPot, max( c_dEpsilonZero,
				(double)Pot.m_vecdStdevs[ iCondition ] ) ) );
			AdjustProbabilities( sDataset.m_dZ, dPCluster, dPPot ); }
		else {
			vector<float>	vecdGene;

			vecdGene.resize( sDataset.GetConditions( ) );
			for( iCondition = 0; iCondition < sDataset.GetConditions( ); ++iCondition )
				vecdGene[ iCondition ] = PCL.Get( iGene, sDataset.GetCondition( iCondition ) );
			dPCluster = max( c_dEpsilonZero, CStatistics::MultivariateNormalPDF( vecdGene,
				sDataset.m_vecdCentroid, sDataset.m_psDataset->m_dSigmaDetSqrt,
				sDataset.m_psDataset->m_MatSigmaInv ) );
			dPPot = max( c_dEpsilonZero, CStatistics::MultivariateNormalPDF( vecdGene,
				Pot.m_vecsDatasets[ veciDatasets[ iDataset ] ].m_vecdCentroid,
				sDataset.m_psDataset->m_dSigmaDetSqrt, sDataset.m_psDataset->m_MatSigmaInv ) );
			AdjustProbabilities( sDataset.m_dZ, dPCluster, dPPot ); }
		dPIn *= dPCluster;
		dPOut *= dPPot;
		if( ( dPIn < DBL_MIN ) || ( dPOut < DBL_MIN ) ) {
			dLogPIn += (float)log( dPIn );
			dLogPOut += (float)log( dPOut );
			dPIn = dPOut = 1; } }
	dLogPIn += (float)log( dPIn );
	dLogPOut += (float)log( dPOut );

	return true; }

bool CCoalesceClusterImpl::CalculateProbabilityMotifs( const CCoalesceGeneScores& GeneScores, size_t iGene,
	const CCoalesceGroupHistograms& HistsCluster, const CCoalesceGroupHistograms& HistsPot,
	bool fClustered, size_t iMinimum, float& dLogPIn, float& dLogPOut ) const {
	set<SMotifMatch>::const_iterator	iterMotif;
	size_t								iType, iCluster, iPot;
	unsigned short						sCluster, sPot;
	double								dPCluster, dPPot;
	float								dGene;
	const float*						adValues;
	long double							dPIn, dPOut;

	dLogPIn = dLogPOut = 0;
	if( m_setsMotifs.empty( ) )
		return true;

// TODO: this is the slowest part of the entire algorithm
	dPIn = dPOut = 1;
	for( iterMotif = m_setsMotifs.begin( ); iterMotif != m_setsMotifs.end( ); ++iterMotif ) {
		if( ( ( iType = GeneScores.GetType( iterMotif->m_strType ) ) == -1 ) ||
			!( adValues = GeneScores.Get( iType, iterMotif->m_eSubsequence, iGene ) ) ||
			( HistsCluster.GetType( iterMotif->m_strType ) == -1 ) ||
			( HistsPot.GetType( iterMotif->m_strType ) == -1 ) )
			continue;
		dGene = adValues[ iterMotif->m_iMotif ];
		{
			const CCoalesceHistogramSet<>&	HistSetCluster	= HistsCluster.Get( iterMotif->m_strType,
																iterMotif->m_eSubsequence );
			const CCoalesceHistogramSet<>&	HistSetPot		= HistsPot.Get( iterMotif->m_strType,
																iterMotif->m_eSubsequence );

			sCluster = HistSetCluster.Get( iterMotif->m_iMotif, dGene );
			iCluster = HistSetCluster.GetTotal( );
			sPot = HistSetPot.Get( iterMotif->m_iMotif, dGene );
			iPot = HistSetPot.GetTotal( );
			if( ( iCluster < iMinimum ) || ( iPot < iMinimum ) )
				continue;

			if( fClustered ) {
				if( sCluster ) {
					sCluster--;
					iCluster--; }
				else
					g_CatSleipnir.warn( "CCoalesceClusterImpl::CalculateProbabilityMotifs( %d, %d, %g, %g ) no motifs of %d in cluster: %g, %s",
						iGene, fClustered, dLogPIn, dLogPOut, iCluster, dGene,
						iterMotif->Save( NULL ).c_str( ) ); }
			else if( sPot ) {
				sPot--;
				iPot--; }
			else
				g_CatSleipnir.warn( "CCoalesceClusterImpl::CalculateProbabilityMotifs( %d, %d, %g, %g ) no motifs of %d in pot: %g, %s",
					iGene, fClustered, dLogPIn, dLogPOut, iPot, dGene, iterMotif->Save( NULL ).c_str( ) );
			dPCluster = (double)( sCluster + 1 ) / ( iCluster + HistSetCluster.GetEdges( ) );
			dPPot = (double)( sPot + 1 ) / ( iPot + HistSetPot.GetEdges( ) );
		}
		AdjustProbabilities( iterMotif->m_dZ, dPCluster, dPPot );
		dPIn *= dPCluster;
		dPOut *= dPPot;
		if( ( dPIn < DBL_MIN ) || ( dPOut < DBL_MIN ) ) {
			dLogPIn += (float)log( dPIn );
			dLogPOut += (float)log( dPOut );
			dPIn = dPOut = 1; } }
	dLogPIn += (float)log( dPIn );
	dLogPOut += (float)log( dPOut );

	return true; }

bool CCoalesceCluster::Save( const std::string& strDirectory, size_t iID, const CPCL& PCL,
	const CCoalesceMotifLibrary* pMotifs ) const {
	ofstream							ofsm;
	char*								szTmp;
	CPCL								PCLCopy;
	size_t								iGeneFrom, iGeneTo, iExpFrom, iExpTo, iLength;
	string								strBase;
	set<SMotifMatch>::const_iterator	iterMotif;
	set<size_t>							setiConditions;

	GetConditions( setiConditions );
	szTmp = new char[ iLength = ( strDirectory.length( ) + 16 ) ];
	sprintf_s( szTmp, iLength - 1, "%s/c%04d_XXXXXX", strDirectory.empty( ) ? "." : strDirectory.c_str( ),
		iID );
	_mktemp_s( szTmp, iLength - 1 );
	szTmp[ iLength - 1 ] = 0;
	strBase = szTmp;
	delete[] szTmp;

	ofsm.open( ( strBase + c_szMotifs ).c_str( ) );
	if( !ofsm.is_open( ) )
		return false;
	for( iterMotif = m_setsMotifs.begin( ); iterMotif != m_setsMotifs.end( ); ++iterMotif )
		ofsm << iterMotif->Save( pMotifs ) << endl;
	ofsm.close( );

	PCLCopy.Open( PCL );
// Reorder header
	for( iExpTo = iExpFrom = 0; iExpFrom < PCL.GetExperiments( ); ++iExpFrom )
		if( setiConditions.find( iExpFrom ) != setiConditions.end( ) )
			PCLCopy.SetExperiment( iExpTo++, c_cStar + PCL.GetExperiment( iExpFrom ) );
	for( iExpFrom = 0; iExpFrom < PCL.GetExperiments( ); ++iExpFrom )
		if( setiConditions.find( iExpFrom ) == setiConditions.end( ) )
			PCLCopy.SetExperiment( iExpTo++, PCL.GetExperiment( iExpFrom ) );
// Reorder genes
	for( iGeneTo = iGeneFrom = 0; iGeneFrom < PCL.GetGenes( ); ++iGeneFrom )
		if( IsGene( iGeneFrom ) && !SaveCopy( PCL, setiConditions, iGeneFrom, PCLCopy, iGeneTo++, false ) )
			return false;
	for( iGeneFrom = 0; iGeneFrom < PCL.GetGenes( ); ++iGeneFrom )
		if( !IsGene( iGeneFrom ) )
			PCLCopy.MaskGene( iGeneTo++ );

	ofsm.clear( );
	ofsm.open( ( strBase + CPCL::GetExtension( ) ).c_str( ) );
	if( !ofsm.is_open( ) )
		return false;
	PCLCopy.Save( ofsm );
	ofsm.close( );

	return true; }

bool CCoalesceClusterImpl::SaveCopy( const CPCL& PCLFrom, const set<size_t>& setiConditions, size_t iGeneFrom,
	CPCL& PCLTo, size_t iGeneTo, bool fClustered ) const {
	size_t	i, iExpTo, iExpFrom;

	PCLTo.SetGene( iGeneTo, PCLFrom.GetGene( iGeneFrom ) );
	for( i = 1; i < PCLFrom.GetFeatures( ); ++i )
		PCLTo.SetFeature( iGeneTo, i, (string)( ( fClustered && ( i == 1 ) ) ? "*" : "" ) +
			PCLFrom.GetFeature( iGeneFrom, i ) );
	for( iExpTo = iExpFrom = 0; iExpFrom < PCLFrom.GetExperiments( ); ++iExpFrom )
		if( setiConditions.find( iExpFrom ) != setiConditions.end( ) )
			PCLTo.Set( iGeneTo, iExpTo++, PCLFrom.Get( iGeneFrom, iExpFrom ) );
	for( iExpFrom = 0; iExpFrom < PCLFrom.GetExperiments( ); ++iExpFrom )
		if( setiConditions.find( iExpFrom ) == setiConditions.end( ) )
			PCLTo.Set( iGeneTo, iExpTo++, PCLFrom.Get( iGeneFrom, iExpFrom ) );

	return true; }

void CCoalesceCluster::Save( std::ostream& ostm, size_t iID, const CPCL& PCL,
	const CCoalesceMotifLibrary* pMotifs, float dCutoffPWMs, float dPenaltyGap, float dPenaltyMismatch,
	bool fNoRCs ) const {
	set<size_t>::const_iterator			iterID;
	set<SMotifMatch>::const_iterator	iterMotif;
	size_t								i;
	string								strMotif;

	ostm << "Cluster\t" << iID << endl;
	ostm << "Genes";
	for( iterID = GetGenes( ).begin( ); iterID != GetGenes( ).end( ); ++iterID )
		ostm << '\t' << PCL.GetGene( *iterID );
	ostm << endl << "Conditions";
	for( iterID = m_setiDatasets.begin( ); iterID != m_setiDatasets.end( ); ++iterID ) 
		for( i = 0; i < GetConditions( *iterID ); ++i )
			ostm << '\t' << PCL.GetExperiment( GetCondition( *iterID, i ) );
	ostm << endl << "Motifs" << endl;
	for( iterMotif = GetMotifs( ).begin( ); iterMotif != GetMotifs( ).end( ); ++iterMotif )
		if( !( strMotif = iterMotif->Save( pMotifs, true, dCutoffPWMs, dPenaltyGap, dPenaltyMismatch,
			fNoRCs ) ).empty( ) )
			ostm << strMotif << endl; }

size_t CCoalesceCluster::Open( const string& strPCL, size_t iSkip, const CPCL& PCL,
	CCoalesceMotifLibrary* pMotifs ) {
	CPCL		PCLCluster;
	size_t		i, j, iGene;
	string		strMotifs;
	ifstream	ifsm;

	Clear( );
	if( !PCLCluster.Open( strPCL.c_str( ), iSkip ) )
		return -1;

	for( i = 0; i < PCLCluster.GetExperiments( ); ++i ) {
		if( PCLCluster.GetExperiment( i )[ 0 ] != c_cStar )
			break;
		for( j = 0; j < PCL.GetExperiments( ); ++j )
			if( PCL.GetExperiment( j ) == ( PCLCluster.GetExperiment( i ).c_str( ) + 1 ) ) {
				m_setiDatasets.insert( j );
				break; } }
	for( i = 0; i < PCLCluster.GetGenes( ); ++i ) {
		if( ( iGene = PCL.GetGene( PCLCluster.GetGene( i ) ) ) == -1 ) {
			g_CatSleipnir.error( "CCoalesceCluster::Open( %s, %i ) unrecognized gene: %s", strPCL.c_str( ),
				iSkip, PCLCluster.GetGene( i ).c_str( ) );
			return -1; }
		m_setiGenes.insert( iGene ); }

	if( !pMotifs || ( ( i = strPCL.rfind( CPCL::GetExtension( ) ) ) !=
		( strPCL.length( ) - strlen( CPCL::GetExtension( ) ) ) ) )
		return PCLCluster.GetExperiments( );
	strMotifs = strPCL.substr( 0, i ) + c_szMotifs;
	ifsm.open( strMotifs.c_str( ) );
	if( !ifsm.is_open( ) )
		return PCLCluster.GetExperiments( );
	while( !ifsm.eof( ) && ( ifsm.peek( ) != -1 ) ) {
		SMotifMatch	sMotif;

		if( !sMotif.Open( ifsm, *pMotifs ) ) {
			g_CatSleipnir.warn( "CCoalesceCluster::Open( %s, %d ) could not open: %s", strPCL.c_str( ),
				iSkip, strMotifs.c_str( ) );
			return -1; }
		m_setsMotifs.insert( sMotif ); }

	return PCLCluster.GetExperiments( ); }

bool CCoalesceCluster::Open( const CHierarchy& Hierarchy, const vector<CCoalesceCluster>& vecClusters,
	const vector<string>& vecstrClusters, float dFraction, float dCutoff, size_t iCutoff,
	CCoalesceMotifLibrary* pMotifs ) {
	map<size_t, size_t>						mapiiGenes, mapiiDatasets;
	size_t									i, iClusters;
	map<size_t, size_t>::const_iterator		iterItem;
	vector<map<string, set<SMotifMatch> > >	vecmapstrsetsMotifs;

	vecmapstrsetsMotifs.resize( CCoalesceSequencerBase::ESubsequenceEnd );
	g_CatSleipnir.notice( "CCoalesceCluster::Open( %g ) merging clusters:", dFraction );
	if( !( iClusters = CCoalesceClusterImpl::Open( Hierarchy, vecClusters, vecstrClusters, mapiiGenes,
		mapiiDatasets, vecmapstrsetsMotifs ) ) ) {
		g_CatSleipnir.error( "CCoalesceCluster::Open( %g ) no clusters found", dFraction );
		return false; }

	Clear( );
	for( iterItem = mapiiGenes.begin( ); iterItem != mapiiGenes.end( ); ++iterItem )
		if( ( (float)iterItem->second / iClusters ) >= dFraction )
			m_setiGenes.insert( iterItem->first );
	for( iterItem = mapiiDatasets.begin( ); iterItem != mapiiDatasets.end( ); ++iterItem )
//		if( ( (float)iterItem->second / iClusters ) >= dFraction )
			m_setiDatasets.insert( iterItem->first );
	if( !pMotifs )
		return true;

	for( i = 0; i < vecmapstrsetsMotifs.size( ); ++i ) {
		const map<string, set<SMotifMatch> >&			mapstrsetsMotifs	= vecmapstrsetsMotifs[ i ];
		map<string, set<SMotifMatch> >::const_iterator	iterMotifs;

		for( iterMotifs = mapstrsetsMotifs.begin( ); iterMotifs != mapstrsetsMotifs.end( ); ++iterMotifs ) {
			const set<SMotifMatch>&	setsMotifs	= iterMotifs->second;

			if( !( ( setsMotifs.size( ) < iCutoff ) ?
				CCoalesceClusterImpl::OpenMotifs( setsMotifs, *pMotifs, dCutoff ) :
				CCoalesceClusterImpl::OpenMotifsHeuristic( setsMotifs, *pMotifs, dCutoff, iCutoff ) ) )
				return false; } }

	return true; }

bool CCoalesceClusterImpl::OpenMotifsHeuristic( const set<SMotifMatch>& setsMotifs,
	CCoalesceMotifLibrary& Motifs, float dCutoff, size_t iCutoff ) {
	vector<SMotifMatch>	vecsMotifs;
	bool				fDone;
	size_t				i, iMotifs;
	set<SMotifMatch>	setsMerged;

	iMotifs = setsMotifs.size( );
	g_CatSleipnir.notice( "CCoalesceClusterImpl::OpenMotifsHeuristic( %g ) resolving %d motifs", dCutoff,
		iMotifs );

	vecsMotifs.resize( iMotifs );
	copy( setsMotifs.begin( ), setsMotifs.end( ), vecsMotifs.begin( ) );
	do {
		fDone = true;
		sort( vecsMotifs.begin( ), vecsMotifs.end( ) );
		for( i = 0; ( i + 1 ) < vecsMotifs.size( ); ++i ) {
			SMotifMatch&	sOne	= vecsMotifs[ i ];
			SMotifMatch&	sTwo	= vecsMotifs[ i + 1 ];

			if( ( sOne.m_iMotif == -1 ) || ( sTwo.m_iMotif == -1 ) )
				break;
			if( Motifs.Align( sOne.m_iMotif, sTwo.m_iMotif, dCutoff ) > dCutoff )
				continue;
			if( sTwo.Open( sOne, sTwo, Motifs ) == -1 )
				return false;
			if( Motifs.GetPST( sTwo.m_iMotif )->Integrate( ) > iCutoff )
				Motifs.Simplify( sTwo.m_iMotif );
			fDone = false;
			iMotifs--;
			sOne.m_iMotif = -1; } }
	while( !fDone );

	for( i = 0; i < vecsMotifs.size( ); ++i )
		if( vecsMotifs[ i ].m_iMotif != -1 )
			setsMerged.insert( vecsMotifs[ i ] );

	return OpenMotifs( setsMerged, Motifs, dCutoff ); }

bool CCoalesceClusterImpl::OpenMotifs( const set<SMotifMatch>& setsMotifs, CCoalesceMotifLibrary& Motifs,
	float dCutoff ) {
	vector<SMotifMatch>	vecsMotifs;
	CDistanceMatrix		MatSimilarity;
	CHierarchy*			pHierMotifs;
	size_t				i, j;
	bool				fRet;

	g_CatSleipnir.notice( "CCoalesceClusterImpl::OpenMotifs( %g ) resolving %d motifs", dCutoff,
		setsMotifs.size( ) );

	vecsMotifs.resize( setsMotifs.size( ) );
	copy( setsMotifs.begin( ), setsMotifs.end( ), vecsMotifs.begin( ) );
	MatSimilarity.Initialize( vecsMotifs.size( ) );
	for( i = 0; i < vecsMotifs.size( ); ++i )
		for( j = ( i + 1 ); j < vecsMotifs.size( ); ++j )
			MatSimilarity.Set( i, j, -Motifs.Align( vecsMotifs[ i ].m_iMotif,
				vecsMotifs[ j ].m_iMotif, dCutoff ) );
	if( !( pHierMotifs = CClustHierarchical::Cluster( MatSimilarity ) ) )
		return false;
	fRet = CCoalesceClusterImpl::OpenMotifs( Motifs, *pHierMotifs, vecsMotifs, dCutoff, m_setsMotifs );
	pHierMotifs->Destroy( );
	return fRet; }

size_t CCoalesceClusterImpl::Open( const CHierarchy& Hier, const vector<CCoalesceCluster>& vecClusters,
	const vector<string>& vecstrClusters, map<size_t, size_t>& mapiiGenes, map<size_t, size_t>& mapiiDatasets,
	TVecMapStrSetSMotifs& vecmapstrsetsMotifs ) {
	set<size_t>::const_iterator			iterFrom;
	map<size_t, size_t>::iterator		iterTo;
	set<SMotifMatch>::const_iterator	iterMotif;

	if( !Hier.IsGene( ) )
		return ( Open( Hier.Get( false ), vecClusters, vecstrClusters, mapiiGenes, mapiiDatasets,
			vecmapstrsetsMotifs ) + Open( Hier.Get( true ), vecClusters, vecstrClusters, mapiiGenes,
			mapiiDatasets, vecmapstrsetsMotifs ) );

	const CCoalesceCluster&	Cluster	= vecClusters[ Hier.GetID( ) ];

	g_CatSleipnir.notice( "CCoalesceClusterImpl::Open( ) cluster %s",
		vecstrClusters[ Hier.GetID( ) ].c_str( ) );
	for( iterFrom = Cluster.GetGenes( ).begin( ); iterFrom != Cluster.GetGenes( ).end( ); ++iterFrom )
		if( ( iterTo = mapiiGenes.find( *iterFrom ) ) == mapiiGenes.end( ) )
			mapiiGenes[ *iterFrom ] = 1;
		else
			iterTo->second++;
	for( iterFrom = Cluster.m_setiDatasets.begin( ); iterFrom != Cluster.m_setiDatasets.end( ); ++iterFrom )
		if( ( iterTo = mapiiDatasets.find( *iterFrom ) ) == mapiiDatasets.end( ) )
			mapiiDatasets[ *iterFrom ] = 1;
		else
			iterTo->second++;
	for( iterMotif = Cluster.m_setsMotifs.begin( ); iterMotif != Cluster.m_setsMotifs.end( ); ++iterMotif )
		vecmapstrsetsMotifs[ iterMotif->m_eSubsequence ][ iterMotif->m_strType ].insert( *iterMotif );

	return 1; }

bool CCoalesceClusterImpl::OpenMotifs( CCoalesceMotifLibrary& Motifs, const CHierarchy& Hier,
	const vector<SMotifMatch>& vecsMotifs, float dCutoff, set<SMotifMatch>& setsMotifs ) {

	if( Hier.IsGene( ) || ( -Hier.GetSimilarity( ) < dCutoff ) ) {
		SMotifMatch	sMotif;
		size_t		iCount;

		sMotif.m_dZ = 0;
		if( !sMotif.Open( Hier, vecsMotifs, Motifs, iCount = 0 ) )
			return false;
		sMotif.m_dZ /= iCount;
		setsMotifs.insert( sMotif );
		return true; }

	return ( OpenMotifs( Motifs, Hier.Get( false ), vecsMotifs, dCutoff, setsMotifs ) &&
		OpenMotifs( Motifs, Hier.Get( true ), vecsMotifs, dCutoff, setsMotifs ) ); }

float CCoalesceCluster::GetSimilarity( const CCoalesceCluster& Cluster, size_t iGenes,
	size_t iDatasets ) const {
	size_t						iOverlapGenes;//, iOverlapDatasets;
	set<size_t>::const_iterator	iterItem;
	float						dRet, dGenes;//, dDatasets;

	for( iOverlapGenes = 0,iterItem = GetGenes( ).begin( ); iterItem != GetGenes( ).end( ); ++iterItem )
		if( Cluster.IsGene( *iterItem ) )
			iOverlapGenes++;
	dRet = dGenes = (float)iOverlapGenes / min( GetGenes( ).size( ), Cluster.GetGenes( ).size( ) );
/*
	for( iOverlapDatasets = 0,iterItem = m_setiDatasets.begin( ); iterItem != m_setiDatasets.end( ); ++iterItem )
		if( Cluster.IsDataset( *iterItem ) )
			iOverlapDatasets++;
	dGenes = sqrt( (float)iOverlapGenes / min( GetGenes( ).size( ), Cluster.GetGenes( ).size( ) ) );
	dDatasets = sqrt( (float)iOverlapDatasets / min( m_setiDatasets.size( ),
		Cluster.m_setiDatasets.size( ) ) );
	dRet = exp( ( ( iGenes * log( dGenes ) ) + ( iDatasets * log( dDatasets ) ) ) / ( iGenes + iDatasets ) );
*/
	g_CatSleipnir.debug( "CCoalesceCluster::GetSimilarity( %d, %d ) genes: %d, %d, %d = %g (%g)", iGenes,
		iDatasets, iOverlapGenes, GetGenes( ).size( ), Cluster.GetGenes( ).size( ), dGenes, dRet );
//	g_CatSleipnir.debug( "CCoalesceCluster::GetSimilarity( %d, %d ) datasets: %d, %d, %d = %g (%g)", iGenes,
//		iDatasets, iOverlapDatasets, m_setiDatasets.size( ), Cluster.m_setiDatasets.size( ), dDatasets, dRet );
	return dRet; }

void CCoalesceCluster::Snapshot( const CCoalesceGeneScores& GeneScores,
	CCoalesceGroupHistograms& Histograms ) {

	Histograms.SetTotal( GeneScores, GetGenes( ) );
	m_setiHistory.insert( GetHash( ) );
	CCoalesceClusterImpl::Snapshot( m_setiDatasets, m_veciPrevDatasets );
	CCoalesceClusterImpl::Snapshot( m_setsMotifs, m_vecsPrevMotifs );
	CCoalesceClusterImpl::Snapshot( GetGenes( ), m_veciPrevGenes ); }

}

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
#include "pcl.h"
#include "measure.h"
#include "statistics.h"
#include "fasta.h"
#include "motif.h"

namespace Sleipnir {

// CCoalesceMotifLibrary

const char	CCoalesceMotifLibrary::c_acBases[]	= "ACGT";

// CCoalesceGeneScores

bool CCoalesceGeneScores::Add( CCoalesceMotifLibrary& Motifs, const SFASTASequence& sSequence ) {
	size_t	i, iType;

	iType = AddType( sSequence.m_strType );
	for( i = 0; i < sSequence.m_vecstrSequences.size( ); ++i )
		if( !Add( Motifs, sSequence.m_vecstrSequences[ i ], iType, sSequence.m_fIntronFirst && !( i % 2 ) ) )
			return false;

	return true; }

bool CCoalesceGeneScores::Add( CCoalesceMotifLibrary& Motifs, const string& strSequence, size_t iType,
	bool fIntron ) {
	size_t			i;
	uint32_t		iMotifThem, iMotifUs;
	ESubsequence	eSubsequence;

// BUGBUG: make me span intron/exon boundaries
	for( i = 0; ( i + Motifs.GetK( ) ) <= strSequence.size( ); ++i ) {
		const string&	strKMer	= strSequence.substr( i, Motifs.GetK( ) );

		if( CCoalesceMotifLibrary::IsIgnorableKMer( strKMer ) )
			continue;
		if( ( iMotifThem = Motifs.GetID( strKMer ) ) == -1 ) {
			g_CatSleipnir.error( "CCoalesceHistograms::Add( %s, %d, %d ) unrecognized kmer: %s",
				strSequence.c_str( ), iType, fIntron, strKMer.c_str( ) );
			return false; }

		eSubsequence = fIntron ? ESubsequenceIntrons : ESubsequenceExons;
		iMotifUs = AddMotif( iType, eSubsequence, iMotifThem );
		CCoalesceSequencer::Get( iType, eSubsequence )[ iMotifUs ]++; }

	return true; }

// CCoalesceGroupHistograms

bool CCoalesceGroupHistograms::Add( const CCoalesceGeneScores& GeneScores ) {
	size_t		iTypeThem, iTypeUs, iSubsequence;
	uint32_t	iMotif;

	for( iTypeThem = 0; iTypeThem < GeneScores.GetTypes( ); ++iTypeThem ) {
		iTypeUs = AddType( GeneScores.GetType( iTypeThem ) );
		for( iSubsequence = ( ESubsequenceTotal + 1 );
			iSubsequence < GeneScores.GetSubsequences( iTypeThem ); ++iSubsequence ) {
			CCoalesceHistogramSet<>&	HistogramsTotal	= Get( iTypeUs, ESubsequenceTotal );
			CCoalesceHistogramSet<>&	Histograms		= Get( iTypeUs, (ESubsequence)iSubsequence );

			if( GeneScores.IsEmpty( iTypeThem, iSubsequence ) )
				continue;
			if( !HistogramsTotal.GetMembers( ) )
				HistogramsTotal.Initialize( GetMotifs( ), m_iBins );
			if( !Histograms.GetMembers( ) )
				Histograms.Initialize( GetMotifs( ), m_iBins );
			for( iMotif = 0; iMotif < GetMotifs( ); ++iMotif ) {
				unsigned short	s	= GeneScores.Get( iTypeThem, iSubsequence, iMotif );

				HistogramsTotal.Add( iMotif, s );
				Histograms.Add( iMotif, s ); } } }

	return true; }

// SMotifMatch

string SMotifMatch::Save( const CCoalesceMotifLibrary* pMotifs ) const {
	ostringstream	ossm;

	ossm << m_strType << '\t' << CCoalesceSequencerBase::GetSubsequence( m_eSubsequence ) << endl;
	if( pMotifs )
		ossm << pMotifs->GetMotif( m_iMotif );
	else
		ossm << m_iMotif;

	return ossm.str( ); }

// CCoalesceCluster

bool CCoalesceCluster::Initialize( const CPCL& PCL, CCoalesceCluster& Pot, float dPValue ) {
	size_t	i;

	m_setiConditions.clear( );
	m_setiGenes.clear( );
	m_setsMotifs.clear( );
	for( i = 0; i < PCL.GetExperiments( ); ++i )
		m_setiConditions.insert( i );

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

void CCoalesceCluster::CalculateHistograms( const vector<CCoalesceGeneScores>& vecGeneScores,
	CCoalesceGroupHistograms& Histograms ) const {
	set<size_t>::const_iterator	iterGene;

// BUGBUG: make me incremental
	Histograms.Clear( );
	if( vecGeneScores.empty( ) )
		return;
	for( iterGene = m_setiGenes.begin( ); iterGene != m_setiGenes.end( ); ++iterGene )
		Histograms.Add( vecGeneScores[ *iterGene ] ); }

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
	for( iterGene = m_setiGenes.begin( ); iterGene != m_setiGenes.end( ); ++iterGene )
		for( i = 0; i < m_veciCounts.size( ); ++i )
			if( !CMeta::IsNaN( d = PCL.Get( *iterGene, i ) ) ) {
				m_veciCounts[ i ]++;
				m_vecdCentroid[ i ] += d; }
	for( i = 0; i < m_veciCounts.size( ); ++i )
		if( m_veciCounts[ i ] ) {
			m_vecdCentroid[ i ] /= m_veciCounts[ i ]; } }

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

bool CCoalesceCluster::SelectMotifs( const CCoalesceGroupHistograms& HistsCluster,
	const CCoalesceGroupHistograms& HistsPot, float dPValue, const CCoalesceMotifLibrary* pMotifs ) {
	uint32_t	iMotif;

	m_setsMotifs.clear( );
	for( iMotif = 0; iMotif < HistsCluster.GetMotifs( ); ++iMotif )
		AddSignificant( iMotif, pMotifs, HistsCluster, HistsPot, dPValue );

	return true; }

void CCoalesceClusterImpl::AddSignificant( uint32_t iMotif, const CCoalesceMotifLibrary* pMotifs,
	const CCoalesceGroupHistograms& HistsCluster, const CCoalesceGroupHistograms& HistsPot, float dPValue ) {
	size_t									iType;
	double									dP;
	CCoalesceSequencerBase::ESubsequence	eSubsequence;

	for( iType = 0; iType < HistsCluster.GetTypes( ); ++iType )
		for( eSubsequence = CCoalesceSequencerBase::ESubsequenceBegin;
			(size_t)eSubsequence < HistsCluster.GetSubsequences( iType );
			eSubsequence = (CCoalesceSequencerBase::ESubsequence)( (size_t)eSubsequence + 1 ) ) {
			const CCoalesceHistogramSet<>&	HistSetCluster	= HistsCluster.Get( iType, eSubsequence );
			const CCoalesceHistogramSet<>&	HistSetPot		= HistsPot.Get( iType, eSubsequence );

			if( ( HistSetCluster.GetMembers( ) <= iMotif ) ||
				( HistSetPot.GetMembers( ) <= iMotif ) )
				continue;
// BUGBUG: multiple hypothesis correction?
			dP = HistSetCluster.PoissonTest( iMotif, HistSetPot );
// KS and Chi2 tests are too sensitive
			if( dP < dPValue ) {
				SMotifMatch	sMotif( iMotif, HistsCluster.GetType( iType ), eSubsequence );

				if( g_CatSleipnir.isInfoEnabled( ) ) {
					g_CatSleipnir.info( "CCoalesceClusterImpl::AddSignificant( %d, %g ) adding at %g:\n%s",
						iMotif, dPValue, dP, sMotif.Save( pMotifs ).c_str( ) );
					g_CatSleipnir.info( "Cluster	%s", HistSetCluster.Save( iMotif ).c_str( ) );
					g_CatSleipnir.info( "Pot	%s", HistSetPot.Save( iMotif ).c_str( ) ); }
				m_setsMotifs.insert( SMotifMatch( iMotif, HistsCluster.GetType( iType ), eSubsequence ) ); } } }

bool CCoalesceCluster::SelectGenes( const CPCL& PCL, const vector<CCoalesceGeneScores>& vecGeneScores,
	const CCoalesceGroupHistograms& HistsCluster, const CCoalesceGroupHistograms& HistsPot,
	CCoalesceCluster& Pot, float dPValue, const CCoalesceMotifLibrary* pMotifs ) {
	static const CCoalesceGeneScores	GeneScoresDummy;
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
		if( IsSignificant( iGene, PCL, pMotifs, vecGeneScores.empty( ) ? GeneScoresDummy :
			vecGeneScores[ iGene ], HistsCluster, HistsPot, Pot, dPValue ) )
			Add( iGene );
		else
			Pot.Add( iGene );

	return true; }

bool CCoalesceClusterImpl::IsSignificant( size_t iGene, const CPCL& PCL, const CCoalesceMotifLibrary* pMotifs,
	const CCoalesceGeneScores& GeneScores, const CCoalesceGroupHistograms& HistsCluster,
	const CCoalesceGroupHistograms& HistsPot, const CCoalesceCluster& Pot, float dProbability ) const {
	long double	dPNotGivenMotifs, dPGivenExpression;
	bool		fClustered;

	fClustered = IsGene( iGene );
	dPGivenExpression = CalculateProbabilityExpression( iGene, PCL, Pot, fClustered );
	if( dPGivenExpression < dProbability )
		return false;
	dPNotGivenMotifs = CalculateProbabilityMotifs( GeneScores, HistsCluster, HistsPot, fClustered );
	if( dPNotGivenMotifs > dProbability ) {
		g_CatSleipnir.debug( "CCoalesceClusterImpl::IsSignificant( %s ) rejected, exp. p=%g, seq. p=%g",
			PCL.GetGene( iGene ).c_str( ), dPGivenExpression, dPNotGivenMotifs );
		if( g_CatSleipnir.isDebugEnabled( ) ) {
			set<SMotifMatch>::const_iterator	iterMotif;

			for( iterMotif = m_setsMotifs.begin( ); iterMotif != m_setsMotifs.end( ); ++iterMotif )
				g_CatSleipnir.debug( "%d	%s", GeneScores.Get( GeneScores.GetType( iterMotif->m_strType ),
					iterMotif->m_eSubsequence, iterMotif->m_iMotif ), iterMotif->Save( pMotifs ).c_str( ) ); }
		return false; }

	return true; }

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

float CCoalesceClusterImpl::CalculateProbabilityMotifs( const CCoalesceGeneScores& GeneScores,
	const CCoalesceGroupHistograms& HistsCluster, const CCoalesceGroupHistograms& HistsPot,
	bool fClustered ) const {
	set<SMotifMatch>::const_iterator	iterMotif;
	size_t								i, iType, iCluster, iPot;
	unsigned short						sGene, sCluster, sPot;
	float								dRet, dCur, dPCluster, dPPot;

	if( m_setsMotifs.empty( ) )
		return 0;

	m_vecdPIn.resize( GeneScores.GetTypes( ) * CCoalesceSequencerBase::ESubsequenceEnd );
	fill( m_vecdPIn.begin( ), m_vecdPIn.end( ), -1 );
	m_vecdPOut.resize( m_vecdPIn.size( ) );
	fill( m_vecdPOut.begin( ), m_vecdPOut.end( ), -1 );
	for( iterMotif = m_setsMotifs.begin( ); iterMotif != m_setsMotifs.end( ); ++iterMotif ) {
		if( ( ( iType = GeneScores.GetType( iterMotif->m_strType ) ) == -1 ) ||
			( GeneScores.IsEmpty( iType, iterMotif->m_eSubsequence ) ) ||
			( HistsCluster.GetType( iterMotif->m_strType ) == -1 ) ||
			( HistsPot.GetType( iterMotif->m_strType ) == -1 ) )
			continue;
		sGene = GeneScores.Get( iType, iterMotif->m_eSubsequence, iterMotif->m_iMotif );
		{
			const CCoalesceHistogramSet<>&	HistSetCluster	= HistsCluster.Get( iterMotif->m_strType,
																iterMotif->m_eSubsequence );
			const CCoalesceHistogramSet<>&	HistSetPot		= HistsPot.Get( iterMotif->m_strType,
																iterMotif->m_eSubsequence );


			sCluster = HistSetCluster.Get( iterMotif->m_iMotif, sGene );
			iCluster = HistSetCluster.GetTotal( iterMotif->m_iMotif );
			sPot = HistSetPot.Get( iterMotif->m_iMotif, sGene );
			iPot = HistSetPot.GetTotal( iterMotif->m_iMotif );

			if( fClustered ) {
				sCluster--;
				iCluster--; }
			else {
				sPot--;
				iPot--; }
			dPCluster = (float)( sCluster + 1 ) / ( iCluster + HistSetCluster.GetEdges( ) );
			dPPot = (float)( sPot + 1 ) / ( iPot + HistSetPot.GetEdges( ) );
		}
		i = ( GeneScores.GetType( iterMotif->m_strType ) * CCoalesceSequencerBase::ESubsequenceEnd ) +
			iterMotif->m_eSubsequence;
		if( m_vecdPIn[ i ] < 0 )
			m_vecdPIn[ i ] = dPCluster;
		else
			m_vecdPIn[ i ] *= dPCluster;
		if( m_vecdPOut[ i ] < 0 )
			m_vecdPOut[ i ] = dPPot;
		else
			m_vecdPOut[ i ] *= dPPot; }

// BUGBUG: how to best handle this - multiple hypothesis, and, or, something?
	dRet = -1;
	for( i = 0; i < m_vecdPIn.size( ); ++i ) {
		if( ( m_vecdPIn[ i ] < 0 ) || ( m_vecdPOut[ i ] < 0 ) )
			continue;
		dCur = (float)( m_vecdPOut[ i ] / ( m_vecdPIn[ i ] + m_vecdPOut[ i ] ) );
		if( ( dRet < 0 ) || ( dCur < dRet ) )
			dRet = dCur; }

	return max( dRet, 0 ); }

bool CCoalesceCluster::Save( const string& strDirectory, size_t iID, const CPCL& PCL,
	const CCoalesceMotifLibrary* pMotifs ) const {
	ofstream							ofsm;
	char*								szTmp;
	CPCL								PCLCopy;
	size_t								iGeneFrom, iGeneTo, iExpFrom, iExpTo, iLength;
	string								strBase;
	set<SMotifMatch>::const_iterator	iterMotif;

	szTmp = new char[ iLength = ( strDirectory.length( ) + 16 ) ];
	sprintf_s( szTmp, iLength - 1, "%s/c%04d_XXXXXX", strDirectory.empty( ) ? "." : strDirectory.c_str( ),
		iID );
	_mktemp_s( szTmp, iLength - 1 );
	szTmp[ iLength - 1 ] = 0;
	strBase = szTmp;
	delete[] szTmp;

	ofsm.open( ( strBase + "_motifs.txt" ).c_str( ) );
	if( !ofsm.is_open( ) )
		return false;
	for( iterMotif = m_setsMotifs.begin( ); iterMotif != m_setsMotifs.end( ); ++iterMotif )
		ofsm << iterMotif->Save( pMotifs ) << endl;
	ofsm.close( );

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

	ofsm.clear( );
	ofsm.open( ( strBase + ".txt" ).c_str( ) );
	if( !ofsm.is_open( ) )
		return false;
	PCLCopy.Save( ofsm );
	ofsm.close( );

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

void CCoalesceCluster::Save( ostream& ostm, size_t iID, const CPCL& PCL,
	const CCoalesceMotifLibrary* pMotifs ) const {
	set<size_t>::const_iterator			iterID;
	set<SMotifMatch>::const_iterator	iterMotif;

	ostm << "Cluster\t" << iID << endl;
	ostm << "Genes";
	for( iterID = GetGenes( ).begin( ); iterID != GetGenes( ).end( ); ++iterID )
		ostm << '\t' << PCL.GetGene( *iterID );
	ostm << endl << "Conditions";
	for( iterID = GetConditions( ).begin( ); iterID != GetConditions( ).end( ); ++iterID )
		ostm << '\t' << *iterID;
	ostm << endl << "Motifs" << endl;
	for( iterMotif = GetMotifs( ).begin( ); iterMotif != GetMotifs( ).end( ); ++iterMotif )
		ostm << iterMotif->Save( pMotifs ) << endl; }

// CCoalesce

CCoalesceImpl::~CCoalesceImpl( ) {

	Clear( ); }

void CCoalesceImpl::Clear( ) {

	if( m_fMotifs && m_pMotifs )
		delete m_pMotifs;
	m_pMotifs = NULL; }

size_t CCoalesceImpl::GetMotifCount( ) const {

	return ( m_pMotifs ? m_pMotifs->GetMotifs( ) : 0 ); }

bool CCoalesce::Cluster( const CPCL& PCL, const CFASTA& FASTA, vector<CCoalesceCluster>& vecClusters ) {
	CPCL						PCLCopy;
	size_t						i, j;
	vector<size_t>				veciPCL2FASTA;
	vector<CCoalesceGeneScores>	vecGeneScores;

	PCLCopy.Open( PCL );
	if( FASTA.GetGenes( ) ) {
		if( !m_pMotifs ) {
			m_fMotifs = true;
			m_pMotifs = new CCoalesceMotifLibrary( GetK( ) ); }
		veciPCL2FASTA.resize( PCL.GetGenes( ) );
		vecGeneScores.resize( PCL.GetGenes( ) );
		for( i = 0; i < veciPCL2FASTA.size( ); ++i )
			if( ( veciPCL2FASTA[ i ] = FASTA.GetGene( PCL.GetGene( i ) ) ) != -1 ) {
				vector<SFASTASequence>	vecsSequences;

				FASTA.Get( veciPCL2FASTA[ i ], vecsSequences );
				for( j = 0; j < vecsSequences.size( ); ++j )
					if( !vecGeneScores[ i ].Add( *m_pMotifs, vecsSequences[ j ] ) )
						return false; } }
	else
		Clear( );
//DWORD iTmp, iHists, iConds, iMotifs, iGenes;
//iHists = iConds = iMotifs = iGenes = 0;
	while( true ) {
		CCoalesceCluster			Cluster, Pot;
		CCoalesceGroupHistograms	HistsCluster( GetMotifCount( ), GetBins( ) );
		CCoalesceGroupHistograms	HistsPot( GetMotifCount( ), GetBins( ) );

		for( i = 0; i < PCLCopy.GetGenes( ); ++i )
			Pot.Add( i );
		if( !Cluster.Initialize( PCLCopy, Pot, GetPValueCorrelation( ) ) )
			break;
		g_CatSleipnir.notice( "CCoalesce::Cluster( ) initialized %d genes with: %s, %s",
			Cluster.GetGenes( ).size( ), PCLCopy.GetGene( *Cluster.GetGenes( ).begin( ) ).c_str( ),
			PCLCopy.GetGene( *(++Cluster.GetGenes( ).begin( )) ).c_str( ) );
		while( !( Cluster.IsConverged( ) || Cluster.IsEmpty( ) ) ) {
//iTmp = GetTickCount( );
			Cluster.CalculateHistograms( vecGeneScores, HistsCluster );
			Pot.CalculateHistograms( vecGeneScores, HistsPot );
			Cluster.Snapshot( );
//iHists += GetTickCount( ) - iTmp;
//iTmp = GetTickCount( );
//Cluster.SelectConditions( PCLCopy, Pot, GetPValueCondition( ) );
//iConds += GetTickCount( ) - iTmp;
//iTmp = GetTickCount( );
//Cluster.SelectMotifs( HistsCluster, HistsPot, GetPValueMotif( ), m_pMotifs );
//iMotifs += GetTickCount( ) - iTmp;
//iTmp = GetTickCount( );
//Cluster.SelectGenes( PCLCopy, vecGeneScores, HistsCluster, HistsPot, Pot, GetProbabilityGene( ), m_pMotifs );
//iGenes += GetTickCount( ) - iTmp;
/*
			if( !( Cluster.SelectConditions( PCLCopy, Pot, GetPValueCondition( ) ) &&
				Cluster.SelectMotifs( HistsCluster, HistsPot, GetPValueMotif( ), m_pMotifs ) &&
				Cluster.SelectGenes( PCLCopy, vecGeneScores, HistsCluster, HistsPot, Pot,
				GetProbabilityGene( ), m_pMotifs ) ) )
				return false;
*/
			g_CatSleipnir.info( "CCoalesce::Cluster( ) processed %d genes, %d conditions, %d motifs",
				Cluster.GetGenes( ).size( ), Cluster.GetConditions( ).size( ), Cluster.GetMotifs( ).size( ) );
//cerr << "TIMER	" << iHists << '\t' << iConds << '\t' << iMotifs << '\t' << iGenes << endl;
			if( IsOutputIntermediate( ) )
				Cluster.Save( GetDirectoryIntermediate( ), vecClusters.size( ), PCLCopy, m_pMotifs ); }
		if( Cluster.IsConverged( ) ) {
			vecClusters.push_back( Cluster );
			Cluster.Subtract( PCLCopy ); } }

	return true; }

}

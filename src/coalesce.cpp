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
#include "pst.h"

namespace Sleipnir {

// CCoalesceMotifLibrary

// Order independent, but complements must be adjacent
const char	CCoalesceMotifLibraryImpl::c_acBases[]	= "ATCG";

CCoalesceMotifLibraryImpl::~CCoalesceMotifLibraryImpl( ) {
	size_t	i;

	for( i = 0; i < m_vecpPSTs.size( ); ++i )
		delete m_vecpPSTs[ i ]; }

float CCoalesceMotifLibrary::GetMatch( const string& strSequence, uint32_t iMotif ) const {
	size_t		i;
	float		dRet;
	const CPST*	pPST;

// BUGBUG: this could in theory be implemented
	if( ( GetType( iMotif ) != ETypePST ) || !( pPST = GetPST( iMotif ) ) )
		return CMeta::GetNaN( );

	dRet = 0;
	for( i = 1; i < pPST->GetDepth( ); ++i )
		dRet += pPST->GetMatch( strSequence, pPST->GetDepth( ) - i );
	for( i = 0; i < strSequence.length( ); ++i )
		dRet += pPST->GetMatch( strSequence.substr( i ) );

	return dRet; }

std::string CCoalesceMotifLibraryImpl::GetMotif( uint32_t iMotif ) const {
	std::string	strKMer;

// kmer
	if( iMotif < GetKMers( ) )
		return ID2KMer( iMotif, m_iK );
// reverse complement
	if( iMotif < GetBasePSTs( ) ) {
		strKMer = GetRCOne( iMotif );
		return ( strKMer + c_cSeparator + GetReverseComplement( strKMer ) ); }
// pst
	return GetPST( iMotif )->GetMotif( ); }

CPST* CCoalesceMotifLibraryImpl::CreatePST( uint32_t& iMotif ) {
	CPST*	pRet;

	iMotif = GetBasePSTs( ) + GetPSTs( );
	m_vecpPSTs.push_back( pRet = new CPST( strlen( c_acBases ) ) );
	return pRet; }

uint32_t CCoalesceMotifLibraryImpl::MergeKMers( const std::string& strOne, const std::string& strTwo,
	float dCutoff ) {
	int			iOffset;
	float		dScore;
	uint32_t	iRet;
	CPST*		pPST;

	if( ( dScore = Align( strOne, strTwo, dCutoff, iOffset ) ) > dCutoff )
		return -1;

	pPST = CreatePST( iRet );
	pPST->Add( strOne, strTwo, iOffset );
	if( g_CatSleipnir.isInfoEnabled( ) )
		g_CatSleipnir.info( "CCoalesceMotifLibraryImpl::MergeKMers( %s, %s, %g ) merged at %g to %s",
			strOne.c_str( ), strTwo.c_str( ), dCutoff, dScore, pPST->GetMotif( ).c_str( ) );
	return iRet; }

uint32_t CCoalesceMotifLibraryImpl::MergeKMerRC( const string& strKMer, uint32_t iRC, float dCutoff ) {
	string		strOne, strTwo;
	float		dOne, dTwo;
	int			iOne, iTwo;
	uint32_t	iRet;
	CPST*		pPST;

	strOne = GetRCOne( iRC );
	strTwo = GetReverseComplement( strOne );
	dOne = Align( strKMer, strOne, dCutoff, iOne );
	dTwo = Align( strKMer, strTwo, dCutoff, iTwo );
	if( ( dOne > dCutoff ) && ( dTwo > dCutoff ) )
		return -1;

	pPST = CreatePST( iRet );
	if( dOne < dTwo ) {
		pPST->Add( strKMer, strOne, iOne );
		pPST->Add( strTwo ); }
	else {
		pPST->Add( strKMer, strTwo, iTwo );
		pPST->Add( strOne ); }
	if( g_CatSleipnir.isInfoEnabled( ) )
		g_CatSleipnir.info( "CCoalesceMotifLibraryImpl::MergeKMerRC( %s, %s, %g ) merged at %g to %s",
			strKMer.c_str( ), GetMotif( iRC ).c_str( ), dCutoff, min( dOne, dTwo ),
			pPST->GetMotif( ).c_str( ) );
	return iRet; }

struct SCrossRCs {
	string	m_strOne;
	string	m_strTwo;
	float	m_dScore;
	int		m_iOffset;
};

uint32_t CCoalesceMotifLibraryImpl::MergeRCs( uint32_t iOne, uint32_t iTwo, float dCutoff ) {
	SCrossRCs	asCrosses[ 4 ];
	uint32_t	iRet;
	CPST*		pPST;
	size_t		i, iMin;
	float		dMin;

	asCrosses[ 0 ].m_strOne = asCrosses[ 1 ].m_strOne = GetRCOne( iOne );
	asCrosses[ 0 ].m_strTwo = asCrosses[ 2 ].m_strTwo = GetRCOne( iTwo );
	asCrosses[ 1 ].m_strTwo = asCrosses[ 3 ].m_strTwo = GetReverseComplement( asCrosses[ 0 ].m_strTwo );
	asCrosses[ 2 ].m_strOne = asCrosses[ 3 ].m_strOne = GetReverseComplement( asCrosses[ 0 ].m_strOne );
	dMin = FLT_MAX;
	for( iMin = i = 0; i < ARRAYSIZE(asCrosses); ++i ) {
		asCrosses[ i ].m_dScore = Align( asCrosses[ i ].m_strOne, asCrosses[ i ].m_strTwo, dCutoff,
			asCrosses[ i ].m_iOffset );
		if( asCrosses[ i ].m_dScore < dMin ) {
			dMin = asCrosses[ i ].m_dScore;
			iMin = i; } }
	if( dMin > dCutoff )
		return -1;

	pPST = CreatePST( iRet );
	pPST->Add( asCrosses[ iMin ].m_strOne, asCrosses[ iMin ].m_strTwo, asCrosses[ iMin ].m_iOffset );
	pPST->Add( asCrosses[ ( iMin + 2 ) % ARRAYSIZE(asCrosses) ].m_strOne );
	pPST->Add( asCrosses[ ( iMin + 1 ) % ARRAYSIZE(asCrosses) ].m_strTwo );
	if( g_CatSleipnir.isInfoEnabled( ) )
		g_CatSleipnir.info( "CCoalesceMotifLibraryImpl::MergeRCs( %s, %s, %g ) merged at %g to %s",
			GetMotif( iOne ).c_str( ), GetMotif( iTwo ).c_str( ), dCutoff, dMin, pPST->GetMotif( ).c_str( ) );
	return iRet; }

uint32_t CCoalesceMotifLibraryImpl::MergeKMerPST( const string& strKMer, const CPST& PSTIn, float dCutoff ) {
	int			iOffset;
	float		dScore;
	uint32_t	iRet;
	CPST*		pPSTOut;

	if( ( dScore = PSTIn.Align( strKMer, m_dPenaltyGap, m_dPenaltyMismatch, dCutoff, iOffset ) ) > dCutoff )
		return -1;

	pPSTOut = CreatePST( iRet );
	pPSTOut->Add( strKMer, PSTIn, iOffset );
	if( g_CatSleipnir.isInfoEnabled( ) )
		g_CatSleipnir.info( "CCoalesceMotifLibraryImpl::MergeKMerPST( %s, %s, %g ) merged at %g to %s",
			strKMer.c_str( ), PSTIn.GetMotif( ).c_str( ), dCutoff, dScore, pPSTOut->GetMotif( ).c_str( ) );
	return iRet; }

uint32_t CCoalesceMotifLibraryImpl::MergeRCPST( uint32_t iRC, const CPST& PSTIn, float dCutoff ) {
	int			iOne, iTwo;
	uint32_t	iRet;
	CPST*		pPSTOut;
	string		strOne, strTwo;
	float		dOne, dTwo;

	strOne = GetRCOne( iRC );
	strTwo = GetReverseComplement( strOne );
	dOne = PSTIn.Align( strOne, m_dPenaltyGap, m_dPenaltyMismatch, dCutoff, iOne );
	dTwo = PSTIn.Align( strTwo, m_dPenaltyGap, m_dPenaltyMismatch, dCutoff, iTwo );
	if( ( dOne > dCutoff ) && ( dTwo > dCutoff ) )
		return -1;

	pPSTOut = CreatePST( iRet );
	if( dOne < dTwo ) {
		pPSTOut->Add( strOne, PSTIn, iOne );
		pPSTOut->Add( strTwo ); }
	else {
		pPSTOut->Add( strTwo, PSTIn, iTwo );
		pPSTOut->Add( strOne ); }
	if( g_CatSleipnir.isInfoEnabled( ) )
		g_CatSleipnir.info( "CCoalesceMotifLibraryImpl::MergeRCPST( %s, %s, %g ) merged at %g to %s",
			GetMotif( iRC ).c_str( ), PSTIn.GetMotif( ).c_str( ), dCutoff, min( dOne, dTwo ),
			pPSTOut->GetMotif( ).c_str( ) );
	return iRet; }

uint32_t CCoalesceMotifLibraryImpl::MergePSTs( const CPST& PSTOne, const CPST& PSTTwo, float dCutoff ) {
	int			iOffset;
	uint32_t	iRet;
	CPST*		pPSTOut;
	float		dScore;

	if( ( dScore = PSTOne.Align( PSTTwo, m_dPenaltyGap, m_dPenaltyMismatch, dCutoff, iOffset ) ) > dCutoff )
		return -1;

	pPSTOut = CreatePST( iRet );
	if( iOffset < 0 ) {
		pPSTOut->Add( PSTTwo );
		pPSTOut->Add( PSTOne, -iOffset ); }
	else {
		pPSTOut->Add( PSTOne );
		pPSTOut->Add( PSTTwo, iOffset ); }
	if( g_CatSleipnir.isInfoEnabled( ) )
		g_CatSleipnir.info( "CCoalesceMotifLibraryImpl::MergePSTs( %s, %s, %g ) merged at %g to %s",
			PSTOne.GetMotif( ).c_str( ), PSTTwo.GetMotif( ).c_str( ), dCutoff, dScore,
			pPSTOut->GetMotif( ).c_str( ) );
	return iRet; }

// CCoalesceGeneScores

bool CCoalesceGeneScores::Add( CCoalesceMotifLibrary& Motifs, const SFASTASequence& sSequence,
	vector<vector<unsigned short> >& vecvecsCounts, vector<size_t>& veciLengths ) {
	size_t			i, iType, iSubsequence, iLength;
	ESubsequence	eSubsequence;
	uint32_t		iMotifUs, iMotifThem;
	unsigned short	s;

	vecvecsCounts.resize( ESubsequenceEnd );
	for( i = 0; i < vecvecsCounts.size( ); ++i )
		fill( vecvecsCounts[ i ].begin( ), vecvecsCounts[ i ].end( ), 0 );
	veciLengths.resize( ESubsequenceEnd );
	fill( veciLengths.begin( ), veciLengths.end( ), 0 );

	iType = AddType( sSequence.m_strType );
	for( i = 0; i < sSequence.m_vecstrSequences.size( ); ++i )
		if( !Add( Motifs, sSequence.m_vecstrSequences[ i ], iType, sSequence.m_fIntronFirst && !( i % 2 ),
			vecvecsCounts, veciLengths ) )
			return false;
	for( iSubsequence = ESubsequenceBegin; iSubsequence < vecvecsCounts.size( ); ++iSubsequence ) {
		const vector<unsigned short>&	vecsCounts	= vecvecsCounts[ iSubsequence ];

		if( iLength = veciLengths[ iSubsequence ] ) {
			eSubsequence = (ESubsequence)iSubsequence;
			for( iMotifThem = 0; iMotifThem < vecsCounts.size( ); ++iMotifThem )
				if( s = vecsCounts[ iMotifThem ] ) {
					iMotifUs = AddMotif( iType, eSubsequence, iMotifThem );
					CCoalesceSequencer<TVecD>::Get( iType, eSubsequence )[ iMotifUs ] = (float)s /
						iLength; } } }

	return true; }

bool CCoalesceGeneScores::Add( CCoalesceMotifLibrary& Motifs, const string& strSequence, size_t iType,
	bool fIntron, vector<vector<unsigned short> >& vecvecsCounts, vector<size_t>& veciLengths ) {
	size_t			i, j;
	ESubsequence	eSubsequence;

	eSubsequence = fIntron ? ESubsequenceIntrons : ESubsequenceExons;
	veciLengths[ ESubsequenceTotal ] += strSequence.length( );
	veciLengths[ eSubsequence ] += strSequence.length( );

// BUGBUG: make me span intron/exon boundaries
	for( i = 0; ( i + Motifs.GetK( ) ) <= strSequence.size( ); ++i ) {
		string				strKMer	= strSequence.substr( i, Motifs.GetK( ) );
		vector<uint32_t>	veciMotifs;

		if( !Motifs.GetMatches( strKMer, veciMotifs ) ) {
			g_CatSleipnir.error( "CCoalesceHistograms::Add( %s, %d, %d ) unrecognized kmer: %s",
				strSequence.c_str( ), iType, fIntron, strKMer.c_str( ) );
			return false; }
		for( j = 0; j < veciMotifs.size( ); ++j )
			Add( eSubsequence, veciMotifs[ j ], Motifs.GetMotifs( ), vecvecsCounts ); }

	return true; }

bool CCoalesceGeneScores::Add( CCoalesceMotifLibrary& Motifs, const SFASTASequence& sSequence,
	uint32_t iMotifThem, vector<float>& vecdScores, vector<size_t>& veciLengths ) {
	size_t			i, iType, iSubsequence, iLength;
	ESubsequence	eSubsequence;
	uint32_t		iMotifUs;
	float			dScore;

	vecdScores.resize( ESubsequenceEnd );
	fill( vecdScores.begin( ), vecdScores.end( ), 0.0f );
	veciLengths.resize( ESubsequenceEnd );
	fill( veciLengths.begin( ), veciLengths.end( ), 0 );

	iType = AddType( sSequence.m_strType );
	for( i = 0; i < sSequence.m_vecstrSequences.size( ); ++i )
		if( !Add( Motifs, sSequence.m_vecstrSequences[ i ], iType, sSequence.m_fIntronFirst && !( i % 2 ),
			iMotifThem, vecdScores, veciLengths ) )
			return false;
	for( iSubsequence = ESubsequenceBegin; iSubsequence < vecdScores.size( ); ++iSubsequence )
		if( ( iLength = veciLengths[ iSubsequence ] ) && ( dScore = vecdScores[ iSubsequence ] ) ) {
			eSubsequence = (ESubsequence)iSubsequence;
			iMotifUs = AddMotif( iType, eSubsequence, iMotifThem );
			CCoalesceSequencer<TVecD>::Get( iType, eSubsequence )[ iMotifUs ] = dScore / iLength; }

	return true; }

bool CCoalesceGeneScores::Add( CCoalesceMotifLibrary& Motifs, const string& strSequence, size_t iType,
	bool fIntron, uint32_t iMotif, vector<float>& vecdScores, vector<size_t>& veciLengths ) {
	ESubsequence	eSubsequence;
	float			dScore;

	eSubsequence = fIntron ? ESubsequenceIntrons : ESubsequenceExons;
	veciLengths[ ESubsequenceTotal ] += strSequence.length( );
	veciLengths[ eSubsequence ] += strSequence.length( );

// BUGBUG: make me span intron/exon boundaries
	dScore = Motifs.GetMatch( strSequence, iMotif );
	vecdScores[ ESubsequenceTotal ] += dScore;
	vecdScores[ eSubsequence ] += dScore;

	return true; }

void CCoalesceGeneScores::Subtract( const SMotifMatch& sMotif ) {
	size_t	iType;

	if( ( iType = GetType( sMotif.m_strType ) ) == -1 )
		return;
	{
		TVecD&	vecdScores	= Get( iType, sMotif.m_eSubsequence );

		if( !vecdScores.empty( ) ) {
			const TMapII&			mapiiMotifs	= m_MotifMaps.Get( iType, sMotif.m_eSubsequence );
			TMapII::const_iterator	iterMotif;

			if( ( iterMotif = mapiiMotifs.find( sMotif.m_iMotif ) ) != mapiiMotifs.end( ) )
				vecdScores[ iterMotif->second ] -= sMotif.m_dAverage; }
	} }

// CCoalesceGroupHistograms

bool CCoalesceGroupHistograms::Add( const CCoalesceGeneScores& GeneScores, bool fSubtract,
	uint32_t iMotifThem ) {
	size_t			iTypeThem, iTypeUs, iSubsequence;
	uint32_t		iMotifUs;
	unsigned short	sDelta;
	float			dValue;

	sDelta = fSubtract ? -1 : 1;
	for( iTypeThem = 0; iTypeThem < GeneScores.GetTypes( ); ++iTypeThem ) {
		iTypeUs = AddType( GeneScores.GetType( iTypeThem ) );
		for( iSubsequence = ESubsequenceBegin; iSubsequence < GeneScores.GetSubsequences( iTypeThem );
			++iSubsequence ) {
			CCoalesceHistogramSet<>&	Histograms	= Get( iTypeUs, (ESubsequence)iSubsequence );

			if( !GeneScores.GetMotifs( iTypeThem, iSubsequence ) )
				continue;
			if( !Histograms.GetMembers( ) )
				Histograms.Initialize( GetMotifs( ), m_iBins, m_dStep );
			if( iMotifThem == -1 ) {
				for( iMotifUs = 0; iMotifUs < GeneScores.GetMotifs( iTypeThem, iSubsequence ); ++iMotifUs )
					if( dValue = GeneScores.GetLocal( iTypeThem, iSubsequence, iMotifUs ) ) {
						if( !Histograms.Add( GeneScores.GetMotif( iTypeThem, iSubsequence, iMotifUs ), dValue,
							sDelta ) )
							return false; } }
			else if( dValue = GeneScores.GetGlobal( iTypeThem, iSubsequence, iMotifThem ) ) {
				if( !Histograms.Add( iMotifThem, dValue, sDelta ) )
					return false; } } }

	return true; }

void CCoalesceGroupHistograms::Save( ostream& ostm, const CCoalesceMotifLibrary* pMotifs ) const {
	size_t		iType, iSubsequence;
	uint32_t	iMotif;

	for( iType = 0; iType < GetTypes( ); ++iType )
		for( iSubsequence = 0; iSubsequence < GetSubsequences( iType ); ++iSubsequence ) {
			const CCoalesceHistogramSet<>&	Histograms	= Get( iType, (ESubsequence)iSubsequence );

			ostm << GetType( iType ) << '\t' << GetSubsequence( (ESubsequence)iSubsequence ) << endl;
			for( iMotif = 0; iMotif < ( pMotifs ? pMotifs->GetMotifs( ) : GetMotifs( ) ); ++iMotif ) {
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

// SMotifMatch

string SMotifMatch::Save( const CCoalesceMotifLibrary* pMotifs ) const {
	ostringstream	ossm;

	ossm << m_strType << '\t' << CCoalesceSequencerBase::GetSubsequence( m_eSubsequence ) << '\t' << m_dZ <<
		endl;
	if( pMotifs )
		ossm << pMotifs->GetMotif( m_iMotif );
	else
		ossm << m_iMotif;

	return ossm.str( ); }

// CCoalesceCluster

bool CCoalesceCluster::Initialize( const CPCL& PCL, CCoalesceCluster& Pot,
	set<pair<size_t, size_t> >& setpriiSeeds, float dFraction, float dPValue ) {
	size_t	i;

	m_setiConditions.clear( );
	m_setiGenes.clear( );
	m_setsMotifs.clear( );
	for( i = 0; i < PCL.GetExperiments( ); ++i )
		m_setiConditions.insert( i );
	m_vecdPriors.resize( PCL.GetGenes( ) );
	fill( m_vecdPriors.begin( ), m_vecdPriors.end( ), 1.0f );

	return ( AddSeedPair( PCL, Pot, setpriiSeeds, dFraction, dPValue ) &&
		AddCorrelatedGenes( PCL, Pot, dPValue ) ); }

void CCoalesceClusterImpl::Add( size_t iGene, CCoalesceCluster& Pot ) {

	m_setiGenes.insert( iGene );
	Pot.m_setiGenes.erase( iGene ); }

bool CCoalesceClusterImpl::AddSeedPair( const CPCL& PCL, CCoalesceCluster& Pot,
	set<pair<size_t, size_t> >& setpriiSeeds, float dFraction, float dPValue ) {
	size_t					i, j, iOne, iTwo, iN;
	double					dR, dCur, dMin;
	pair<size_t, size_t>	priiSeed;

	if( PCL.GetGenes( ) < 2 ) {
		g_CatSleipnir.error( "CCoalesceClusterImpl::AddSeedPair( %g ) found no genes", dPValue );
		return false; }
	dMin = DBL_MAX;
	for( i = 0; i < PCL.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < PCL.GetGenes( ); ++j ) {
			if( ( ( (float)rand( ) / RAND_MAX ) > dFraction ) ||
				( ( dR = CMeasurePearson::Pearson( PCL.Get( i ), PCL.GetExperiments( ), PCL.Get( j ),
				PCL.GetExperiments( ), IMeasure::EMapNone, NULL, NULL, &iN ) ) < 0 ) )
			continue;
		if( ( dCur = CStatistics::PValuePearson( dR, iN ) ) < dMin ) {
			priiSeed.first = i;
			priiSeed.second = j;
			if( setpriiSeeds.find( priiSeed ) == setpriiSeeds.end( ) ) {
				dMin = dCur;
				iOne = i;
				iTwo = j; } } }
	if( ( dMin * PCL.GetGenes( ) * PCL.GetGenes( ) * dFraction / 2 ) < dPValue ) {
		priiSeed.first = iOne;
		priiSeed.second = iTwo;
		setpriiSeeds.insert( priiSeed );
		Add( iOne, Pot );
		Add( iTwo, Pot );
		return true; }

	g_CatSleipnir.notice( "CCoalesceClusterImpl::AddSeedPair( %g, %g ) no seed pair found", dFraction,
		dPValue );
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

void CCoalesceCluster::CalculateHistograms( const vector<CCoalesceGeneScores>& vecGeneScores,
	CCoalesceGroupHistograms& HistogramsCluster, CCoalesceGroupHistograms* pHistogramsPot ) const {
	set<size_t>::const_iterator	iterGene;
	size_t						i;

	if( vecGeneScores.empty( ) )
		return;
	for( iterGene = m_setiGenes.begin( ); iterGene != m_setiGenes.end( ); ++iterGene )
		if( !binary_search( m_veciPrevGenes.begin( ), m_veciPrevGenes.end( ), *iterGene ) ) {
			const CCoalesceGeneScores&	GeneScores	= vecGeneScores[ *iterGene ];

			HistogramsCluster.Add( GeneScores, false );
			if( pHistogramsPot )
				pHistogramsPot->Add( GeneScores, true ); }
	for( i = 0; i < m_veciPrevGenes.size( ); ++i )
		if( m_setiGenes.find( m_veciPrevGenes[ i ] ) == m_setiGenes.end( ) ) {
			const CCoalesceGeneScores&	GeneScores	= vecGeneScores[ m_veciPrevGenes[ i ] ];

			HistogramsCluster.Add( GeneScores, true );
			if( pHistogramsPot )
				pHistogramsPot->Add( GeneScores, false ); } }

void CCoalesceCluster::Subtract( CPCL& PCL ) const {
	set<size_t>::const_iterator	iterGene, iterCondition;

	for( iterGene = m_setiGenes.begin( ); iterGene != m_setiGenes.end( ); ++iterGene )
		for( iterCondition = m_setiConditions.begin( ); iterCondition != m_setiConditions.end( );
			++iterCondition )
			PCL.Get( *iterGene, *iterCondition ) -= m_vecdCentroid[ *iterCondition ]; }

void CCoalesceCluster::Subtract( vector<CCoalesceGeneScores>& vecGeneScores ) const {
	set<size_t>::const_iterator			iterGene;
	set<SMotifMatch>::const_iterator	iterMotif;

	for( iterGene = m_setiGenes.begin( ); iterGene != m_setiGenes.end( ); ++iterGene )
		for( iterMotif = m_setsMotifs.begin( ); iterMotif != m_setsMotifs.end( ); ++iterMotif )
			vecGeneScores[ *iterGene ].Subtract( *iterMotif ); }

void CCoalesceClusterImpl::CalculateCentroid( const CPCL& PCL ) {
	set<size_t>::const_iterator	iterGene;
	size_t						i;
	float						d;

	m_veciCounts.resize( PCL.GetExperiments( ) );
	fill( m_veciCounts.begin( ), m_veciCounts.end( ), 0 );
	m_vecdCentroid.resize( PCL.GetExperiments( ) );
	fill( m_vecdCentroid.begin( ), m_vecdCentroid.end( ), 0.0f );
	m_vecdStdevs.resize( PCL.GetExperiments( ) );
	fill( m_vecdStdevs.begin( ), m_vecdStdevs.end( ), 0.0f );
	for( iterGene = m_setiGenes.begin( ); iterGene != m_setiGenes.end( ); ++iterGene )
		for( i = 0; i < m_veciCounts.size( ); ++i )
			if( !CMeta::IsNaN( d = PCL.Get( *iterGene, i ) ) ) {
				m_veciCounts[ i ]++;
				m_vecdCentroid[ i ] += d;
				m_vecdStdevs[ i ] += d * d; }
	for( i = 0; i < m_veciCounts.size( ); ++i ) {
		if( m_veciCounts[ i ] ) {
			m_vecdCentroid[ i ] /= m_veciCounts[ i ];
			m_vecdStdevs[ i ] = ( m_vecdStdevs[ i ] / m_veciCounts[ i ] ) - ( m_vecdCentroid[ i ] *
				m_vecdCentroid[ i ] );
			m_vecdStdevs[ i ] = ( m_vecdStdevs[ i ] < 0 ) ? 0 : sqrt( m_vecdStdevs[ i ] ); }
		g_CatSleipnir.info( "CCoalesceClusterImpl::CalculateCentroid( ) condition %d: count %d, mean %g, stdev %g",
			i, m_veciCounts[ i ], m_vecdCentroid[ i ], m_vecdStdevs[ i ] ); } }

bool CCoalesceCluster::SelectConditions( const CPCL& PCL, const vector<CCoalesceImpl::SDataset>& vecsDatasets,
	const CCoalesceCluster& Pot, float dPValue ) {
	vector<float>				vecdCluster, vecdPot, vecdDatasetCluster, vecdDatasetPot;
	vector<size_t>				veciDatasetCluster, veciDatasetPot;
	size_t						iDataset, iCondition;
	set<size_t>::const_iterator	iterGene;
	double						dP, dZ;
	float						d;

	m_vecdZs.resize( PCL.GetExperiments( ) );
	m_setiConditions.clear( );
	for( iDataset = 0; iDataset < vecsDatasets.size( ); ++iDataset ) {
		const CCoalesceImpl::SDataset&	sDataset	= vecsDatasets[ iDataset ];

		if( sDataset.GetConditions( ) == 1 ) {
			iCondition = sDataset.GetCondition( 0 );
			vecdCluster.clear( );
			vecdCluster.reserve( m_setiGenes.size( ) );
			vecdPot.clear( );
			vecdPot.reserve( Pot.m_setiGenes.size( ) );
			for( iterGene = m_setiGenes.begin( ); iterGene != m_setiGenes.end( ); ++iterGene )
				if( !CMeta::IsNaN( d = PCL.Get( *iterGene, iCondition ) ) )
					vecdCluster.push_back( d );
			for( iterGene = Pot.m_setiGenes.begin( ); iterGene != Pot.m_setiGenes.end( ); ++iterGene )
				if( !CMeta::IsNaN( d = PCL.Get( *iterGene, iCondition ) ) )
					vecdPot.push_back( d );
			dZ = CStatistics::CohensD( vecdCluster, vecdPot );
			dP = CStatistics::ZTest( dZ, vecdCluster.size( ) ) * vecsDatasets.size( );
			if( dP < dPValue ) {
				g_CatSleipnir.info( "CCoalesceCluster::SelectConditions( %g ) selected condition %d at %g, z=%g",
					dPValue, iCondition, dP, dZ );
				m_setiConditions.insert( iCondition );
				m_vecdZs[ iCondition ] = (float)dZ; } }
		else {
			vecdDatasetCluster.resize( sDataset.GetConditions( ) );
			fill( vecdDatasetCluster.begin( ), vecdDatasetCluster.end( ), 0.0f );
			vecdDatasetPot.resize( sDataset.GetConditions( ) );
			fill( vecdDatasetPot.begin( ), vecdDatasetPot.end( ), 0.0f );
			veciDatasetCluster.resize( sDataset.GetConditions( ) );
			fill( veciDatasetCluster.begin( ), veciDatasetCluster.end( ), 0 );
			veciDatasetPot.resize( sDataset.GetConditions( ) );
			fill( veciDatasetPot.begin( ), veciDatasetPot.end( ), 0 );
			for( iterGene = m_setiGenes.begin( ); iterGene != m_setiGenes.end( ); ++iterGene )
				for( iCondition = 0; iCondition < sDataset.GetConditions( ); ++iCondition )
					if( !CMeta::IsNaN( d = PCL.Get( *iterGene, sDataset.GetCondition( iCondition ) ) ) ) {
						veciDatasetCluster[ iCondition ]++;
						vecdDatasetCluster[ iCondition ] += d; }
			for( iterGene = Pot.m_setiGenes.begin( ); iterGene != Pot.m_setiGenes.end( ); ++iterGene )
				for( iCondition = 0; iCondition < sDataset.GetConditions( ); ++iCondition )
					if( !CMeta::IsNaN( d = PCL.Get( *iterGene, sDataset.GetCondition( iCondition ) ) ) ) {
						veciDatasetPot[ iCondition ]++;
						vecdDatasetPot[ iCondition ] += d; }
			for( iCondition = 0; iCondition < vecdDatasetCluster.size( ); ++iCondition ) {
				if( veciDatasetCluster[ iCondition ] )
					vecdDatasetCluster[ iCondition ] /= veciDatasetCluster[ iCondition ];
				if( veciDatasetPot[ iCondition ] )
					vecdDatasetPot[ iCondition ] /= veciDatasetPot[ iCondition ]; }
			dP = CStatistics::MultivariateNormalCDF( vecdDatasetCluster, vecdDatasetPot,
				sDataset.m_MatCovChol );
			if( dP > 0.5 )
				dP = 1 - dP;
			if( dP < dPValue ) {
				dZ = 0;
				for( iCondition = 0; iCondition < sDataset.GetConditions( ); ++iCondition )
					if( sDataset.m_vecdStdevs[ iCondition ] ) {
						d = ( vecdDatasetCluster[ iCondition ] - vecdDatasetPot[ iCondition ] ) /
							sDataset.m_vecdStdevs[ iCondition ];
						dZ += d * d; }
				dZ = sqrt( dZ );
				g_CatSleipnir.info( "CCoalesceCluster::SelectConditions( %g ) selected dataset %d at %g, z=%g",
					dPValue, iDataset, dP, dZ );
				for( iCondition = 0; iCondition < sDataset.GetConditions( ); ++iCondition ) {
					g_CatSleipnir.info( "CCoalesceCluster::SelectConditions( %g ) contains condition %d",
						dPValue, sDataset.GetCondition( iCondition ) );
					m_setiConditions.insert( sDataset.GetCondition( iCondition ) );
					m_vecdZs[ sDataset.GetCondition( iCondition ) ] = (float)dZ; } } } }

	return true; }

bool CCoalesceCluster::SelectMotifs( const vector<CCoalesceGeneScores>& vecGeneScores,
	const CCoalesceGroupHistograms& HistsCluster, const CCoalesceGroupHistograms& HistsPot, float dPValue,
	const CCoalesceMotifLibrary* pMotifs ) {
	uint32_t	iMotif;

	m_setsMotifs.clear( );
	if( !pMotifs )
		return true;
	for( iMotif = 0; iMotif < pMotifs->GetMotifs( ); ++iMotif )
		if( !AddSignificant( *pMotifs, iMotif, HistsCluster, HistsPot, dPValue ) )
			return false;

	return true; }

bool CCoalesceClusterImpl::AddSignificant( const CCoalesceMotifLibrary& Motifs, uint32_t iMotif,
	const CCoalesceGroupHistograms& HistsCluster, const CCoalesceGroupHistograms& HistsPot, float dPValue ) {
	size_t									iTypeCluster, iTypePot;
	double									dP, dAverage, dZ;
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
			dP = HistSetCluster.CohensD( iMotif, HistSetPot, dAverage, dZ );
			if( dP < dPValue ) {
				SMotifMatch	sMotif( iMotif, strTypeCluster, eSubsequence, (float)dZ, (float)dAverage );

				if( g_CatSleipnir.isInfoEnabled( ) ) {
					g_CatSleipnir.info( "CCoalesceClusterImpl::AddSignificant( %d, %g ) adding at %g:\n%s",
						iMotif, dPValue, dP, sMotif.Save( &Motifs ).c_str( ) );
					g_CatSleipnir.info( "Cluster	%s", HistSetCluster.Save( iMotif ).c_str( ) );
					g_CatSleipnir.info( "Pot	%s", HistSetPot.Save( iMotif ).c_str( ) ); }
				m_setsMotifs.insert( sMotif ); }
			else if( g_CatSleipnir.isDebugEnabled( ) ) {
				g_CatSleipnir.debug( "CCoalesceClusterImpl::AddSignificant( %d, %g ) failed at %g:\n%s",
					iMotif, dPValue, dP, SMotifMatch( iMotif, strTypeCluster, eSubsequence, (float)dZ,
					(float)dAverage ).Save( &Motifs ).c_str( ) );
				g_CatSleipnir.debug( "Cluster	%s", HistSetCluster.Save( iMotif ).c_str( ) );
				g_CatSleipnir.debug( "Pot	%s", HistSetPot.Save( iMotif ).c_str( ) ); } } }

	return true; }

bool CCoalesceCluster::SelectGenes( const CPCL& PCL, const vector<CCoalesceGeneScores>& vecGeneScores,
	const CCoalesceGroupHistograms& HistsCluster, const CCoalesceGroupHistograms& HistsPot,
	CCoalesceCluster& Pot, float dProbability, const CCoalesceMotifLibrary* pMotifs ) {
	static const CCoalesceGeneScores	GeneScoresDummy;
	size_t								iGene;

	CalculateCentroid( PCL );
	Pot.CalculateCentroid( PCL );
	m_setiGenes.clear( );
	Pot.m_setiGenes.clear( );
	for( iGene = 0; iGene < PCL.GetGenes( ); ++iGene )
		if( IsSignificant( iGene, PCL, pMotifs, vecGeneScores.empty( ) ? GeneScoresDummy :
			vecGeneScores[ iGene ], HistsCluster, HistsPot, Pot, dProbability ) ) {
			Add( iGene );
			m_vecdPriors[ iGene ] = 1; }
		else {
			Pot.Add( iGene );
			m_vecdPriors[ iGene ] = dProbability; }

	return true; }

bool CCoalesceClusterImpl::IsSignificant( size_t iGene, const CPCL& PCL, const CCoalesceMotifLibrary* pMotifs,
	const CCoalesceGeneScores& GeneScores, const CCoalesceGroupHistograms& HistsCluster,
	const CCoalesceGroupHistograms& HistsPot, const CCoalesceCluster& Pot, float dProbability ) const {
	long double	dP, dPMotifsGivenIn, dPMotifsGivenOut, dPExpressionGivenIn, dPExpressionGivenOut;
	bool		fClustered;

	fClustered = IsGene( iGene );
// We want P(g in C|S,E) =
// P(S,E|g in C)P(g in C)/P(S,E) =
// P(S|g in C)P(E|g in C)P(g in C)/(P(S,E|g in C) + P(S,E|g notin C)) =
// P(S|g in C)P(E|g in C)P(g in C)/(P(S|g in C)P(E|g in C) + P(S|g notin C)P(E|g notin C))
	if( !( CalculateProbabilityExpression( iGene, PCL, Pot, fClustered, dPExpressionGivenIn,
		dPExpressionGivenOut ) && CalculateProbabilityMotifs( GeneScores, HistsCluster, HistsPot, fClustered,
		dPMotifsGivenIn, dPMotifsGivenOut ) ) )
		return false;
	dP = dPExpressionGivenIn * dPMotifsGivenIn;
	dP = dP * m_vecdPriors[ iGene ] / ( dP + ( dPExpressionGivenOut * dPMotifsGivenOut ) );

	g_CatSleipnir.debug( "CCoalesceClusterImpl::IsSignificant( %s ) is %Lg, exp. p=%Lg vs. %Lg, seq. p=%Lg vs %Lg",
		PCL.GetGene( iGene ).c_str( ), dP, dPExpressionGivenIn, dPExpressionGivenOut, dPMotifsGivenIn,
		dPMotifsGivenOut );
	if( g_CatSleipnir.isDebugEnabled( ) ) {
		set<SMotifMatch>::const_iterator	iterMotif;
		size_t								iType;

		for( iterMotif = m_setsMotifs.begin( ); iterMotif != m_setsMotifs.end( ); ++iterMotif )
			if( ( iType = GeneScores.GetType( iterMotif->m_strType ) ) != -1 )
				g_CatSleipnir.debug( "%g	%s", GeneScores.GetGlobal( iType, iterMotif->m_eSubsequence,
					iterMotif->m_iMotif ), iterMotif->Save( pMotifs ).c_str( ) ); }

	return ( dP > dProbability ); }

bool CCoalesceClusterImpl::CalculateProbabilityExpression( size_t iGene, const CPCL& PCL,
	const CCoalesceCluster& Pot, bool fClustered, long double& dPIn, long double& dPOut ) const {
	set<size_t>::const_iterator	iterCondition;
	float						dGene, dCluster, dPot;
	double						dPCluster, dPPot;
	size_t						iPot;

	iPot = PCL.GetGenes( ) - m_veciPrevGenes.size( );
	dPIn = dPOut = 1;
	for( iterCondition = m_setiConditions.begin( ); iterCondition != m_setiConditions.end( );
		++iterCondition ) {
		if( !m_vecdStdevs[ *iterCondition ] || CMeta::IsNaN( dGene = PCL.Get( iGene, *iterCondition ) ) )
			continue;
		dCluster = m_vecdCentroid[ *iterCondition ];
		dPot = Pot.m_vecdCentroid[ *iterCondition ];
		if( fClustered )
			dCluster = ( ( dCluster * m_veciPrevGenes.size( ) ) - dGene ) / ( m_veciPrevGenes.size( ) - 1 );
		else
			dPot = ( ( dPot * iPot ) - dGene ) / ( iPot - 1 );

		dPCluster = CStatistics::NormalPDF( dGene, dCluster, m_vecdStdevs[ *iterCondition ] );
		dPPot = CStatistics::NormalPDF( dGene, dPot, m_vecdStdevs[ *iterCondition ] );
		AdjustProbabilities( m_vecdZs[ *iterCondition ], dPCluster, dPPot );
		dPIn *= dPCluster;
		dPOut *= dPPot; }

	return true; }

bool CCoalesceClusterImpl::CalculateProbabilityMotifs( const CCoalesceGeneScores& GeneScores,
	const CCoalesceGroupHistograms& HistsCluster, const CCoalesceGroupHistograms& HistsPot,
	bool fClustered, long double& dPIn, long double& dPOut ) const {
	set<SMotifMatch>::const_iterator	iterMotif;
	size_t								iType, iCluster, iPot;
	unsigned short						sCluster, sPot;
	double								dPCluster, dPPot;
	float								dGene;

	dPIn = dPOut = 0.5;
	if( m_setsMotifs.empty( ) )
		return true;

// TODO: this is the slowest part of the entire algorithm
	for( iterMotif = m_setsMotifs.begin( ); iterMotif != m_setsMotifs.end( ); ++iterMotif ) {
		if( ( ( iType = GeneScores.GetType( iterMotif->m_strType ) ) == -1 ) ||
			( !GeneScores.GetMotifs( iType, iterMotif->m_eSubsequence ) ) ||
			( HistsCluster.GetType( iterMotif->m_strType ) == -1 ) ||
			( HistsPot.GetType( iterMotif->m_strType ) == -1 ) )
			continue;
		dGene = GeneScores.GetGlobal( iType, iterMotif->m_eSubsequence, iterMotif->m_iMotif );
		{
			const CCoalesceHistogramSet<>&	HistSetCluster	= HistsCluster.Get( iterMotif->m_strType,
																iterMotif->m_eSubsequence );
			const CCoalesceHistogramSet<>&	HistSetPot		= HistsPot.Get( iterMotif->m_strType,
																iterMotif->m_eSubsequence );

			sCluster = HistSetCluster.Get( iterMotif->m_iMotif, dGene );
			iCluster = HistSetCluster.GetTotal( );
			sPot = HistSetPot.Get( iterMotif->m_iMotif, dGene );
			iPot = HistSetPot.GetTotal( );

			if( fClustered ) {
				sCluster--;
				iCluster--; }
			else {
				sPot--;
				iPot--; }
			dPCluster = (double)( sCluster + 1 ) / ( iCluster + HistSetCluster.GetEdges( ) );
			dPPot = (double)( sPot + 1 ) / ( iPot + HistSetPot.GetEdges( ) );
		}
		AdjustProbabilities( iterMotif->m_dZ, dPCluster, dPPot );
		dPIn *= dPCluster;
		dPOut *= dPPot; }

	return true; }

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
		ostm << '\t' << PCL.GetExperiment( *iterID );
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

void CCoalesceImpl::Save( const vector<CCoalesceGeneScores>& vecGeneScores ) const {
	ofstream	ofsm;
	size_t		i;
	uint32_t	iSize;

	ofsm.open( m_strSequenceCache.c_str( ), ios_base::binary );
	if( !ofsm.is_open( ) )
		return;
	iSize = vecGeneScores.size( );
	ofsm.write( (const char*)&iSize, sizeof(iSize) );
	for( i = 0; i < vecGeneScores.size( ); ++i )
		vecGeneScores[ i ].Save( ofsm ); }

void CCoalesceImpl::Open( vector<CCoalesceGeneScores>& vecGeneScores ) const {
	ifstream	ifsm;
	size_t		i;
	uint32_t	iSize;

	ifsm.open( m_strSequenceCache.c_str( ), ios_base::binary );
	if( !ifsm.is_open( ) )
		return;
	ifsm.read( (char*)&iSize, sizeof(iSize) );
	vecGeneScores.resize( iSize );
	for( i = 0; i < vecGeneScores.size( ); ++i )
		vecGeneScores[ i ].Open( ifsm ); }

bool CCoalesceImpl::CombineMotifs( const CFASTA& FASTA, const vector<size_t>& veciPCL2FASTA,
	const CCoalesceCluster& Cluster, vector<CCoalesceGeneScores>& vecGeneScores,
	CCoalesceGroupHistograms& HistsCluster, CCoalesceGroupHistograms& HistsPot ) const {
	set<SMotifMatch>::const_iterator	iterMotifOne, iterMotifTwo;
	uint32_t							iMotif;
	size_t								i, j;
	vector<float>						vecdScores;
	vector<size_t>						veciLengths;

	if( vecGeneScores.empty( ) || !m_pMotifs || ( Cluster.GetMotifs( ).size( ) > m_iSizeMaximum ) )
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

			for( i = 0; i < veciPCL2FASTA.size( ); ++i )
				if( veciPCL2FASTA[ i ] != -1 ) {
					vector<SFASTASequence>	vecsSequences;

					if( FASTA.Get( veciPCL2FASTA[ i ], vecsSequences ) ) {
						for( j = 0; j < vecsSequences.size( ); ++j )
							if( !vecGeneScores[ i ].Add( *m_pMotifs, vecsSequences[ j ], iMotif, vecdScores,
								veciLengths ) )
								return false; }
						if( Cluster.IsGene( i ) )
							HistsCluster.Add( vecGeneScores[ i ], false, iMotif );
						else
							HistsPot.Add( vecGeneScores[ i ], false, iMotif ); } }

	return true; }

bool CCoalesceImpl::InitializeDatasets( const CPCL& PCL, vector<SDataset>& vecsDatasets ) {
	size_t			i, j, k, iDataset, iCondition;
	bool			fMultivariate;
	vector<bool>	vecfDataset;
	float			dOne, dTwo;

	fMultivariate = false;
	vecfDataset.resize( PCL.GetExperiments( ) );
	for( i = 0; i < m_vecsDatasets.size( ); ++i ) {
		const SDataset&	sDataset	= m_vecsDatasets[ i ];
		vector<size_t>	veciTo;

		veciTo.reserve( sDataset.GetConditions( ) );
		for( j = 0; j < sDataset.GetConditions( ); ++j )
			if( ( iCondition = sDataset.GetCondition( j ) ) < PCL.GetExperiments( ) ) {
				vecfDataset[ iCondition ] = true;
				veciTo.push_back( iCondition ); }
		if( !veciTo.empty( ) ) {
			vecsDatasets.push_back( SDataset( veciTo ) );
			if( veciTo.size( ) > 1 )
				fMultivariate = true; } }
	for( i = 0; i < vecfDataset.size( ); ++i )
		if( !vecfDataset[ i ] )
			vecsDatasets.push_back( SDataset( i ) );
	if( !fMultivariate )
		return true;

	for( iDataset = 0; iDataset < vecsDatasets.size( ); ++iDataset ) {
		SDataset&		sDataset		= vecsDatasets[ iDataset ];
		CDataMatrix&	MatCovariance	= sDataset.m_MatCovChol;
		vector<float>	vecdAves;

		if( sDataset.GetConditions( ) == 1 )
			continue;
		sDataset.m_vecdStdevs.resize( sDataset.GetConditions( ) );
		MatCovariance.Initialize( sDataset.GetConditions( ), sDataset.GetConditions( ) );
		MatCovariance.Clear( );
		vecdAves.resize( sDataset.GetConditions( ) );
		for( i = 0; i < PCL.GetGenes( ); ++i )
			for( j = 0; j < vecdAves.size( ); ++j )
				if( !CMeta::IsNaN( dOne = PCL.Get( i, sDataset.GetCondition( j ) ) ) )
					vecdAves[ j ] += dOne;
		for( i = 0; i < vecdAves.size( ); ++i )
			vecdAves[ i ] /= PCL.GetGenes( );
		for( i = 0; i < PCL.GetGenes( ); ++i )
			for( j = 0; j < sDataset.GetConditions( ); ++j ) {
				if( CMeta::IsNaN( dOne = PCL.Get( i, sDataset.GetCondition( j ) ) ) )
					continue;
				dOne -= vecdAves[ j ];
				for( k = j; k < sDataset.GetConditions( ); ++k )
					if( !CMeta::IsNaN( dTwo = PCL.Get( i, sDataset.GetCondition( k ) ) ) )
						MatCovariance.Get( j, k ) += dOne * ( dTwo - vecdAves[ k ] ); }
		for( i = 0; i < MatCovariance.GetRows( ); ++i ) {
			for( j = i; j < MatCovariance.GetColumns( ); ++j )
				MatCovariance.Get( i, j ) /= PCL.GetGenes( );
			sDataset.m_vecdStdevs[ i ] = sqrt( MatCovariance.Get( i, i ) ); }
		CStatistics::CholeskyDecomposition( MatCovariance ); }

	return true; }

bool CCoalesceImpl::InitializeGeneScores( const CPCL& PCL, const CFASTA& FASTA, CPCL& PCLCopy,
	vector<size_t>& veciPCL2FASTA, vector<CCoalesceGeneScores>& vecGeneScores ) {
	size_t	i, j;

	PCLCopy.Open( PCL );
	PCLCopy.Normalize( CPCL::ENormalizeColumn );
	if( FASTA.GetGenes( ) ) {
		veciPCL2FASTA.resize( PCL.GetGenes( ) );
		for( i = 0; i < veciPCL2FASTA.size( ); ++i )
			veciPCL2FASTA[ i ] = FASTA.GetGene( PCL.GetGene( i ) ); }
	if( !m_strSequenceCache.empty( ) )
		Open( vecGeneScores );
	if( vecGeneScores.empty( ) && FASTA.GetGenes( ) ) {
		vector<vector<unsigned short> >	vecvecsCounts;
		vector<size_t>					veciLengths;

		if( !m_pMotifs ) {
			m_fMotifs = true;
			m_pMotifs = new CCoalesceMotifLibrary( m_iK ); }
		vecGeneScores.resize( PCL.GetGenes( ) );
		for( i = 0; i < veciPCL2FASTA.size( ); ++i )
			if( veciPCL2FASTA[ i ] != -1 ) {
				vector<SFASTASequence>	vecsSequences;

				if( FASTA.Get( veciPCL2FASTA[ i ], vecsSequences ) )
					for( j = 0; j < vecsSequences.size( ); ++j )
						if( !vecGeneScores[ i ].Add( *m_pMotifs, vecsSequences[ j ], vecvecsCounts,
							veciLengths ) )
							return false; }
		if( !m_strSequenceCache.empty( ) )
			Save( vecGeneScores ); }
	if( vecGeneScores.empty( ) )
		Clear( );
	else {
		m_pMotifs->SetPenaltyGap( m_dPenaltyGap );
		m_pMotifs->SetPenaltyMismatch( m_dPenaltyMismatch ); }

	return true; }

bool CCoalesce::Cluster( const CPCL& PCL, const CFASTA& FASTA, vector<CCoalesceCluster>& vecClusters ) {
	CPCL						PCLCopy;
	size_t						i;
	vector<size_t>				veciPCL2FASTA;
	vector<CCoalesceGeneScores>	vecGeneScores;
	vector<SDataset>			vecsDatasets;
	set<pair<size_t, size_t> >	setpriiSeeds;

	if( !( InitializeDatasets( PCL, vecsDatasets ) &&
		InitializeGeneScores( PCL, FASTA, PCLCopy, veciPCL2FASTA, vecGeneScores ) ) )
		return false;
	while( true ) {
		CCoalesceCluster			Cluster, Pot;
		CCoalesceGroupHistograms	HistsCluster( GetMotifCount( ), GetBins( ), 1.0f / GetBasesPerMatch( ) );
		CCoalesceGroupHistograms	HistsPot( GetMotifCount( ), GetBins( ), 1.0f / GetBasesPerMatch( ) );

		for( i = 0; i < PCLCopy.GetGenes( ); ++i )
			Pot.Add( i );
		Pot.CalculateHistograms( vecGeneScores, HistsPot, NULL );
		if( !Cluster.Initialize( PCLCopy, Pot, setpriiSeeds, GetFractionCorrelation( ),
			GetPValueCorrelation( ) ) )
			break;
		g_CatSleipnir.notice( "CCoalesce::Cluster( ) initialized %d genes", Cluster.GetGenes( ).size( ) );
		if( Cluster.GetGenes( ).size( ) < GetSizeMinimum( ) )
			continue;
		while( !( Cluster.IsConverged( ) || Cluster.IsEmpty( ) ) ) {
			Cluster.CalculateHistograms( vecGeneScores, HistsCluster, &HistsPot );
			Cluster.Snapshot( vecGeneScores, HistsCluster );
			Pot.Snapshot( vecGeneScores, HistsPot );
			if( !Cluster.SelectConditions( PCLCopy, vecsDatasets, Pot, GetPValueCondition( ) ) )
				return false;
			if( Cluster.IsEmpty( ) ) {
				g_CatSleipnir.info( "CCoalesce::Cluster( ) selected no conditions" );
				break; }
			if( ( Cluster.GetGenes( ).size( ) >= GetSizeMinimum( ) ) &&
				!( CombineMotifs( FASTA, veciPCL2FASTA, Cluster, vecGeneScores, HistsCluster, HistsPot ) &&
				Cluster.SelectMotifs( vecGeneScores, HistsCluster, HistsPot, GetPValueMotif( ),
				GetMotifs( ) ) ) )
				return false;
			if( !Cluster.SelectGenes( PCLCopy, vecGeneScores, HistsCluster, HistsPot, Pot,
				GetProbabilityGene( ), GetMotifs( ) ) )
				return false;
			g_CatSleipnir.info( "CCoalesce::Cluster( ) processed %d genes, %d conditions, %d motifs",
				Cluster.GetGenes( ).size( ), Cluster.GetConditions( ).size( ), Cluster.GetMotifs( ).size( ) ); }
		if( Cluster.IsConverged( ) && ( Cluster.GetGenes( ).size( ) >= GetSizeMinimum( ) ) ) {
			g_CatSleipnir.info( "CCoalesce::Cluster( ) finalizing cluster" );
			setpriiSeeds.clear( );
			if( IsOutputIntermediate( ) )
				Cluster.Save( GetDirectoryIntermediate( ), vecClusters.size( ), PCLCopy, GetMotifs( ) );
			vecClusters.push_back( Cluster );
			Cluster.Subtract( PCLCopy );
			Cluster.Subtract( vecGeneScores ); }
		else
			g_CatSleipnir.info( "CCoalesce::Cluster( ) discarding cluster" ); }

	return true; }

}
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
		const string&		strKMer	= strSequence.substr( i, Motifs.GetK( ) );
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
				vecdScores[ iterMotif->second ] = max( 0.0f, vecdScores[ iterMotif->second ] -
					sMotif.m_dAverage ); }
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
#if 0
		static size_t	s_iOne	= 0;
		static size_t	s_iTwo	= 0;
		if( ++s_iTwo >= PCL.GetGenes( ) ) {
			s_iTwo = 0;
			s_iOne++; }
		iOne = s_iOne;
		iTwo = s_iTwo;
#else
		iOne = rand( ) % PCL.GetGenes( );
		do
			iTwo = rand( ) % PCL.GetGenes( );
		while( iOne == iTwo );
#endif

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
	for( iGene = 0; iGene < PCL.GetGenes( ); ++iGene )
		if( !IsGene( iGene ) &&
			( ( dR = CMeasurePearson::Pearson( &m_vecdCentroid.front( ), PCL.GetExperiments( ),
			PCL.Get( iGene ), PCL.GetExperiments( ), IMeasure::EMapNone, NULL, NULL ) ) > 0 ) &&
			( ( CStatistics::PValuePearson( dR, PCL.GetExperiments( ) ) * ( PCL.GetGenes( ) -
			m_setiGenes.size( ) ) ) < dPValue ) )
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
	m_vecdCentroid.resize( m_veciCounts.size( ) );
	fill( m_vecdCentroid.begin( ), m_vecdCentroid.end( ), 0.0f );
	for( iterGene = m_setiGenes.begin( ); iterGene != m_setiGenes.end( ); ++iterGene )
		for( i = 0; i < m_veciCounts.size( ); ++i )
			if( !CMeta::IsNaN( d = PCL.Get( *iterGene, i ) ) ) {
				m_veciCounts[ i ]++;
				m_vecdCentroid[ i ] += d; }
	for( i = 0; i < m_veciCounts.size( ); ++i )
		if( m_veciCounts[ i ] )
			m_vecdCentroid[ i ] /= m_veciCounts[ i ]; }

bool CCoalesceCluster::SelectConditions( const CPCL& PCL, const CCoalesceCluster& Pot, float dPValue ) {
	vector<float>				vecdCluster, vecdPot;
	size_t						iGene, iCondition;
	set<size_t>::const_iterator	iterGene;
	double						dP, dZ;

	dPValue /= PCL.GetExperiments( );
	m_vecdZs.resize( PCL.GetExperiments( ) );
	m_setiConditions.clear( );
	vecdCluster.resize( m_setiGenes.size( ) );
	vecdPot.resize( Pot.m_setiGenes.size( ) );
	for( iCondition = 0; iCondition < PCL.GetExperiments( ); ++iCondition ) {
		for( iterGene = m_setiGenes.begin( ),iGene = 0; iterGene != m_setiGenes.end( ); ++iterGene,++iGene )
			vecdCluster[ iGene ] = PCL.Get( *iterGene, iCondition );
		for( iterGene = Pot.m_setiGenes.begin( ),iGene = 0; iterGene != Pot.m_setiGenes.end( );
			++iterGene,++iGene )
			vecdPot[ iGene ] = PCL.Get( *iterGene, iCondition );
		dZ = CStatistics::CohensD( vecdCluster, vecdPot );
		dP = CStatistics::ZTest( dZ, vecdCluster.size( ) );
		if( dP < dPValue ) {
			g_CatSleipnir.info( "CCoalesceCluster::SelectConditions( %g ) selected %d at %g, z=%g",
				dPValue, iCondition, dP, dZ );
			m_setiConditions.insert( iCondition );
			m_vecdZs[ iCondition ] = (float)dZ; } }

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
// I have no justification at all for this multiple hypothesis correction, but it seems to work
			dP = HistSetCluster.CohensD( iMotif, HistSetPot, dAverage, dZ ) *
				sqrt( (double)HistSetCluster.GetTotal( ) );
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
	long double	dP, dPInGivenMotifs, dPOutGivenMotifs, dPInGivenExpression, dPOutGivenExpression;
	bool		fClustered;

	fClustered = IsGene( iGene );
	if( !( CalculateProbabilityExpression( iGene, PCL, Pot, fClustered, dPInGivenExpression,
		dPOutGivenExpression ) && CalculateProbabilityMotifs( GeneScores, HistsCluster, HistsPot, fClustered,
		dPInGivenMotifs, dPOutGivenMotifs ) ) )
		return false;
	dP = dPInGivenExpression * dPInGivenMotifs;
	dP = dP / ( dP + ( dPOutGivenExpression * dPOutGivenMotifs ) );

	g_CatSleipnir.debug( "CCoalesceClusterImpl::IsSignificant( %s ) is %Lg, exp. p=%Lg vs. %Lg, seq. p=%Lg vs %Lg",
		PCL.GetGene( iGene ).c_str( ), dP, dPInGivenExpression, dPOutGivenExpression, dPInGivenMotifs,
		dPOutGivenMotifs );
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
		dGene = PCL.Get( iGene, *iterCondition );
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

	dPIn = dPOut = 1;
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

bool CCoalesce::Cluster( const CPCL& PCL, const CFASTA& FASTA, vector<CCoalesceCluster>& vecClusters ) {
	CPCL						PCLCopy;
	size_t						i, j;
	vector<size_t>				veciPCL2FASTA;
	vector<CCoalesceGeneScores>	vecGeneScores;

	PCLCopy.Open( PCL );
	PCLCopy.Normalize( CPCL::ENormalizeColumn );
	if( FASTA.GetGenes( ) ) {
		veciPCL2FASTA.resize( PCL.GetGenes( ) );
		for( i = 0; i < veciPCL2FASTA.size( ); ++i )
			veciPCL2FASTA[ i ] = FASTA.GetGene( PCL.GetGene( i ) ); }
	if( !GetSequenceCache( ).empty( ) )
		Open( vecGeneScores );
	if( vecGeneScores.empty( ) && FASTA.GetGenes( ) ) {
		vector<vector<unsigned short> >	vecvecsCounts;
		vector<size_t>					veciLengths;

		if( !GetMotifs( ) ) {
			m_fMotifs = true;
			m_pMotifs = new CCoalesceMotifLibrary( GetK( ) ); }
		vecGeneScores.resize( PCL.GetGenes( ) );
		for( i = 0; i < veciPCL2FASTA.size( ); ++i )
			if( veciPCL2FASTA[ i ] != -1 ) {
				vector<SFASTASequence>	vecsSequences;

				if( FASTA.Get( veciPCL2FASTA[ i ], vecsSequences ) )
					for( j = 0; j < vecsSequences.size( ); ++j )
						if( !vecGeneScores[ i ].Add( *m_pMotifs, vecsSequences[ j ], vecvecsCounts,
							veciLengths ) )
							return false; }
		if( !GetSequenceCache( ).empty( ) )
			Save( vecGeneScores ); }
	if( vecGeneScores.empty( ) )
		Clear( );
	else {
		m_pMotifs->SetPenaltyGap( GetPenaltyGap( ) );
		m_pMotifs->SetPenaltyMismatch( GetPenaltyMismatch( ) ); }
	while( true ) {
		CCoalesceCluster			Cluster, Pot;
		CCoalesceGroupHistograms	HistsCluster( GetMotifCount( ), GetBins( ), 1.0f / GetBasesPerMatch( ) );
		CCoalesceGroupHistograms	HistsPot( GetMotifCount( ), GetBins( ), 1.0f / GetBasesPerMatch( ) );

		for( i = 0; i < PCLCopy.GetGenes( ); ++i )
			Pot.Add( i );
		Pot.CalculateHistograms( vecGeneScores, HistsPot, NULL );
		if( !Cluster.Initialize( PCLCopy, Pot, GetPValueCorrelation( ) ) )
			break;
		g_CatSleipnir.notice( "CCoalesce::Cluster( ) initialized %d genes", Cluster.GetGenes( ).size( ) );
		while( !( Cluster.IsConverged( ) || Cluster.IsEmpty( ) ) ) {
			Cluster.CalculateHistograms( vecGeneScores, HistsCluster, &HistsPot );
			Cluster.Snapshot( vecGeneScores, HistsCluster );
			Pot.Snapshot( vecGeneScores, HistsPot );
			if( !Cluster.SelectConditions( PCLCopy, Pot, GetPValueCondition( ) ) )
				return false;
			if( Cluster.IsEmpty( ) )
				break;
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
		g_CatSleipnir.info( "CCoalesce::Cluster( ) finalized cluster" );
		if( Cluster.IsConverged( ) && ( Cluster.GetGenes( ).size( ) >= GetSizeMinimum( ) ) ) {
			if( IsOutputIntermediate( ) )
				Cluster.Save( GetDirectoryIntermediate( ), vecClusters.size( ), PCLCopy, GetMotifs( ) );
			vecClusters.push_back( Cluster );
			Cluster.Subtract( PCLCopy );
			Cluster.Subtract( vecGeneScores ); } }

	return true; }

}

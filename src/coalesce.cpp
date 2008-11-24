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
#include "genome.h"
#include "clusthierarchical.h"

namespace Sleipnir {

const char	CCoalesceClusterImpl::c_szMotifs[]			= "_motifs.txt";
const char*	CCoalesceSequencerBase::c_aszSubsequences[]	= {"total", "introns", "exons", NULL};

// SCoalesceModifiers

void SCoalesceModifiers::Initialize( const CPCL& PCL ) {
	size_t	i, j;

	m_vecveciPCL2Wiggles.resize( m_vecpWiggles.size( ) );
	for( i = 0; i < m_vecveciPCL2Wiggles.size( ); ++i ) {
		m_vecveciPCL2Wiggles[ i ].resize( PCL.GetGenes( ) );
		for( j = 0; j < m_vecveciPCL2Wiggles[ i ].size( ); ++j )
			m_vecveciPCL2Wiggles[ i ][ j ] = m_vecpWiggles[ i ]->GetGene( PCL.GetGene( j ) ); } }

// SCoalesceModifierCache

void SCoalesceModifierCache::Get( size_t iPCL ) {
	size_t	i;

	m_vecvecsWiggles.resize( m_Modifiers.m_vecpWiggles.size( ) );
	for( i = 0; i < m_vecvecsWiggles.size( ); ++i ) {
		m_vecvecsWiggles[ i ].clear( );
		m_Modifiers.m_vecpWiggles[ i ]->Get( m_Modifiers.m_vecveciPCL2Wiggles[ i ][ iPCL ],
			m_vecvecsWiggles[ i ] ); } }

void SCoalesceModifierCache::SetType( const std::string& strType ) {
	size_t	i, j;

	m_veciWiggleTypes.resize( m_vecvecsWiggles.size( ) );
	for( i = 0; i < m_vecvecsWiggles.size( ); ++i ) {
		for( j = 0; j < m_vecvecsWiggles[ i ].size( ); ++j )
			if( strType == m_vecvecsWiggles[ i ][ j ].m_strType )
				break;
		m_veciWiggleTypes[ i ] = ( j < m_vecvecsWiggles[ i ].size( ) ) ? j : -1; } }

void SCoalesceModifierCache::InitializeWeight( size_t iK, size_t iOffset ) {
	size_t	i, j;

	m_dWeight = 0;
	for( i = 0; i < m_vecvecsWiggles.size( ); ++i )
		if( ( j = m_veciWiggleTypes[ i ] ) != -1 ) {
			const vector<float>&	vecdWiggle	= m_vecvecsWiggles[ i ][ j ].m_vecdValues;

			for( j = 0; ( j < iK ) && ( ( iOffset + j ) < vecdWiggle.size( ) ); ++j )
				m_dWeight += vecdWiggle[ iOffset + j ];
			for( ; j < iK; ++j )
				m_dWeight += 1; } }

void SCoalesceModifierCache::AddWeight( size_t iK, size_t iOffset, size_t iDelta ) {
	size_t	i, j;
	float	dSub, dAdd;

	for( i = 0; i < m_vecvecsWiggles.size( ); ++i )
		if( ( j = m_veciWiggleTypes[ i ] ) != -1 ) {
			const std::vector<float>&	vecdWiggle	= m_vecvecsWiggles[ i ][ j ].m_vecdValues;

			dSub = ( ( iOffset + iDelta ) < vecdWiggle.size( ) ) ? vecdWiggle[ iOffset + iDelta ] : 1;
			if( iK ) {
				j = iOffset + iDelta + iK;
				dAdd = ( j < vecdWiggle.size( ) ) ? vecdWiggle[ j ] : 1; }
			else {
				dAdd = dSub;
				dSub = 0; }
			m_dWeight -= dSub;
			m_dWeight += dAdd; } }

// SCoalesceDataset

bool SCoalesceDataset::CalculateCovariance( const CPCL& PCL ) {
	size_t			i, j, k;
	vector<bool>	vecfDataset;
	float			dOne, dTwo;
	vector<float>	vecdAves;
	CDataMatrix		MatSigma;
	vector<size_t>	veciIndices;
	bool			fEven;

	if( GetConditions( ) == 1 )
		return false;

	m_vecdStdevs.resize( GetConditions( ) );
	MatSigma.Initialize( GetConditions( ), GetConditions( ) );
	MatSigma.Clear( );
	vecdAves.resize( GetConditions( ) );
	for( i = 0; i < PCL.GetGenes( ); ++i )
		for( j = 0; j < vecdAves.size( ); ++j )
			if( !CMeta::IsNaN( dOne = PCL.Get( i, GetCondition( j ) ) ) )
				vecdAves[ j ] += dOne;
	for( i = 0; i < vecdAves.size( ); ++i )
		vecdAves[ i ] /= PCL.GetGenes( );
	for( i = 0; i < PCL.GetGenes( ); ++i )
		for( j = 0; j < GetConditions( ); ++j ) {
			if( CMeta::IsNaN( dOne = PCL.Get( i, GetCondition( j ) ) ) )
				continue;
			dOne -= vecdAves[ j ];
			for( k = j; k < GetConditions( ); ++k )
				if( !CMeta::IsNaN( dTwo = PCL.Get( i, GetCondition( k ) ) ) )
					MatSigma.Get( j, k ) += dOne * ( dTwo - vecdAves[ k ] ); }
	for( i = 0; i < MatSigma.GetRows( ); ++i ) {
		for( j = i; j < MatSigma.GetColumns( ); ++j )
			MatSigma.Set( j, i, MatSigma.Get( i, j ) /= PCL.GetGenes( ) );
		m_vecdStdevs[ i ] = sqrt( MatSigma.Get( i, i ) ); }

	m_MatSigmaChol.Open( MatSigma );
	CStatistics::CholeskyDecomposition( m_MatSigmaChol );

	CStatistics::MatrixLUDecompose( MatSigma, veciIndices, fEven );
	CStatistics::MatrixLUInvert( MatSigma, veciIndices, m_MatSigmaInv );
	m_dSigmaDetSqrt = sqrt( CStatistics::MatrixLUDeterminant( MatSigma, fEven ) );

	return true; }

// CCoalesceMotifLibrary

// Order independent, but complements must be adjacent
const char	CCoalesceMotifLibraryImpl::c_acBases[]	= "ATCG";

CCoalesceMotifLibraryImpl::~CCoalesceMotifLibraryImpl( ) {
	size_t	i;

	for( i = 0; i < m_vecpPSTs.size( ); ++i )
		delete m_vecpPSTs[ i ]; }

/*!
 * \brief
 * Calculates the length-normalized match strength of the given motif against the appropriate number of
 * characters in the input sequence at the requested offset.
 * 
 * \param strSequence
 * Sequence against which motif is matched.
 * 
 * \param iMotif
 * ID of motif to be matched.
 * 
 * \param iOffset
 * Zero-based offset within strSequence at which the match is performed.
 * 
 * \param sModifiers
 * A modifier cache containing any prior weights to be incorporated into the match.
 * 
 * \returns
 * Length-normalized strength of motif match against the given sequence and offset.
 * 
 * Write detailed description for GetMatch here.
 * 
 * \remarks
 * iMotif must represent a valid motif for the current library, and iOffset must fall within strSequence,
 * although motifs extending from a valid iOffset past the end of the sequence will be handled appropriately.
 * Only PST motifs are currently supported, as there should never be any need to match non-PST motifs at
 * runtime, but support for kmers and RCs could be added in a straightforward manner.
 * 
 * \see
 * GetMatches
 */
float CCoalesceMotifLibrary::GetMatch( const std::string& strSequence, uint32_t iMotif, size_t iOffset,
	SCoalesceModifierCache& sModifiers ) const {
	size_t		i;
	float		dRet;
	const CPST*	pPST;

// TODO: this could in theory be implemented
	if( ( GetType( iMotif ) != ETypePST ) || !( pPST = GetPST( iMotif ) ) ) {
		g_CatSleipnir.error( "CCoalesceMotifLibrary::GetMatch( %s, %d, %d ) attempted to match a non-PST motif",
			strSequence.c_str( ), iMotif, iOffset );
		return CMeta::GetNaN( ); }

	sModifiers.InitializeWeight( 0, 0 );
	dRet = 0;
	for( i = 1; i < pPST->GetDepth( ); ++i ) {
		sModifiers.AddWeight( 0, iOffset, i - 1 );
		dRet += pPST->GetMatch( strSequence, pPST->GetDepth( ) - i ) * sModifiers.GetWeight( ) / i; }
	sModifiers.AddWeight( 0, iOffset, i - 1 );
	for( i = 0; i < strSequence.length( ); ++i ) {
		dRet += pPST->GetMatch( strSequence.substr( i ) ) * sModifiers.GetWeight( ) / pPST->GetDepth( );
		sModifiers.AddWeight( pPST->GetDepth( ), iOffset, i ); }
	if( dRet < 0 ) {
		g_CatSleipnir.error( "CCoalesceMotifLibrary::GetMatch( %s, %d, %d ) found negative score: %g",
			strSequence.c_str( ), iMotif, iOffset, dRet );
		return CMeta::GetNaN( ); }

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

uint32_t CCoalesceMotifLibrary::Open( const std::string& strMotif ) {
	uint32_t	iMotif;
	CPST*		pPST;

	if( strMotif.length( ) == GetK( ) && !IsIgnorableKMer( strMotif ) )
		return KMer2ID( strMotif );
	if( ( strMotif.length( ) == ( ( 2 * GetK( ) ) + 1 ) ) &&
		!IsIgnorableKMer( strMotif.substr( 0, GetK( ) ) ) &&
		( strMotif.substr( 0, GetK( ) ) == GetReverseComplement( strMotif.substr( GetK( ) + 1 ) ) ) &&
		( ( iMotif = KMer2ID( strMotif.substr( 0, GetK( ) ) ) ) != -1 ) )
		return ( GetBaseRCs( ) + m_veciKMer2RC[ iMotif ] );

	pPST = CreatePST( iMotif );
	if( !pPST->Open( strMotif ) ) {
		delete pPST;
		m_vecpPSTs.pop_back( );
		return -1; }

	return iMotif; }

float CCoalesceMotifLibraryImpl::AlignKMers( const std::string& strOne, const std::string& strTwo,
	float dCutoff ) const {
	int	iOffset;

	return Align( strOne, strTwo, dCutoff, iOffset ); }

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

float CCoalesceMotifLibraryImpl::AlignKMerRC( const std::string& strKMer, uint32_t iRC, float dCutoff ) const {
	string	strOne, strTwo;
	float	dOne, dTwo;
	int		iOne, iTwo;

	strOne = GetRCOne( iRC );
	strTwo = GetReverseComplement( strOne );
	dOne = Align( strKMer, strOne, dCutoff, iOne );
	dTwo = Align( strKMer, strTwo, dCutoff, iTwo );

	return min( dOne, dTwo ); }

uint32_t CCoalesceMotifLibraryImpl::MergeKMerRC( const std::string& strKMer, uint32_t iRC, float dCutoff ) {
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

float CCoalesceMotifLibraryImpl::AlignRCs( uint32_t iOne, uint32_t iTwo, float dCutoff ) const {
	SCrossRCs	asCrosses[ 4 ];
	size_t		i;
	float		dMin;

	asCrosses[ 0 ].m_strOne = asCrosses[ 1 ].m_strOne = GetRCOne( iOne );
	asCrosses[ 0 ].m_strTwo = asCrosses[ 2 ].m_strTwo = GetRCOne( iTwo );
	asCrosses[ 1 ].m_strTwo = asCrosses[ 3 ].m_strTwo = GetReverseComplement( asCrosses[ 0 ].m_strTwo );
	asCrosses[ 2 ].m_strOne = asCrosses[ 3 ].m_strOne = GetReverseComplement( asCrosses[ 0 ].m_strOne );
	dMin = FLT_MAX;
	for( i = 0; i < ARRAYSIZE(asCrosses); ++i ) {
		asCrosses[ i ].m_dScore = Align( asCrosses[ i ].m_strOne, asCrosses[ i ].m_strTwo, dCutoff,
			asCrosses[ i ].m_iOffset );
		if( asCrosses[ i ].m_dScore < dMin )
			dMin = asCrosses[ i ].m_dScore; }

	return dMin; }

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

float CCoalesceMotifLibraryImpl::AlignKMerPST( const std::string& strKMer, const CPST& PSTIn,
	float dCutoff ) const {
	int	iOffset;

	return PSTIn.Align( strKMer, m_dPenaltyGap, m_dPenaltyMismatch, dCutoff, iOffset ); }

uint32_t CCoalesceMotifLibraryImpl::MergeKMerPST( const std::string& strKMer, const CPST& PSTIn,
	float dCutoff ) {
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

float CCoalesceMotifLibraryImpl::AlignRCPST( uint32_t iRC, const CPST& PSTIn, float dCutoff ) const {
	int		iOne, iTwo;
	string	strOne, strTwo;
	float	dOne, dTwo;

	strOne = GetRCOne( iRC );
	strTwo = GetReverseComplement( strOne );
	dOne = PSTIn.Align( strOne, m_dPenaltyGap, m_dPenaltyMismatch, dCutoff, iOne );
	dTwo = PSTIn.Align( strTwo, m_dPenaltyGap, m_dPenaltyMismatch, dCutoff, iTwo );

	return min( dOne, dTwo ); }

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

float CCoalesceMotifLibraryImpl::AlignPSTs( const CPST& PSTOne, const CPST& PSTTwo, float dCutoff ) const {
	int	iOffset;

	return PSTOne.Align( PSTTwo, m_dPenaltyGap, m_dPenaltyMismatch, dCutoff, iOffset ); }

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
			Add( eSubsequence, veciMotifs[ j ], Motifs.GetMotifs( ), vecvecdCounts, sModifiers.GetWeight( ) /
				Motifs.GetK( ) );
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
	float**	aadScores;
	float*	adScores;

	if( ( ( iType = GetType( sMotif.m_strType ) ) != -1 ) &&
		( aadScores = Get( iType, sMotif.m_eSubsequence ) ) && ( adScores = aadScores[ iGene ] ) )
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
					vecdEdges[ i ] = ( i - 1 ) * ( m_dStep / 2 );
				Histograms.Initialize( GetMotifs( ), vecdEdges ); }
// unlock
			if( iMotifThem == -1 ) {
				for( iMotifUs = 0; iMotifUs < GeneScores.GetMotifs( ); ++iMotifUs )
					if( ( dValue = adScores[ iMotifUs ] ) && !Histograms.Add( iMotifUs, dValue, sDelta ) )
							return false; }
			else if( dValue = adScores[ iMotifThem ] )
				if( !Histograms.Add( iMotifThem, dValue, sDelta ) )
					return false; } }

	return true; }

void CCoalesceGroupHistograms::Save( std::ostream& ostm, const CCoalesceMotifLibrary* pMotifs ) const {
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

bool SMotifMatch::Open( istream& istm, CCoalesceMotifLibrary& Motifs ) {
	char			acLine[ 1024 ];
	vector<string>	vecstrLine;

	istm.getline( acLine, ARRAYSIZE(acLine) - 1 );
	CMeta::Tokenize( acLine, vecstrLine );
	if( vecstrLine.size( ) != 3 ) {
		g_CatSleipnir.error( "SMotifMatch::Open( ) invalid line: %s", acLine );
		return false; }
	m_strType = vecstrLine[ 0 ];
	if( ( m_eSubsequence = CCoalesceSequencerBase::GetSubsequence( vecstrLine[ 1 ] ) ) ==
		CCoalesceSequencerBase::ESubsequenceEnd ) {
		g_CatSleipnir.error( "SMotifMatch::Open( ) invalid subsequence: %s", vecstrLine[ 1 ].c_str( ) );
		return false; }
	m_dZ = (float)atof( vecstrLine[ 2 ].c_str( ) );
	istm.getline( acLine, ARRAYSIZE(acLine) - 1 );
	if( ( m_iMotif = Motifs.Open( acLine ) ) == -1 ) {
		g_CatSleipnir.error( "SMotifMatch::Open( ) invalid motif: %s", acLine );
		return false; }

	return true; }

uint32_t SMotifMatch::Open( const CHierarchy& Hier, const vector<SMotifMatch>& vecsMotifs,
	CCoalesceMotifLibrary& Motifs, size_t& iCount ) {
	uint32_t	iLeft, iRight;

	if( Hier.IsGene( ) ) {
		const SMotifMatch&	sMotif	= vecsMotifs[ Hier.GetID( ) ];

		m_eSubsequence = sMotif.m_eSubsequence;
		m_strType = sMotif.m_strType;
		m_dZ += sMotif.m_dZ;
		iCount++;
		return ( m_iMotif = Motifs.Merge( sMotif.m_iMotif ) ); }

	return ( ( ( ( iLeft = Open( Hier.Get( false ), vecsMotifs, Motifs, iCount ) ) == -1 ) ||
		( ( iRight = Open( Hier.Get( true ), vecsMotifs, Motifs, iCount ) ) == -1 ) ) ? -1 :
		( m_iMotif = Motifs.Merge( iLeft, iRight, FLT_MAX ) ) ); }

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
	for( iterGene = m_setiGenes.begin( ); iterGene != m_setiGenes.end( ); ++iterGene )
		if( !binary_search( m_veciPrevGenes.begin( ), m_veciPrevGenes.end( ), *iterGene ) ) {
			HistogramsCluster.Add( GeneScores, *iterGene, false );
			if( pHistogramsPot )
				pHistogramsPot->Add( GeneScores, *iterGene, true ); }
	for( i = 0; i < m_veciPrevGenes.size( ); ++i )
		if( m_setiGenes.find( m_veciPrevGenes[ i ] ) == m_setiGenes.end( ) ) {
			HistogramsCluster.Add( GeneScores, m_veciPrevGenes[ i ], true );
			if( pHistogramsPot )
				pHistogramsPot->Add( GeneScores, m_veciPrevGenes[ i ], false ); } }

void CCoalesceCluster::Subtract( CPCL& PCL ) const {
	set<size_t>::const_iterator	iterGene, iterDataset;
	size_t						i;
	float						d;

	for( iterGene = m_setiGenes.begin( ); iterGene != m_setiGenes.end( ); ++iterGene )
		for( iterDataset = m_setiDatasets.begin( ); iterDataset != m_setiDatasets.end( ); ++iterDataset )
			for( i = 0; i < GetConditions( *iterDataset ); ++i )
				if( !CMeta::IsNaN( d = m_vecdCentroid[ GetCondition( *iterDataset, i ) ] ) )
					PCL.Get( *iterGene, GetCondition( *iterDataset, i ) ) -= d; }

void CCoalesceCluster::Subtract( CCoalesceGeneScores& GeneScores ) const {
	set<size_t>::const_iterator			iterGene;
	set<SMotifMatch>::const_iterator	iterMotif;

	for( iterGene = m_setiGenes.begin( ); iterGene != m_setiGenes.end( ); ++iterGene )
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
				dZ = CStatistics::CohensD( adCluster, adCluster + iCluster, adPot, adPot + iPot );
				dP = CStatistics::ZTest( dZ, iCluster ) * psData->m_pvecsDatasets->size( ); }
			if( dP < psData->m_dPValue ) {
				g_CatSleipnir.info( "CCoalesceClusterImpl::ThreadSelectCondition( %g ) selected condition %d at %g, z=%g",
					psData->m_dPValue, iCondition, dP, dZ );
				(*psData->m_pvecfSignificant)[ iDataset ] = true;
				sDataset.m_dZ = (float)dZ; } }
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
				sDataset.m_psDataset->m_MatSigmaChol );
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
				sDataset.m_dZ = (float)dZ; } } }
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
	veciCluster.resize( m_setiGenes.size( ) );
	copy( m_setiGenes.begin( ), m_setiGenes.end( ), veciCluster.begin( ) );
	veciPot.resize( Pot.m_setiGenes.size( ) );
	copy( Pot.m_setiGenes.begin( ), Pot.m_setiGenes.end( ), veciPot.begin( ) );
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
				vecsMotifs.push_back( sMotif ); }
			else if( g_CatSleipnir.isDebugEnabled( ) ) {
				g_CatSleipnir.debug( "CCoalesceClusterImpl::AddSignificant( %d, %g ) failed at %g:\n%s",
					iMotif, dPValue, dP, SMotifMatch( iMotif, strTypeCluster, eSubsequence, (float)dZ,
					(float)dAverage ).Save( &Motifs ).c_str( ) );
				g_CatSleipnir.debug( "Cluster	%s", HistSetCluster.Save( iMotif ).c_str( ) );
				g_CatSleipnir.debug( "Pot	%s", HistSetPot.Save( iMotif ).c_str( ) ); } } }

	return true; }

void* CCoalesceClusterImpl::ThreadCentroid( void* pData ) {
	SThreadCentroid*	psData;

	psData = (SThreadCentroid*)pData;
	psData->m_pCluster->CalculateCentroid( psData->m_PCL );
	psData->m_pCluster->m_setiGenes.clear( );

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
	for( i = 0; i < vecfSignificant.size( ); ++i )
		if( vecfSignificant[ i ] ) {
			Add( i );
			m_vecdPriors[ i ] = 1; }
		else {
			Pot.Add( i );
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
	double				dPCluster, dPPot, dStdCluster, dStdPot;
	size_t				iDataset, iPot, iCondition;
	long double			dPIn, dPOut;

	iPot = PCL.GetGenes( ) - m_veciPrevGenes.size( );
	dLogPIn = dLogPOut = 0;
	dPIn = dPOut = 1;
	for( iDataset = 0; iDataset < veciDatasets.size( ); ++iDataset ) {
		const SDataset&	sDataset	= m_vecsDatasets[ veciDatasets[ iDataset ] ];

		if( sDataset.GetConditions( ) == 1 ) {
			iCondition = sDataset.GetCondition( 0 );
			if( CMeta::IsNaN( dGene = PCL.Get( iGene, iCondition ) ) ||
				CMeta::IsNaN( dCluster = m_vecdCentroid[ iCondition ] ) ||
				CMeta::IsNaN( dPot = Pot.m_vecdCentroid[ iCondition ] ) ||
				( ( dStdCluster = m_vecdStdevs[ iCondition ] ) < c_dEpsilonZero ) ||
				( ( dStdPot = Pot.m_vecdStdevs[ iCondition ] ) < c_dEpsilonZero ) )
				continue;
			
			if( fClustered )
				dCluster = ( ( dCluster * m_veciPrevGenes.size( ) ) - dGene ) /
					( m_veciPrevGenes.size( ) - 1 );
			else
				dPot = ( ( dPot * iPot ) - dGene ) / ( iPot - 1 );
			dPCluster = max( c_dEpsilonZero, CStatistics::NormalPDF( dGene, dCluster, dStdCluster ) );
			dPPot = max( c_dEpsilonZero, CStatistics::NormalPDF( dGene, dPot, dStdPot ) );
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

	GetConditions( setiConditions, PCL.GetExperiments( ) );
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
		if( IsGene( iGeneFrom ) && !SaveCopy( PCL, iGeneFrom, PCLCopy, iGeneTo++, false ) ) // true ) )
			return false;
	for( iGeneFrom = 0; iGeneFrom < PCL.GetGenes( ); ++iGeneFrom )
		if( !IsGene( iGeneFrom ) )
			PCLCopy.MaskGene( iGeneTo++ );
//		if( !IsGene( iGeneFrom ) && !SaveCopy( PCL, iGeneFrom, PCLCopy, iGeneTo++, false ) )
//			return false;

	ofsm.clear( );
	ofsm.open( ( strBase + CPCL::GetExtension( ) ).c_str( ) );
	if( !ofsm.is_open( ) )
		return false;
	PCLCopy.Save( ofsm );
	ofsm.close( );

	return true; }

bool CCoalesceClusterImpl::SaveCopy( const CPCL& PCLFrom, size_t iGeneFrom, CPCL& PCLTo, size_t iGeneTo,
	bool fClustered ) const {
	size_t		i, iExpTo, iExpFrom;
	set<size_t>	setiConditions;

	GetConditions( setiConditions, PCLFrom.GetExperiments( ) );
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
	const CCoalesceMotifLibrary* pMotifs ) const {
	set<size_t>::const_iterator			iterID;
	set<SMotifMatch>::const_iterator	iterMotif;
	size_t								i;

	ostm << "Cluster\t" << iID << endl;
	ostm << "Genes";
	for( iterID = GetGenes( ).begin( ); iterID != GetGenes( ).end( ); ++iterID )
		ostm << '\t' << PCL.GetGene( *iterID );
	ostm << endl << "Conditions";
	for( iterID = m_setiDatasets.begin( ); iterID != m_setiDatasets.end( ); ++iterID ) 
		for( i = 0; i < GetConditions( *iterID, PCL.GetExperiments( ) ); ++i )
			ostm << '\t' << PCL.GetExperiment( GetCondition( *iterID, i ) );
	ostm << endl << "Motifs" << endl;
	for( iterMotif = GetMotifs( ).begin( ); iterMotif != GetMotifs( ).end( ); ++iterMotif )
		ostm << iterMotif->Save( pMotifs ) << endl; }

size_t CCoalesceCluster::Open( const string& strPCL, size_t iSkip, const CPCL& PCL,
	CCoalesceMotifLibrary* pMotifs ) {
	CPCL		PCLCluster;
	size_t		i, iGene;
	string		strMotifs;
	ifstream	ifsm;

	Clear( );
	if( !PCLCluster.Open( strPCL.c_str( ), iSkip ) )
		return -1;

	for( i = 0; i < PCLCluster.GetExperiments( ); ++i ) {
		if( PCLCluster.GetExperiment( i )[ 0 ] != c_cStar )
			break;
		m_setiDatasets.insert( i ); }
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
	const vector<string>& vecstrClusters, float dFraction, float dCutoff, CCoalesceMotifLibrary* pMotifs ) {
	map<size_t, size_t>						mapiiGenes, mapiiDatasets;
	size_t									i, j, k, iClusters;
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
		if( ( (float)iterItem->second / iClusters ) > dFraction )
			m_setiGenes.insert( iterItem->first );
	for( iterItem = mapiiDatasets.begin( ); iterItem != mapiiDatasets.end( ); ++iterItem )
		if( ( (float)iterItem->second / iClusters ) > dFraction )
			m_setiDatasets.insert( iterItem->first );
	if( !pMotifs )
		return true;

	for( i = 0; i < vecmapstrsetsMotifs.size( ); ++i ) {
		const map<string, set<SMotifMatch> >&			mapstrsetsMotifs	= vecmapstrsetsMotifs[ i ];
		map<string, set<SMotifMatch> >::const_iterator	iterMotifs;

		for( iterMotifs = mapstrsetsMotifs.begin( ); iterMotifs != mapstrsetsMotifs.end( ); ++iterMotifs ) {
			const set<SMotifMatch>&	setsMotifs	= iterMotifs->second;
			vector<SMotifMatch>		vecsMotifs;
			CDistanceMatrix			MatSimilarity;
			CHierarchy*				pHierMotifs;

			vecsMotifs.resize( setsMotifs.size( ) );
			copy( setsMotifs.begin( ), setsMotifs.end( ), vecsMotifs.begin( ) );
			MatSimilarity.Initialize( vecsMotifs.size( ) );
			for( j = 0; j < vecsMotifs.size( ); ++j )
				for( k = ( j + 1 ); k < vecsMotifs.size( ); ++k )
					MatSimilarity.Set( j, k, -pMotifs->Align( vecsMotifs[ j ].m_iMotif,
						vecsMotifs[ k ].m_iMotif, dCutoff ) );
			if( !( ( pHierMotifs = CClustHierarchical::Cluster( MatSimilarity ) ) &&
				CCoalesceClusterImpl::Open( *pMotifs, *pHierMotifs, vecsMotifs, dCutoff, m_setsMotifs ) ) )
				return false;
			pHierMotifs->Destroy( ); } }

	return true; }

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
	for( iterFrom = Cluster.m_setiGenes.begin( ); iterFrom != Cluster.m_setiGenes.end( ); ++iterFrom )
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

bool CCoalesceClusterImpl::Open( CCoalesceMotifLibrary& Motifs, const CHierarchy& Hier,
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

	return ( Open( Motifs, Hier.Get( false ), vecsMotifs, dCutoff, setsMotifs ) &&
		Open( Motifs, Hier.Get( true ), vecsMotifs, dCutoff, setsMotifs ) ); }

float CCoalesceCluster::GetSimilarity( const CCoalesceCluster& Cluster, size_t iGenes,
	size_t iDatasets ) const {
	size_t						iOverlapGenes, iOverlapDatasets;
	set<size_t>::const_iterator	iterItem;

	for( iOverlapGenes = 0,iterItem = m_setiGenes.begin( ); iterItem != m_setiGenes.end( ); ++iterItem )
		if( Cluster.IsGene( *iterItem ) )
			iOverlapGenes++;
	for( iOverlapDatasets = 0,iterItem = m_setiDatasets.begin( ); iterItem != m_setiDatasets.end( ); ++iterItem )
		if( Cluster.IsDataset( *iterItem ) )
			iOverlapDatasets++;
	return (float)( ( 2 - CStatistics::HypergeometricCDF( iOverlapGenes, m_setiGenes.size( ),
		Cluster.m_setiGenes.size( ), iGenes ) - CStatistics::HypergeometricCDF( iOverlapDatasets,
		m_setiDatasets.size( ), Cluster.m_setiDatasets.size( ), iDatasets ) ) / 2 ); }

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
	if( g_CatSleipnir.isInfoEnabled( ) ) {
		size_t			iSequences;
		ostringstream	ossm;

		for( iSequences = i = 0; i < veciPCL2FASTA.size( ); ++i )
			if( veciPCL2FASTA[ i ] != -1 )
				iSequences++;
		g_CatSleipnir.info( "CCoalesce::Cluster( ) running with %d genes, %d conditions, and %d sequences",
			PCL.GetGenes( ), PCL.GetExperiments( ), iSequences );
		for( i = 0; i < PCL.GetExperiments( ); ++i )
			ossm << ( i ? "\t" : "" ) << PCL.GetExperiment( i );
		g_CatSleipnir.info( ossm.str( ) );
		g_CatSleipnir.info( "k %d, P gene %g, p condition %g, p motif %g, p correlation %g", GetK( ),
			GetProbabilityGene( ), GetPValueCondition( ), GetPValueMotif( ), GetPValueCorrelation( ) );
		g_CatSleipnir.info( "p merge %g, cutoff merge %g, penalty gap %g, penalty mismatch %g",
			GetPValueMerge( ), GetCutoffMerge( ), GetMotifs( ) ? GetMotifs( )->GetPenaltyGap( ) : 0,
			GetMotifs( ) ? GetMotifs( )->GetPenaltyMismatch( ) : 0 );
		g_CatSleipnir.info( "correlation pairs %d, bases %d, min size %d, max size %d",
			GetNumberCorrelation( ), GetBasesPerMatch( ), GetSizeMinimum( ), GetSizeMaximum( ) ); }
	for( dFailure = 1; dFailure > c_dEpsilon; dFailure *= GetPValueCorrelation( ) ) {
		CCoalesceCluster			Cluster, Pot;
		CCoalesceGroupHistograms	HistsCluster( GetMotifCount( ), GetBins( ), 1.0f / GetBasesPerMatch( ) );
		CCoalesceGroupHistograms	HistsPot( GetMotifCount( ), GetBins( ), 1.0f / GetBasesPerMatch( ) );

		for( i = 0; i < PCLCopy.GetGenes( ); ++i )
			Pot.Add( i );
		Pot.CalculateHistograms( GeneScores, HistsPot, NULL );
		if( !Cluster.Initialize( PCLCopy, Pot, m_vecsDatasets, setpriiSeeds, GetNumberCorrelation( ),
			GetPValueCorrelation( ), GetThreads( ) ) )
			continue;
		g_CatSleipnir.notice( "CCoalesce::Cluster( ) initialized %d genes", Cluster.GetGenes( ).size( ) );
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

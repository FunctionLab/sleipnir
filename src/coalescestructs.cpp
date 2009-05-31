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
#include "coalescestructsi.h"
#include "coalescemotifs.h"
#include "fasta.h"
#include "pcl.h"
#include "statistics.h"
#include "clusthierarchical.h"

namespace Sleipnir {

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
	vector<size_t>	veciIndices, veciCounts;
	bool			fEven;

	if( GetConditions( ) == 1 )
		return false;

	m_vecdStdevs.resize( GetConditions( ) );
	MatSigma.Initialize( GetConditions( ), GetConditions( ) );
	MatSigma.Clear( );
	vecdAves.resize( GetConditions( ) );
	veciCounts.resize( GetConditions( ) );
	for( i = 0; i < PCL.GetGenes( ); ++i )
		for( j = 0; j < vecdAves.size( ); ++j )
			if( !CMeta::IsNaN( dOne = PCL.Get( i, GetCondition( j ) ) ) ) {
				veciCounts[ j ]++;
				vecdAves[ j ] += dOne; }
	for( i = 0; i < vecdAves.size( ); ++i )
		if( j = veciCounts[ i ] )
			vecdAves[ i ] /= j;
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

// SMotifMatch

bool SMotifMatch::Open( std::istream& istm, CCoalesceMotifLibrary& Motifs ) {
	string			strLine;
	vector<string>	vecstrLine;
	size_t			i;

	m_vecprstrdKnown.clear( );
	strLine.resize( CFile::GetBufferSize( ) );
	istm.getline( &strLine[ 0 ], strLine.size( ) - 1 );
	CMeta::Tokenize( strLine.c_str( ), vecstrLine );
	if( vecstrLine.size( ) < 3 ) {
		g_CatSleipnir.error( "SMotifMatch::Open( ) invalid line: %s", strLine.c_str( ) );
		return false; }
	m_strType = vecstrLine[ 0 ];
	if( ( m_eSubsequence = CCoalesceSequencerBase::GetSubsequence( vecstrLine[ 1 ] ) ) ==
		CCoalesceSequencerBase::ESubsequenceEnd ) {
		g_CatSleipnir.error( "SMotifMatch::Open( ) invalid subsequence: %s", vecstrLine[ 1 ].c_str( ) );
		return false; }
	m_dZ = (float)atof( vecstrLine[ 2 ].c_str( ) );

	if( vecstrLine.size( ) > 3 ) {
		vector<string>	vecstrKnowns;

		CMeta::Tokenize( vecstrLine[ 3 ].c_str( ), vecstrKnowns, "|" );
		for( i = 0; i < vecstrKnowns.size( ); ++i ) {
			vector<string>	vecstrKnown;

			CMeta::Tokenize( vecstrKnowns[ i ].c_str( ), vecstrKnown, ":" );
			if( vecstrKnown.size( ) != 2 ) {
				g_CatSleipnir.error( "SMotifMatch::Open( ) invalid known: %s", vecstrKnowns[ i ].c_str( ) );
				return false; }
			m_vecprstrdKnown.push_back( pair<string, float>( vecstrKnown[ 0 ],
				(float)atof( vecstrKnown[ 1 ].c_str( ) ) ) ); } }

	istm.getline( &strLine[ 0 ], strLine.size( ) - 1 );
	if( ( m_iMotif = Motifs.Open( strLine.c_str( ) ) ) == -1 ) {
		g_CatSleipnir.error( "SMotifMatch::Open( ) invalid motif: %s", strLine.c_str( ) );
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
		return ( m_iMotif = sMotif.m_iMotif ); }

	return ( ( ( ( iLeft = Open( Hier.Get( false ), vecsMotifs, Motifs, iCount ) ) == -1 ) ||
		( ( iRight = Open( Hier.Get( true ), vecsMotifs, Motifs, iCount ) ) == -1 ) ) ? -1 :
		( m_iMotif = Motifs.Merge( iLeft, iRight, FLT_MAX, true ) ) ); }

uint32_t SMotifMatch::Open( const SMotifMatch& sOne, const SMotifMatch& sTwo,
	CCoalesceMotifLibrary& Motifs ) {

	m_eSubsequence = sOne.m_eSubsequence;
	m_strType = sOne.m_strType;
	m_dZ = ( sOne.m_dZ + sTwo.m_dZ ) / 2;
	return ( m_iMotif = Motifs.Merge( sOne.m_iMotif, sTwo.m_iMotif, FLT_MAX, true ) ); }

string SMotifMatch::Save( const CCoalesceMotifLibrary* pMotifs, bool fPWM, float dCutoffPWMs,
	float dPenaltyGap, float dPenaltyMismatch, bool fNoRCs ) const {
	ostringstream	ossm;
	string			strPWM;
	size_t			i;

	ossm << m_strType << '\t' << CCoalesceSequencerBase::GetSubsequence( m_eSubsequence ) << '\t' << m_dZ;
	for( i = 0; i < m_vecprstrdKnown.size( ); ++i ) {
		const pair<string, float>&	prstrdKnown	= m_vecprstrdKnown[ i ];

		ossm << ( i ? "|" : "\t" ) << prstrdKnown.first << ':' << prstrdKnown.second; }
	ossm << endl;
	if( pMotifs ) {
		ossm << pMotifs->GetMotif( m_iMotif );
		if( fPWM ) {
			if( ( strPWM = pMotifs->GetPWM( m_iMotif, dCutoffPWMs, dPenaltyGap, dPenaltyMismatch,
				fNoRCs ) ).empty( ) )
				return "";
			ossm << endl << strPWM; } }
	else
		ossm << m_iMotif;

	return ossm.str( ); }

bool SMotifMatch::Label( const CCoalesceMotifLibrary& Motifs, EType eMatchType, float dPenaltyGap,
	float dPenaltyMismatch, float dPValue ) {

	m_vecprstrdKnown.clear( );
	return Motifs.GetKnown( m_iMotif, eMatchType, dPenaltyGap, dPenaltyMismatch, m_vecprstrdKnown, dPValue ); }

}

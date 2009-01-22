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
#include "coalescemotifs.h"
#include "coalescestructsi.h"
#include "pst.h"

namespace Sleipnir {

const char	CCoalesceMotifLibraryImpl::c_szBases[]			= "ACGT";
const char	CCoalesceMotifLibraryImpl::c_szComplements[]	= "TGCA";

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
		dRet += pPST->GetMatch( strSequence, pPST->GetDepth( ) - i ) * sModifiers.GetWeight( i ) / i; }
	sModifiers.AddWeight( 0, iOffset, i - 1 );
	for( i = 0; i < strSequence.length( ); ++i ) {
		dRet += pPST->GetMatch( strSequence.substr( i ) ) * sModifiers.GetWeight( pPST->GetDepth( ) ) /
			pPST->GetDepth( );
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
	m_vecpPSTs.push_back( pRet = new CPST( strlen( c_szBases ) ) );
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
	float dCutoff, bool fAllowDuplicates ) {
	int			iOffset;
	float		dScore;
	uint32_t	iRet;
	CPST*		pPST;

	if( ( ( dScore = Align( strOne, strTwo, dCutoff, iOffset ) ) > dCutoff ) ||
		!( fAllowDuplicates || dScore ) )
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

uint32_t CCoalesceMotifLibraryImpl::MergeKMerRC( uint32_t iKMer, uint32_t iRC, float dCutoff,
	bool fAllowDuplicates ) {
	string		strKMer, strOne, strTwo;
	float		dOne, dTwo, dMin;
	int			iOne, iTwo;
	uint32_t	iRet;
	CPST*		pPST;

	if( m_veciKMer2RC[ iKMer ] == iRC )
		return -1;
	strKMer = GetMotif( iKMer );
	strOne = GetRCOne( iRC );
	strTwo = GetReverseComplement( strOne );
	dOne = Align( strKMer, strOne, dCutoff, iOne );
	dTwo = Align( strKMer, strTwo, dCutoff, iTwo );
	dMin = min( dOne, dTwo );
	if( ( dMin > dCutoff ) || !( fAllowDuplicates || dMin ) )
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

uint32_t CCoalesceMotifLibraryImpl::MergeRCs( uint32_t iOne, uint32_t iTwo, float dCutoff,
	bool fAllowDuplicates ) {
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
	if( ( dMin > dCutoff ) || !( fAllowDuplicates || dMin ) )
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
	float dCutoff, bool fAllowDuplicates ) {
	int			iOffset;
	float		dScore;
	uint32_t	iRet;
	CPST*		pPSTOut;

	if( ( ( dScore = PSTIn.Align( strKMer, m_dPenaltyGap, m_dPenaltyMismatch, dCutoff, iOffset ) ) >
		dCutoff ) || !( fAllowDuplicates || dScore ) )
		return -1;

	pPSTOut = CreatePST( iRet );
	pPSTOut->Add( strKMer, PSTIn, iOffset );
	if( g_CatSleipnir.isInfoEnabled( ) ) {
		ostringstream	ossm;

		ossm << "CCoalesceMotifLibraryImpl::MergeKMerPST( " << strKMer << ", " << PSTIn.GetMotif( ) <<
			", " << dCutoff << " ) merged at " << dScore << " to " << pPSTOut->GetMotif( );
		g_CatSleipnir.info( ossm.str( ) ); }
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

uint32_t CCoalesceMotifLibraryImpl::MergeRCPST( uint32_t iRC, const CPST& PSTIn, float dCutoff,
	bool fAllowDuplicates ) {
	int			iOne, iTwo;
	uint32_t	iRet;
	CPST*		pPSTOut;
	string		strOne, strTwo;
	float		dOne, dTwo, dMin;

	strOne = GetRCOne( iRC );
	strTwo = GetReverseComplement( strOne );
	dOne = PSTIn.Align( strOne, m_dPenaltyGap, m_dPenaltyMismatch, dCutoff, iOne );
	dTwo = PSTIn.Align( strTwo, m_dPenaltyGap, m_dPenaltyMismatch, dCutoff, iTwo );
	dMin = min( dOne, dTwo );
	if( ( dMin > dCutoff ) || !( fAllowDuplicates || dMin ) )
		return -1;

	pPSTOut = CreatePST( iRet );
	if( dOne < dTwo ) {
		pPSTOut->Add( strOne, PSTIn, iOne );
		pPSTOut->Add( strTwo ); }
	else {
		pPSTOut->Add( strTwo, PSTIn, iTwo );
		pPSTOut->Add( strOne ); }
	if( g_CatSleipnir.isInfoEnabled( ) ) {
		ostringstream	ossm;

		ossm << "CCoalesceMotifLibraryImpl::MergeRCPST( " << GetMotif( iRC ) << ", " << PSTIn.GetMotif( ) <<
			", " << dCutoff << " ) merged at " << min( dOne, dTwo ) << " to " << pPSTOut->GetMotif( );
		g_CatSleipnir.info( ossm.str( ) ); }
	return iRet; }

float CCoalesceMotifLibraryImpl::AlignPSTs( const CPST& PSTOne, const CPST& PSTTwo, float dCutoff ) const {
	int	iOffset;

	return PSTOne.Align( PSTTwo, m_dPenaltyGap, m_dPenaltyMismatch, dCutoff, iOffset ); }

uint32_t CCoalesceMotifLibraryImpl::MergePSTs( const CPST& PSTOne, const CPST& PSTTwo, float dCutoff,
	bool fAllowDuplicates ) {
	int			iOffset;
	uint32_t	iRet;
	CPST*		pPSTOut;
	float		dScore;

	if( ( ( dScore = PSTOne.Align( PSTTwo, m_dPenaltyGap, m_dPenaltyMismatch, dCutoff, iOffset ) ) >
		dCutoff ) || !( fAllowDuplicates || dScore ) )
		return -1;

	pPSTOut = CreatePST( iRet );
	if( iOffset < 0 ) {
		pPSTOut->Add( PSTOne );
		pPSTOut->Add( PSTTwo, -iOffset ); }
	else {
		pPSTOut->Add( PSTTwo );
		pPSTOut->Add( PSTOne, iOffset ); }
	if( g_CatSleipnir.isInfoEnabled( ) ) {
		ostringstream	ossm;

		ossm << "CCoalesceMotifLibraryImpl::MergePSTs( " << PSTOne.GetMotif( ) << ", " <<
			PSTTwo.GetMotif( ) << ", " << dCutoff << " ) merged at " << dScore << " to " <<
			pPSTOut->GetMotif( );
		g_CatSleipnir.info( ossm.str( ) ); }

	return iRet; }

uint32_t CCoalesceMotifLibraryImpl::RemoveRCs( const CPST& PST, float dPenaltyGap, float dPenaltyMismatch ) {
	CPST*		pPST;
	uint32_t	iRet;

	if( !( pPST = CreatePST( iRet ) ) )
		return -1;
	PST.RemoveRCs( dPenaltyGap, dPenaltyMismatch, *pPST );
	return iRet; }

string CCoalesceMotifLibrary::GetPWM( uint32_t iMotif, float dCutoffPWMs, float dPenaltyGap,
	float dPenaltyMismatch, bool fNoRCs ) const {
	CFullMatrix<uint16_t>	MatPWM;
	string					strMotif;
	size_t					i, j;
	ostringstream			ossm;
	float					d;

	if( fNoRCs )
		iMotif = ((CCoalesceMotifLibrary*)this)->RemoveRCs( iMotif, dPenaltyGap, dPenaltyMismatch );
	switch( GetType( iMotif ) ) {
		case ETypeKMer:
			if( !CCoalesceMotifLibraryImpl::GetPWM( GetMotif( iMotif ), MatPWM ) )
				return "";
			break;

		case ETypeRC:
			if( !( CCoalesceMotifLibraryImpl::GetPWM( strMotif = GetRCOne( iMotif ), MatPWM ) &&
				CCoalesceMotifLibraryImpl::GetPWM( GetReverseComplement( strMotif ), MatPWM ) ) )
				return "";
			break;

		case ETypePST:
			GetPST( iMotif )->GetPWM( MatPWM, c_szBases );
			break; }

	if( dCutoffPWMs ) {
		d = GetInformation( MatPWM );
		if( d < dCutoffPWMs ) {
			if( g_CatSleipnir.isInfoEnabled( ) ) {
				ostringstream	ossm;

				ossm << "CCoalesceMotifLibrary::GetPWM( " << iMotif << ", " << dCutoffPWMs << ", " <<
					fNoRCs << " ) rejected (" << d << "):" << endl;
				MatPWM.Save( ossm, false );
				g_CatSleipnir.info( ossm.str( ) ); }
			return ""; }
		if( g_CatSleipnir.isDebugEnabled( ) ) {
			ostringstream	ossm;

			ossm << "CCoalesceMotifLibrary::GetPWM( " << iMotif << ", " << dCutoffPWMs << ", " <<
				fNoRCs << " ) got information (" << d << "):" << endl;
			MatPWM.Save( ossm, false );
			g_CatSleipnir.debug( ossm.str( ) ); } }
	for( i = 0; i < MatPWM.GetRows( ); ++i ) {
		for( j = 0; j < MatPWM.GetColumns( ); ++j )
			ossm << ( j ? "\t" : "" ) << MatPWM.Get( i, j );
		ossm << endl; }

	return ossm.str( ); }

bool CCoalesceMotifLibrary::Simplify( uint32_t iMotif ) const {

	return ( ( GetType( iMotif ) == ETypePST ) ? CCoalesceMotifLibraryImpl::GetPST( iMotif )->Simplify( ) :
		false ); }

}

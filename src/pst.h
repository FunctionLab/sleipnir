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
#ifndef PST_H
#define PST_H

#include "psti.h"

namespace Sleipnir {

class CPST : protected CPSTImpl {
public:
	CPST( size_t iArity ) : CPSTImpl(iArity) { }

	void Add( const std::string& strSequence, const CPST& PST, int iOffset ) {

		if( iOffset < 0 ) {
			Add( strSequence );
			Add( PST, -iOffset ); }
		else {
			Add( PST );
			Add( strSequence, iOffset ); } }

	void Add( const std::string& strOne, const std::string& strTwo, int iOffset ) {
		size_t	i, j;

		if( iOffset < 0 )
			return Add( strTwo, strOne, -iOffset );
		i = CPSTImpl::Add( strOne, 0, m_sRoot );
		j = CPSTImpl::Add( strTwo, iOffset, m_sRoot );
		if( ( i = max( i, j ) ) > m_iDepth )
			m_iDepth = i; }

	void Add( const std::string& strSequence, size_t iOffset = 0 ) {
		size_t	i;

		if( ( i = CPSTImpl::Add( strSequence, iOffset, m_sRoot ) ) > GetDepth( ) )
			m_iDepth = i; }

	void Add( const CPST& PST, size_t iOffset = 0 ) {
		size_t	i;

		if( ( i = CPSTImpl::Add( PST.m_sRoot, iOffset, m_sRoot ) ) > GetDepth( ) )
			m_iDepth = i; }

	float GetMatch( const std::string& strTarget, size_t iOffset = 0 ) const {
		float	dPMatch;
		size_t	iMatched;

		iMatched = 0;
		return ( ( ( dPMatch = CPSTImpl::GetMatch( strTarget, m_sRoot, iOffset, iMatched ) ) && iMatched ) ?
			( dPMatch * pow( 1.0f / m_iArity, (float)( GetDepth( ) - iMatched ) ) ) : 0 ); }

	size_t GetDepth( ) const {

		return m_iDepth; }

	std::string GetMotif( ) const {

		return CPSTImpl::GetMotif( m_sRoot ); }

	float Align( const std::string& strSequence, float dPenaltyGap, float dPenaltyMismatch, float dCutoff,
		int& iOffset ) const {
		float	dRet;

		dRet = CPSTImpl::Align( strSequence, strSequence.length( ), dPenaltyGap, dPenaltyMismatch, dCutoff,
			iOffset );
		iOffset *= -1;
		return dRet; }

	float Align( const CPST& PST, float dPenaltyGap, float dPenaltyMismatch, float dCutoff,
		int& iOffset ) const {

		return CPSTImpl::Align( PST.m_sRoot, PST.GetDepth( ), dPenaltyGap, dPenaltyMismatch, dCutoff,
			iOffset ); }
};

}

#endif // PST_H

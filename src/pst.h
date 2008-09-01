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

	void Add( const std::string& strOne, const std::string& strTwo, size_t iOffset ) {
		size_t	i;

		if( ( i = CPSTImpl::Add( strOne, strTwo, iOffset, m_sRoot ) ) > GetDepth( ) )
			m_iDepth = i; }

	float GetMatch( const std::string& strTarget ) const {
		long double	dPMatch, dPMismatch;

		if( !( dPMatch = CPSTImpl::GetMatch( strTarget, m_sRoot ) ) )
			return 0;

		dPMismatch = pow( 1 - ( 1.0 / m_iArity ), (int)min( strTarget.length( ), GetDepth( ) ) );
		return (float)( dPMatch / ( dPMatch + dPMismatch ) ); }

	size_t GetDepth( ) const {

		return m_iDepth; }

	std::string GetMotif( ) const {

		return CPSTImpl::GetMotif( m_sRoot ); }
};

}

#endif // PST_H

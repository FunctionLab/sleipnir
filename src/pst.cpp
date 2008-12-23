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
#include "pst.h"
#include "coalescemotifs.h"

namespace Sleipnir {

void CPSTImpl::RemoveRCs( const map<unsigned char, unsigned char>& mapccComplements, const SNode& sNode,
	size_t iOffset, string& strSeq, vector<SRC>& vecsOut ) {
	size_t	i;

	strSeq.push_back( sNode.m_cCharacter );
	if( sNode.m_vecsChildren.size( ) == 0 )
		vecsOut.push_back( ( CCoalesceMotifLibrary::GetPurines( strSeq ) < 0.5 ) ?
			SRC( CCoalesceMotifLibrary::GetReverseComplement( strSeq ), iOffset ) : SRC( strSeq, 0 ) );
	else {
		iOffset--;
		for( i = 0; i < sNode.m_vecsChildren.size( ); ++i )
			RemoveRCs( mapccComplements, sNode.m_vecsChildren[ i ], iOffset, strSeq, vecsOut ); }
	strSeq.resize( strSeq.size( ) - 1 ); }

}

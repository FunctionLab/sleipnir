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

struct SSortRCs {
	bool operator()( const string& strOne, const string& strTwo ) const {

		return ( strOne.length( ) > strTwo.length( ) ); }
};

void CPST::RemoveRCs( const map<unsigned char, unsigned char>& mapccComplements, float dPenaltyGap,
	float dPenaltyMismatch, CPST& PSTOut ) const {
	size_t			i;
	string			str;
	vector<string>	vecstrAdd;
	float			dOne, dTwo;
	int				iOne, iTwo;

	for( i = 0; i < m_sRoot.m_vecsChildren.size( ); ++i )
		CPSTImpl::RemoveRCs( mapccComplements, m_sRoot.m_vecsChildren[ i ], str, vecstrAdd );
	if( vecstrAdd.empty( ) )
		return;

	sort( vecstrAdd.begin( ), vecstrAdd.end( ), SSortRCs( ) );
	PSTOut.Add( vecstrAdd[ 0 ], 0 );
	for( i = 1; i < vecstrAdd.size( ); ++i ) {
		dOne = PSTOut.Align( vecstrAdd[ i ], dPenaltyGap, dPenaltyMismatch, FLT_MAX, iOne );
		dTwo = PSTOut.Align( str = CCoalesceMotifLibrary::GetReverseComplement( vecstrAdd[ i ] ),
			dPenaltyGap, dPenaltyMismatch, FLT_MAX, iTwo );
		if( dTwo < dOne )
			PSTOut.Add( str, iTwo );
		else
			PSTOut.Add( vecstrAdd[ i ], iOne ); } }

}

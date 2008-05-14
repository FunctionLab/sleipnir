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
#include "bayesnet.h"
#include "dataset.h"

namespace Sleipnir {

const char	CBayesNetImpl::c_szFR[]		= "FR";
const char	CBayesNetImpl::c_szZero[]	= "zero";

void CBayesNetImpl::EncodeData( const IDataset* pData, TMapData& mapData ) {
	size_t				i, j;
	string				strCur;
	TMapData::iterator	iterDatum;

	for( i = 0; i < pData->GetGenes( ); ++i )
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j ) {
			if( !pData->IsExample( i, j ) || ( pData->GetDiscrete( i, j, 0 ) == -1 ) )
				continue;
			strCur = EncodeDatum( pData, i, j );
			if( ( iterDatum = mapData.find( strCur ) ) == mapData.end( ) )
				mapData[ strCur ] = 1;
			else
				iterDatum->second++; } }

string CBayesNetImpl::EncodeDatum( const IDataset* pData, size_t iOne, size_t iTwo ) {
	string	strRet;
	size_t	i, iCur;

	for( i = 0; i < pData->GetExperiments( ); ++i )
		strRet += ( ( iCur = pData->GetDiscrete( iOne, iTwo, i ) ) == -1 ) ? c_cMissing :
			(char)( c_cBase + ( iCur & 0xFF ) );

	return strRet; }

string CBayesNetImpl::EncodeDatum( const CPCLPair& PCL, size_t iGene, const vector<size_t>& veciMap ) {
	string	strRet;
	size_t	i, iCur;

	for( i = 0; i < veciMap.size( ); ++i )
		strRet += ( ( veciMap[ i ] == -1 ) || ( ( iCur = PCL.Quantize( PCL.Get( iGene, veciMap[ i ] ),
			veciMap[ i ] ) ) == -1 ) ) ? c_cMissing : (char)( c_cBase + ( iCur & 0xFF ) );

	return strRet; }

void CBayesNetImpl::DecodeDatum( const std::string& strDatum, std::vector<size_t>& veciDatum ) {
	size_t	i;

	for( i = 0; i < strDatum.length( ); ++i )
		veciDatum[ i ] = ( strDatum[ i ] == c_cMissing ) ? -1 : ( strDatum[ i ] - c_cBase ); }

bool CBayesNetImpl::IsAnswer( const std::string& strDatum ) {

	return ( strDatum[ 0 ] != c_cMissing ); }

CBayesNetImpl::CBayesNetImpl( bool fGroup ) : m_fGroup(fGroup) { }

}

#include "stdafx.h"
#include "bayesnet.h"
#include "dataset.h"

namespace libBioUtils {

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

void CBayesNetImpl::DecodeDatum( const string& strDatum, vector<size_t>& veciDatum ) {
	size_t	i;

	for( i = 0; i < strDatum.length( ); ++i )
		veciDatum[ i ] = ( strDatum[ i ] == c_cMissing ) ? -1 : ( strDatum[ i ] - c_cBase ); }

bool CBayesNetImpl::IsAnswer( const string& strDatum ) {

	return ( strDatum[ 0 ] != c_cMissing ); }

CBayesNetImpl::CBayesNetImpl( bool fGroup ) : m_fGroup(fGroup) { }

}

#include "stdafx.h"
#include "bayesnet.h"
#include "dataset.h"

namespace libBioUtils {

string CBayesNetImpl::EncodeDatum( const IDataset* pData, size_t iOne, size_t iTwo ) {
	string	strRet;
	size_t	i, iCur;

	for( i = 0; i < pData->GetExperiments( ); ++i )
		strRet += ( ( iCur = pData->GetDiscrete( iOne, iTwo, i ) ) == -1 ) ? '_' :
			(char)( 'A' + ( iCur & 0xFF ) );

	return strRet; }

void CBayesNetImpl::DecodeDatum( const string& strDatum, vector<size_t>& veciDatum ) {
	size_t	i;

	for( i = 0; i < strDatum.length( ); ++i )
		veciDatum[ i ] = ( strDatum[ i ] == '_' ) ? -1 : ( strDatum[ i ] - 'A' ); }

CBayesNetImpl::CBayesNetImpl( bool fGroup ) : m_fGroup(fGroup) { }

}

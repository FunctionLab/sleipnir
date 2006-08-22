#include "stdafx.h"
#include "mathb.h"

namespace libBioUtils {

vector<double>	CMathImpl::s_vecdLogFact;

size_t CMath::GCD( size_t iA, size_t iB ) {
	size_t	i;

	while( iB ) {
		i = iA;
		iA = iB;
		iB = i % iB; }

	return iA; }

double CMath::Min( double dA, double dB ) {

	return ( ( dA < dB ) ? dA : dB ); }

size_t CMath::Round( double d ) {
	size_t	i;

	i = (size_t)d;
	return ( ( ( d - i ) >= 0.5 ) ? ( i + 1 ) : i ); }

double CMath::LogFact( size_t iN ) {
	size_t	i;

	if( s_vecdLogFact.empty( ) ) {
		s_vecdLogFact.resize( c_iFactCutoff );
		s_vecdLogFact[ 0 ] = s_vecdLogFact[ 1 ] = 0;
		for( i = 2; i < s_vecdLogFact.size( ); ++i )
			s_vecdLogFact[ i ] = s_vecdLogFact[ i - 1 ] + log( (double)i ); }

	return ( ( iN >= c_iFactCutoff ) ? LogFactStirling( iN ) : s_vecdLogFact[ iN ] ); }

double CMathImpl::LogFactStirling( size_t iN ) {

	if( iN < 3 )
		return ( ( iN == 2 ) ? log( 2.0 ) : 0 );

	return ( log( 2.5066282746271 ) + ( ( iN + 0.5 ) * log( (double)iN ) ) - iN ); }

double CMathImpl::LogFactRec( size_t iN ) {

	if( iN < 3 )
		return ( ( iN == 2 ) ? log( 2.0 ) : 0 );

	return ( log( (double)iN ) + LogFactRec( iN - 1 ) ); }

}

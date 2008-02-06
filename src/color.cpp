#include "stdafx.h"
#include "color.h"
#include "mathb.h"
#include "meta.h"

namespace libBioUtils {

const CColor	CColor::c_Black		= CColor( 0x00, 0x00, 0x00 );
const CColor	CColor::c_Cyan		= CColor( 0x00, 0xFF, 0xFF );
const CColor	CColor::c_Green		= CColor( 0x00, 0xFF, 0x00 );
const CColor	CColor::c_Red		= CColor( 0xFF, 0x00, 0x00 );
const CColor	CColor::c_White		= CColor( 0xFF, 0xFF, 0xFF );
const CColor	CColor::c_Yellow	= CColor( 0xFF, 0xFF, 0x00 );

CColor CColor::Interpolate( float dValue, const CColor& Min, const CColor& Med,
	const CColor& Max ) {
	float			dOther;
	const CColor*	pOther;

	dOther = 2 * dValue;
	if( dValue < 0.5 ) {
		pOther = &Min;
		dOther = 1 - dOther; }
	else {
		pOther = &Max;
		dOther -= 1; }

	return ( ( Med * ( 1 - dOther ) ) + ( *pOther * dOther ) ); }

CColor::CColor( const unsigned char* ac ) {
	size_t	i;

	for( i = 0; i < c_iChannels; ++i )
		m_acRGB[ i ] = ac[ i ]; }

CColor::CColor( unsigned char cR, unsigned char cG, unsigned char cB ) {

	m_acRGB[ 0 ] = cR;
	m_acRGB[ 1 ] = cG;
	m_acRGB[ 2 ] = cB; }

CColor CColor::operator+( const CColor& Color ) const {
	size_t			ai[ c_iChannels ];
	unsigned char	ac[ c_iChannels ];
	size_t			i;

	for( i = 0; i < c_iChannels; ++i ) {
		if( ( ai[ i ] = m_acRGB[ i ] + Color.m_acRGB[ i ] ) > UCHAR_MAX )
			ai[ i ] = UCHAR_MAX;
		ac[ i ] = ai[ i ]; }

	return CColor( ac ); }

CColor CColor::operator*( float dValue ) const {
	size_t			ai[ c_iChannels ];
	size_t			i;
	unsigned char	ac[ c_iChannels ];

	for( i = 0; i < c_iChannels; ++i ) {
		if( ( ai[ i ] = CMath::Round( m_acRGB[ i ] * dValue ) ) > UCHAR_MAX )
			ai[ i ] = UCHAR_MAX;
		ac[ i ] = ai[ i ]; }

	return CColor( ac ); }

CColor& CColor::operator=( const CColor& Color ) {
	size_t	i;

	for( i = 0; i < c_iChannels; ++i )
		m_acRGB[ i ] = Color.m_acRGB[ i ];

	return *this; }

string CColor::ToRGB( ) const {
	char	ac[ 7 ];

	sprintf_s( ac, ARRAYSIZE(ac), "%02x%02x%02x", m_acRGB[ 0 ], m_acRGB[ 1 ], m_acRGB[ 2 ] );

	return ac; }

void CColorImpl::ToHSV( float& dHue, float& dSat, float& dVal ) const {
	float	dMin, dMax, dDelta;
	float	adRGB[ c_iChannels ];
	size_t	i;

	for( i = 0; i < c_iChannels; ++i )
		adRGB[ i ] = (float)m_acRGB[ i ] / UCHAR_MAX;
	dMin = dMax = adRGB[ 0 ];
	for( i = 1; i < c_iChannels; ++i ) {
		if( adRGB[ i ] < dMin )
			dMin = adRGB[ i ];
		if( adRGB[ i ] > dMax )
			dMax = adRGB[ i ]; }

	if( !( dVal = dMax ) ) {
		dSat = 0;
		dHue = -1;
		return; }		
	dDelta = dMax - dMin;
	dSat = dDelta / dMax;

	if( adRGB[ 0 ] == dMax )
		dHue = ( adRGB[ 1 ] - adRGB[ 2 ] ) / dDelta;
	else if( adRGB[ 1 ] == dMax )
		dHue = 2 + ( ( adRGB[ 2 ] - adRGB[ 0 ] ) / dDelta );
	else
		dHue = 4 + ( ( adRGB[ 0 ] - adRGB[ 1 ] ) / dDelta );
	if( ( dHue *= 60 ) < 0 )
		dHue += 360;
	dHue /= 360; }

bool CColor::IsDark( ) const {
	float	dHue, dSat, dVal;

	ToHSV( dHue, dSat, dVal );
	
	return ( ( dVal < 0.5 ) || ( dHue < 0.05 ) || ( dHue > 0.666 ) ); }

}

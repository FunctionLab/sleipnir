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
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, and Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "stdafx.h"
#include "color.h"
#include "mathb.h"
#include "meta.h"

namespace Sleipnir {

const CColor	CColor::c_Black		= CColor( 0x00, 0x00, 0x00 );
const CColor	CColor::c_Cyan		= CColor( 0x00, 0xFF, 0xFF );
const CColor	CColor::c_Green		= CColor( 0x00, 0xFF, 0x00 );
const CColor	CColor::c_Red		= CColor( 0xFF, 0x00, 0x00 );
const CColor	CColor::c_White		= CColor( 0xFF, 0xFF, 0xFF );
const CColor	CColor::c_Yellow	= CColor( 0xFF, 0xFF, 0x00 );

/*!
 * \brief
 * Given a fractional value, interpolate between three colors.
 * 
 * \param dValue
 * Fractional value between 0 and 1.
 * 
 * \param ColorMinimum
 * Minimum color (corresponding to value zero).
 * 
 * \param ColorMedium
 * Midpoint color (corresponding to value 0.5).
 * 
 * \param ColorMaximum
 * Maximum color (corresponding to value one).
 * 
 * \returns
 * A new color blending the minimum/midpoint or midpoint/maximum.
 * 
 * Blend continuously between three colors given a value between 0 and 1.  Given 0, the minimum color is
 * returned.  Given 0.5, the midpoint color is returned.  Between 0 and 0.5, the minimum and midpoint
 * colors are mixed; between 0.5 and 1, the midpoint and maximum colors are mixed.  Given 1, the maximum
 * color is returned.
 */
CColor CColor::Interpolate( float dValue, const CColor& ColorMinimum, const CColor& ColorMedium,
	const CColor& ColorMaximum ) {
	float			dOther;
	const CColor*	pOther;

	if( dValue < 0 )
		dValue = 0;
	else if( dValue > 1 )
		dValue = 1;
	dOther = 2 * dValue;
	if( dValue < 0.5 ) {
		pOther = &ColorMinimum;
		dOther = 1 - dOther; }
	else {
		pOther = &ColorMaximum;
		dOther -= 1; }

	return ( ( ColorMedium * ( 1 - dOther ) ) + ( *pOther * dOther ) ); }

/*!
 * \brief
 * Construct a new color with the given red, green, and blue values.
 * 
 * \param abRGB
 * Array of three bytes representing red, green, and blue values between 0 and 255.
 */
CColor::CColor( const unsigned char* abRGB ) {
	size_t	i;

	for( i = 0; i < c_iChannels; ++i )
		m_abRGB[ i ] = abRGB[ i ]; }

/*!
 * \brief
 * Construct a new color with the given red, green, and blue values.
 * 
 * \param bRed
 * Red value between 0 and 255.
 * 
 * \param bGreen
 * Green value between 0 and 255.
 * 
 * \param bBlue
 * Blue value between 0 and 255.
 */
CColor::CColor( unsigned char bRed, unsigned char bGreen, unsigned char bBlue ) {

	m_abRGB[ 0 ] = bRed;
	m_abRGB[ 1 ] = bGreen;
	m_abRGB[ 2 ] = bBlue; }

/*!
 * \brief
 * Blend two colors equally.
 * 
 * \param Color
 * Color to which the current color is added.
 * 
 * \returns
 * New color blending the current and given colors.
 */
CColor CColor::operator+( const CColor& Color ) const {
	size_t			ai[ c_iChannels ];
	unsigned char	ac[ c_iChannels ];
	size_t			i;

	for( i = 0; i < c_iChannels; ++i ) {
		if( ( ai[ i ] = m_abRGB[ i ] + Color.m_abRGB[ i ] ) > UCHAR_MAX )
			ai[ i ] = UCHAR_MAX;
		ac[ i ] = ai[ i ]; }

	return CColor( ac ); }

/*!
 * \brief
 * Scale the current color by the requested amount.
 * 
 * \param dValue
 * Scalar by which the channels of the current color are multiplied.
 * 
 * \returns
 * New color scaling the current color by the requested amount.
 */
CColor CColor::operator*( float dValue ) const {
	size_t			ai[ c_iChannels ];
	size_t			i;
	unsigned char	ac[ c_iChannels ];

	for( i = 0; i < c_iChannels; ++i ) {
		if( ( ai[ i ] = CMath::Round( m_abRGB[ i ] * dValue ) ) > UCHAR_MAX )
			ai[ i ] = UCHAR_MAX;
		ac[ i ] = ai[ i ]; }

	return CColor( ac ); }

/*!
 * \brief
 * Copies the given color into the current color.
 * 
 * \param Color
 * Color whose channel values are copied into the current color.
 * 
 * \returns
 * Reference to the current color with values copied from the given color.
 */
CColor& CColor::operator=( const CColor& Color ) {
	size_t	i;

	for( i = 0; i < c_iChannels; ++i )
		m_abRGB[ i ] = Color.m_abRGB[ i ];

	return *this; }

/*!
 * \brief
 * Constructs a standard hexadecimal string encoding the current color.
 * 
 * \returns
 * Hexadecimal string encoding the current colors RGB values.
 * 
 * \remarks
 * String does not include an initial # mark.
 */
string CColor::ToRGB( ) const {
	char	ac[ 7 ];

	sprintf_s( ac, ARRAYSIZE(ac), "%02x%02x%02x", m_abRGB[ 0 ], m_abRGB[ 1 ], m_abRGB[ 2 ] );

	return ac; }

void CColorImpl::ToHSV( float& dHue, float& dSat, float& dVal ) const {
	float	dMin, dMax, dDelta;
	float	adRGB[ c_iChannels ];
	size_t	i;

	for( i = 0; i < c_iChannels; ++i )
		adRGB[ i ] = (float)m_abRGB[ i ] / UCHAR_MAX;
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

/*!
 * \brief
 * Returns true if the current color is dark enough to warrant a light foreground.
 * 
 * \returns
 * True if the current color is perceptibly dark.
 * 
 * If this method returns true, it is generally safe to assume that the current color is dark enough to
 * warrant a white foreground; if it returns false, a black foreground is appropriate.  This calculation
 * is performed heuristically by testing whether the color's value is low or its hue is extreme in HSV
 * space.
 */
bool CColor::IsDark( ) const {
	float	dHue, dSat, dVal;

	ToHSV( dHue, dSat, dVal );
	
	return ( ( dVal < 0.5 ) || ( dHue < 0.05 ) || ( dHue > 0.666 ) ); }

}

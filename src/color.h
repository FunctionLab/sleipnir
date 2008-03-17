#ifndef COLOR_H
#define COLOR_H

#include "colori.h"

namespace Sleipnir {

/*!
 * \brief
 * Simple representation of a color triple in RGB space.
 */
class CColor : CColorImpl {
public:
	/*!
	 * \brief
	 * Black
	 */
	static const CColor	c_Black;
	/*!
	 * \brief
	 * Cyan
	 */
	static const CColor	c_Cyan;
	/*!
	 * \brief
	 * Green
	 */
	static const CColor	c_Green;
	/*!
	 * \brief
	 * Red
	 */
	static const CColor	c_Red;
	/*!
	 * \brief
	 * White
	 */
	static const CColor	c_White;
	/*!
	 * \brief
	 * Yellow
	 */
	static const CColor	c_Yellow;

	static CColor Interpolate( float dValue, const CColor& ColorMinimum, const CColor& ColorMedium,
		const CColor& ColorMaximum );

	CColor( const unsigned char* abRGB );
	CColor( unsigned char bRed, unsigned char bGreen, unsigned char bBlue );

	CColor operator+( const CColor& Color ) const;
	CColor operator*( float dValue ) const;
	CColor& operator=( const CColor& Color );
	std::string ToRGB( ) const;
	bool IsDark( ) const;
};

}

#endif // COLOR_H

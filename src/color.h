#ifndef COLOR_H
#define COLOR_H

#include "colori.h"

namespace Sleipnir {

class CColor : CColorImpl {
public:
	static const CColor	c_Black;
	static const CColor	c_Cyan;
	static const CColor	c_Green;
	static const CColor	c_Red;
	static const CColor	c_White;
	static const CColor	c_Yellow;

	static CColor Interpolate( float, const CColor&, const CColor&, const CColor& );

	CColor( const unsigned char* );
	CColor( unsigned char, unsigned char, unsigned char );

	CColor operator+( const CColor& ) const;
	CColor operator*( float ) const;
	CColor& operator=( const CColor& );
	std::string ToRGB( ) const;
	bool IsDark( ) const;
};

}

#endif // COLOR_H

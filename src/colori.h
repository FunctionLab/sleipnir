#ifndef COLORI_H
#define COLORI_H

namespace libBioUtils {

class CColorImpl {
protected:
	static const size_t	c_iChannels	= 3;

	void ToHSV( float&, float&, float& ) const;

	unsigned char	m_acRGB[ c_iChannels ];
};

}

#endif // COLORI_H

#ifndef HALFMATRIXI_H
#define HALFMATRIXI_H

namespace Sleipnir {

class CHalfMatrixBase {
protected:
	CHalfMatrixBase( ) : m_iSize(0) { }

	static void HalfIndex( size_t& iX, size_t& iY ) {
		size_t	i;

		if( iX > iY ) {
			i = iX;
			iX = iY;
			iY = i - iY - 1; }
		else
			iY -= iX + 1; }

	size_t	m_iSize;
};

}

#endif // HALFMATRIXI_H

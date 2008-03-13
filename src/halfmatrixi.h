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

	void UnHalfIndex( size_t iIndex, size_t& iX, size_t& iY ) const {
		size_t	i, iRow;

		iRow = m_iSize - 1;
		for( i = 0; ( i + 1 ) < m_iSize; ++i,--iRow ) {
			if( iIndex < iRow )
				break;
			iIndex -= iRow; }

		iX = i;
		iY = iIndex; }

	size_t	m_iSize;
};

}

#endif // HALFMATRIXI_H

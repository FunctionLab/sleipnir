#ifndef HALFMATRIXI_H
#define HALFMATRIXI_H

#include "pstdint.h"

namespace libBioUtils {

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

class CCompactMatrixImpl : protected CHalfMatrixBase {
protected:
	CCompactMatrixImpl( );
	~CCompactMatrixImpl( );

	size_t* GetWord( size_t, size_t, unsigned char& ) const;
	size_t CountWords( ) const;

	bool			m_fMemory;
	uint32_t		m_iSize;
	unsigned char	m_cBits;
	size_t*			m_aiData;
};

}

#endif // HALFMATRIXI_H

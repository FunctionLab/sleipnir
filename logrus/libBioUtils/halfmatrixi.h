#ifndef HALFMATRIXI_H
#define HALFMATRIXI_H

#ifdef _MSC_VER
typedef UINT16	uint16_t;
typedef UINT32	uint32_t;
typedef UINT64	uint64_t;
#endif // _MSC_VER

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

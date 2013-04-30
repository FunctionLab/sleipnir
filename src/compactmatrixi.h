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
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#ifndef COMPACTMATRIXI_H
#define COMPACTMATRIXI_H

//may need to enable the following lines for Cygwin compilation
//#ifndef SIZE_MAX
//#define SIZE_MAX (4294967295U) //assumes 32bit GCC
//#endif

#include "halfmatrixi.h"

namespace Sleipnir {

class CCompactMatrixBase {
protected:
	CCompactMatrixBase( ) : m_cBits(0), m_aiData(NULL), m_fMemory(true) { }

	~CCompactMatrixBase( ) {

		if( m_fMemory && m_aiData )
			delete[] m_aiData; }

	void Initialize( unsigned char, bool );

	unsigned char Get( size_t iX, size_t iY ) const {
		size_t*			pi;
		unsigned char	cRet, cShift;

		if( !( m_cBits && m_aiData ) )
			return 0;
		pi = GetWord( iX, iY, cShift );
		// Bits starting at Shift, mask m_cBits long
		cRet = (unsigned char)( ( *pi >> cShift ) & ( SIZE_MAX >> ( ( 8 * sizeof(*m_aiData) ) - m_cBits ) ) );
		// If we overflow a boundary...
		if( ( cShift + m_cBits ) > ( 8 * sizeof(*m_aiData) ) )
		// Bits starting at 0, mask (m_cBits-Shift) long, shift left by bits we got
			cRet |= ( *( pi + 1 ) & ( SIZE_MAX >> ( ( 16 * sizeof(*m_aiData) ) - m_cBits -
				cShift ) ) ) << ( ( 8 * sizeof(*m_aiData) ) - cShift );

		return cRet;
	}

	void Set( size_t iX, size_t iY, unsigned char cValue ) {
		unsigned char	cShift;
		size_t			iMask;
		size_t*			pi;

		if( !( m_cBits && m_aiData ) )
			return;
		pi = GetWord( iX, iY, cShift );
		iMask = ( SIZE_MAX >> ( ( 8 * sizeof(*m_aiData) ) - m_cBits ) ) << cShift;
		*pi = ( *pi & ~iMask ) | ( ( (size_t)cValue << cShift ) & iMask );
		if( ( cShift + m_cBits ) > ( 8 * sizeof(*m_aiData) ) ) {
			pi++;
			iMask = SIZE_MAX >> ( ( 16 * sizeof(*m_aiData) ) - m_cBits - cShift );
			*pi = ( *pi & ~iMask ) |
				( ( cValue >> ( ( 8 * sizeof(*m_aiData) ) - cShift ) ) & iMask );
		}
	}

	virtual size_t CountWords( ) const = 0;
	virtual size_t* GetWord( size_t, size_t, unsigned char& ) const = 0;

	bool			m_fMemory;
	unsigned char	m_cBits;
	size_t*			m_aiData;
};

class CCompactMatrixImpl : protected CHalfMatrixBase, protected CCompactMatrixBase {
protected:
	CCompactMatrixImpl( ) : m_iSize(0) { }

	size_t* GetWord( size_t iX, size_t iY, unsigned char& cShift ) const {
		size_t	iIndex;

		// Closed form for sum(m_iSize - i - 1, i=0..(iX-1)) + iY
		iIndex = ( iX * ( m_iSize - 1 ) ) - ( iX * ( iX - 1 ) / 2 ) + iY;
		iIndex *= m_cBits;
		cShift = (unsigned char)( iIndex % ( 8 * sizeof(*m_aiData) ) );
		iIndex /= 8 * sizeof(*m_aiData);

		return &m_aiData[ iIndex ]; }

	size_t CountWords( ) const {
		size_t	iRet;

		return ( ( m_cBits && ( iRet = m_iSize * ( m_iSize - 1 ) / 2 ) ) ?
			( ( ( ( iRet * m_cBits ) - 1 ) / ( 8 * sizeof(*m_aiData) ) ) + 1 ) : 0 ); }

	uint32_t	m_iSize;
};

class CCompactFullMatrixImpl : protected CCompactMatrixBase {
protected:
	CCompactFullMatrixImpl( ) : m_iRows(0), m_iColumns(0) { }

	size_t CountWords( ) const {
		size_t	iRet;
		iRet = m_iRows * m_iColumns;
		return ( ( ( ( iRet * m_cBits ) - 1 ) / ( 8 * sizeof(*m_aiData) ) ) + 1 );
	}

	size_t* GetWord( size_t iY, size_t iX, unsigned char& cShift ) const {
		size_t	iIndex;

		iIndex = ( ( iY * m_iColumns ) + iX ) * m_cBits;
		cShift = (unsigned char)( iIndex % ( 8 * sizeof(*m_aiData) ) );
		iIndex /= 8 * sizeof(*m_aiData);

		return &m_aiData[ iIndex ]; }

	uint32_t	m_iRows;
	uint32_t	m_iColumns;
};

}

#endif // COMPACTMATRIXI_H

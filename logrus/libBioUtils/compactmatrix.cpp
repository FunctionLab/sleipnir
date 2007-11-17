#include "stdafx.h"
#include "compactmatrix.h"

namespace libBioUtils {

///////////////////////////////////////////////////////////////////////////////
// CCompactMatrixBase
///////////////////////////////////////////////////////////////////////////////

void CCompactMatrixBase::Initialize( unsigned char cValues, bool fClear ) {
	size_t	iWords;

	if( m_fMemory && m_aiData )
		delete[] m_aiData;
	m_fMemory = true;
	for( m_cBits = 0,--cValues; cValues; ++m_cBits,cValues >>= 1 );
	m_aiData = new size_t[ iWords = CountWords( ) ];
	if( fClear )
		memset( m_aiData, 0, iWords * sizeof(*m_aiData) ); }

///////////////////////////////////////////////////////////////////////////////
// CCompactMatrixImpl
///////////////////////////////////////////////////////////////////////////////

bool CCompactMatrix::Open( istream& istm ) {
	uint32_t	iBytes;
	size_t		iWords;

	if( m_fMemory && m_aiData )
		delete[] m_aiData;
	m_fMemory = true;
	istm.read( (char*)&m_iSize, sizeof(m_iSize) );
	istm.read( (char*)&m_cBits, sizeof(m_cBits) );
	istm.read( (char*)&iBytes, sizeof(iBytes) );
	m_aiData = new size_t[ iWords = CountWords( ) ];
	istm.read( (char*)m_aiData, iWords *= sizeof(*m_aiData) );
	istm.seekg( iBytes - iWords, ios_base::cur );

	return true; }

const unsigned char* CCompactMatrix::Open( const unsigned char* pbData ) {
	uint32_t	iBytes;

	if( m_fMemory && m_aiData )
		delete[] m_aiData;
	m_fMemory = false;
	m_iSize = *(uint32_t*)pbData;
	pbData += sizeof(m_iSize);
	m_cBits = *(unsigned char*)pbData;
	pbData += sizeof(m_cBits);
	iBytes = *(uint32_t*)pbData;
	pbData += sizeof(iBytes);
	m_aiData = (size_t*)pbData;
	pbData += iBytes;

	return pbData; }

void CCompactMatrix::Save( ostream& ostm ) const {
	static const char	abPad[]	= "\0\0\0\0\0\0\0\0";
	uint32_t		iPad, iBytes;
	size_t			iWords;

	ostm.write( (char*)&m_iSize, sizeof(m_iSize) );
	ostm.write( (char*)&m_cBits, sizeof(m_cBits) );
	iBytes = ( iWords = CountWords( ) ) * sizeof(*m_aiData);
	// Pad for 64-bit memmapping compatibility
	iBytes += ( iPad = ( iBytes % 8 ) );
	ostm.write( (char*)&iBytes, sizeof(iBytes) );
	ostm.write( (char*)m_aiData, iWords * sizeof(*m_aiData) );
	ostm.write( abPad, iPad ); }

}

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
#include "stdafx.h"
#include "compactmatrix.h"

namespace Sleipnir {

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
	//printf("size %d\n", iWords);
	if( fClear )
		memset( m_aiData, 0, iWords * sizeof(*m_aiData) ); }

///////////////////////////////////////////////////////////////////////////////
// CCompactMatrixImpl
///////////////////////////////////////////////////////////////////////////////

/*!
 * \brief
 * Load a compact half matrix from the given binary stream.
 * 
 * \param istm
 * Stream from which matrix is loaded.
 * 
 * \returns
 * True if matrix was loaded successfully.
 * 
 * \remarks
 * istm must be binary and contain a compact half matrix stored by CCompactMatrix::Save.
 */
bool CCompactMatrix::Open( std::istream& istm ) {
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

/*!
 * \brief
 * Initialize a new compact half matrix backed by the given bytes.
 * 
 * \param pbData
 * Bytes corresponding to a saved compact half matrix.
 * 
 * \returns
 * Pointer to the first byte of pbData not used by the compact half matrix.  This can be used to store
 * several compact matrices adjacently in memory and open them in serial.
 * 
 * \remarks
 * pbData must correspond to a compact half matrix stored by CCompactMatrix::Save.  This is useful for
 * memory mapping a stored file, since the opened matrix will store a pointer to the given data rather
 * than copying it.
 * 
 * \see
 * CCompactMatrix::Open
 */
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

/*!
 * \brief
 * Save a compact half matrix to the given binary stream.
 * 
 * \param ostm
 * Stream to which matrix is saved.
 * 
 * \see
 * CCompactMatrix::Open
 */
void CCompactMatrix::Save( std::ostream& ostm ) const {
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

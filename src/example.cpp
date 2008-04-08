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
#include "examplei.h"
#include "meta.h"

namespace Sleipnir {

CExampleImpl::CExampleImpl( ) {

	m_auData = NULL; }

CExampleImpl::~CExampleImpl( ) {

	Reset( ); }

void CExampleImpl::Reset( ) {

	if( !m_auData )
		return;

	delete[] m_auData;
	m_auData = NULL; }

void CExampleImpl::Set( size_t iLoc, float dValue, const CDataPair& Datum, size_t iMax ) {
	size_t	i;

	if( !m_auData ) {
		if( !iMax )
			return;
		m_auData = new UDatum[ iMax ];
		for( i = 0; i < iMax; ++i )
			m_auData[ i ].m_d = CMeta::GetNaN( ); }

	if( Datum.IsContinuous( ) || CMeta::IsNaN( dValue ) )
		m_auData[ iLoc ].m_d = dValue;
	else
		m_auData[ iLoc ].m_i = Datum.Quantize( dValue ); }

bool CExampleImpl::Equals( const CExampleImpl& Example, size_t iSize ) const {
	size_t	i;

	for( i = 0; i < iSize; ++i )
		if( m_auData[ i ].m_i != Example.m_auData[ i ].m_i )
			return false;

	return true; }

bool CExampleImpl::IsEvidence( size_t iMax ) const {
	size_t	i;

	if( !m_auData )
		return false;

	for( i = 1; i < iMax; ++i )
		if( !CMeta::IsNaN( m_auData[ i ].m_d ) )
			break;

	return ( i < iMax ); }

}

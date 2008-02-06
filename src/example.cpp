#include "stdafx.h"
#include "examplei.h"
#include "meta.h"

namespace libBioUtils {

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

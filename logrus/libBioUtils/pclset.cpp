#include "stdafx.h"
#include "pclset.h"
#include "pcl.h"
#include "meta.h"

namespace libBioUtils {

CPCLSetImpl::CPCLSetImpl( ) : m_aPCLs(NULL), m_iPCLs(0) { }

CPCLSetImpl::~CPCLSetImpl( ) {

	Reset( ); }

void CPCLSetImpl::Reset( ) {

	if( m_aPCLs ) {
		delete[] m_aPCLs;
		m_aPCLs = NULL; }
	m_iPCLs = NULL;
	m_vecstrGenes.clear( );
	m_Genes.Reset( ); }

bool CPCLSet::Open( const vector<string>& vecstrData, size_t iSkip ) {
	size_t						i, j;
	ifstream					ifsm;
	set<string>					setstrGenes;
	set<string>::const_iterator	iterGene;

	Reset( );
	m_aPCLs = new CPCL[ m_iPCLs = vecstrData.size( ) ];
	for( i = 0; i < m_iPCLs; ++i ) {
		ifsm.clear( );
		ifsm.open( vecstrData[ i ].c_str( ) );
		if( !m_aPCLs[ i ].Open( ifsm, iSkip ) )
			return false;
		ifsm.close( );
		for( j = 0; j < m_aPCLs[ i ].GetGenes( ); ++j )
			setstrGenes.insert( m_aPCLs[ i ].GetGene( j ) ); }
	m_vecstrGenes.resize( setstrGenes.size( ) );
	for( iterGene = setstrGenes.begin( ),i = 0; iterGene != setstrGenes.end( );
		++iterGene,++i )
		m_vecstrGenes[ i ] = *iterGene;

	m_Genes.Initialize( m_iPCLs, m_vecstrGenes.size( ) );
	for( i = 0; i < m_iPCLs; ++i )
		for( j = 0; j < m_vecstrGenes.size( ); ++j )
			m_Genes.Set( i, j, m_aPCLs[ i ].GetGene( m_vecstrGenes[ j ] ) );

	return true; }

size_t CPCLSet::GetGenes( ) const {

	return m_vecstrGenes.size( ); }

size_t CPCLSet::GetPCLs( ) const {

	return m_iPCLs; }

size_t CPCLSet::GetExperiments( size_t iPCL ) const {

	return m_aPCLs[ iPCL ].GetExperiments( ); }

float CPCLSet::Get( size_t iPCL, size_t iGene, size_t iExp ) const {
	size_t	iMap;

	if( ( iMap = m_Genes.Get( iPCL, iGene ) ) == -1 )
		return CMeta::GetNaN( );

	return m_aPCLs[ iPCL ].Get( iMap, iExp ); }

size_t CPCLSet::GetGene( const string& strGene ) const {
	TMapStrI::const_iterator	iterGene;

	return ( ( ( iterGene = m_mapGenes.find( strGene ) ) == m_mapGenes.end( ) ) ? -1 :
		iterGene->second ); }

const string& CPCLSet::GetGene( size_t iGene ) const {

	return m_vecstrGenes[ iGene ]; }

const vector<string>& CPCLSet::GetGeneNames( ) const {

	return m_vecstrGenes; }

}

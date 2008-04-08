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
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, and Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "stdafx.h"
#include "pclset.h"

namespace Sleipnir {

CPCLSetImpl::CPCLSetImpl( ) : m_aPCLs(NULL), m_iPCLs(0) { }

CPCLSetImpl::~CPCLSetImpl( ) {

	Reset( ); }

void CPCLSetImpl::Reset( ) {

	if( m_aPCLs ) {
		delete[] m_aPCLs;
		m_aPCLs = NULL; }
	m_iPCLs = 0;
	m_vecstrGenes.clear( );
	m_Genes.Reset( ); }

/*!
 * \brief
 * Create a new PCL set by opening and aligning the given PCL filenames.
 * 
 * \param vecstrFiles
 * Vector of PCL filenames to be opened.
 * 
 * \param iSkip
 * Number of feature columns to skip between the gene IDs and experimental data columns in each file.
 * 
 * \param eNormalize
 * Way in which each PCL should be normalized.
 * 
 * \returns
 * True if the PCL set was opened successfully.
 * 
 * \remarks
 * The same number of skip columns must be used for each PCL in the set.
 */
bool CPCLSet::Open( const std::vector<std::string>& vecstrFiles, size_t iSkip, CPCL::ENormalize eNormalize ) {
	size_t						i, j;
	ifstream					ifsm;
	set<string>					setstrGenes;
	set<string>::const_iterator	iterGene;

	Reset( );
	m_aPCLs = new CPCL[ m_iPCLs = vecstrFiles.size( ) ];
	for( i = 0; i < m_iPCLs; ++i ) {
		ifsm.clear( );
		ifsm.open( vecstrFiles[ i ].c_str( ) );
		if( !m_aPCLs[ i ].Open( ifsm, iSkip ) )
			return false;
		ifsm.close( );
		m_aPCLs[ i ].Normalize( eNormalize );
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

}

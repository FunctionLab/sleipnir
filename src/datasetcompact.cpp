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
#include "dataset.h"
#include "bayesnet.h"
#include "genome.h"
#include "compactmatrix.h"

namespace Sleipnir {

CDatasetCompactImpl::CDatasetCompactImpl( ) : m_iData(0), m_aData(NULL) {

	m_fContinuous = false; }

CDatasetCompactImpl::~CDatasetCompactImpl( ) {

	if( m_aData )
		delete[] m_aData; }

/*!
 * \brief
 * Construct a dataset corresponding to the given files.
 * 
 * \param vecstrDataFiles
 * Vector of file paths to load.
 * 
 * \param fMemmap
 * If true, memory map data files while they are being discretized rather than loading them into memory.
 * 
 * \returns
 * True if dataset was constructed successfully.
 * 
 * Creates a dataset with nodes corresponding to the given data files.
 */
bool CDatasetCompact::Open( const std::vector<std::string>& vecstrDataFiles, bool fMemmap ) {
	size_t	i;

	if( !OpenGenes( vecstrDataFiles ) )
		return false;
	if( m_aData )
		delete[] m_aData;
	m_aData = new CCompactMatrix[ m_iData = (uint32_t)vecstrDataFiles.size( ) ];

	for( i = 0; i < vecstrDataFiles.size( ); ++i ) {
		CDataPair	Datum;

		if( !( Datum.Open( vecstrDataFiles[ i ].c_str( ), false, fMemmap ) &&
			CDatasetCompactImpl::Open( Datum, i ) ) )
			return false; }

	return true; }

struct SIsGene {
	const CGenes&	m_Genes;
	bool			m_fIn;

	SIsGene( const CGenes& Genes, bool fIn ) : m_Genes(Genes), m_fIn(fIn) { }

	bool operator()( const string& strGene ) {

		return ( m_fIn == m_Genes.IsGene( strGene ) ); }
};

/*!
 * \brief
 * Construct a dataset corresponding to the given Bayes net using the provided answer file and data
 * files from the given directory.
 * 
 * \param Answers
 * Pre-loaded answer file which will become the first node of the dataset.
 * 
 * \param szDataDirectory
 * Directory from which data files are loaded.
 * 
 * \param pBayesNet
 * Bayes nets whose nodes will correspond to files in the dataset.
 * 
 * \param fEverything
 * If true, load all data; if false, load only data for gene pairs with values in the given answer file.
 * 
 * \returns
 * True if dataset was constructed successfully.
 * 
 * Creates a dataset with nodes corresponding to the given Bayes net structure; the given answer file
 * is always inserted as the first (0th) data file, and thus corresponds to the first node in the Bayes
 * net (generally the class node predicting functional relationships).  Nodes for which a corresponding
 * data file (i.e. one with the same name followed by an appropriate CDat extension) cannot be located
 * are marked as hidden.
 * 
 * \see
 * CDataset::Open
 */
bool CDatasetCompact::Open( const CDataPair& Answers, const char* szDataDirectory, const IBayesNet* pBayesNet,
	bool fEverything ) {
	CGenome	Genome;
	CGenes	GenesIn( Genome ), GenesEx( Genome );

	return Open( Answers, szDataDirectory, pBayesNet, GenesIn, GenesEx, fEverything ); }

/*!
 * \brief
 * Construct a dataset corresponding to the given Bayes net using the provided answer file and data
 * files from the given directory.
 * 
 * \param Answers
 * Pre-loaded answer file which will become the first node of the dataset.
 * 
 * \param szDataDirectory
 * Directory from which data files are loaded.
 * 
 * \param pBayesNet
 * Bayes nets whose nodes will correspond to files in the dataset.
 * 
 * \param GenesInclude
 * Data is filtered using FilterGenes with CDat::EFilterInclude and the given gene set (unless empty).
 * 
 * \param GenesExclude
 * Data is filtered using FilterGenes with CDat::EFilterExclude and the given gene set (unless empty).
 * 
 * \param fEverything
 * If true, load all data; if false, load only data for gene pairs with values in the given answer file.
 * 
 * \returns
 * True if dataset was constructed successfully.
 * 
 * Creates a dataset with nodes corresponding to the given Bayes net structure; the given answer file
 * is always inserted as the first (0th) data file, and thus corresponds to the first node in the Bayes
 * net (generally the class node predicting functional relationships).  Nodes for which a corresponding
 * data file (i.e. one with the same name followed by an appropriate CDat extension) cannot be located
 * are marked as hidden.
 * 
 * \see
 * CDataset::Open
 */
bool CDatasetCompact::Open( const CDataPair& Answers, const char* szDataDirectory, const IBayesNet* pBayesNet,
	const CGenes& GenesInclude, const CGenes& GenesExclude, bool fEverything ) {
	size_t			i;
	vector<string>	vecstrData, vecstrNodes;
	set<string>		setstrGenes;

	if( pBayesNet->IsContinuous( ) )
		return false;

	pBayesNet->GetNodes( vecstrNodes );
	m_iData = 1 + (uint32_t)OpenMax( szDataDirectory, vecstrNodes, true, vecstrData, fEverything ?
		&setstrGenes : NULL );
	m_veccQuants.resize( m_iData );
	if( m_aData )
		delete[] m_aData;
	m_aData = new CCompactMatrix[ m_iData ];

	if( fEverything ) {
		m_vecstrGenes.resize( setstrGenes.size( ) );
		copy( setstrGenes.begin( ), setstrGenes.end( ), m_vecstrGenes.begin( ) ); }
	else {
		m_vecstrGenes.resize( Answers.GetGenes( ) );
		for( i = 0; i < m_vecstrGenes.size( ); ++i )
			m_vecstrGenes[ i ] = Answers.GetGene( i ); }
	if( GenesInclude.GetGenes( ) )
		remove_if( m_vecstrGenes.begin( ), m_vecstrGenes.end( ), SIsGene( GenesInclude, false ) );
	if( GenesExclude.GetGenes( ) )
		remove_if( m_vecstrGenes.begin( ), m_vecstrGenes.end( ), SIsGene( GenesExclude, true ) );

	if( !CDatasetCompactImpl::Open( Answers, 0 ) )
		return false;
	for( i = 0; i < vecstrData.size( ); ++i ) {
		CDataPair	Datum;

		if( !( Datum.Open( vecstrData[ i ].c_str( ), false ) &&
			CDatasetCompactImpl::Open( Datum, i + 1 ) ) )
			return false; }

/*
	for( i = 0; i < m_vecstrGenes.size( ); ++i )
		for( j = ( i + 1 ); j < m_vecstrGenes.size( ); ++j ) {
			for( k = 1; k < m_iData; ++k )
				if( m_aData[ k ].Get( i, j ) )
					break;
			if( k >= m_iData )
				m_aData[ 0 ].Set( i, j, 0 ); }
*/

	return true; }

/*!
 * \brief
 * Construct a dataset corresponding to the given Bayes net using the provided answer file and data files.
 * 
 * \param Answers
 * Pre-loaded answer file which will become the first node of the dataset.
 * 
 * \param vecstrDataFiles
 * Vector of file paths to load.
 * 
 * \param fEverything
 * If true, load all data; if false, load only data for gene pairs with values in the given answer file.
 * 
 * \param fMemmap
 * If true, memory map data files while they are being discretized rather than loading them into memory.
 * 
 * \param iSkip
 * If any of the given files is a PCL, the number of columns to skip between the ID and experiments.
 * 
 * \param fZScore
 * If true and any of the given files is a PCL, z-score similarity measures after pairwise calculation.
 * 
 * \returns
 * True if dataset was constructed successfully.
 * 
 * Creates a dataset with nodes corresponding to the given answer and data files.  The given data file
 * names are each loaded using CDat::Open.  The answer file is always inserted as the first (0th) data file.
 * 
 * \remarks
 * The same number of skip columns and z-score setting will be used for all PCLs.
 */
bool CDatasetCompact::Open( const CDataPair& Answers, const std::vector<std::string>& vecstrDataFiles,
	bool fEverything, bool fMemmap, size_t iSkip, bool fZScore ) {
	size_t	i, j, k;

	if( Answers.GetGenes( ) && Answers.IsContinuous( ) )
		return false;

	m_veciMapping.resize( m_iData = 1 + vecstrDataFiles.size( ) );
	for( i = 0; i < m_veciMapping.size( ); ++i )
		m_veciMapping[ i ] = i;
	m_veccQuants.resize( m_iData );
	if( m_aData )
		delete[] m_aData;
	m_aData = new CCompactMatrix[ m_iData ];

	if( fEverything ) {
		set<string>	setstrGenes;

		for( i = 0; i < Answers.GetGenes( ); ++i )
			setstrGenes.insert( Answers.GetGene( i ) );
		for( i = 0; i < vecstrDataFiles.size( ); ++i ) {
			CDat	Dat;

			if( !Dat.OpenGenes( vecstrDataFiles[ i ].c_str( ), iSkip ) )
					return false;
			for( j = 0; j < Dat.GetGenes( ); ++j )
				setstrGenes.insert( Dat.GetGene( j ) ); }
		m_vecstrGenes.resize( setstrGenes.size( ) );
		copy( setstrGenes.begin( ), setstrGenes.end( ), m_vecstrGenes.begin( ) ); }
	else {
		m_vecstrGenes.resize( Answers.GetGenes( ) );
		for( i = 0; i < m_vecstrGenes.size( ); ++i )
			m_vecstrGenes[ i ] = Answers.GetGene( i ); }

	if( !CDatasetCompactImpl::Open( Answers, 0 ) )
		return false;
	for( i = 0; i < vecstrDataFiles.size( ); ++i ) {
		CDataPair	Datum;

		if( !( Datum.Open( vecstrDataFiles[ i ].c_str( ), false, fMemmap, iSkip, fZScore ) &&
			CDatasetCompactImpl::Open( Datum, i + 1 ) ) )
			return false; }

	if( !fEverything && ( m_iData > 1 ) )
		for( i = 0; i < m_vecstrGenes.size( ); ++i )
			for( j = ( i + 1 ); j < m_vecstrGenes.size( ); ++j ) {
				for( k = 1; k < m_iData; ++k )
					if( m_aData[ k ].Get( i, j ) )
						break;
				if( k >= m_iData )
					m_aData[ 0 ].Set( i, j, 0 ); }

	return true; }

bool CDatasetCompactImpl::Open( const CDataPair& Datum, size_t iExp ) {
	vector<size_t>	veciGenes;
	size_t			i, j, iOne, iTwo;
	float			d;
	CCompactMatrix&	Target	= m_aData[ iExp ];

	m_veccQuants[ iExp ] = Datum.IsContinuous( ) ? -1 : Datum.GetValues( );
	Target.Initialize( m_vecstrGenes.size( ), (unsigned char)( Datum.GetValues( ) + 1 ),
		true );
	veciGenes.resize( m_vecstrGenes.size( ) );
	for( i = 0; i < veciGenes.size( ); ++i )
		veciGenes[ i ] = Datum.GetGene( m_vecstrGenes[ i ] );

	for( i = 0; i < veciGenes.size( ); ++i )
		if( ( iOne = veciGenes[ i ] ) != -1 )
			for( j = ( i + 1 ); j < veciGenes.size( ); ++j )
				if( ( ( iTwo = veciGenes[ j ] ) != -1 ) &&
					!CMeta::IsNaN( d = Datum.Get( iOne, iTwo ) ) )
					Target.Set( i, j, (unsigned char)( Datum.Quantize( d ) + 1 ) );

	return true; }

bool CDatasetCompactImpl::Open( const char* szDataDir, const IBayesNet* pBayesNet,
	const CGenes* pGenesIn, const CGenes* pGenesEx ) {
	size_t						i;
	vector<string>				vecstrData, vecstrNodes;
	set<string>					setstrGenes;
	set<string>::const_iterator	iterGenes;

	if( pBayesNet->IsContinuous( ) )
		return false;

	pBayesNet->GetNodes( vecstrNodes );
	m_iData = (uint32_t)OpenMax( szDataDir, vecstrNodes, false, vecstrData, &setstrGenes );
	m_veccQuants.resize( m_iData );
	if( pGenesIn )
		for( i = 0; i < pGenesIn->GetGenes( ); ++i )
			setstrGenes.insert( pGenesIn->GetGene( i ).GetName( ) );
	if( pGenesEx )
		for( i = 0; i < pGenesEx->GetGenes( ); ++i )
			setstrGenes.erase( pGenesEx->GetGene( i ).GetName( ) );
	m_vecstrGenes.resize( setstrGenes.size( ) );
	for( i = 0,iterGenes = setstrGenes.begin( ); iterGenes != setstrGenes.end( );
		++iterGenes )
		m_vecstrGenes[ i++ ] = *iterGenes;

	if( m_aData )
		delete[] m_aData;
	m_aData = new CCompactMatrix[ m_iData ];

	for( i = 0; i < vecstrData.size( ); ++i ) {
		CDataPair	Datum;

		if( !( Datum.Open( vecstrData[ i ].c_str( ), false ) &&
			CDatasetCompactImpl::Open( Datum, i ) ) )
			return false; }

	return true; }

/*!
 * \brief
 * Remove values from the dataset based on the given gene file and filter type.
 * 
 * \param szGenes
 * File from which gene names are loaded, one per line.
 * 
 * \param eFilter
 * Way in which to use the given genes to remove values.
 * 
 * Remove values and genes (by removing all incident edges) from the dataset based on one of several
 * algorithms.  For details, see CDat::EFilter.
 * 
 * \remarks
 * Generally implemented using Remove; clears the filtered data.
 * 
 * \see
 * CDat::FilterGenes
 */
bool CDatasetCompact::FilterGenes( const char* szGenes, CDat::EFilter eFilter ) {
	ifstream	ifsm;
	CGenome		Genome;
	CGenes		Genes( Genome );

	ifsm.open( szGenes );
	if( !( ifsm.is_open( ) && Genes.Open( ifsm ) ) )
		return false;
	FilterGenes( Genes, eFilter );

	return true; }

/*!
 * \brief
 * Removes all data for gene pairs lacking a value in the answer (0th) data file.
 * 
 * \see
 * Remove
 */
void CDatasetCompact::FilterAnswers( ) {
	size_t	i, j;

	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( IsExample( i, j ) && ( GetDiscrete( i, j, 0 ) == -1 ) )
				Remove( i, j ); }

size_t CDatasetCompactImpl::GetDiscrete( size_t iX, size_t iY, size_t iNode ) const {
	size_t	iMap;

	if( ( iMap = m_veciMapping[ iNode ] ) == -1 )
		return -1;

	return ( m_aData[ iMap ].Get( iX, iY ) - 1 ); }

bool CDatasetCompactImpl::IsExample( size_t iX, size_t iY ) const {
	size_t	i;

	for( i = 0; i < m_iData; ++i )
		if( m_aData[ i ].Get( iX, iY ) )
			return true;

	return false; }

void CDatasetCompactImpl::Remove( size_t iX, size_t iY ) {
	size_t	i;

	for( i = 0; i < m_iData; ++i )
		m_aData[ i ].Set( iX, iY, 0 ); }

bool CDatasetCompactImpl::Open( const unsigned char* pbData ) {
	size_t	i;

	if( m_aData )
		delete[] m_aData;

	if( !( pbData = CDataImpl::OpenBinary( pbData ) ) )
		return false;
	m_iData = *(uint32_t*)pbData;
	pbData += sizeof(m_iData);
	m_aData = new CCompactMatrix[ m_iData ];
	for( i = 0; i < m_iData; ++i )
		if( !( pbData = m_aData[ i ].Open( pbData ) ) )
			return false;

	return true; }

/*!
 * \brief
 * Load a binary DAD dataset from the given binary stream.
 * 
 * \param istm
 * Stream from which dataset is loaded.
 * 
 * \returns
 * True if dataset was loaded successfully.
 * 
 * \remarks
 * Should be generated by Save; only used with binary DADs.
 */
bool CDatasetCompact::Open( std::istream& istm ) {
	size_t	i;

	if( m_aData )
		delete[] m_aData;

	if( !CDataImpl::OpenBinary( istm ) )
		return false;
	istm.read( (char*)&m_iData, sizeof(m_iData) );
	m_aData = new CCompactMatrix[ m_iData ];
	for( i = 0; i < m_iData; ++i )
		if( !m_aData[ i ].Open( istm ) )
			return false;

	return true; }

void CDatasetCompactImpl::SaveBinary( std::ostream& ostm ) const {
	size_t	i;

	CDataImpl::SaveBinary( ostm );
	ostm.write( (char*)&m_iData, sizeof(m_iData) );
	for( i = 0; i < m_iData; ++i )
		m_aData[ i ].Save( ostm ); }

void CDatasetCompactImpl::SaveText( std::ostream& ostm ) const {
	size_t	i, j, k, iVal;

	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( IsExample( i, j ) ) {
				ostm << GetGene( i ) << '\t' << GetGene( j );
				for( k = 0; k < GetExperiments( ); ++k ) {
					ostm << '\t';
					if( ( iVal = GetDiscrete( i, j, k ) ) == -1 )
						ostm << "-1";
					else
						ostm << iVal; }
				ostm << endl; } }

/*!
 * \brief
 * Constructs a dataset corresponding to the given Bayes net using the provided answer file and data matrices
 * generated from the given PCLs and similarity measure.
 * 
 * \param GenesInclude
 * Data is filtered using FilterGenes with CDat::EFilterInclude and the given gene set (unless empty).
 * 
 * \param GenesExclude
 * Data is filtered using FilterGenes with CDat::EFilterExclude and the given gene set (unless empty).
 * 
 * \param Answers
 * Pre-loaded answer file which will become the first node of the dataset.
 * 
 * \param vecstrPCLs
 * Vector of PCL file paths from which pairwise scores are calculated.
 * 
 * \param iSkip
 * The number of columns to skip between the ID and experiments in each PCL.
 * 
 * \param pMeasure
 * Similarity measure used to calculate pairwise scores between genes.
 * 
 * \param vecdBinEdges
 * Vector of values corresponding to discretization bin edges (the last of which is ignored) for all PCLs.
 * 
 * \returns
 * True if the dataset was constructed successfully.
 * 
 * Constructs a dataset by loading each PCL, converting it to pairwise scores using the given similarity
 * measure, and discretizing these scores using the given bin edges.  The given answer file and the
 * resulting matrices are collected together in order in the dataset.
 * 
 * \remarks
 * The same number of skip columns and discretization bin edges will be used for all PCLs.
 * 
 * \see
 * CPCL
 */
bool CDatasetCompact::Open( const CGenes& GenesInclude, const CGenes& GenesExclude, const CDataPair& Answers,
	const std::vector<std::string>& vecstrPCLs, size_t iSkip, const IMeasure* pMeasure,
	const std::vector<float>& vecdBinEdges ) {
	size_t					i, j, iPCL;
	set<string>				setstrGenes;
	set<string>::iterator	iterGene;

	g_CatSleipnir.notice( "CDatasetCompact::Open( %d ) opening PCL files",
		iSkip );

	m_veciMapping.resize( m_iData = 1 + (uint32_t)vecstrPCLs.size( ) );
	for( i = 0; i < m_veciMapping.size( ); ++i )
		m_veciMapping[ i ] = i;
	m_veccQuants.resize( m_iData );
	m_veccQuants[ 0 ] = Answers.GetValues( );
	for( i = 1; i < m_veccQuants.size( ); ++i )
		m_veccQuants[ i ] = (unsigned char)vecdBinEdges.size( );

	for( i = 0; i < Answers.GetGenes( ); ++i )
		setstrGenes.insert( Answers.GetGene( i ) );
	for( iPCL = 0; iPCL < vecstrPCLs.size( ); ++iPCL ) {
		ifstream	ifsm;

		ifsm.open( vecstrPCLs[ iPCL ].c_str( ) );
		if( !CDataImpl::OpenGenes( ifsm, false, true, setstrGenes ) ) {
			g_CatSleipnir.error( "CDatasetCompact::Open( %d ) could not open: %s", iSkip,
				vecstrPCLs[ iPCL ].c_str( ) );
			return false; } }
	if( GenesInclude.GetGenes( ) ) {
		for( iterGene = setstrGenes.begin( ); iterGene != setstrGenes.end( ); ++iterGene )
			if( !GenesInclude.IsGene( *iterGene ) )
				setstrGenes.erase( iterGene );
		for( i = 0; i < GenesInclude.GetGenes( ); ++i )
			setstrGenes.insert( GenesInclude.GetGene( i ).GetName( ) ); }
	if( GenesExclude.GetGenes( ) )
		for( i = 0; i < GenesExclude.GetGenes( ); ++i )
			setstrGenes.erase( GenesExclude.GetGene( i ).GetName( ) );
	m_vecstrGenes.resize( setstrGenes.size( ) );
	copy( setstrGenes.begin( ), setstrGenes.end( ), m_vecstrGenes.begin( ) );

	if( m_aData )
		delete[] m_aData;
	m_aData = new CCompactMatrix[ m_iData ];
	if( !CDatasetCompactImpl::Open( Answers, 0 ) )
		return false;

	for( iPCL = 0; iPCL < vecstrPCLs.size( ); ++iPCL ) {
		CPCL			PCL;
		ifstream		ifsm;
		CDistanceMatrix	Dist;
		CDataPair		Datum;
		vector<size_t>	veciGenes;
		vector<string>	vecstrGenes;
		size_t			iGenes, iOne, iTwo;
		const float*	adOne;

		g_CatSleipnir.notice( "CDatasetCompact::Open( %d ) opening: %s", iSkip, vecstrPCLs[ iPCL ].c_str( ) );
		ifsm.open( vecstrPCLs[ iPCL ].c_str( ) );
		if( !PCL.Open( ifsm, iSkip ) ) {
			g_CatSleipnir.error( "CDatasetCompact::Open( %d ) could not open: %s", iSkip, vecstrPCLs[ iPCL ].c_str( ) );
			return 1; }
		if( pMeasure->IsRank( ) )
			PCL.RankTransform( );

		veciGenes.resize( PCL.GetGenes( ) );
		if( GenesInclude.GetGenes( ) || GenesExclude.GetGenes( ) )
			for( i = 0; i < PCL.GetGenes( ); ++i ) {
				const string&	strGene	= PCL.GetGene( i );

				if( GenesExclude.GetGenes( ) && GenesExclude.IsGene( strGene ) )
					veciGenes[ i ] = -1;
				else if( GenesInclude.GetGenes( ) )
					veciGenes[ i ] = (unsigned int)( GenesInclude.IsGene( strGene ) ? iGenes++ : -1 );
				else
					veciGenes[ i ] = (unsigned int)iGenes++;
				if( veciGenes[ i ] != -1 )
					vecstrGenes.push_back( strGene ); }
		else {
			vecstrGenes.resize( PCL.GetGenes( ) );
			copy( PCL.GetGeneNames( ).begin( ), PCL.GetGeneNames( ).end( ), vecstrGenes.begin( ) );
			for( i = 0; i < veciGenes.size( ); ++i )
				veciGenes[ i ] = i; }
		Dist.Initialize( vecstrGenes.size( ) );
		for( i = 0; i < Dist.GetSize( ); ++i )
			for( j = ( i + 1 ); j < Dist.GetSize( ); ++j )
				Dist.Set( i, j, CMeta::GetNaN( ) );
		for( i = 0; i < PCL.GetGenes( ); ++i ) {
			if( ( iOne = veciGenes[ i ] ) == -1 )
				continue;
			adOne = PCL.Get( i );
			for( j = ( i + 1 ); j < PCL.GetGenes( ); ++j )
				if( ( iTwo = veciGenes[ j ] ) != -1 )
					Dist.Set( iOne, iTwo, (float)pMeasure->Measure( adOne, PCL.GetExperiments( ), PCL.Get( j ),
						PCL.GetExperiments( ) ) ); }

		Datum.Open( vecstrGenes, Dist );
		Datum.Normalize( false );
		Datum.SetQuants( vecdBinEdges );
		if( !CDatasetCompactImpl::Open( Datum, iPCL + 1 ) )
			return false; }

	return true; }

/*!
 * \brief
 * Randomizes the contents all except the answer (0th) data file.
 * 
 * \remarks
 * The first (0th) data node in the dataset is assumed to be an answer file and is left unchanged.
 * 
 * \see
 * CCompactMatrix::Randomize
 */
void CDatasetCompact::Randomize( ) {
	size_t	i;

	if( !m_aData )
		return;

	for( i = 1; i < m_iData; ++i )
		m_aData[ i ].Randomize( ); }

CDatasetCompactMap::CDatasetCompactMap( ) : m_pbData(NULL), m_hndlMap(0) { }

CDatasetCompactMap::~CDatasetCompactMap( ) {

	CMeta::Unmap( m_pbData, m_hndlMap, m_iData ); }

/*!
 * \brief
 * Opens a binary DAD file using memory mapping and unmasks it.
 * 
 * \param szFile
 * File from which dataset is opened.
 * 
 * \returns
 * True if dataset was opened successfully.
 */
bool CDatasetCompactMap::Open( const char* szFile ) {
	size_t	i, j;

	CMeta::MapRead( m_pbData, m_hndlMap, m_iData, szFile );
	if( !CDatasetCompactImpl::Open( m_pbData ) ) {
		CMeta::Unmap( m_pbData, m_hndlMap, m_iData );
		return false; }

	m_Mask.Initialize( GetGenes( ) );
	for( i = 0; i < m_Mask.GetSize( ); ++i )
		for( j = ( i + 1 ); j < m_Mask.GetSize( ); ++j )
			m_Mask.Set( i, j, CDatasetCompact::IsExample( i, j ) );
	return true; }

}

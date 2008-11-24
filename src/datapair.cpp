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
#include "datapair.h"
#include "meta.h"
#include "genome.h"

namespace Sleipnir {

const char	CPairImpl::c_szQuantExt[]	= ".quant";

bool CPairImpl::Open( const char* szDatafile, std::ifstream& ifsm ) {
	string		strToken;
	const char*	pc;

	strToken = szDatafile;
	strToken += c_szQuantExt;
	ifsm.open( strToken.c_str( ) );
	if( !ifsm.is_open( ) ) {
		if( !( pc = strrchr( szDatafile, '.' ) ) ) {
			g_CatSleipnir.error( "CPairImpl::Open( %s ) could not replace extension for quant file",
				szDatafile );
			return false; }
		strToken = szDatafile;
		strToken.resize( pc - szDatafile );
		strToken += c_szQuantExt;
		ifsm.clear( );
		ifsm.open( strToken.c_str( ) );
		if( !ifsm.is_open( ) ) {
			g_CatSleipnir.error( "CPairImpl::Open( %s ) could not open quant file: %s", szDatafile,
				strToken.c_str( ) );
			return false; } }

	return true; }

bool CPairImpl::Open( const char* szLine, vector<float>& vecdQuant ) {
	vector<string>	vecstrQuant;
	size_t			i;

	CMeta::Tokenize( szLine, vecstrQuant, CMeta::c_szWS, true );
	vecdQuant.resize( vecstrQuant.size( ) );
	for( i = 0; i < vecdQuant.size( ); ++i )
		vecdQuant[ i ] = (float)atof( vecstrQuant[ i ].c_str( ) );

	return true; }

/*!
 * \brief
 * Construct an unbinned CDat from the given ontology slim.
 * 
 * \param Slim
 * Set of ontology terms from which to generate a QUANT-less data pair.
 * 
 * \returns
 * True if data pair was generated successfully.
 * 
 * \remarks
 * Quantize will behave inconsistently if the data pair is not assigned bin edges through some other means.
 * 
 * \see
 * CDat::Open
 */
bool CDataPair::Open( const CSlim& Slim ) {

	Reset( false );
	return CDat::Open( Slim ); }

/*!
 * \brief
 * Open the given data file as a CDat and load discretization bin edges from an accompanying QUANT file.
 * 
 * \param szDatafile
 * Filename from which CDat is loaded.
 * 
 * \param fContinuous
 * If true, do not load an associated QUANT file and only open the underlying CDat.
 * 
 * \param fMemmap
 * If true, memory map file rather than allocating memory and copying its contents.
 * 
 * \param iSkip
 * If the given file is a PCL, the number of columns to skip between the ID and experiments.
 * 
 * \param fZScore
 * If true and the given file is a PCL, z-score similarity measures after pairwise calculation.
 * 
 * \returns
 * True if data pair was successfully opened.
 * 
 * \see
 * CDat::Open
 */
bool CDataPair::Open( const char* szDatafile, bool fContinuous, bool fMemmap, size_t iSkip,
	bool fZScore ) {

	g_CatSleipnir.notice( "CDataPair::Open( %s, %d )", szDatafile, fContinuous );

	Reset( fContinuous );
	if( !CDat::Open( szDatafile, fMemmap, iSkip, fZScore ) )
		return false;
	return ( m_fContinuous ? true : OpenQuants( szDatafile ) ); }

/*!
 * \brief
 * Open only the QUANT file associated with the given data file name.
 * 
 * \param szDatafile
 * CDat filename for which the accompanying QUANT file should be loaded.
 * 
 * \returns
 * True if bin edges were loaded successfully.
 * 
 * \remarks
 * Get calls to the underlying CDat will behave inconsistently unless data is loaded through some other means;
 * this method will only load the data pair's bin edges (which allows Quantize calls to be made).
 */
bool CDataPair::OpenQuants( const char* szDatafile ) {
	static const size_t	c_iBuf	= 8192;
	char		szBuf[ c_iBuf ];
	ifstream	ifsm;

	if( !CPairImpl::Open( szDatafile, ifsm ) )
		return false;
	ifsm.getline( szBuf, c_iBuf - 1 );
	ifsm.close( );
	return CPairImpl::Open( szBuf, m_vecdQuant ); }

/*!
 * \brief
 * Return the discretized form of the given value using the data pair's current bin edges.
 * 
 * \param dValue
 * Continuous value to be discretized.
 * 
 * \returns
 * Discretized version of the given value, less than GetValues; -1 if the given value is not finite.
 * 
 * Discretizes a given continuous value using the data pair's bin edges.  Standard usage is:
 * \code
 * DP.Quantize( DP.Get( i, j ) );
 * \endcode
 * 
 * \see
 * SetQuants | CMeta::Quantize
 */
size_t CDataPair::Quantize( float dValue ) const {

	return CMeta::Quantize( dValue, m_vecdQuant ); }

void CDataPairImpl::Reset( bool fContinuous ) {

	m_vecdQuant.clear( );
	m_fContinuous = fContinuous; }

/*!
 * \brief
 * Set the data pair's bin edges.
 * 
 * \param adBinEdges
 * Array of values corresponding to discretization bin edges (the last of which is ignored).
 * 
 * \param iBins
 * Number of discretization bins.
 * 
 * \see
 * GetValues | Quantize
 */
void CDataPair::SetQuants( const float* adBinEdges, size_t iBins ) {

	Reset( false );
	m_vecdQuant.resize( iBins );
	copy( adBinEdges, adBinEdges + iBins, m_vecdQuant.begin( ) ); }

/*!
 * \brief
 * Set the data pair's bin edges.
 * 
 * \param vecdBinEdges
 * Vector of values corresponding to discretization bin edges (the last of which is ignored).
 * 
 * \see
 * GetValues | Quantize
 */
void CDataPair::SetQuants( const vector<float>& vecdBinEdges ) {

	Reset( false );
	m_vecdQuant.resize( vecdBinEdges.size( ) );
	copy( vecdBinEdges.begin( ), vecdBinEdges.end( ), m_vecdQuant.begin( ) ); }

/*!
 * \brief
 * Open the given data file as a PCL and load discretization bin edges from an accompanying QUANT file.
 * 
 * \param szDatafile
 * Filename from which PCL is loaded.
 * 
 * \param iSkip
 * Number of columns to skip between the ID and experiments.
 * 
 * \returns
 * True if PCL pair was generated successfully.
 * 
 * \see
 * CPCL::Open
 */
bool CPCLPair::Open( const char* szDatafile, size_t iSkip ) {
	static const size_t	c_iBuf	= 8192;
	char		szBuf[ c_iBuf ];
	ifstream	ifsm;
	size_t		i;

	g_CatSleipnir.notice( "CPCLPair::Open( %s )", szDatafile );

	ifsm.open( szDatafile );
	if( !CPCL::Open( ifsm, iSkip ) )
		return false;
	ifsm.close( );

	ifsm.clear( );
	if( !CPairImpl::Open( szDatafile, ifsm ) )
		return false;

	m_vecvecdQuants.resize( GetExperiments( ) );
	for( i = 0; i < m_vecvecdQuants.size( ); ++i ) {
		if( ifsm.eof( ) ) {
			g_CatSleipnir.error( "CPCLPair::Open( %s, %d ) invalid quant file", szDatafile, iSkip );
			return false; }
		ifsm.getline( szBuf, c_iBuf - 1 );
		if( !CPairImpl::Open( szBuf, m_vecvecdQuants[ i ] ) )
			return false; }

	return true; }

/*!
 * \brief
 * Return the discretized form of the given value using the PCL pair's current bin edges.
 * 
 * \param dValue
 * Continuous value to be discretized.
 * 
 * \param iExperiment
 * Experiment index whose bin edges should be used for discretization.
 * 
 * \returns
 * Discretized version of the given value, less than the number of bins; -1 if the given value is not finite.
 * 
 * \see
 * CDataPair::Quantize
 */
size_t CPCLPair::Quantize( float dValue, size_t iExperiment ) const {

	return CMeta::Quantize( dValue, m_vecvecdQuants[ iExperiment ] ); }

bool CDatFilterImpl::Attach( const CDataPair* pDat, const CDatFilter* pFilter, const CGenes* pGenes,
	CDat::EFilter eFilter, const CDat* pAnswers ) {
	size_t	i;

	if( !( ( eFilter == CDat::EFilterInclude ) || ( eFilter == CDat::EFilterExclude ) ||
		( eFilter == CDat::EFilterTerm ) || ( eFilter == CDat::EFilterEdge ) ) )
		return false;

	m_pFilter = pFilter;
	m_pDat = pDat;
	m_pAnswers = pAnswers;
	m_eFilter = eFilter;
	if( m_pAnswers ) {
		m_veciAnswers.resize( GetGenes( ) );
		for( i = 0; i < m_veciAnswers.size( ); ++i )
			m_veciAnswers[ i ] = m_pAnswers->GetGene( GetGene( i ) ); }
	else {
		m_veciAnswers.clear( );
		if( m_eFilter == CDat::EFilterTerm )
			m_eFilter = CDat::EFilterEdge; }

	if( pGenes ) {
		m_vecfGenes.resize( GetGenes( ) );
		for( i = 0; i < m_vecfGenes.size( ); ++i )
			m_vecfGenes[ i ] = pGenes->IsGene( GetGene( i ) ); }
	else
		m_vecfGenes.clear( );

	return true; }

size_t CDatFilterImpl::GetGenes( ) const {

	return ( m_pFilter ? m_pFilter->GetGenes( ) : ( m_pDat ? m_pDat->GetGenes( ) : -1 ) ); }

string CDatFilterImpl::GetGene( size_t iGene ) const {

	return ( m_pFilter ? m_pFilter->GetGene( iGene ) : ( m_pDat ? m_pDat->GetGene( iGene ) : "" ) ); }

/*!
 * \brief
 * Associates the data filter with the given CDat, gene set, and filter type.
 * 
 * \param Dat
 * CDat to be associated with the overlaying mask.
 * 
 * \param Genes
 * Gene set used to filter the data.
 * 
 * \param eFilter
 * Way in which to use the given genes to remove gene pairs.
 * 
 * \param pAnswers
 * If non-null, answer set to be used for filter types requiring answers (e.g. CDat::EFilterTerm).
 * 
 * \returns
 * True if filter was attached successfully.
 * 
 * \remarks
 * No calculation occurs during Attach; the gene set and filter type are stored, and calculations are
 * performed dynamically in IsExample.  The gene set is not copied, so the given object should not be
 * destroyed until after the filter.  If a filter types requiring an answer file is given without an
 * accompanying answer file, the filter won't crash, but results might be a little odd.
 */
bool CDatFilter::Attach( const CDataPair& Dat, const CGenes& Genes, CDat::EFilter eFilter,
	const CDat* pAnswers ) {

	return CDatFilterImpl::Attach( &Dat, NULL, &Genes, eFilter, pAnswers ); }

/*!
 * \brief
 * Associates (overlays0 the data filter with the given pre-existing filter, gene set, and filter type.
 * 
 * \param Dat
 * Filter with which this mask will associate (and overlay).
 * 
 * \param Genes
 * Gene set used to filter the data.
 * 
 * \param eFilter
 * Way in which to use the given genes to remove gene pairs.
 * 
 * \param pAnswers
 * If non-null, answer set to be used for filter types requiring answers (e.g. CDat::EFilterTerm).
 * 
 * \returns
 * True if filter was attached successfully.
 * 
 * \remarks
 * This is a slightly ugly way to allow additive filters without requiring a virtual interface.  The
 * CDat class is so fundamental, and CDat::Get calls occur so frequently, that a virtual function call can
 * actually have a noticeable impact on performance.
 */
bool CDatFilter::Attach( const CDatFilter& Dat, const CGenes& Genes, CDat::EFilter eFilter,
	const CDat* pAnswers ) {

	return CDatFilterImpl::Attach( NULL, &Dat, &Genes, eFilter, pAnswers ); }

}

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
#ifndef DATAPAIR_H
#define DATAPAIR_H

#include "datapairi.h"

namespace Sleipnir {

class CSlim;

/*!
 * \brief
 * Encapsulates a CDat paired with a quantization file.
 * 
 * A data pair consists of a CDat (often on disk in DAB format) paired with quantization information.  This
 * information is generally stored in a QUANT file with the same name and location as the CDat.  For example,
 * a DAB file named <tt>data.dab</tt> and a QUANT file named <tt>data.quant</tt> might reside in the same
 * directory; these would be loaded together as a CDataPair.
 * 
 * A QUANT file consists of a single line of text containing tab-delimited increasing numbers.  These numbers
 * represent bin edges for discretizing the CDat associated with the QUANT.  The number of bins is equal to
 * the number of numbers in the QUANT, meaning that the last number will be ignored.  Upper bin edges are
 * inclusive, lower bin edges are exclusive.  This means that for a QUANT file containing:
 * \code
 * -0.1	0.3	0.6
 * \endcode
 * the associated CDat will be discretized into three values:
 * - 0, corresponding to values less than or equal to -0.1.
 * - 1, corresponding to values greater than -0.1 but less than or equal to 0.3.
 * - 2, corresponding to values greater than 0.3.
 * 
 * \see
 * CMeta::Quantize
 */
class CDataPair : public CDataPairImpl {
public:
	bool Open( const char* szDatafile, bool fContinuous, bool fMemmap = false, size_t iSkip = 2,
		bool fZScore = false );
	bool Open( const CSlim& Slim );
	bool OpenQuants( const char* szDatafile );
	void SetQuants( const float* adBinEdges, size_t iBins );
	void SetQuants( const std::vector<float>& vecdBinEdges );
	size_t Quantize( float dValue ) const;

	/*!
	 * \brief
	 * Returns the number of discrete values taken by this data pair.
	 * 
	 * \returns
	 * Number of discrete values taken by this data pair.
	 * 
	 * \remarks
	 * Equivalent to number of bins in the data pair and number of bin edges in the QUANT file.
	 * 
	 * \see
	 * SetQuants | Quantize
	 */
	unsigned char GetValues( ) const {

		return (unsigned char)m_vecdQuant.size( ); }

	/*!
	 * \brief
	 * Returns true if the data pair has no associated discretization information.
	 * 
	 * \returns
	 * True if the data pair has no associated discretization information.
	 * 
	 * \remarks
	 * Generally only useful with continuous Bayes nets, which themselves aren't that useful.
	 */
	bool IsContinuous( ) const {

		return m_fContinuous; }

	/*!
	 * \brief
	 * Construct a data pair from the given known gene relationships and gene sets and with no discretization
	 * information.
	 * 
	 * \param DatKnown
	 * Known pairwise scores, either positive or negative as indicated.
	 * 
	 * \param vecpOther
	 * Gene sets, either positive or nonnegative as indicated (possibly empty).
	 * 
	 * \param Genome
	 * Genome containing all genes of interest.
	 * 
	 * \param fKnownNegatives
	 * If true, DatKnown contains known negative gene pairs (0 scores); if false, it contains known related
	 * gene pairs (1 scores).  In the former case, positives are generated from pairs coannotated to the
	 * given gene sets; in the latter, negatives are generated from pairs not coannotated to the given gene
	 * sets.
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
	bool Open( const CDat& DatKnown, const std::vector<CGenes*>& vecpOther, const CGenome& Genome,
		bool fKnownNegatives ) {

		return CDat::Open( DatKnown, vecpOther, Genome, fKnownNegatives ); }

	/*!
	 * \brief
	 * Construct a new data pair with the given gene names and values and with no discretization information.
	 * 
	 * \param vecstrGenes
	 * Gene names and size to associate with the data pair.
	 * 
	 * \param MatScores
	 * Values to associate with the data pair.
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
	bool Open( const std::vector<std::string>& vecstrGenes, const CDistanceMatrix& MatScores ) {

		return CDat::Open( vecstrGenes, MatScores ); }
};

/*!
 * \brief
 * Encapsulates a CPCL paired with a quantization file.
 * 
 * A PCL pair consists of a CPCL paired with quantization information.  This information is generally stored
 * in a QUANT file with the same name and location as the PCL.  For example, a PCL file named
 * <tt>data.pcl</tt> and a QUANT file named <tt>data.quant</tt> might reside in the same directory; these
 * would be loaded together as a CPCLPair.  The discretization information from the QUANT can be used to
 * convert continuous values in the PCL into discrete values, e.g. for use with a Bayes net.
 * 
 * Unlike a CDataPair, a PCL's QUANT file should contain one line per experiment.  Each line is the equivalent
 * of one standard QUANT file, i.e. it contains tab delimited bin edges in increasing order, the largest of
 * which is ignored.  This allows the values for individual experiments to be discretized differently if so
 * desired.
 * 
 * \remarks
 * The number of lines in the QUANT file must equal the number of experiments in the PCL file.
 * 
 * \see
 * IBayesNet::Evaluate
 */
class CPCLPair : public CPCLPairImpl {
public:
	bool Open( const char* szDatafile, size_t iSkip );
	size_t Quantize( float dValue, size_t iExperiment ) const;
};

/*!
 * \brief
 * Augments a CDat with a dynamically calculated gene set filter.
 * 
 * A filter wraps an underlying CDat with a dynamically calculated filter using a gene set and
 * CDat::EFilter type.  A filtered gene pair will act like missing data; unfiltered gene pairs will be
 * retrieved from the underlying CDat.  This allows data to be temporarily hidden without modifying the
 * underlying (potentially memory mapped) CDat.
 * 
 * \remarks
 * Permanent modifications such as CDat::Set should be performed only on the underlying dataset.  Yes, an
 * interface would make this a lot cleaner, but it also makes it a lot slower (losing the ability to inline
 * calls to Get actually has a non-trivial impact on runtime).
 * 
 * \see
 * CDat::FilterGenes | CDataFilter
 */
class CDatFilter : public CDatFilterImpl {
public:

	bool Attach( const CDataPair& Dat, const CGenes& Genes, CDat::EFilter eFilter,
		const CDat* pAnswers = NULL );
	bool Attach( const CDatFilter& Dat, const CGenes& Genes, CDat::EFilter eFilter,
		const CDat* pAnswers = NULL );

	bool Attach( const CDataPair& Dat ) {

		return CDatFilterImpl::Attach( &Dat, NULL, NULL, CDat::EFilterInclude, NULL ); }

	size_t GetValues( ) const {

		return ( m_pFilter ? m_pFilter->GetValues( ) : ( m_pDat ? m_pDat->GetValues( ) : -1 ) ); }

	size_t GetGene( const std::string& strGene ) const {

		return ( m_pFilter ? m_pFilter->GetGene( strGene ) : ( m_pDat ? m_pDat->GetGene( strGene ) : -1 ) ); }

	std::string GetGene( size_t iGene ) const {

		return CDatFilterImpl::GetGene( iGene ); }

	size_t Quantize( float dValue ) const {

		return ( m_pFilter ? m_pFilter->Quantize( dValue ) : ( m_pDat ? m_pDat->Quantize( dValue ) : -1 ) ); }

	float& Get( size_t iY, size_t iX ) const {
		static float	c_dNaN	= CMeta::GetNaN( );

		if( !( m_pDat || m_pFilter ) )
			return c_dNaN;
		if( m_vecfGenes.empty( ) )
			return ( m_pFilter ? m_pFilter->Get( iY, iX ) : ( m_pDat ? m_pDat->Get( iY, iX ) : c_dNaN ) );

		switch( m_eFilter ) {
			case CDat::EFilterInclude:
				if( !( m_vecfGenes[ iX ] && m_vecfGenes[ iY ] ) )
					return c_dNaN;
				break;

			case CDat::EFilterExclude:
				if( m_vecfGenes[ iX ] || m_vecfGenes[ iY ] )
					return c_dNaN;
				break;

			case CDat::EFilterEdge:
				if( !( m_vecfGenes[ iX ] || m_vecfGenes[ iY ] ) )
					return c_dNaN;
				break;

			case CDat::EFilterTerm:
				float	d;
				size_t	iOne, iTwo;

				if( !m_pAnswers )
					return c_dNaN;
				d = ( ( ( iOne = m_veciAnswers[ iX ] ) != -1 ) && ( ( iTwo = m_veciAnswers[ iY ] ) != -1 ) ) ?
					m_pAnswers->Get( iTwo, iOne ) : CMeta::GetNaN( );
				if( !( m_vecfGenes[ iX ] || m_vecfGenes[ iY ] ) ||
					( ( m_vecfGenes[ iX ] != m_vecfGenes[ iY ] ) && !CMeta::IsNaN( d ) && ( d > 0 ) ) )
					return c_dNaN;
				break; }

		return ( m_pFilter ? m_pFilter->Get( iY, iX ) : m_pDat->Get( iY, iX ) ); }
};

}

#endif // DATAPAIR_H

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
#ifndef DAT_H
#define DAT_H

#include <iostream>
#include <string>
#include <vector>

#include "dati.h"

namespace Sleipnir {

class CGenes;
class CGenome;

/*!
 * \brief
 * Stores a continuously valued half matrix paired with a list of names for matrix elements.
 * 
 * Conceptually, a CDat stores a list of weighted pairs; this is equivalent to a weighted undirected
 * graph with labeled nodes, or a symmetric matrix with labels for each matrix element.  CDat entries are
 * stored as continuous values, although they can be discretized in various ways.  CDats can be constructed
 * in several ways, read from disk, persisted to disk in multiple file formats, or calculated from
 * existing gene sets, microarray data, or gold standards.  In practice, a CDat is simply a continuously
 * valued symmetric matrix (in which zero or more values may be missing) paired with a list of element
 * names (assumed to be genes), but this data structure is sufficiently flexible to represent nearly any
 * biological dataset.
 * 
 * CDats can be loaded (by Open) and/or stored (by Save) from/to disk in the following formats:
 * - DAT.  A tab-delimited text file in which each line contains two identifiers and a score:
 * \code
 * GENE1	GENE2	SCORE1
 * GENE1	GENE3	SCORE2
 * GENE2	GENE3	SCORE3
 * \endcode
 * Element pair order is irrelevant, missing values are allowed, and duplicates can be optionally ignored.
 * The DAT format is most suitable for human readability and manipulation by scripting languages; it is
 * much larger and slower to process than the other formats, however.
 * - DAB.  A binary file containing an integer size, a list of null-terminated element identifiers
 * (generally gene names), and the CDat's values in row-major order.  Missing values are stored as NaNs.
 * Should be generated by Save; usually the smallest and most rapidly parsed format, and the only one
 * amenable to memory mapping.
 * - DAS.  A sparse binary file containing an integer size, a list of null-terminated element identifiers
 * (generally gene names), and the CDat's non-missing values in row-major order, pairing column indices with
 * values.  Should be generated by Save.  Note that this sounds like it should save space for sparse
 * CDats, but because of the overhead of storing column indices, the matrix has to be awfully sparse before
 * it actually does.
 * - PCL.  A standard PCL file, which is loaded and converted to pairwise similarity scores using
 * z-transformed Pearson correlation as calculated by CMeasurePearNorm.  Can be converted once and cached
 * in memory or calculated on-the-fly; the former consumes more memory, the latter is (often) slower.
 * 
 * \see
 * CDataPair | CHalfMatrix
 */
class CDat : protected CDatImpl {
public:
	/*!
	 * \brief
	 * Ways in which nodes/edges can be removed to filter a CDat.
	 * 
	 * \see
	 * FilterGenes
	 */
	enum EFilter {
		/*!
		 * \brief
		 * Remove any edge including a node outside the given set.
		 */
		EFilterInclude		= 0,
		/*!
		 * \brief
		 * Remove any positive edge including a node outside the given set.
		 */
		EFilterTerm			= EFilterInclude + 1,
		/*!
		 * \brief
		 * Remove any edge including a node in the given set.
		 */
		EFilterExclude		= EFilterTerm + 1,
		/*!
		 * \brief
		 * Perform a bioPIXIE query using the given set and remove any edge not in the resulting subgraph.
		 */
		EFilterPixie		= EFilterExclude + 1,
		/*!
		 * \brief
		 * Remove any edge not including a node in the given set.
		 */
		EFilterEdge			= EFilterPixie + 1,
		/*!
		 * \brief
		 * Perform a HEFalMp query using the given set and remove any edge not in the resulting subgraph.
		 */
		EFilterHefalmp		= EFilterEdge + 1
	};

	/*!
	 * \brief
	 * Ways in which a CDat can be persisted to/from disk.
	 * 
	 * \see
	 * Open | Save
	 */
	enum EFormat {
		/*!
		 * \brief
		 * Binary format listing null-terminated element name strings followed by floating point values.
		 */
		EFormatBinary	= 0,
		/*!
		 * \brief
		 * Text format listing element name pairs followed by numerical value strings.
		 */
		EFormatText		= EFormatBinary + 1,
		/*!
		 * \brief
		 * PCL file from which pairwise scores are calculated using some similarity measure.
		 */
		EFormatPCL		= EFormatText + 1,
		/*!
		 * \brief
		 * Binary format listing null-terminated element name strings followed by index/value pairs.
		 */
		EFormatSparse	= EFormatPCL + 1
	};

	/*!
	 * \brief
	 * Ways in which a CDat can have its edge values normalized.
	 * 
	 * \see
	 * Normalize
	 */
	enum ENormalize {
		ENormalizeNone		= 0,
		/*!
		 * \brief
		 * Linearly transform the minimum score to 0 and the maximum to 1.
		 */
		ENormalizeMinMax	= ENormalizeNone + 1,
		/*!
		 * \brief
		 * Z-score all edges (subtract mean, divide by standard deviation).
		 */
		ENormalizeZScore	= ENormalizeMinMax + 1,
		/*!
		 * \brief
		 * Sigmoid transform scores to the range [0, 1].
		 */
		ENormalizeSigmoid	= ENormalizeZScore + 1
	};

	bool Open( const char* szFile, bool fMemmap = false, size_t iSkip = 2, bool fZScore = false,
		bool fDuplicates = false );
	bool Open( std::istream& istm, EFormat eFormat = EFormatBinary, float dDefault = HUGE_VAL,
		bool fDuplicates = false, size_t iSkip = 2, bool fZScore = false );
	bool Open( const CSlim& Slim );
	bool Open( const CSlim& SlimPositives, const CSlim& SlimNonnegatives );
	bool Open( const std::vector<std::string>& vecstrGenes, bool fClear = true, const char* szFile = NULL );
	bool Open( const std::vector<std::string>& vecstrGenes, const CDistanceMatrix& MatValues );
	bool Open( const std::vector<CGenes*>& vecpPositives, const std::vector<CGenes*>& vecpNonnegatives,
		float dPValue, const CGenome& Genome );
	bool Open( const CDat& DatKnown, const std::vector<CGenes*>& vecpOther, const CGenome& Genome,
		bool fKnownNegatives );
	bool Open( const CPCL& PCL, const IMeasure* pMeasure, bool fMeasureMemory );
	bool Open( const CDat& Dat );

	bool OpenGenes( std::istream& istm, bool fBinary, bool fPCL = false );
	bool OpenGenes( const char* szFile, size_t iSkip = 2 );
	void Save( std::ostream& ostm, EFormat eFormat = EFormatBinary ) const;
	void Save( const char* szFile ) const;
	void SaveDOT( std::ostream& ostm, float dCutoff = HUGE_VAL, const CGenome* pGenome = NULL,
		bool fUnlabeled = false, bool fHashes = true, const std::vector<float>* pvecdColors = NULL,
		const std::vector<float>* pvecdBorders = NULL ) const;
	void SaveGDF( std::ostream& ostm, float dCutoff = HUGE_VAL ) const;
	void SaveNET( std::ostream& ostm, float dCutoff = HUGE_VAL ) const;
	void SaveMATISSE( std::ostream& ostm, float dCutoff = HUGE_VAL, const CGenome* pGenome = NULL ) const;
	void Invert( );
	void Rank( );
	bool FilterGenes( const char* szGenes, EFilter eFilter, size_t iLimit = -1 );
	void FilterGenes( const CGenes& Genes, EFilter eFilter, size_t iLimit = -1,
		float dEdgeAggressiveness = 0.5 );

	/*!
	 * \brief
	 * Normalize each finite value in the CDat by a specific function.
	 * 
	 * \param eNormalize
	 * Method by which scores are normalized.
	 * 
	 * \remarks
	 * Values are left unchanged if ( dMax == dMin ) or ( dStd == 0 ).
	 * 
	 * \see
	 * ENormalize | Invert
	 */
	void Normalize( ENormalize eNormalize ) {

		switch( eNormalize ) {
			case ENormalizeMinMax:
				NormalizeMinmax( );
				break;

			case ENormalizeZScore:
				NormalizeStdev( );
				break;

			default:
				NormalizeSigmoid( ); } }

	/*!
	 * \brief
	 * Return the index of the given gene name, or -1 if it is not included in the CDat.
	 * 
	 * \param strGene
	 * Gene name to retrieve.
	 * 
	 * \returns
	 * Index of the requested gene name, or -1 if it is not in the CDat.
	 * 
	 * \see
	 * GetGeneNames
	 */
	size_t GetGene( const std::string& strGene ) const {

		return CDatImpl::GetGene( strGene ); }

	/*!
	 * \brief
	 * Return the value at the requested CDat position.
	 * 
	 * \param iY
	 * CDat row.
	 * 
	 * \param iX
	 * CDat column.
	 * 
	 * \returns
	 * Value at the requested CDat position.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given row and column must be smaller than
	 * GetGenes.  As a symmetric matrix, the value at position XY will always equal the value at position YX.
	 * 
	 * \see
	 * Set
	 */
	float& Get( size_t iY, size_t iX ) const {

		return CDatImpl::Get( iY, iX ); }

	/*!
	 * \brief
	 * Returns the number of elements (genes) in the CDat.
	 * 
	 * \returns
	 * Number of elements (genes) in the CDat.
	 * 
	 * \remarks
	 * Since a symmetric matrix must be square, the number of rows equals the number of columns and is thus
	 * referred to as the number of elements (genes).
	 */
	size_t GetGenes( ) const {

		return CDatImpl::GetGenes( ); }

	/*!
	 * \brief
	 * Returns the symmetric matrix containing the CDat's values.
	 * 
	 * \returns
	 * Symmetric matrix containing the CDat's values.
	 */
	const CDistanceMatrix& Get( ) const {

		return m_Data; }

	/*!
	 * \brief
	 * Returns the symmetric matrix containing the CDat's values.
	 * 
	 * \returns
	 * Symmetric matrix containing the CDat's values.
	 */
	CDistanceMatrix& Get( ) {

		return m_Data; }

	/*!
	 * \brief
	 * Set the value at the requested CDat position.
	 * 
	 * \param iY
	 * CDat row.
	 * 
	 * \param iX
	 * CDat column.
	 * 
	 * \param dValue
	 * Value to store.
	 * 
	 * \returns
	 * True if the value was stored successfully.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given row and column must be smaller than
	 * GetGenes.
	 * 
	 * \see
	 * Get
	 */
	bool Set( size_t iY, size_t iX, float dValue ) {

		return CDatImpl::Set( iY, iX, dValue ); }

	/*!
	 * \brief
	 * Returns the gene name at the given CDat position.
	 * 
	 * \param iGene
	 * Index of gene name to return.
	 * 
	 * \returns
	 * Gene name at the requested index.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given index must be smaller than GetGenes.
	 */
	std::string GetGene( size_t iGene ) const {

		return CDatImpl::GetGene( iGene ); }

	/*!
	 * \brief
	 * Returns the vector of gene names associated with this CDat.
	 * 
	 * \returns
	 * Vector of this CDat's gene names.
	 * 
	 * \remarks
	 * Returned vector size will be identical to GetGenes.
	 */
	const std::vector<std::string>& GetGeneNames( ) const {

		return CDatImpl::GetGeneNames( ); }

	/*!
	 * \brief
	 * Set an entire row of CDat values efficiently.
	 * 
	 * \param iY
	 * CDat row.
	 * 
	 * \param adValues
	 * Values to store.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given row must be smaller than GetGenes, and the
	 * given array must be non-null and have length exactly (size - iY - 1).
	 * 
	 * \see
	 * Get
	 */
	void Set( size_t iY, const float* adValues ) {

		m_Data.Set( iY, adValues ); }

	/*!
	 * \brief
	 * Get an entire row of CDat values efficiently.
	 * 
	 * \param iY
	 * CDat row.
	 * 
	 * \returns
	 * Retrieved values.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given row must be smaller than GetGenes and the
	 * returned array will have length exactly (size - iY - 1).
	 * 
	 * \see
	 * Set
	 */
	const float* Get( size_t iY ) const {

		return m_Data.Get( iY ); }

	/*!
	 * \brief
	 * Get an entire row of CDat values efficiently.
	 * 
	 * \param iY
	 * CDat row.
	 * 
	 * \returns
	 * Retrieved values.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given row must be smaller than GetGenes and the
	 * returned array will have length exactly (size - iY - 1).
	 * 
	 * \see
	 * Set
	 */
	float* Get( size_t iY ) {

		return m_Data.Get( iY ); }

	/*!
	 * \brief
	 * Set the gene name at the given index.
	 * 
	 * \param iGene
	 * Index of gene name to modify.
	 * 
	 * \param strGene
	 * Gene name to store at the requested index.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given index must be smaller than GetGenes.
	 * 
	 * \see
	 * GetGene
	 */
	void SetGene( size_t iGene, const std::string& strGene ) {

		if( m_pPCL )
			m_pPCL->SetGene( iGene, strGene );
		else
			m_vecstrGenes[ iGene ] = strGene; }

	/*!
	 * \brief
	 * Randomizes the CDat's values by iterated swapping.
	 */
	void Randomize( ) {
		size_t	i, j, iOne, iTwo;
		float	dOne, dTwo;

		for( i = 0; i < GetGenes( ); ++i )
			for( j = ( i + 1 ); j < GetGenes( ); ++j ) {
				if( CMeta::IsNaN( dOne = Get( i, j ) ) )
					continue;
				while( true ) {
					iOne = rand( ) % GetGenes( );
					iTwo = rand( ) % GetGenes( );
					if( iOne > iTwo )
						std::swap( iOne, iTwo );
					if( ( ( iOne != i ) || ( iTwo != j ) ) && !CMeta::IsNaN( dTwo = Get( iOne, iTwo ) ) )
						break; }
				Set( i, j, dTwo );
				Set( iOne, iTwo, dOne ); } }
};

}

#endif // DAT_H

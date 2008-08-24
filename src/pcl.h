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
#ifndef PCL_H
#define PCL_H

#include <algorithm>
#include <fstream>
#include <string>

#include "pcli.h"

namespace Sleipnir {

class CDat;
class IMeasure;

/*!
 * \brief
 * Encapsulates a PCL (or, in a pinch, CDT) formatted microarray data file.
 * 
 * A canonical PCL is a tab-delimited text format used to store microarray datasets.  Each row represents a
 * single gene, each column represents a single microarray, and each cell stores the gene's expression
 * value for that condition.  There are usually two header rows, one naming each column and one listing the
 * microarray conditions' weights.  There are usually three header columns: unique gene identifiers in the
 * first column, human-readable gene name synonyms in the second column, and gene weights in the third.
 * Thus, a simple PCL might be:
 * \code
 * GID	NAME	GWEIGHT	Condition 1	Condition 2	Condition 3
 * EWEIGHT			1	1	1
 * YKL032C	IXR1	1	-0.02	0.03	-0.1
 * YMR027W	HRT2	1	-0.018	0.013	-0.016
 * YLR220W	CCC1	1	-0.0050	-0.01	0.015</pre>
 * \endcode
 * 
 * In practice, the format of this file is very flexible: the first row's always a header, the first
 * column's always unique gene IDs, and most other columns are data, but everything else changes.  There
 * can be an arbitrary number of feature columns between the gene IDs and the first experimental column
 * (Sleipnir refers to these as "skip" columns).  The EWEIGHT row may or may not be present.  Missing
 * values might be notated by empty cells or the string NA.  CPCL will correctly parse all of these cases.
 * 
 * So in practice, a PCL is a full matrix of data in which each row represents a gene and each column an
 * experiment.  Associated with each gene are zero or more features, possibly including the name and weight;
 * associated with each experiment is a name (drawn from the header row).  The EWEIGHT row is discarded.
 * More generally, a PCL can be used to hold any data in which vectors of values are associated with a
 * set of genes.
 * 
 * \see
 * CDat | CDataMatrix
 */
class CPCL : protected CPCLImpl {
public:
	/*!
	 * \brief
	 * Ways in which a PCL's data entries can be normalized.
	 */
	enum ENormalize {
		/*!
		 * \brief
		 * Perform no normalization and leave PCL values unchanged.
		 */
		ENormalizeNone,
		/*!
		 * \brief
		 * Subtract the global average from every value and divide by the global standard deviation.
		 */
		ENormalizeZScore,
		/*!
		 * \brief
		 * Subtract the row average from every value and divide by the row standard deviation.
		 */
		ENormalizeRow,
		/*!
		 * \brief
		 * Subtract the global minimum from every value and divide by the global maximum (transforming
		 * all values to the range [0, 1]).
		 */
		ENormalizeMinMax
	};

	static int Distance( const char* szFile, size_t iSkip, const char* szSimilarityMeasure, bool fNormalize,
		bool fZScore, bool fAutocorrelate, const char* szGeneFile, float dCutoff, size_t iLimit, CPCL& PCL,
		CDat& Dat );

	/*!
	 * \brief
	 * Return the default number of skip columns between the gene IDs and experimental values.
	 * 
	 * \returns
	 * Default number of skip columns.
	 */
	static size_t GetSkip( ) {

		return c_iSkip; }

	/*!
	 * \brief
	 * Create a new PCL object with or without a header row.
	 * 
	 * \param fHeader
	 * If true, associate a standard PCL header row with this object.
	 * 
	 * \remarks
	 * Header setting influences both Open and Save (by way of SaveHeader).
	 */
	CPCL( bool fHeader = true ) : CPCLImpl( fHeader ) { }

	void Open( const CPCL& PCL );
	bool Open( std::istream& istm, size_t iSkip );
	void Open( const std::vector<size_t>& veciGenes, const std::vector<std::string>& vecstrGenes,
		const std::vector<std::string>& vecstrExperiments );
	void Open( const std::vector<std::string>& vecstrGenes, const std::vector<std::string>& vecstrExperiments,
		const std::vector<std::string>& vecstrFeatures );
	void OpenBinary( std::istream& istm );
	void Save( std::ostream& ostm, const std::vector<size_t>* pveciGenes = NULL ) const;
	void SaveBinary( std::ostream& ostm ) const;
	void SaveGene( std::ostream& ostm, size_t iGene, size_t iOriginal = -1 ) const;
	void SaveHeader( std::ostream& ostm, bool fCDT = false ) const;
	bool SortGenes( const std::vector<size_t>& veciOrder );
	void RankTransform( );
	bool AddGenes( const std::vector<std::string>& vecstrGenes );
	void Normalize( ENormalize eNormalize = ENormalizeRow );
	void Impute( size_t iNeighbors, float dMinimumPresent, const CDat& DatSimilarity );
	void Impute( size_t iNeighbors, float dMinimumPresent, const IMeasure* pMeasure, bool fPrecompute = true );

	/*!
	 * \brief
	 * Load a PCL from the given text file.
	 * 
	 * \param szFile
	 * File from which PCL file is loaded.
	 * 
	 * \param iSkip
	 * Number of feature columns to skip between the gene IDs and first experimental column.
	 * 
	 * \returns
	 * True if the PCL was opened successfully.
	 * 
	 * \see
	 * Save
	 */
	bool Open( const char* szFile, size_t iSkip ) {
		std::ifstream	ifsm;

		if( szFile )
			ifsm.open( szFile );
		return Open( szFile ? ifsm : cin, iSkip ); }

	/*!
	 * \brief
	 * Load a PCL from the given text stream using the default number of skip columns.
	 * 
	 * \param istm
	 * Stream from which PCL file is loaded.
	 * 
	 * \returns
	 * True if the PCL was opened successfully.
	 * 
	 * \see
	 * Save
	 */
	bool Open( std::istream& istm ) {

		return Open( istm, GetSkip( ) ); }

	/*!
	 * \brief
	 * Empties the PCL and deallocates all associated memory.
	 */
	void Reset( ) {

		CPCLImpl::Reset( ); }

	/*!
	 * \brief
	 * Sets all PCL entries to zero (without changing size or memory allocation).
	 */
	void Clear( ) {

		m_Data.Clear( ); }

	/*!
	 * \brief
	 * Return the number of features (skip columns) in the PCL.
	 * 
	 * \returns
	 * Number of feature columns in the PCL.
	 */
	size_t GetFeatures( ) const {

		return m_vecstrFeatures.size( ); }

	/*!
	 * \brief
	 * Return the label (header) of the requested feature.
	 * 
	 * \param iFeature
	 * Index of the feature label to retrieve.
	 * 
	 * \returns
	 * Feature header at the requested index.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given value must be smaller than GetFeatures.
	 */
	const std::string& GetFeature( size_t iFeature ) const {

		return m_vecstrFeatures[ iFeature ]; }

	/*!
	 * \brief
	 * Return the requested gene's value for the given feature.
	 * 
	 * \param iGene
	 * Gene for which feature value should be returned.
	 * 
	 * \param iFeature
	 * Index of feature value to be retrieved.
	 * 
	 * \returns
	 * Feature value for the given gene and feature indices.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given values must be smaller than GetGenes
	 * and GetFeatures.
	 */
	const std::string& GetFeature( size_t iGene, size_t iFeature ) const {

		return m_vecvecstrFeatures[ iFeature - 1 ][ iGene ]; }

	/*!
	 * \brief
	 * Set the requested gene's value for the given feature.
	 * 
	 * \param iGene
	 * Gene for which feature value should be set.
	 * 
	 * \param iFeature
	 * Index of feature value to be set.
	 * 
	 * \param strValue
	 * Feature value to set for the given gene and feature indices.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given values must be smaller than GetGenes
	 * and GetFeatures.
	 */
	void SetFeature( size_t iGene, size_t iFeature, const std::string& strValue ) {

		m_vecvecstrFeatures[ iFeature - 1 ][ iGene ] = strValue; }

	/*!
	 * \brief
	 * Returns the value at the requested PCL position.
	 * 
	 * \param iGene
	 * Gene row.
	 * 
	 * \param iExperiment
	 * Experiment column.
	 * 
	 * \returns
	 * Value at the requested PCL position.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given gene and experiment must be smaller than
	 * GetGenes and GetExperiments.
	 * 
	 * \see
	 * Set
	 */
	float& Get( size_t iGene, size_t iExperiment ) const {

		return m_Data.Get( iGene, iExperiment ); }

	/*!
	 * \brief
	 * Return a single gene's row from the PCL.
	 * 
	 * \param iGene
	 * PCL row.
	 * 
	 * \returns
	 * Requested PCL row.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given index must be smaller than GetGenes.
	 * The returned array will contain a number of elements equal to GetExperiments.
	 * 
	 * \see
	 * Set
	 */
	float* Get( size_t iGene ) const {

		return m_Data.Get( iGene ); }

	/*!
	 * \brief
	 * Return the PCL's underlying data matrix.
	 * 
	 * \returns
	 * Data matrix by which the PCL is backed.
	 */
	const CDataMatrix& Get( ) const {

		return m_Data; }

	/*!
	 * \brief
	 * Return the number of gene rows in the PCL.
	 * 
	 * \returns
	 * Number of gene rows in the PCL.
	 * 
	 * \see
	 * GetExperiments
	 */
	size_t GetGenes( ) const {

		return m_vecstrGenes.size( ); }

	/*!
	 * \brief
	 * Returns the vector of gene names associated with this PCL.
	 * 
	 * \returns
	 * Vector of this PCL's gene names.
	 * 
	 * \remarks
	 * Returned vector size will be identical to GetGenes.
	 */
	const std::vector<std::string>& GetGeneNames( ) const {

		return m_vecstrGenes; }

	/*!
	 * \brief
	 * Return the number of experiment columns in the PCL.
	 * 
	 * \returns
	 * Number of experiment columns in the PCL.
	 * 
	 * \see
	 * GetGenes
	 */
	size_t GetExperiments( ) const {

		return m_vecstrExperiments.size( ); }

	/*!
	 * \brief
	 * Returns the gene name at the given PCL row.
	 * 
	 * \param iGene
	 * Index of gene name to return.
	 * 
	 * \returns
	 * Gene name at the requested row.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given index must be smaller than GetGenes.
	 */
	const std::string& GetGene( size_t iGene ) const {

		return m_vecstrGenes[ iGene ]; }

	/*!
	 * \brief
	 * Set the gene name at the given index.
	 * 
	 * \param iGene
	 * Row of gene name to modify.
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

		m_vecstrGenes[ iGene ] = strGene; }

	/*!
	 * \brief
	 * Returns the experiment label at the given PCL column.
	 * 
	 * \param iExperiment
	 * Index of experiment label to return.
	 * 
	 * \returns
	 * Experiment label at the requested column.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given index must be smaller than GetExperiments.
	 */
	const std::string& GetExperiment( size_t iExperiment ) const {

		return m_vecstrExperiments[ iExperiment ]; }

	/*!
	 * \brief
	 * Set the experiment label at the given PCL column.
	 * 
	 * \param iExperiment
	 * Index of experiment label to set.
	 * 
	 * \param strExperiment
	 * Experiment label to which the requested column is set.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given index must be smaller than GetExperiments.
	 */
	void SetExperiment( size_t iExperiment, const std::string& strExperiment ) {

		m_vecstrExperiments[ iExperiment ] = strExperiment; }

	/*!
	 * \brief
	 * Set the mask value for a given gene row.
	 * 
	 * \param iGene
	 * Index of gene row to be masked or unmasked.
	 * 
	 * \param fMask
	 * If true, mask given gene; if false, unmask it.
	 * 
	 * PCLs allow individual gene rows to be masked.  This generally excludes those rows from processing,
	 * particularly from any saved output or from use in any machine learning algorithm consuming a PCL.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given value must be smaller than GetGenes.
	 * 
	 * \see
	 * IsMasked
	 */
	void MaskGene( size_t iGene, bool fMask = true ) {

		if( fMask )
			m_setiGenes.insert( iGene );
		else
			m_setiGenes.erase( iGene ); }

	/*!
	 * \brief
	 * Return true if the requested gene index is masked.
	 * 
	 * \param iGene
	 * Gene row to test for masking.
	 * 
	 * \returns
	 * True if the requested gene row is masked.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given value must be smaller than GetGenes.
	 * 
	 * \see
	 * MaskGene
	 */
	bool IsMasked( size_t iGene ) const {

		return ( m_setiGenes.find( iGene ) != m_setiGenes.end( ) ); }

	/*!
	 * \brief
	 * Set the value at the requested PCL position.
	 * 
	 * \param iGene
	 * PCL row.
	 * 
	 * \param iExperiment
	 * PCL column.
	 * 
	 * \param dValue
	 * Value to store.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given row and column must be smaller than
	 * GetGenes and GetExperiments.
	 * 
	 * \see
	 * Get
	 */
	void Set( size_t iGene, size_t iExperiment, float dValue ) {

		m_Data.Set( iGene, iExperiment, dValue ); }

	/*!
	 * \brief
	 * Return the index of the given gene name, or -1 if it is not included in the PCL.
	 * 
	 * \param strGene
	 * Gene name to retrieve.
	 * 
	 * \returns
	 * Index of the requested gene name, or -1 if it is not in the PCL.
	 * 
	 * \see
	 * GetGeneNames
	 */
	size_t GetGene( const std::string& strGene ) const {
		size_t	i;

		for( i = 0; i < m_vecstrGenes.size( ); ++i )
			if( m_vecstrGenes[ i ] == strGene )
				return i;

		return -1; }

	/*!
	 * \brief
	 * Return the requested gene's value for the given feature.
	 * 
	 * \param iGene
	 * Gene for which feature value should be returned.
	 * 
	 * \param szFeature
	 * Label of feature for which value should be retrieved.
	 * 
	 * \returns
	 * Feature value for the given gene and feature, or an empty string if the feature label is
	 * unrecognized.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given value must be smaller than GetGenes.
	 */
	std::string GetFeature( size_t iGene, const char* szFeature ) const {
		size_t	 i;

		for( i = 0; i < m_vecstrFeatures.size( ); ++i )
			if( m_vecstrFeatures[ i ] == szFeature )
				return GetFeature( iGene, i );

		return ""; }

	/*!
	 * \brief
	 * Randomizes the values in each row of the PCL.
	 * 
	 * Randomly shuffles each gene's vector of values.
	 */
	void Randomize( ) {
		size_t	i;

		for( i = 0; i < GetGenes( ); ++i )
			std::random_shuffle( Get( i ), Get( i ) + GetExperiments( ) ); }
};

}

#endif // PCL_H

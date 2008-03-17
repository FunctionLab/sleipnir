#ifndef PCLSET_H
#define PCLSET_H

#include "pclseti.h"
#include "meta.h"
#include "pcl.h"

namespace Sleipnir {

/*!
 * \brief
 * A PCL set manages a collection of CPCL objects and aligns their gene indices.
 * 
 * A PCL set is a set of PCLs whose gene indices have been aligned; that is, the same gene will appear at
 * any index in all contained PCLs.  Genes not actually present in individual PCLs will be given missing
 * values for data and empty strings for feature values.
 * 
 * \see
 * IDataset
 */
class CPCLSet : CPCLSetImpl {
public:
	bool Open( const std::vector<std::string>& vecstrFiles, size_t iSkip = 2,
		CPCL::ENormalize eNormalize = CPCL::ENormalizeNone );

	/*!
	 * \brief
	 * Return the number of genes in the PCL set.
	 * 
	 * \returns
	 * Number of genes in the PCL set.
	 * 
	 * \remarks
	 * This will be the size of the union of all genes in the contained PCLs.
	 */
	size_t GetGenes( ) const {

		return m_vecstrGenes.size( ); }

	/*!
	 * \brief
	 * Return the number of PCLs in the set.
	 * 
	 * \returns
	 * Number of PCLs in the set.
	 */
	size_t GetPCLs( ) const {

		return m_iPCLs; }

	/*!
	 * \brief
	 * Returns the value at the requested PCL position.
	 * 
	 * \param iPCL
	 * PCL index.
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
	 * For efficiency, no bounds checking is performed.  The given values must be smaller than GetPCLs,
	 * GetGenes, and that PCL's GetExperiments.  Note that the gene index is within the PCL set, not the
	 * individual PCL; the experiment index is within the individual PCL.
	 */
	float Get( size_t iPCL, size_t iGene, size_t iExperiment ) const {
		size_t	iMap;

		if( ( iMap = m_Genes.Get( iPCL, iGene ) ) == -1 )
			return CMeta::GetNaN( );

		return m_aPCLs[ iPCL ].Get( iMap, iExperiment ); }

	/*!
	 * \brief
	 * Return a single gene's row from the given PCL.
	 * 
	 * \param iPCL
	 * PCL index.
	 * 
	 * \param iGene
	 * Gene row.
	 * 
	 * \returns
	 * Requested PCL row.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given values must be smaller than GetPCLs
	 * and GetGenes.  Note that the gene index is within the PCL set, not the individual PCL.
	 */
	const float* Get( size_t iPCL, size_t iGene ) const {
		size_t	iMap;

		if( ( iMap = m_Genes.Get( iPCL, iGene ) ) == -1 )
			return NULL;

		return m_aPCLs[ iPCL ].Get( iMap ); }

	/*!
	 * \brief
	 * Return the index of the given gene name, or -1 if it is not included in the PCL set.
	 * 
	 * \param strGene
	 * Gene name to retrieve.
	 * 
	 * \returns
	 * Index of the requested gene name, or -1 if it is not in the PCL set.
	 * 
	 * \see
	 * GetGeneNames
	 */
	size_t GetGene( const std::string& strGene ) const {
		TMapStrI::const_iterator	iterGene;

		return ( ( ( iterGene = m_mapGenes.find( strGene ) ) == m_mapGenes.end( ) ) ? -1 :
			iterGene->second ); }

	/*!
	 * \brief
	 * Returns the gene name at the given index.
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
	const std::string& GetGene( size_t iGene ) const {

		return m_vecstrGenes[ iGene ]; }

	/*!
	 * \brief
	 * Returns the vector of gene names associated with this PCL set.
	 * 
	 * \returns
	 * Vector of this PCL set's gene names.
	 * 
	 * \remarks
	 * Returned vector size will be identical to GetGenes.
	 */
	const std::vector<std::string>& GetGeneNames( ) const {

		return m_vecstrGenes; }

	/*!
	 * \brief
	 * Return the PCL at the given index.
	 * 
	 * \param iPCL
	 * Index of PCL to retrieve from the set.
	 * 
	 * \returns
	 * PCL at the requested index.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given index must be smaller than GetPCLs.  Note
	 * that no gene alignment is done within the original PCL object.
	 */
	const CPCL& Get( size_t iPCL ) const {

		return m_aPCLs[ iPCL ]; }
};

}

#endif // PCLSET_H

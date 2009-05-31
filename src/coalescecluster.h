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
#ifndef COALESCECLUSTER_H
#define COALESCECLUSTER_H

#include "coalesceclusteri.h"

namespace Sleipnir {

/*!
 * \brief
 * Manages a single converging regulatory module for CCoalesce.
 * 
 * A COALESCE cluster represents a regulatory module consisting of genes, conditions where they're coexpressed,
 * and putative regulatory motifs.  The cluster object is responsible for both the data necessary to
 * represent such a module (read and written using Open and Save) and the methods used to generate the
 * module (primarily Initialize, SelectConditions, SelectMotifs, and SelecGenes).  A variety of optimizations
 * make the algorithm both more efficient and more complex, including discretization of various scores and
 * change tracking between iterations.
 * 
 * \remarks
 * A cluster can be generated during COALESCE execution, but it can also be read and written independently.
 * Clusters generally need to be associated with the CPCL object from which they were first created, and
 * may also interact with a CCoalesceMotifLibrary if they contain predicted motifs.
 * 
 * \see
 * CCoalesce
 */
class CCoalesceCluster : public CCoalesceClusterImpl {
public:
	bool Initialize( const CPCL& PCL, CCoalesceCluster& Pot, const std::vector<SCoalesceDataset>& vecsDatasets,
		std::set<std::pair<size_t, size_t> >& setpriiSeeds, size_t iPairs, float dPValue, size_t iThreads );
	void Subtract( CPCL& PCL, const CCoalesceCluster& Pot ) const;
	void Subtract( CCoalesceGeneScores& GeneScores ) const;
	bool SelectConditions( const CPCL& PCL, const CCoalesceCluster& Pot, size_t iThreads, float dPValue,
		float dZScore );
	bool SelectMotifs( const CCoalesceGroupHistograms& HistsCluster, const CCoalesceGroupHistograms& HistsPot,
		float dPValue, float dZScore, size_t iMaxMotifs, size_t iThreads,
		const CCoalesceMotifLibrary* pMotifs = NULL );
	bool SelectGenes( const CPCL& PCL, const CCoalesceGeneScores& GeneScores,
		const CCoalesceGroupHistograms& HistsCluster, const CCoalesceGroupHistograms& HistsPot,
		size_t iMinimum, size_t iThreads, CCoalesceCluster& Pot, float dProbability,
		const CCoalesceMotifLibrary* pMotifs = NULL );
	void CalculateHistograms( const CCoalesceGeneScores& GeneScores,
		CCoalesceGroupHistograms& HistogramsCluster, CCoalesceGroupHistograms* pHistogramsPot ) const;
	size_t Open( const std::string& strPCL, size_t iSkip, const CPCL& PCL,
		CCoalesceMotifLibrary* pMotifs = NULL );
	size_t Open( std::istream& istm, const CPCL& PCL, CCoalesceMotifLibrary* pMotifs = NULL );
	bool Open( const CHierarchy& Hierarchy, const std::vector<CCoalesceCluster>& vecClusters,
		const std::vector<std::string>& vecstrClusters, float dFraction, float dCutoff, size_t iCutoff,
		CCoalesceMotifLibrary* pMotifs = NULL );
	bool Save( const std::string& strDirectory, size_t iID, const CPCL& PCL,
		const CCoalesceMotifLibrary* pMotifs = NULL ) const;
	void Save( std::ostream& ostm, size_t iID, const CPCL& PCL, const CCoalesceMotifLibrary* pMotifs = NULL,
		float dCutoffPWMs = 0, float dPenaltyGap = 0, float dPenaltyMismatch = 0, bool fNoRCs = false ) const;
	float GetSimilarity( const CCoalesceCluster& Cluster, size_t iGenes, size_t iDatasets ) const;
	void Snapshot( const CCoalesceGeneScores& GeneScores, CCoalesceGroupHistograms& Histograms );
	bool LabelMotifs( const CCoalesceMotifLibrary& Motifs, SMotifMatch::EType eMatchType, float dPenaltyGap,
		float dPenaltyMismatch, float dPValue );

	/*!
	 * \brief
	 * Returns true if the cluster has converged to a previously visited state.
	 * 
	 * \returns
	 * True if the cluster is in a previously visited state.
	 * 
	 * \remarks
	 * States are a combination of genes, conditions, and motifs, which are hashed into a single
	 * integer.  This has the potential to cause false collisions, which will at worst prematurely
	 * terminate cluster refinement.
	 */
	bool IsConverged( ) {

		return ( m_setiHistory.find( GetHash( ) ) != m_setiHistory.end( ) ); }

	/*!
	 * \brief
	 * Returns true if the cluster contains no genes or no conditions.
	 * 
	 * \returns
	 * True if the cluster contains no genes or no conditions.
	 */
	bool IsEmpty( ) const {

		return ( m_setiGenes.empty( ) || m_setiDatasets.empty( ) ); }

	/*!
	 * \brief
	 * Replaces the cluster's genes with the requested number of initial genes.
	 * 
	 * \param iGenes
	 * Number of genes to be included in the cluster.
	 * 
	 * \remarks
	 * Sets the cluster's gene set to 0 through iGenes - 1.
	 */
	void SetGenes( size_t iGenes ) {
		size_t	i;

		m_setiGenes.clear( );
		for( i = 0; i < iGenes; ++i )
			m_setiGenes.insert( i ); }

	/*!
	 * \brief
	 * Returns the set of genes significant in this cluster.
	 * 
	 * \returns
	 * Set of genes included in this cluster.
	 */
	const std::set<size_t>& GetGenes( ) const {

		return CCoalesceClusterImpl::GetGenes( ); }

	/*!
	 * \brief
	 * Returns the set of datasets significant in this cluster.
	 * 
	 * \returns
	 * Set of datasets included in this cluster.
	 */
	const std::set<size_t>& GetDatasets( ) const {

		return m_setiDatasets; }

	/*!
	 * \brief
	 * Returns the set of motifs significant in this cluster.
	 * 
	 * \returns
	 * Set of motifs included in this cluster.
	 */
	const std::set<SMotifMatch>& GetMotifs( ) const {

		return m_setsMotifs; }

	/*!
	 * \brief
	 * Returns true if the given gene is significant in this cluster.
	 * 
	 * \param iGene
	 * Gene index to be tested for inclusion.
	 * 
	 * \returns
	 * True if the given gene is included in this cluster.
	 */
	bool IsGene( size_t iGene ) const {

		return CCoalesceClusterImpl::IsGene( iGene ); }

	/*!
	 * \brief
	 * Returns true if the given dataset is significant in this cluster.
	 * 
	 * \param iDataset
	 * Dataset index to be tested for inclusion.
	 * 
	 * \returns
	 * True if the given dataset is included in this cluster.
	 */
	bool IsDataset( size_t iDataset ) const {

		return ( m_setiDatasets.find( iDataset ) != m_setiDatasets.end( ) ); }

	/*!
	 * \brief
	 * Removes the given genes from the cluster.
	 * 
	 * \param veciGenes
	 * Vector of gene IDs to be removed.
	 */
	void RemoveGenes( const std::vector<size_t>& veciGenes ) {
		size_t	i;

		for( i = 0; i < veciGenes.size( ); ++i )
			m_setiGenes.erase( veciGenes[ i ] ); }
};

}

#endif // COALESCECLUSTER_H

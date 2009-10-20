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
#ifndef COALESCE_H
#define COALESCE_H

#include "coalescei.h"
#include "coalescecluster.h"
#include "coalescemotifs.h"

namespace Sleipnir {

/*!
 * \brief
 * Performs regulatory module prediction (gene expression biclustering plus de novo sequence motif discovery)
 * using the COALESCE algorithm of Huttenhower et al. 2009.
 * 
 * The COALESCE algorithm consumes gene expression data and, optionally, DNA sequences, to predict regulatory
 * modules.  These consist of expression biclusters (subsets of genes and conditions) and putative regulatory
 * motifs.  COALESCE predicts modules in a serial manner, seeding each module with a small number of
 * correlated genes.  It then iterates between feature selection and Bayesian integration of the selected
 * features to determine which genes should be in the module.  Feature selection chooses expression conditions
 * in which the cluster's genes are differentially expressed (i.e. significantly different than the genomic
 * background) and sequence motifs over- or under-enriched in sequences associated with the cluster's genes
 * (also relative to genomic background).  Bayesian integration assumes that these features are independent
 * (although prior knowledge of non-independent datasets can be provided and used to incorporate covariance
 * information) and calculates the probability with which each gene in the genome is included in the developing
 * module.  These two steps (feature selection and Bayesian integration) are iterated until the module has
 * converged, at which point its average values (expression and motif frequencies) are subtracted from its
 * genes' data, and COALESCE continues with the next module.  A variety of options and data can be used to
 * modify this procedure, both at the level of the algorithm itself (e.g. the probability threshhold at which
 * genes are included in a module) and at the level of implementation optimizations (e.g. the granularity
 * with which motif frequencies are discretized).
 * 
 * \remarks
 * CCoalesce is tightly coupled to CCoalesceCluster, where many of the details of the COALESCE algorithm
 * are implemented.  CCoalesce serves mainly to store configuration state information for clustering, to
 * initialize associated data structures, and to provide the outermost skeleton of the algorithm.
 * 
 * \see
 * CCoalesceCluster
 */
class CCoalesce : CCoalesceImpl {
public:
	bool Cluster( const CPCL& PCL, const CFASTA& FASTA, std::vector<CCoalesceCluster>& vecClusters );

	void SetSeed( const CPCL& PCL );

	/*!
	 * \brief
	 * Sets the correlation p-value threshhold for genes to be included in a cluster during initialization.
	 * 
	 * \param dPValue
	 * Correlation p-value threshhold for gene inclusion during module initialization.
	 * 
	 * \see
	 * GetPValueCorrelation
	 */
	void SetPValueCorrelation( float dPValue ) {

		m_dPValueCorrelation = dPValue; }

	/*!
	 * \brief
	 * Returns the correlation p-value threshhold for genes to be included in a cluster during initialization.
	 * 
	 * \returns
	 * P-value threshhold for gene inclusion during module initialization.
	 * 
	 * \see
	 * SetPValueCorrelation
	 */
	float GetPValueCorrelation( ) const {

		return m_dPValueCorrelation; }

	/*!
	 * \brief
	 * Sets the number of discretization bins used for calculating motif frequency histograms.
	 * 
	 * \param iBins
	 * Number of bins used to discretize motif frequencies.
	 * 
	 * \see
	 * GetBins
	 */
	void SetBins( size_t iBins ) {

		m_iBins = iBins; }

	/*!
	 * \brief
	 * Returns the number of discretization bins used for calculating motif frequency histograms.
	 * 
	 * \returns
	 * Number of bins used to discretize motif frequencies.
	 * 
	 * \see
	 * SetBins
	 */
	size_t GetBins( ) const {

		return m_iBins; }

	/*!
	 * \brief
	 * Returns the z-score effect size threshhold for including significant expression conditions in a cluster.
	 * 
	 * \returns
	 * Z-score threshhold for inclusion of expression conditions in a cluster.
	 * 
	 * \see
	 * SetZScoreCondition
	 */
	float GetZScoreCondition( ) const {

		return m_dZScoreCondition; }

	/*!
	 * \brief
	 * Sets the z-score effect size threshhold for including significant expression conditions in a cluster.
	 * 
	 * \param dZScore
	 * Z-score threshhold for inclusion of expression conditions in a cluster.
	 * 
	 * \see
	 * GetZScoreCondition
	 */
	void SetZScoreCondition( float dZScore ) {

		m_dZScoreCondition = dZScore; }

	/*!
	 * \brief
	 * Returns the p-value threshhold for including significant expression conditions in a cluster.
	 * 
	 * \returns
	 * P-value threshhold for inclusion of expression conditions in a cluster.
	 * 
	 * \see
	 * SetPValueCondition
	 */
	float GetPValueCondition( ) const {

		return m_dPValueCondition; }

	/*!
	 * \brief
	 * Sets the p-value threshhold for including significant expression conditions in a cluster.
	 * 
	 * \param dPValue
	 * P-value threshhold for inclusion of expression conditions in a cluster.
	 * 
	 * \see
	 * GetPValueCondition
	 */
	void SetPValueCondition( float dPValue ) {

		m_dPValueCondition = dPValue; }

	/*!
	 * \brief
	 * Returns the z-score effect size threshhold for including significant sequence motifs in a cluster.
	 * 
	 * \returns
	 * Z-score threshhold for inclusion of motifs in a cluster.
	 * 
	 * \see
	 * SetZScoreMotif
	 */
	float GetZScoreMotif( ) const {

		return m_dZScoreMotif; }

	/*!
	 * \brief
	 * Sets the z-score effect size threshhold for including significant sequence motifs in a cluster.
	 * 
	 * \param dZScore
	 * Z-score threshhold for inclusion of motifs in a cluster.
	 * 
	 * \see
	 * GetZScoreMotif
	 */
	void SetZScoreMotif( float dZScore ) {

		m_dZScoreMotif = dZScore; }

	/*!
	 * \brief
	 * Returns the p-value threshhold for including significant sequence motifs in a cluster.
	 * 
	 * \returns
	 * P-value threshhold for inclusion of motifs in a cluster.
	 * 
	 * \see
	 * SetPValueMotif
	 */
	float GetPValueMotif( ) const {

		return m_dPValueMotif; }

	/*!
	 * \brief
	 * Sets the p-value threshhold for including significant sequence motifs in a cluster.
	 * 
	 * \param dPValue
	 * P-value threshhold for inclusion of motifs in a cluster.
	 * 
	 * \see
	 * GetPValueMotif
	 */
	void SetPValueMotif( float dPValue ) {

		m_dPValueMotif = dPValue; }

	/*!
	 * \brief
	 * Returns the probability threshhold for including genes in a cluster.
	 * 
	 * \returns
	 * Probability threshhold for inclusion of genes in a cluster.
	 * 
	 * \see
	 * SetProbabilityGene
	 */
	float GetProbabilityGene( ) const {

		return m_dProbabilityGene; }

	/*!
	 * \brief
	 * Sets the probability threshhold for including genes in a cluster.
	 * 
	 * \param dProbability
	 * Probability threshhold for inclusion of genes in a cluster.
	 * 
	 * \see
	 * GetProbabilityGene
	 */
	void SetProbabilityGene( float dProbability ) {

		m_dProbabilityGene = dProbability; }

	/*!
	 * \brief
	 * Returns true if a module output directory has been set.
	 * 
	 * \returns
	 * True if a module output directory has been set; false if modules are to be output only to standard out.
	 * 
	 * \see
	 * GetDirectoryIntermediate | SetDirectoryIntermediate
	 */
	bool IsDirectoryIntermediate( ) const {

		return !GetDirectoryIntermediate( ).empty( ); }

	/*!
	 * \brief
	 * Returns the output directory for predicted modules.
	 * 
	 * \returns
	 * Output directory in which predicted modules are saved.
	 * 
	 * \see
	 * SetDirectoryIntermediate
	 */
	const std::string& GetDirectoryIntermediate( ) const {

		return m_strDirectoryIntermediate; }

	/*!
	 * \brief
	 * Sets the output directory for predicted modules.
	 * 
	 * \param strDirectoryIntermediate
	 * Output directory in which predicted modules are saved.
	 * 
	 * \see
	 * GetDirectoryIntermediate
	 */
	void SetDirectoryIntermediate( const std::string& strDirectoryIntermediate ) {

		m_strDirectoryIntermediate = strDirectoryIntermediate; }

	/*!
	 * \brief
	 * Sets the motif library used to manage gene sequences and motifs.
	 * 
	 * \param Motifs
	 * Motif library used to manage gene sequences and motifs during clustering.
	 * 
	 * \see
	 * GetMotifs
	 */
	void SetMotifs( CCoalesceMotifLibrary& Motifs ) {

		if( m_fMotifs && m_pMotifs && ( m_pMotifs != &Motifs ) )
			delete m_pMotifs;
		m_pMotifs = &Motifs; }

	/*!
	 * \brief
	 * Returns the motif library used to manage gene sequences and motifs.
	 * 
	 * \returns
	 * Motif library used to manage gene sequences and motifs during clustering; null if none has been set.
	 * 
	 * \see
	 * SetMotifs
	 */
	const CCoalesceMotifLibrary* GetMotifs( ) const {

		return m_pMotifs; }

	/*!
	 * \brief
	 * Returns the length of k-mer motifs.
	 * 
	 * \returns
	 * K-mer length of predicted motifs; also used as building blocks for more complex motifs.
	 * 
	 * \see
	 * SetK
	 */
	size_t GetK( ) const {

		return m_iK; }

	/*!
	 * \brief
	 * Sets the length of k-mer motifs.
	 * 
	 * \param iK
	 * K-mer length of predicted motifs; also used as building blocks for more complex motifs.
	 * 
	 * \see
	 * GetK
	 */
	void SetK( size_t iK ) {

		m_iK = iK; }

	/*!
	 * \brief
	 * Returns the granularity in base pairs with which motif frequency histograms are calculated.
	 * 
	 * \returns
	 * Number of base pairs per match used to calculated motif frequency histograms.
	 * 
	 * \see
	 * SetBasesPerMatch
	 */
	size_t GetBasesPerMatch( ) const {

		return m_iBasesPerMatch; }

	/*!
	 * \brief
	 * Sets the granularity in base pairs with which motif frequency histograms are calculated.
	 * 
	 * \param iBasesPerMatch
	 * Number of base pairs per match used to calculated motif frequency histograms.
	 * 
	 * \remarks
	 * Each bin in a motif frequency histogram will be of width 1 / iBasesPerMatch.
	 * 
	 * \see
	 * GetBasesPerMatch
	 */
	void SetBasesPerMatch( size_t iBasesPerMatch ) {

		m_iBasesPerMatch = iBasesPerMatch; }

	/*!
	 * \brief
	 * Returns the p-value threshhold at which motifs are merged to build PSTs.
	 * 
	 * \returns
	 * P-value threshhold at which motifs are merged to build PSTs.
	 * 
	 * \see
	 * SetPValueMerge
	 */
	float GetPValueMerge( ) const {

		return m_dPValueMerge; }

	/*!
	 * \brief
	 * Sets the p-value threshhold at which motifs are merged to build PSTs.
	 * 
	 * \param dPValue
	 * P-value threshhold at which motifs are merged to build PSTs.
	 * 
	 * \see
	 * GetPValueMerge
	 */
	void SetPValueMerge( float dPValue ) {

		m_dPValueMerge = dPValue; }

	/*!
	 * \brief
	 * Returns the edit distance threshhold at which motifs are merged to build PSTs.
	 * 
	 * \returns
	 * Edit distance threshhold at which motifs are merged to build PSTs.
	 * 
	 * \see
	 * SetCutoffMerge
	 */
	float GetCutoffMerge( ) const {

		return m_dCutoffMerge; }

	/*!
	 * \brief
	 * Sets the edit distance threshhold at which motifs are merged to build PSTs.
	 * 
	 * \param dCutoff
	 * Edit distance threshhold at which motifs are merged to build PSTs.
	 * 
	 * \see
	 * GetCutoffMerge
	 */
	void SetCutoffMerge( float dCutoff ) {

		m_dCutoffMerge = dCutoff; }

	/*!
	 * \brief
	 * Returns the minimum number of genes that must be present in a successful module.
	 * 
	 * \returns
	 * Minimum number of genes present in a successful module.
	 * 
	 * \see
	 * SetSizeMinimum
	 */
	size_t GetSizeMinimum( ) const {

		return m_iSizeMinimum; }

	/*!
	 * \brief
	 * Sets the minimum number of genes that must be present in a successful module.
	 * 
	 * \param iSizeGenes
	 * Minimum number of genes present in a successful module.
	 * 
	 * \see
	 * GetSizeMinimum
	 */
	void SetSizeMinimum( size_t iSizeGenes ) {

		m_iSizeMinimum = iSizeGenes; }

	/*!
	 * \brief
	 * Returns the maximum number of motifs that may be associated with a converging module.
	 * 
	 * \returns
	 * Maximum number of motifs associated with a converging module.
	 * 
	 * \see
	 * SetSizeMaximum
	 */
	size_t GetSizeMaximum( ) const {

		return m_iSizeMaximum; }

	/*!
	 * \brief
	 * Sets the maximum number of motifs that may be associated with a converging module.
	 * 
	 * \param iSizeMotifs
	 * Maximum number of motifs associated with a converging module.
	 * 
	 * \remarks
	 * Additional motifs may be associated with a module during a final pass after convergence.
	 * 
	 * \see
	 * GetSizeMaximum
	 */
	void SetSizeMaximum( size_t iSizeMotifs ) {

		m_iSizeMaximum = iSizeMotifs; }

	/*!
	 * \brief
	 * Returns the maximum number of motifs that are considered for merging into PSTs during module convergence.
	 * 
	 * \returns
	 * Maximum number of motifs considered for PSTs construction during module convergence.
	 * 
	 * \see
	 * SetSizeMerge
	 */
	size_t GetSizeMerge( ) const {

		return m_iSizeMerge; }

	/*!
	 * \brief
	 * Sets the maximum number of motifs that are considered for merging into PSTs during module convergence.
	 * 
	 * \param iSizeMerge
	 * Maximum number of motifs considered for PSTs construction during module convergence.
	 * 
	 * \remarks
	 * Additional motifs may be merged during module postprocessing.
	 * 
	 * \see
	 * GetSizeMerge
	 */
	void SetSizeMerge( size_t iSizeMerge ) {

		m_iSizeMerge = iSizeMerge; }

	/*!
	 * \brief
	 * Removes all currently set dataset blocks.
	 * 
	 * \see
	 * AddDataset
	 */
	void ClearDatasets( ) {

		m_vecsDatasets.clear( ); }

	/*!
	 * \brief
	 * Adds a block of conditions known to form a non-independent dataset.
	 * 
	 * \param setiDataset
	 * Set of condition indices forming a dataset.
	 * 
	 * \returns
	 * True if the dataset was added successfully; false if it already existed in the current configuration.
	 * 
	 * Adds a dataset block to subsequent executions of COALESCE.  A dataset block consists of two or more
	 * expression conditions known to be non-independent, e.g. multiple conditions belonging to the same
	 * time course.  Such dataset blocks are treated as units for inclusion in/exclusion from predicted
	 * modules, and their covariance is determined and incorporated into significance calculations for
	 * differential expression.
	 * 
	 * \remarks
	 * Condition indices must correspond to columns in a PCL file subsequently provided to a call to Cluster.
	 * 
	 * \see
	 * ClearDatasets
	 */
	bool AddDataset( const std::set<size_t>& setiDataset ) {
		size_t								i;
		std::set<size_t>::const_iterator	iterExperiment;

		if( setiDataset.empty( ) )
			return true;
		for( iterExperiment = setiDataset.begin( ); iterExperiment != setiDataset.end( ); ++iterExperiment )
			for( i = 0; i < m_vecsDatasets.size( ); ++i )
				if( m_vecsDatasets[ i ].IsCondition( *iterExperiment ) )
					return false;

		m_vecsDatasets.push_back( SCoalesceDataset( setiDataset ) );
		return true; }

	/*!
	 * \brief
	 * Sets the maximum number of gene pairs subsampled for seed pair discovery during module initialization.
	 * 
	 * \param iPairs
	 * Maximum number of gene pairs subsampled for module seeding.
	 * 
	 * \see
	 * GetNumberCorrelation
	 */
	void SetNumberCorrelation( size_t iPairs ) {

		m_iNumberCorrelation = iPairs; }

	/*!
	 * \brief
	 * Returns the maximum number of gene pairs subsampled for seed pair discovery during module initialization.
	 * 
	 * \returns
	 * Maximum number of gene pairs subsampled for module seeding.
	 * 
	 * \see
	 * SetNumberCorrelation
	 */
	size_t GetNumberCorrelation( ) const {

		return m_iNumberCorrelation; }

	/*!
	 * \brief
	 * Sets the maximum number of simultaneous threads used for clustering.
	 * 
	 * \param iThreads
	 * Maximum number of simultaneous threads used during clustering.
	 * 
	 * \see
	 * GetThreads
	 */
	void SetThreads( size_t iThreads ) {

		m_iThreads = iThreads; }

	/*!
	 * \brief
	 * Returns the maximum number of simultaneous threads used for clustering.
	 * 
	 * \returns
	 * Maximum number of simultaneous threads used during clustering.
	 * 
	 * \see
	 * SetThreads
	 */
	size_t GetThreads( ) const {

		return m_iThreads; }

	/*!
	 * \brief
	 * Adds a wiggle track of supporting data to be used to weight sequence information.
	 * 
	 * \param FASTA
	 * FASTA file containing peudo-wiggle-track formatted per-base weights for gene sequences.
	 * 
	 * Adds a wiggle track of supporting information used to weight gene sequence positions during
	 * COALESCE clustering.  A wiggle track as used by COALESCE is not precisely in the wiggle track
	 * format as defined by the ENCODE project; instead, it is a FASTA file in which sequence base pairs
	 * have been replaced by per-base-pair scores, one floating point value per line.  In COALESCE, one or
	 * more wiggle tracks can be used to weight the individual base pairs used to determine motif
	 * occurrence and frequencies.  Lower weights (down to zero) will downweight the base pairs at those
	 * positions (and thus the effective frequencies of any motifs that occur there), and higher weights
	 * will upweight them.  In the absence of wiggle tracks, the default weight of all base pairs is one.
	 * 
	 * \remarks
	 * Contents of the provided pseudo-wiggle file must align with the FASTA file of gene sequences
	 * provided to subsequent calls to Cluster.
	 * 
	 * \see
	 * ClearWiggles | CFASTA
	 */
	void AddWiggle( const CFASTA& FASTA ) {

		m_vecpWiggles.push_back( &FASTA ); }

	/*!
	 * \brief
	 * Removes all currently active wiggle tracks.
	 * 
	 * \see
	 * AddWiggle
	 */
	void ClearWiggles( ) {

		m_vecpWiggles.clear( ); }

	/*!
	 * \brief
	 * Adds an output stream to which module information is printed after convergence.
	 * 
	 * \param ostm
	 * Output stream to which each module will be printed after it converges.
	 * 
	 * \remarks
	 * Usually a single output stream, standard output, is sufficient; this is provided for useless
	 * convenience.
	 * 
	 * \see
	 * RemoveOutputIntermediate
	 */
	void AddOutputIntermediate( std::ostream& ostm ) {

		m_vecpostm.push_back( &ostm ); }

	/*!
	 * \brief
	 * Removes an output stream to which module information was printed after convergence.
	 * 
	 * \param ostm
	 * Output stream to which modules were to be printed.
	 * 
	 * \remarks
	 * Removal of an output stream not in the current set will be ignored.
	 * 
	 * \see
	 * AddOutputIntermediate
	 */
	void RemoveOutputIntermediate( std::ostream& ostm ) {
		std::vector<std::ostream*>::iterator	iter;

		if( ( iter = std::find( m_vecpostm.begin( ), m_vecpostm.end( ), &ostm ) ) != m_vecpostm.end( ) )
			m_vecpostm.erase( iter ); }

	/*!
	 * \brief
	 * Removes all currently active intermediate output streams.
	 * 
	 * \see
	 * AddOutputIntermediate | RemoveOutputIntermediate
	 */
	void ClearOutputIntermediate( ) {

		m_vecpostm.clear( ); }

	/*!
	 * \brief
	 * Sets the normalization behavior for automatically detected single channel expression conditions.
	 * 
	 * \param fNormalize
	 * If true, single channel conditions are detected and normalized; otherwise, they are left unchanged.
	 * 
	 * \remarks
	 * Single channel normalization is time-consuming and often degrades performance; it should usually
	 * be left disabled.  However, it can find some interesting clusters given the right input data.
	 * 
	 * \see
	 * GetNormalize
	 */
	void SetNormalize( bool fNormalize ) {

		m_fNormalize = fNormalize; }

	/*!
	 * \brief
	 * Returns true if automatic detection and normalization of single channel expression data is enabled.
	 * 
	 * \returns
	 * True if single channel condition detection and normalization is enabled.
	 * 
	 * \see
	 * SetNormalize
	 */
	bool GetNormalize( ) const {

		return m_fNormalize; }

	void ClearSeed( ) {

		m_vecdSeed.clear( ); }
};

}

#endif // COALESCE_H

#ifndef DATASET_H
#define DATASET_H

#include <assert.h>

#include "dataseti.h"

namespace Sleipnir {

class IBayesNet;

/*!
 * \brief
 * An IDataset abstracts a collection of individual datasets, usually CDats, using various continuous and/or
 * discrete encodings.
 * 
 * An IDataset is intended to manage a collection of individual datasets, usually CDats.  This is often
 * used for integration of many datasets in a model such as a Bayes net or SVM, and as such, IDatasets can
 * be used to learn or evaluate these models.  Although most datasets will be backed by discretized CDats
 * with no hidden data (e.g. CDatasetCompact), the IDataset interface allows:
 * - Alignment of gene indices across all underlying datasets.
 * - Heterogeneous mixtures of discrete and continuous datasets for SVMs or Bayes nets with continuous
 * nodes.
 * - "Hidden" datasets corresponding to unobserved nodes in Bayes nets.
 * Data is organized into "nodes" inside a dataset; usually, each node corresponds to one data file, but
 * certain nodes may be hidden if the dataset is constructed from a model allowing such nodes (e.g. a Bayes
 * net with unobserved nodes).
 * 
 * The IDataset interface merges the gene lists from all contained data files into a single gene list,
 * which it exposes through GetGenes/GetGene/GetGeneNames/etc.  Gene indices are similarly normalized;
 * requesting gene pair i,j will "mean" the same thing in each encapsulated dataset.  Missing values will
 * be filled in as necessary for data files not containing information for the requested pair.  QUANT files
 * associated with non-continuous data files will be loaded automatically.
 * 
 * \remarks
 * The IDataset interface is something of a mess, since it evolved over time from something meant to match
 * data files with Bayes nets to something meant to generically load lots of data.  You're usually best off
 * using CDatasetCompact and/or CDataFilter directly.
 * 
 * \see
 * CDat | CPCLSet | CSVM | IBayesNet
 */
class IDataset {
public:
	/*!
	 * \brief
	 * Returns true if the requested experimental node is hidden (does not correspond to a data file).
	 * 
	 * \param iNode
	 * Experimental node to investigate.
	 * 
	 * \returns
	 * True if the requested experimental node is hidden.
	 * 
	 * Since a dataset can be constructed either directly on a collection of data files or by tying a model
	 * such as a Bayes net to data files, IDataset can determine which model nodes are hidden by testing
	 * whether a data file exists for them.  If no such file exists, the node is hidden and, for example,
	 * can be treated specially during Bayesian learning.
	 * 
	 * \remarks
	 * Datasets constructed directly from data files will never have hidden nodes.
	 */
	virtual bool IsHidden( size_t iNode ) const = 0;
	/*!
	 * \brief
	 * Return the discretized value at the requested position.
	 * 
	 * \param iY
	 * Data row.
	 * 
	 * \param iX
	 * Data column.
	 * 
	 * \param iNode
	 * Experimental node from which to retrieve the requested pair's value.
	 * 
	 * \returns
	 * Discretized value from the requested position and data file using that file's quantization information;
	 * -1 if the value is missing.
	 * 
	 * \remarks
	 * Equivalent to using CDataPair::Quantize and GetContinuous or CDataPair::Get on the encapsulated data
	 * file with the appropriate indices.  Behavior not defined when no discretization information is
	 * available for the requested data node.
	 */
	virtual size_t GetDiscrete( size_t iY, size_t iX, size_t iNode ) const = 0;
	/*!
	 * \brief
	 * Return the continuous value at the requested position.
	 * 
	 * \param iY
	 * Data row.
	 * 
	 * \param iX
	 * Data column.
	 * 
	 * \param iNode
	 * Experimental node from which to retrieve the requested pair's value.
	 * 
	 * \returns
	 * Continuous value from the requested position and data file; not-a-number (NaN) if the value is missing.
	 * 
	 * \remarks
	 * Equivalent to using CDataPair::Get on the encapsulated data file with the appropriate indices.
	 * Behavior not defined when the corresponding data node is inherently discrete.
	 * 
	 * \see
	 * GetDiscrete
	 */
	virtual float GetContinuous( size_t iY, size_t iX, size_t iNode ) const = 0;
	/*!
	 * \brief
	 * Returns the gene name at the requested index.
	 * 
	 * \param iGene
	 * Index of gene name to return.
	 * 
	 * \returns
	 * Gene name at the requested index.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.
	 * 
	 * \see
	 * GetGenes
	 */
	virtual const std::string& GetGene( size_t iGene ) const = 0;
	/*!
	 * \brief
	 * Returns the number of genes in the dataset.
	 * 
	 * \returns
	 * Number of genes in the dataset.
	 * 
	 * \remarks
	 * Equal to the union of all genes in encapsulated data files.
	 * 
	 * \see
	 * GetGene
	 */
	virtual size_t GetGenes( ) const = 0;
	/*!
	 * \brief
	 * Returns true if some data file can be accessed at the requested position.
	 * 
	 * \param iY
	 * Data row.
	 * 
	 * \param iX
	 * Data column.
	 * 
	 * \returns
	 * True if a data file can be accessed at the requested position.
	 * 
	 * A dataset position is a usable example if at least one data file can be accessed at that position;
	 * that is, if some data file provides a non-missing value for that gene pair.  Implementations that
	 * filter pairs in some manner can also prevent particular positions from being usable examples.
	 */
	virtual bool IsExample( size_t iY, size_t iX ) const = 0;
	/*!
	 * \brief
	 * Return a vector of all gene names in the dataset.
	 * 
	 * \returns
	 * Vector of gene names in the dataset.
	 * 
	 * \see
	 * GetGenes | GetGene
	 */
	virtual const std::vector<std::string>& GetGeneNames( ) const = 0;
	/*!
	 * \brief
	 * Return the number of experimental nodes in the dataset.
	 * 
	 * \returns
	 * Number of experimental nodes in the dataset.
	 * 
	 * \remarks
	 * For most datasets (those not containing hidden nodes), this will be equal to the number of encapsulated
	 * data files.
	 */
	virtual size_t GetExperiments( ) const = 0;
	/*!
	 * \brief
	 * Return the index of the given gene name, or -1 if it is not included in the dataset.
	 * 
	 * \param strGene
	 * Gene name to retrieve.
	 * 
	 * \returns
	 * Index of the requested gene name, or -1 if it is not in the dataset.
	 * 
	 * \see
	 * GetGeneNames
	 */
	virtual size_t GetGene( const std::string& strGene ) const = 0;
	/*!
	 * \brief
	 * Return the number of discrete values in the requested experimental node; -1 if the node is hidden
	 * or continuous.
	 * 
	 * \param iNode
	 * Experimental node for which bin number should be returned.
	 * 
	 * \returns
	 * Number of discrete values taken by the given experimental node; -1 if the node is hidden or continuous.
	 * 
	 * \see
	 * GetDiscrete
	 */
	virtual size_t GetBins( size_t iNode ) const = 0;
	/*!
	 * \brief
	 * Remove all data for the given dataset position.
	 * 
	 * \param iY
	 * Data row.
	 * 
	 * \param iX
	 * Data column.
	 * 
	 * Unloads or masks data from all encapsulated files for the requested gene pair.
	 * 
	 * \remarks
	 * For efficiency, bounds checking is not performed; the given row and column should both be less than
	 * GetGenes.  Not supported by all implementations.
	 */
	virtual void Remove( size_t iY, size_t iX ) = 0;
	/*!
	 * \brief
	 * Remove values from the dataset based on the given gene set and filter type.
	 * 
	 * \param Genes
	 * Gene set used to filter the dataset.
	 * 
	 * \param eFilter
	 * Way in which to use the given genes to remove values.
	 * 
	 * Remove values and genes (by removing all incident edges) from the dataset based on one of several
	 * algorithms.  For details, see CDat::EFilter.
	 * 
	 * \remarks
	 * Generally implemented using Remove, so may not be supported by all implementations and may either
	 * mask or unload the filtered data.
	 * 
	 * \see
	 * CDat::FilterGenes
	 */
	virtual void FilterGenes( const CGenes& Genes, CDat::EFilter eFilter ) = 0;
};

/*!
 * \brief
 * A simple implementation of IDataset directly loading unmodified CDats for each non-hidden data node.
 * 
 * \remarks
 * For any purpose not requiring continuous values, CDatasetCompact is more appropriate.  CDistanceMatrix
 * objects are used to store continuous data, CCompactMatrix objects for discrete data.
 */
class CDataset : CDatasetImpl, public IDataset {
public:
	bool Open( const char* szAnswerFile, const std::vector<std::string>& vecstrDataFiles );
	bool Open( const std::vector<std::string>& vecstrDataFiles );
	bool Open( const char* szAnswerFile, const char* szDataDirectory, const IBayesNet* pBayesNet );

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
	 * Bayes net whose nodes will correspond to files in the dataset.
	 * 
	 * \returns
	 * True if dataset was constructed successfully.
	 * 
	 * Creates a dataset with nodes corresponding to the given Bayes net structure; the given answer file
	 * is always inserted as the first (0th) data file, and thus corresponds to the first node in the Bayes
	 * net (generally the class node predicting functional relationships).  Data is loaded continuously or
	 * discretely as indicated by the Bayes net, and nodes for which a corresponding data file (i.e. one
	 * with the same name followed by an appropriate CDat extension) cannot be located are marked as hidden.
	 * 
	 * \remarks
	 * Each data file is loaded more or less as-is; continuous data files will be loaded directly into
	 * memory, and discrete files are pre-discretized and stored in compact matrices.
	 */
	bool Open( const CDataPair& Answers, const char* szDataDirectory, const IBayesNet* pBayesNet ) {

		return CDatasetImpl::Open( &Answers, szDataDirectory, pBayesNet ); }

	/*!
	 * \brief
	 * Construct a dataset corresponding to the given Bayes net using files from the given directory.
	 * 
	 * \param szDataDirectory
	 * Directory from which data files are loaded.
	 * 
	 * \param pBayesNet
	 * Bayes net whose nodes will correspond to files in the dataset.
	 * 
	 * \returns
	 * True if dataset was constructed successfully.
	 * 
	 * Creates a dataset (without an answer file) with nodes corresponding to the given Bayes net structure.
	 * Data is loaded continuously or discretely as indicated by the Bayes net, and nodes for which a
	 * corresponding data file (i.e. one with the same name followed by an appropriate CDat extension)
	 * cannot be located are marked as hidden.
	 * 
	 * \remarks
	 * Each data file is loaded more or less as-is; continuous data files will be loaded directly into
	 * memory, and discrete files are pre-discretized and stored in compact matrices.
	 */
	bool Open( const char* szDataDirectory, const IBayesNet* pBayesNet ) {

		return CDatasetImpl::Open( NULL, szDataDirectory, pBayesNet ); }

	/*!
	 * \brief
	 * Open only the merged gene list from the given data files.
	 * 
	 * \param vecstrDataFiles
	 * Vector of file paths to load.
	 * 
	 * \returns
	 * True if gene lists were loaded successfully.
	 * 
	 * Provides a way to rapidly list the set of all genes present in a given collection of data files
	 * while avoiding the overhead of loading the data itself.
	 * 
	 * \remarks
	 * Attempting to access data in the dataset without also opening the files themselves won't do anything
	 * good.
	 * 
	 * \see
	 * CDat::OpenGenes
	 */
	bool OpenGenes( const std::vector<std::string>& vecstrDataFiles ) {

		return CDataImpl::OpenGenes( vecstrDataFiles ); }

	size_t GetDiscrete( size_t iY, size_t iX, size_t iNode ) const;
	float GetContinuous( size_t iY, size_t iX, size_t iNode ) const;
	bool IsExample( size_t iY, size_t iX ) const;
	void Remove( size_t iY, size_t iX );

	void FilterGenes( const CGenes& Genes, CDat::EFilter eFilter ) {

		CDataImpl::FilterGenes( this, Genes, eFilter ); }

	const std::vector<std::string>& GetGeneNames( ) const {

		return CDataImpl::GetGeneNames( ); }

	bool IsHidden( size_t iNode ) const {

		return CDataImpl::IsHidden( iNode ); }

	const std::string& GetGene( size_t iGene ) const {

		return CDataImpl::GetGene( iGene ); }

	size_t GetGenes( ) const {

		return CDataImpl::GetGenes( ); }

	size_t GetExperiments( ) const {

		return CDataImpl::GetExperiments( ); }

	size_t GetGene( const std::string& strGene ) const {

		return CDataImpl::GetGene( strGene ); }

	size_t GetBins( size_t iNode ) const {

		return CDataImpl::GetBins( iNode ); }
};

/*!
 * \brief
 * An implementation of IDataset optimized for compactly storying discrete data.
 * 
 * A compact dataset represents a collection of pre-discretized data files in a compact form.  This can be
 * stored independently of the original continuous data (usually DAB/QUANT file pairs) for rapid reloading
 * (or memory mapping) without the overhead of repeated discretization.  Such a file is referred to as a
 * DAD file and can be stored in either binary or text (human readable) form.  As text, it is a large
 * tab-delimited table of the form:
 * \code
 * GENE1	GENE2	VALUE1	VALUE2	...	VALUEN
 * GENE1	GENE3	VALUE1	VALUE2	...	VALUEN
 * GENE2	GENE3	VALUE1	VALUE2	...	VALUEN
 * \endcode
 * Like a DAT file, gene pair order is arbitrary, and duplicate gene pairs are not recommended.  Missing
 * values are indicated by blank cells, and all other values should be small integers (i.e. discretized
 * values).
 * 
 * \remarks
 * Stores all loaded data in CCompactMatrix objects.  Attempts to use with continuous Bayes nets or with
 * QUANTless CDats will fail in one way or another.
 * 
 * \see
 * CDat
 */
class CDatasetCompact : protected CDatasetCompactImpl, public IDataset {
public:
	bool Open( const CDataPair& Answers, const char* szDataDirectory, const IBayesNet* pBayesNet,
		bool fEverything = false );
	bool Open( const CDataPair& Answers, const char* szDataDirectory, const IBayesNet* pBayesNet,
		const CGenes& GenesInclude, const CGenes& GenesExclude, bool fEverything = false );
	bool Open( const std::vector<std::string>& vecstrDataFiles, bool fMemmap = false );
	bool Open( std::istream& istm );
	bool Open( const CGenes& GenesInclude, const CGenes& GenesExclude, const CDataPair& Answers,
		const std::vector<std::string>& vecstrPCLs, size_t iSkip, const IMeasure* pMeasure,
		const std::vector<float>& vecdBinEdges );
	bool Open( const CDataPair& Answers, const std::vector<std::string>& vecstrDataFiles,
		bool fEverything = false, bool fMemmap = false, size_t iSkip = 2, bool fZScore = false );
	bool FilterGenes( const char* szGenes, CDat::EFilter eFilter );
	void FilterAnswers( );
	void Randomize( );

	/*!
	 * \brief
	 * Construct a dataset corresponding to the given Bayes net using files from the given directory.
	 * 
	 * \param szDataDirectory
	 * Directory from which data files are loaded.
	 * 
	 * \param pBayesNet
	 * Bayes net whose nodes will correspond to files in the dataset.
	 * 
	 * \returns
	 * True if dataset was constructed successfully.
	 * 
	 * Creates a dataset (without an answer file) with nodes corresponding to the given Bayes net structure.
	 * Nodes for which a corresponding data file (i.e. one with the same name followed by an appropriate
	 * CDat extension) cannot be located are marked as hidden.
	 * 
	 * \remarks
	 * Missing QUANT files or requests for continuous data from the Bayes net will result in an error.
	 */
	bool Open( const char* szDataDirectory, const IBayesNet* pBayesNet ) {

		return CDatasetCompactImpl::Open( szDataDirectory, pBayesNet ); }

	/*!
	 * \brief
	 * Construct a dataset corresponding to the given Bayes net using files from the given directory.
	 * 
	 * \param szDataDirectory
	 * Directory from which data files are loaded.
	 * 
	 * \param pBayesNet
	 * Bayes net whose nodes will correspond to files in the dataset.
	 * 
	 * \param GenesInclude
	 * Data is filtered using FilterGenes with CDat::EFilterInclude and the given gene set (unless empty).
	 * 
	 * \param GenesExclude
	 * Data is filtered using FilterGenes with CDat::EFilterExclude and the given gene set (unless empty).
	 * 
	 * \returns
	 * True if dataset was constructed successfully.
	 * 
	 * Creates a dataset (without an answer file) with nodes corresponding to the given Bayes net structure.
	 * Nodes for which a corresponding data file (i.e. one with the same name followed by an appropriate
	 * CDat extension) cannot be located are marked as hidden.
	 * 
	 * \remarks
	 * Missing QUANT files or requests for continuous data from the Bayes net will result in an error.
	 */
	bool Open( const char* szDataDirectory, const IBayesNet* pBayesNet, const CGenes& GenesInclude,
		const CGenes& GenesExclude ) {

		if( !CDatasetCompactImpl::Open( szDataDirectory, pBayesNet, &GenesInclude, &GenesExclude ) )
			return false;
		CDataImpl::FilterGenes( this, GenesInclude, CDat::EFilterInclude );
		CDataImpl::FilterGenes( this, GenesExclude, CDat::EFilterExclude );

		return true; }

	/*!
	 * \brief
	 * Open only the merged gene list from the given data files.
	 * 
	 * \param vecstrDataFiles
	 * Vector of file paths to load.
	 * 
	 * \returns
	 * True if gene lists were loaded successfully.
	 * 
	 * Provides a way to rapidly list the set of all genes present in a given collection of data files
	 * while avoiding the overhead of loading the data itself.
	 * 
	 * \remarks
	 * Attempting to access data in the dataset without also opening the files themselves won't do anything
	 * good.
	 * 
	 * \see
	 * CDat::OpenGenes
	 */
	bool OpenGenes( const std::vector<std::string>& vecstrDataFiles ) {

		return CDataImpl::OpenGenes( vecstrDataFiles ); }

	/*!
	 * \brief
	 * Save a dataset to the given stream in binary or tabular (human readable) form.
	 * 
	 * \param ostm
	 * Stream into which dataset is saved.
	 * 
	 * \param fBinary
	 * If true, save the dataset as a binary DAD; if false, save it as a text-based tab-delimited DAD.
	 * 
	 * \remarks
	 * If fBinary is true, output stream must be binary.
	 */
	void Save( std::ostream& ostm, bool fBinary ) const {

		fBinary ? SaveBinary( ostm ) : SaveText( ostm ); }

	float GetContinuous( size_t iY, size_t iX, size_t iNode ) const {
		UNUSED_PARAMETER(iY);
		UNUSED_PARAMETER(iX);
		UNUSED_PARAMETER(iNode);

		return CMeta::GetNaN( ); }

	const std::string& GetGene( size_t iGene ) const {

		return CDataImpl::GetGene( iGene ); }

	size_t GetGenes( ) const {

		return CDataImpl::GetGenes( ); }

	bool IsExample( size_t iY, size_t iX ) const {

		return CDatasetCompactImpl::IsExample( iY, iX ); }

	void FilterGenes( const CGenes& Genes, CDat::EFilter eFilter ) {

		CDataImpl::FilterGenes( this, Genes, eFilter ); }

	bool IsHidden( size_t iNode ) const {

		return CDataImpl::IsHidden( iNode ); }

	size_t GetDiscrete( size_t iY, size_t iX, size_t iNode ) const {

		return CDatasetCompactImpl::GetDiscrete( iY, iX, iNode ); }

	const std::vector<std::string>& GetGeneNames( ) const {

		return CDataImpl::GetGeneNames( ); }

	size_t GetExperiments( ) const {

		return CDataImpl::GetExperiments( ); }

	size_t GetGene( const std::string& strGene ) const {

		return CDataImpl::GetGene( strGene ); }

	size_t GetBins( size_t iNode ) const {

		return CDataImpl::GetBins( iNode ); }

	void Remove( size_t iY, size_t iX ) {

		CDatasetCompactImpl::Remove( iY, iX ); }
};

/*!
 * \brief
 * Augments a compact dataset with a mask to dynamically exclude specific gene pairs.
 * 
 * A compact dataset mask memory maps a binary DAD.  This saves memory but prevents the underlying
 * dataset from being directly modified.  To overcome this, the map adds a binary matrix allowing each
 * gene pair to be individually masked; a masked gene pair will return false from IsExample and act like
 * missing data.  Unmasked gene pairs will be retrieved from the underlying dataset.  Can be used in
 * combination with CDataFilter.
 * 
 * \see
 * CDataMask
 */
class CDatasetCompactMap : public CDatasetCompact {
public:
	CDatasetCompactMap( );
	~CDatasetCompactMap( );

	bool Open( const char* szFile );

	void Remove( size_t iY, size_t iX ) {

		m_Mask.Set( iY, iX, false ); }

	bool IsExample( size_t iY, size_t iX ) const {

		return m_Mask.Get( iY, iX ); }

private:
	unsigned char*	m_pbData;
	CBinaryMatrix	m_Mask;
	size_t			m_iData;
	HANDLE			m_hndlMap;
};

/*!
 * \brief
 * Augments a dataset with a mask to dynamically exclude specific gene pairs.
 * 
 * A data mask wraps an underlying dataset with a binary matrix allowing each gene pair to be individually
 * masked; a masked gene pair will return false from IsExample and act like missing data.  Unmasked gene
 * pairs will be retrieved from the underlying dataset.  This allows data to be temporarily hidden without
 * modifying the underlying dataset.  Can be used in combination with CDataFilter.
 * 
 * \see
 * CDatasetCompactMap
 */
class CDataMask : CDataMaskImpl, public IDataset {
public:
	void Attach( const IDataset* pDataset );
	void AttachRandom( const IDataset* pDataset, float dFraction );
	void AttachComplement( const CDataMask& DataMask );

	bool IsExample( size_t iY, size_t iX ) const {

		return m_Mask.Get( iY, iX ); }

	void Remove( size_t iY, size_t iX ) {

		m_Mask.Set( iY, iX, false ); }

	const std::vector<std::string>& GetGeneNames( ) const {

		return CDataOverlayImpl::GetGeneNames( ); }

	size_t GetExperiments( ) const {

		return CDataOverlayImpl::GetExperiments( ); }

	size_t GetGene( const std::string& strGene ) const {

		return CDataOverlayImpl::GetGene( strGene ); }

	size_t GetBins( size_t iNode ) const {

		return CDataOverlayImpl::GetBins( iNode ); }

	size_t GetGenes( ) const {

		return CDataOverlayImpl::GetGenes( ); }

	bool IsHidden( size_t iNode ) const {

		return CDataOverlayImpl::IsHidden( iNode ); }

	size_t GetDiscrete( size_t iY, size_t iX, size_t iNode ) const {

		return CDataOverlayImpl::GetDiscrete( iY, iX, iNode ); }

	float GetContinuous( size_t iY, size_t iX, size_t iNode ) const {

		return CDataOverlayImpl::GetContinuous( iY, iX, iNode ); }

	const std::string& GetGene( size_t iGene ) const {

		return CDataOverlayImpl::GetGene( iGene ); }

	void FilterGenes( const CGenes& Genes, CDat::EFilter eFilter ) {

		CDataImpl::FilterGenes( this, Genes, eFilter ); }
};

/*!
 * \brief
 * Augments a dataset with a dynamically calculated gene set filter.
 * 
 * A data filter wraps an underlying dataset with a dynamically calculated filter using a gene set and
 * CDat::EFilter type.  A filtered gene pair will return false from IsExample and act like missing data.
 * Unfiltered gene pairs will be retrieved from the underlying dataset.  This allows data to be temporarily
 * hidden without modifying the underlying dataset.
 * 
 * \remarks
 * Remove (and thus FilterGenes) are unimplemented for this class, as their place is essentially taken
 * by its Attach function.  Permanent modifications such as these should instead be performed on the
 * underlying dataset.
 * 
 * \see
 * CDat::FilterGenes | CDataMask
 */
class CDataFilter : CDataFilterImpl, public IDataset {
public:
	void Attach( const IDataset* pDataset, const CGenes& Genes, CDat::EFilter eFilter,
		const CDat* pAnswers = NULL );

	bool IsExample( size_t iY, size_t iX ) const;

	void Remove( size_t iY, size_t iX ) {

		assert( !"Unimplemented" ); }

	const std::vector<std::string>& GetGeneNames( ) const {

		return CDataOverlayImpl::GetGeneNames( ); }

	size_t GetExperiments( ) const {

		return CDataOverlayImpl::GetExperiments( ); }

	size_t GetGene( const std::string& strGene ) const {

		return CDataOverlayImpl::GetGene( strGene ); }

	size_t GetBins( size_t iNode ) const {

		return CDataOverlayImpl::GetBins( iNode ); }

	size_t GetGenes( ) const {

		return CDataOverlayImpl::GetGenes( ); }

	bool IsHidden( size_t iNode ) const {

		return CDataOverlayImpl::IsHidden( iNode ); }

	size_t GetDiscrete( size_t iY, size_t iX, size_t iNode ) const {

		return ( IsExample( iY, iX ) ? CDataOverlayImpl::GetDiscrete( iY, iX, iNode ) : -1 ); }

	float GetContinuous( size_t iY, size_t iX, size_t iNode ) const {

		return ( IsExample( iY, iX ) ? CDataOverlayImpl::GetContinuous( iY, iX, iNode ) :
			CMeta::GetNaN( ) ); }

	const std::string& GetGene( size_t iGene ) const {

		return CDataOverlayImpl::GetGene( iGene ); }

	void FilterGenes( const CGenes& Genes, CDat::EFilter eFilter ) {

		CDataImpl::FilterGenes( this, Genes, eFilter ); }
};

/*!
 * \brief
 * An IDataset implementation that loads subsets of continuous data in serial to conserve memory.
 * 
 * A data subset loads data continuously, like CDataset; this can be very expensive in memory usage for
 * large datasets.  To make this more tractable, CDataSubset allows subsets of genes of a specific size
 * to be loaded in serial.  A data subset is initialized with a specific size; each time a subset is
 * opened, all data for this many genes (the size) will be loaded, starting with the gene at the
 * requested offset.  This allows the user to first load and use genes 0 through S-1, then S through 2S-1,
 * and so forth.  Typical usage resembles:
 * \code
 * CDataSubset Data;
 * iGeneSize = 1000;
 * Data.Initialize( "data/dir/path/", pBayesNet, iGeneSize );
 * for( iGeneOffset = 0; iGeneOffset < iTotalGeneCount; iGeneOffset += iGeneSize ) {
 *   Data.Open( iGeneOffset );
 *   // do stuff with Data
 * }
 * \endcode
 */
class CDataSubset : CDataSubsetImpl, public IDataset {
public:
	bool Initialize( const char* szDataDirectory, const IBayesNet* pBayesNet, size_t iGeneSize );
	bool Initialize( const std::vector<std::string>& vecstrDataFiles, size_t iGeneSize );
	bool Open( size_t iGeneOffset );

	bool IsHidden( size_t iNode ) const {

		return CDataImpl::IsHidden( iNode ); }

	size_t GetDiscrete( size_t iY, size_t iX, size_t iNode ) const {
		size_t	iMap;

		return ( ( ( iMap = m_veciMapping[ iNode ] ) == -1 ) ? -1 :
			m_Examples.Get( iY - m_iOffset, iX ).GetDiscrete( iMap ) ); }

	float GetContinuous( size_t iY, size_t iX, size_t iNode ) const {
		size_t	iMap;

		return ( ( ( iMap = m_veciMapping[ iNode ] ) == -1 ) ? CMeta::GetNaN( ) :
			m_Examples.Get( iY - m_iOffset, iX ).GetContinuous( iMap ) ); }

	const std::string& GetGene( size_t iGene ) const {

		return CDataImpl::GetGene( iGene ); }

	size_t GetGenes( ) const {

		return CDataImpl::GetGenes( ); }

	bool IsExample( size_t iY, size_t iX ) const {

		return ( ( ( iY < m_iOffset ) || ( ( iY - m_iOffset ) >= m_iSize ) ) ? false :
			m_Examples.Get( iY - m_iOffset, iX ).IsSet( ) ); }

	const std::vector<std::string>& GetGeneNames( ) const {

		return CDataImpl::GetGeneNames( ); }

	size_t GetExperiments( ) const {

		return CDataImpl::GetExperiments( ); }

	size_t GetGene( const std::string& strGene ) const {

		return CDataImpl::GetGene( strGene ); }

	size_t GetBins( size_t iNode ) const {

		return CDataImpl::GetBins( iNode ); }

	void Remove( size_t iY, size_t iX ) {

		m_Examples.Get( iY - m_iOffset, iX ).Reset( ); }

	void FilterGenes( const CGenes& Genes, CDat::EFilter eFilter ) {

		CDataImpl::FilterGenes( this, Genes, eFilter ); }
};

}

#endif // DATASET_H

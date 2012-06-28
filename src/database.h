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
#ifndef DATABASE_H
#define DATABASE_H

#include "databasei.h"

namespace Sleipnir {

class IBayesNet;

/*!
 * \brief
 * Encapsulates a simple indexless database allowing rapid per-gene access to values from many datasets.
 * 
 * CDatabase stores essentially the same data as IDataset; that is, a list of genes and, for each gene pair,
 * zero or more discrete values drawn from many datasets.  IDataset exposes this information so that values
 * for many genes can be drawn rapidly from a single dataset.  CDatabase exposes the same information so
 * that values for one gene can be drawn rapidly from many datasets.  A database can be constructed from
 * a Bayes net and its accompanying data files (DABs/QUANTs), although this amounts to something like a
 * very large matrix transposition, so care must be taken with blocking for large sets of data.  CDatabase
 * consumes very little memory, as data is organized on disk and read efficiently on an as-needed basis.
 * 
 * \remarks
 * Data is stored in discretized form (and thus drawn from DAB/QUANT pairs) using only four bits per gene
 * pair per dataset.  This of course means that no dataset can be discretized into more than 15 bins (one
 * value is reserved to indicate missing data).  Data is spread across an arbitrary number of database
 * subfiles, each containing all data associated with one or more genes.  For G genes spread across N
 * subfiles, the nth subfile will contain the G/N genes n, n + G/N, n + 2G/N, and so forth.  Each gene has
 * one value associated with it per gene per dataset, and these are stored in row-major order: gene g's
 * values are stored first for dataset 0, genes 0 through G, then dataset 1, genes 0 through G, and so forth.
 * The layout of each subfile on disk is thus:
 * \code
 * 4 byte unsigned int, size of header (non-data)
 * 4 byte unsigned int G, total number of genes
 * 4 byte unsigned int N, number of datasets
 * 4 byte unsigned int X = ~G/N, number of genes in subfile
 * X null-terminated ASCII strings storing gene names
 * X rows for genes g of the form:
 *   4-bit data values for pairs g,0 through g,G in dataset 0
 *   4-bit data values for pairs g,0 through g,G in dataset 1
 *   ...
 *   4-bit data values for pairs g,0 through g,G in dataset N
 * \endcode
 * It is possible for a CDatabase to be larger than the input DABs, even though it reduces 32-bit floating
 * point values to 4-bit discrete values, since A) each data point is stored twice, once for pair g,h and
 * once for pair h,g, and B) a value must be stored for every gene pair in every dataset, even if it's
 * missing.  These two facts guarantee very rapid database lookups from disk, but they can in the worst
 * case roughly double the size of the data on disk relative to input DABs.
 * 
 * \see
 * IBayesNet | CDat | CDataPair
 */
class CDatabase : CDatabaseImpl {
public:

	/*!
	 * \brief
	 * Construct a new database over the given genes from the given datasets and Bayes net.
	 * 
	 * \param vecstrGenes
	 * Gene names (and size) for the new database.
	 * 
	 * \param strInputDirectory
	 * Directory containing DAB datasets to be collected into the new database.
	 * 
	 * \param pBayesNet
	 * Bayes net with which the resulting database will be used (indicates the order in which the datasets are
	 * stored in the database).
	 * 
	 * \param strOutputDirectory
	 * Directory in which the new database files are generated.
	 * 
	 * \param iFiles
	 * Number of files across which the database should be constructed.
	 * 
	 * \returns
	 * True if the database was constructed successfully.
	 * 
	 * Constructs a new database over the indicated G genes by creating the requested number of new files in the
	 * output directory and, for each gene g, collecting its data in row major order: values for pairs g,0 through
	 * g,G in dataset 0, then dataset 1, and so forth.  Datasets are loaded from the indicated directory and
	 * ordered as specified in the given Bayes net; this allows data to be loaded very rapidly from the resulting
	 * database and immediately inserted into the Bayes net (or one with identical structure) for inference.
	 * 
	 * Note that for large collections of data, this can be an extremely time- and memory-intensive process.
	 * This is due to the fact that, to construct a database row for some gene, every single dataset must be
	 * inspected to collect (and discretize) all values for all pairs including that gene.  For organisms with
	 * large genomes, it's impossible to load more than a few datasets into memory simultaneously - but reloading
	 * every dataset for every gene wastes far too much time loading and reloading data from DAB files, even if
	 * they're memory mapped.
	 * 
	 * The solution is to use blocking: for gene block size b and dataset block size B, load the first B datasets.
	 * Then copy out values for the first b genes.  Then load the next B datasets, copy b genes' values, and so
	 * forth.  Then repeat the whole process for the next b genes, and so on until the entire given gene set has
	 * been covered.  Correctly balancing the input block size (number of datasets to load) and output block size
	 * (number of database subfiles, and thus genes, to produce at once) can reduce the generation time of large
	 * databases from months to hours.
	 * 
	 * \remarks
	 * iFiles should not be substantially larger than 1000 due to filesystem and file handle limits.
	 * 
	 * \see
	 * SetBlockIn | SetBlockOut
	 */
	bool Open( const std::vector<std::string>& vecstrGenes, const std::string& strInputDirectory,
		const IBayesNet* pBayesNet, const std::string& strOutputDirectory, size_t iFiles);


	//Qian
	bool Open( const std::vector<std::string>& vecstrGenes, const std::vector<std::string>& vecstrDatasets,
		const std::string& strInputDirectory, const std::string& strOutputDirectory, size_t iFiles);

	//Qian
	CDatabase(bool isNibble) : CDatabaseImpl(isNibble){
	}

	bool GetGene(string, vector<unsigned char>&);

	/*!
	 * \brief
	 * Open an existing database from subfiles in the given directory.
	 * 
	 * \param strInputDirectory
	 * Directory from which the database files are read.
	 * 
	 * \returns
	 * True if the database was opened successfully.
	 */
	bool Open( const std::string& strInputDirectory );


	/*!
	 * \brief
	 * Retrieve data values from all datasets for a given gene pair.
	 * 
	 * \param iOne
	 * First gene index.
	 * 
	 * \param iTwo
	 * Second gene index.
	 * 
	 * \param vecbData
	 * Output vector containing the retrieved values, two per byte.
	 * 
	 * \returns
	 * True if data was retrieved successfully.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed; however, gene indices will be wrapped into the
	 * appropriate range using modulus, so something will be returned without a crash even for bad input.
	 * vecbData is automatically resized to the appropriate length.  Note that each data value is returned
	 * using only four bits, with even numbered datasets in the high order bits.
	 */
	bool Get( size_t iOne, size_t iTwo, std::vector<unsigned char>& vecbData ) const {

		return m_vecpDBs[ iOne % m_vecpDBs.size( ) ]->Get( iOne / m_vecpDBs.size( ), iTwo, vecbData ); }

	/*!
	 * \brief
	 * Retrieve data values from all gene pairs over all datasets for a given gene.
	 * 
	 * \param iGene
	 * Gene index.
	 * 
	 * \param vecbData
	 * Output vector containing the retrieved values, two per byte.
	 * 
	 * \param fReplace
	 * If true, replace values in output vector rather than appending.
	 * 
	 * \returns
	 * True if data was retrieved successfully.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed; however, the given gene index will be wrapped into the
	 * appropriate range using modulus, so something will be returned without a crash even for bad input.
	 * vecbData is automatically resized to the appropriate length.  Note that each data value is returned
	 * using only four bits, with even numbered datasets in the high order bits.  This is equivalent to
	 * calling Get with two gene indices repeatedly for iTwo from 0 to G and concatenating the results.
	 */
	bool Get( size_t iGene, std::vector<unsigned char>& vecbData, bool fReplace = false ) const {

		return m_vecpDBs[ iGene % m_vecpDBs.size( ) ]->Get( iGene / m_vecpDBs.size( ), vecbData, fReplace ); }

	/*!
	 * \brief
	 * Retrieve data values from the indicated gene pairs over all datasets for a given gene.
	 * 
	 * \param iGene
	 * First gene index.
	 * 
	 * \param veciGenes
	 * Second gene indices.
	 * 
	 * \param vecbData
	 * Output vector containing the retrieved values, two per byte.
	 * 
	 * \param fReplace
	 * If true, replace values in output vector rather than appending.
	 * 
	 * \returns
	 * True if data was retrieved successfully.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed; however, the gene indices will be wrapped into the
	 * appropriate range using modulus, so something will be returned without a crash even for bad input.
	 * vecbData is automatically resized to the appropriate length.  Note that each data value is returned
	 * using only four bits, with even numbered datasets in the high order bits.  This is equivalent to
	 * calling Get with two gene indices repeatedly for iTwo in veciGenes and concatenating the results.
	 */
	bool Get( size_t iGene, const std::vector<size_t>& veciGenes,
		std::vector<unsigned char>& vecbData, bool fReplace = false ) const {

		return m_vecpDBs[ iGene % m_vecpDBs.size( ) ]->Get( iGene / m_vecpDBs.size( ), veciGenes, vecbData,
			fReplace ); }

	/*!
	 * \brief
	 * Returns the total number of genes in the database.
	 * 
	 * \returns
	 * Total number of genes in the database.
	 */
	size_t GetGenes( ) const {

		return m_mapstriGenes.size( ); }

	/*!
	 * \brief
	 * Return the index of the given gene name, or -1 if it is not included in the database.
	 * 
	 * \param strGene
	 * Gene name to retrieve.
	 * 
	 * \returns
	 * Index of the requested gene name, or -1 if it is not in the CDat.
	 */
	size_t GetGene( const std::string& strGene ) const {

		return CDatabaseImpl::GetGene( strGene ); }

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
	 * For efficiency, no bounds checking is performed; however, the given gene index will be wrapped into the
	 * appropriate range using modulus, so something will be returned without a crash even for bad input.
	 * 
	 * \see
	 * GetGenes
	 */
	const std::string& GetGene( size_t iGene ) const {
		static const std::string	c_strEmpty	= "";

		return ( m_vecpDBs.empty( ) ? c_strEmpty : m_vecpDBs[ iGene % m_vecpDBs.size( ) ]->GetGene( iGene /
			m_vecpDBs.size( ) ) ); }

	/*!
	 * \brief
	 * Return the number of datasets stored in the database.
	 * 
	 * \returns
	 * Number of datasets in the database.
	 * 
	 * \remarks
	 * Number of datasets is constant across subfiles, so the number in the first subfile (if present) is
	 * actually returned.
	 */
	size_t GetDatasets( ) const {

		return ( m_vecpDBs.empty( ) ? 0 : m_vecpDBs[ 0 ]->GetDatasets( ) ); }


	bool Open(string &, vector<string> &, size_t &, size_t &);

	/*!
	 * \brief
	 * Set memory mapping behavior when opening DAB files.
	 * 
	 * \param fMemmap
	 * Value to store for memory mapping behavior.
	 * 
	 * If memory mapping is true, DAB files will be memory mapped when a database is constructed using Open.
	 * Default is false.
	 */
	void SetMemmap( bool fMemmap ) {

		m_fMemmap = fMemmap; }

	/*!
	 * \brief
	 * Set output block size (number of subfiles to create at once).
	 * 
	 * \param iSize
	 * Number of subfiles to create at once.
	 * 
	 * Set the number of output files (and thus the number of genes) to process simultaneously when creating
	 * a new database from DAB files.  Defaults to -1, indicating that all genes should be processed in one
	 * pass.
	 * 
	 * \remarks
	 * Optimal settings depend on genome size, number of datasets, physical memory of the host machine,
	 * and the input block size.
	 * 
	 * \see
	 * Open | SetBlockIn
	 */
	void SetBlockOut( size_t iSize ) {

		m_iBlockOut = iSize; }

	/*!
	 * \brief
	 * Set input block size (number of datasets to load at once).
	 * 
	 * \param iSize
	 * Number of datasets to load at once.
	 * 
	 * Set the number of input DAB files (and thus the number of datasets) to process simultaneously when
	 * creating a new database.  Defaults to -1, indicating that all datasets should be loaded at once.
	 * 
	 * \remarks
	 * Optimal settings depend on genome size, number of datasets, physical memory of the host machine,
	 * and the output block size.
	 * 
	 * \see
	 * Open | SetBlockOut
	 */
	void SetBlockIn( size_t iSize ) {

		m_iBlockIn = iSize; }

	/*!
	 * \brief
	 * Set buffering behavior when creating new database subfiles.
	 * 
	 * \param fBuffer
	 * Value to store for buffering behavior.
	 * 
	 * If buffering is true, database subfiles will be constructed in memory and written to disk as a single
	 * unit at the end of each block.  Default is false.
	 * 
	 * \remarks
	 * The jury's still out on whether this improves performance or not; it may be OS-dependent, since some
	 * operating systems do a better job of write-buffering than others.
	 * 
	 * \see
	 * Open
	 */
	void SetBuffer( bool fBuffer ) {

		m_fBuffer = fBuffer; }
};

}

#endif // DATABASE_H

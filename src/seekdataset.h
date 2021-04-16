/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a part of SEEK (Search-based exploration of expression compendium)
* which is authored and maintained by: Qian Zhu (qzhu@princeton.edu)
*
* If you use this file, please cite the following publication:
* Qian Zhu, Aaron K Wong, Arjun Krishnan, Miriam R Aure, Alicja Tadych, 
* Ran Zhang, David C Corney, Casey S Greene, Lars A Bongo, 
* Vessela N Kristensen, Moses Charikar, Kai Li & Olga G Troyanskaya
* "Targeted exploration and analysis of large cross-platform human 
* transcriptomic compendia" Nat Methods (2015)
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library for development, or use any other Sleipnir executable
* tools, please also cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#ifndef SEEKDATASET_H
#define SEEKDATASET_H

#include "seekbasic.h"
#include "seekmap.h"
#include "datapair.h"
#include "seekplatform.h"

namespace Sleipnir {

    class CSeekDBSetting {
    public:
        // For the default constructor create a link from the public reference
        //  values to the private actual values
        CSeekDBSetting() : dbDir(m_dbDirectory),
                           prepDir(m_prepDirectory),
                           platDir(m_platformDirectory),
                           sinfoDir(m_sinfoDirectory),
                           gvarDir(m_gvarDirectory),
                           geneMapFile(m_geneMapFile),
                           geneSymbolFile(m_geneSymbolFile),
                           datasetFile(m_dsetFile),
                           dsetSizeFile(m_dsetSizeFile),
                           quantFile(m_quantFile),
                           numDBs(m_numDB) {}

        CSeekDBSetting(const string &gvar, const string &sinfo,
                       const string &plat, const string &prep,
                       const string &db, const string &gene,
                       const string &gene_symbol, const string &quant,
                       const string &dset, const string &dset_size,
                       const utype &numDB) : CSeekDBSetting()
        {
            m_gvarDirectory = gvar;
            m_sinfoDirectory = sinfo;
            m_platformDirectory = plat;
            m_prepDirectory = prep;
            m_dbDirectory = db;
            m_geneMapFile = gene;
            m_geneSymbolFile = gene_symbol;
            m_quantFile = quant;
            m_dsetFile = dset;
            m_numDB = numDB;
            m_dsetSizeFile = dset_size;
        }

        CSeekDBSetting(CSeekDBSetting const *g) : CSeekDBSetting()
        {
            m_gvarDirectory = g->m_gvarDirectory;
            m_sinfoDirectory = g->m_sinfoDirectory;
            m_platformDirectory = g->m_platformDirectory;
            m_prepDirectory = g->m_prepDirectory;
            m_dbDirectory = g->m_dbDirectory;
            m_geneMapFile = g->m_geneMapFile;
            m_geneSymbolFile = g->m_geneSymbolFile;
            m_quantFile = g->m_quantFile;
            m_dsetFile = g->m_dsetFile;
            m_numDB = g->m_numDB;
            m_dsetSizeFile = g->m_dsetSizeFile;
        }

        ~CSeekDBSetting() {
        }

        utype GetNumDB() {
            return m_numDB;
        }

        friend ostream& operator<<(ostream& os, const CSeekDBSetting& dbs) {
            os << "DB:" << dbs.m_dbDirectory << endl;
            os << "Prep:" << dbs.m_prepDirectory << endl;
            os << "Plat:" << dbs.m_platformDirectory << endl;
            os << "Sinfo:" << dbs.m_sinfoDirectory << endl;
            os << "GVar:" << dbs.m_gvarDirectory << endl;
            os << "GeneMapFile:" << dbs.m_geneMapFile << endl;
            os << "GeneSymbolFile:" << dbs.m_geneSymbolFile << endl;
            os << "DatasetsFile:" << dbs.m_dsetFile << endl;
            os << "DatasetSizes:" << dbs.m_dsetSizeFile << endl;
            os << "QuantFile:" << dbs.m_quantFile << endl;
            os << "NumDBFiles:" << dbs.m_numDB << endl;
            return os;
        }

        const string &dbDir;
        const string &prepDir;
        const string &platDir;
        const string &sinfoDir;
        const string &gvarDir;
        const string &geneMapFile;
        const string &geneSymbolFile;
        const string &datasetFile;
        const string &dsetSizeFile;
        const string &quantFile;
        const utype &numDBs;

    private:
        string m_gvarDirectory;
        string m_sinfoDirectory;
        string m_platformDirectory;
        string m_prepDirectory;
        string m_dbDirectory;
        string m_geneMapFile;
        string m_geneSymbolFile;
        string m_quantFile;
        string m_dsetFile;
        string m_dsetSizeFile;
        utype m_numDB;
    };


/*!
 * \brief Representation of a microarray dataset that is used by Seek
 *
 * A \c CSeekDataset encapsulates the following information about the dataset:
 *
 * \li The gene-gene \a correlation matrix
 * \li Each gene's expression variance
 * \li Each gene's average \a correlation
 * \li The genes in the query that are present in the dataset
 * \li The genes in the genome that are present in the dataset
 * \li The platform which the dataset belongs to
 * \li The weight of the dataset that is assigned by the search algorithm
 *
 * This dataset structure is designed to be used by Seek.
 * \remarks
 * The word \a correlation refers to the standardized z-scores of Pearson correlations, which
 * is derived from a 2-step process:
 * \f[f(x,y)=\frac{1}{2}ln\frac{1+p(x,y)}{1-p(x,y)}\f]
 * where \f$p(x,y)\f$ is the Pearson correlation, \f$f(x,y)\f$ is the Fisher's transformed score.
 * \f[z(x,y)=\frac{f(x,y) - \bar{f}}{\sigma_{f}}\f]
 * where \f$z(x,y)\f$ is the z-score, \f$\bar{f}\f$ is the mean, and \f$\sigma_{f}\f$ is
 * the standard deviation.
 * 
 * From here on, \a correlation always refers to the above z-score definition.
 */

    class CSeekDataset {
    public:

        /*!
         * \enum DistanceMeasure
         * \brief Distance measure (see main section for descriptions)
         */
        enum DistanceMeasure {
            CORRELATION = 0, /**< Pearson correlations */
            Z_SCORE = CORRELATION + 1 /**< Z-score of Pearson correlations */
        };

        /*!
         * \brief
         * Constructor
         */
        CSeekDataset();

        /*!
         * \brief
         * Destructor
         */
        ~CSeekDataset();

        /*!
         * \brief
         * Read the \c *.sinfo file
         *
         * \param strFileName The file name
         *
         * The \c *.sinfo file contains the mean and the standard deviation
         * of the global gene-gene Pearson distribution for this dataset.
         */
        bool ReadDatasetAverageStdev(const string &);

        /*!
         * \brief
         * Read the gene average \a correlation file \c *.gavg
         *
         * \param strFileName The file name
         *
         * The \c *.gavg is an array that stores the average \a correlation of each
         * gene.
         */
        bool ReadGeneAverage(const string &);

        /*!
         * \brief
         * Read the gene variance file \c *.gvar
         *
         * \param strFileName The file name
         *
         * The \c *.gvar file is an array that stores the expression variance of each
         * gene.
         */
        bool ReadGeneVariance(const string &);

        /*!
         * \brief
         * Read the gene presence file \c *.gpres
         *
         * \param strFileName The file name
         *
         * The \c *.gpres is a 2-value array that contains the presence (1), absence (0)
         * status of genes.
         */
        bool ReadGenePresence(const string &);

        /*!
         * \brief
         * Initialize the genome presence map
         *
         * Indicates which genes of the genome are present in the dataset.
         */
        bool InitializeGeneMap();

        /*!
         * \brief
         * Initialize the query presence map
         *
         * \param query The query genes
         *
         * Indicates which query genes are present in the dataset.
         */
        bool InitializeQuery(const vector <utype> &);

        /*!
         * \brief
         * Initialize a presence map for a block of queries
         *
         * \param queryBlock A vector of queries
         *
         * Flattens all the queries into one vector that contains only the unique query genes, then
         * constructs a presence map based on this vector.
         */
        bool InitializeQueryBlock(const vector <utype> &);

        /*!
         * \brief
         * Delete the query
         *
         * Resets all query-related data, such as dataset weight, query presence map, etc.
         */
        bool DeleteQuery();

        /*!
         * \brief
         * Delete query block
         *
         * Resets all query-block related data.
         */
        bool DeleteQueryBlock();

        /*!
         * \brief
         * Initialize the gene-gene \a correlation matrix
         *
         * \param rD A two-dimensional array storing the discretized gene-gene \a correlations
         * \param quant The discretization function
         * \param iRows The number of rows for the \a correlation matrix
         * \param iColumns The number of columns for the \a correlation matrix
         * \param bSubtractAvg If true, subtract the \a correlation by the dataset average
         * \param bNormPlatform If true, subtract the \a correlation by the platform average and divide by standard deviation
         * \param logit If true, apply the logit transform on \a correlations
         * \param dist_measure Distance measure: z-score or correlations
         * \param cutoff Apply a hard cutoff on \a correlations
         * \param bRandom If true, shuffle the \a correlation vector
         * \param rand The random generator for the shuffling operation above
         * \remarks
         * The discretized \a correlation in the matrix \c rD is bounded by 0 to 255 (the limit of
         * \c unsigned \c char). The parameter \c quant specifies how a \a correlation is
         * discretized. For example, if the \c quant has 5 bins:
         * \code
         * [0, 1, 2, 3, 4]
         * \endcode
         * Then if a \a correlation is 2.5, the discretized value would be 2.
         */
        bool InitializeDataMatrix(utype **, const vector<float> &,
                                  const utype &, const utype &, const bool= true,
                                  const bool= false, const bool= false,
                                  const enum DistanceMeasure= Z_SCORE,
                                  const float cutoff = -1.0 * CMeta::GetNaN(),
                                  const bool= false, gsl_rng *rand = NULL);

        /*!
         * \brief
         * Copy constructor
         * \param src A given dataset
         */
        bool Copy(CSeekDataset *);

        /*!
         * \brief
         * Get the gene-gene \a correlation matrix
         *
         * \return A two-dimensional array of type \c utype. Note that the
         * \a correlation has been scaled to a integer range from 0 to 640.
         * See CSeekDataset::InitializeDataMatrix.
         *
         */
        utype **GetDataMatrix();

        /*!
         * \brief
         * Get the gene-gene \a correlation matrix
         *
         * \return A two-dimensional array of type \c unsigned \c char**.
         */
        unsigned char **GetMatrix();

        /*!
         * \brief Get the genome presence map
         * \return The genome presence map
         */
        CSeekIntIntMap *GetGeneMap();

        /*!
         * \brief Get the query-block presence map
         * \return The query-block presence map
         */
        CSeekIntIntMap *GetDBMap();

        /*!
         * \brief Get the query presence map
         * \return The query presence map
         */
        CSeekIntIntMap *GetQueryMap();

        /*!
         * \brief Get the query genes
         * \return A vector of queries
         */
        const vector <utype> &GetQuery() const;

        /*!
         * \brief Get the query gene indices
         * \return A vector of query gene indices
         */
        const vector <utype> &GetQueryIndex() const;

        /*!
         * \brief Get the gene expression variance vector
         * \return The variance vector
         */
        float GetGeneVariance(const utype &) const;

        /*!
         * \brief Get the gene average \a correlation vector
         * \return The average \a correlation vector
         */
        float GetGeneAverage(const utype &) const;

        /*!
         * \brief Get the mean of the global gene-gene Pearson distribution
         * \return The mean Pearson for the dataset
         */
        float GetDatasetAverage() const;

        /*!
         * \brief Get the standard deviation of the global gene-gene Pearson distribution
         * \return The standard deviation of the Pearson distribution
         */
        float GetDatasetStdev() const;

        /*!
         * \brief Get the genome size
         * \return The genome size
         */
        utype GetNumGenes() const;

        /*!
         * \brief Initialize the weight of the dataset
         *
         * \param i The number of cross-validations
         *
         * Initializes the total dataset weight, and the score of the
         * individual cross-validation (CV) runs.
         */
        bool InitializeCVWeight(const utype &);

        /*!
         * \brief Set the score for a particular cross-validation
         * \param i The index
         * \param f The validation score
         */
        bool SetCVWeight(const utype &, const float &);

        /*!
         * \brief Get the score for a particular cross-validation
         * \param i The index
         */
        float GetCVWeight(const utype &);

        /*!
         * \brief Get all the cross-validation scores
         * \return A vector of cross-validation scores
         */
        const vector<float> &GetCVWeight() const;

        /*!
         * \brief Get the dataset weight
         * \return The dataset weight
         */
        float GetDatasetSumWeight();

        /*!
         * \brief Set the platform
         * \param cp The platform
         */
        void SetPlatform(CSeekPlatform &);

        /*!
         * \brief Get the platform
         * \return The platform of this dataset
         */
        CSeekPlatform &GetPlatform() const;

    private:
        CSeekPlatform *platform;
        vector<float> geneAverage;
        vector<float> geneVariance;
        vector<char> genePresence;
        CSeekIntIntMap *geneMap;

        /* previously known as sinfo file */
        float m_fDsetAverage;
        float m_fDsetStdev;

        CSeekIntIntMap *dbMap;
        CSeekIntIntMap *queryMap;
        vector <utype> query;
        vector <utype> queryIndex;

        utype iQuerySize;
        utype iNumGenes;
        utype iDBSize;

        vector<float> weight;

        utype **rData;
        // r is gene-gene correlation matrix with gene correlations for one dataset followed by next dataset
        unsigned char **r;

        float sum_weight;
        bool m_bIsNibble;
    };


}
#endif

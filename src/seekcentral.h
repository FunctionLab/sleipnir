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
#ifndef SEEKCENTRAL_H
#define SEEKCENTRAL_H

#include <queue>

#include "seekbasic.h"
#include "seekdataset.h"
#include "seekplatform.h"
#include "seekmap.h"
#include "seekreader.h"
#include "seekquery.h"
#include "seekevaluate.h"
#include "database.h"
#include "datapair.h"
#include "seekweight.h"
#include "seekhelper.h"
#include "seekerror.h"

namespace Sleipnir {

/*!
 * \brief A suite of search algorithms that are supported by Seek
 *
 * The Seek search algorithms perform the coexpression search of the user's
 * query genes in a large compendium of microarray datasets. 
 * The output of the search algorithms is a ranking of genes based on their
 * gene score, which is determined by the overall weighted coexpression
 * to the query genes. 
 *
 * One of the first steps in a search is to weight
 * the datasets in such a way to prioritize informative datasets.
 * Then, with the weights generated, the final gene-score is given by:
 * \f[FS(g, Q)=\alpha\sum_{d \in D}{w_d \cdot s_d(g, Q)}\f]
 * where \f$w_d\f$ is the weight of the dataset, \f$s_d(g, Q)\f$ is the score
 * of \f$g\f$ to the query in the dataset, \f$\alpha\f$ is the normalization 
 * constant.
 * 
 * Currently the following dataset weighting algorithms are supported in Seek.
 * \li The query cross-validated (CV) weighting (CSeekCentral::CV): 
 * This is a weighting based on the query coexpression. The idea is to 
 * measure how well query genes are able to retrieve each other under a 
 * cross-validation setting.
 * To do so, we first divide the query into \a N
 * parts, use 1 part to build a small search instance, and use \f$N-1\f$ parts 
 * for evaluating the instance. The score of each instance \f$i\f$ is given by:
 * \f[s(i)=\sum_{g \in U}{(1-p)p^{rank(g)}}\f]
 * where \f$U\f$ is the genes in \f$N-1\f$ parts, \f$p\f$ is an exponential rate
 * parameter, \f$rank(g)\f$ is the position of \f$g\f$ in the ranking of genes 
 * generated by the search instance.
 *
 * \li Equal weighting (CSeekCentral::EQUAL): the weight is 1 for all datasets.
 *
 * \li User-supplied weight vector (CSeekCentral::USE_WEIGHT).
 * (ie., Seek does not calculate dataset weights)
 *
 * \li User-supplied gene-sets for weighting datasets, and also use cross-validations (CSeekCentral::CV_CUSTOM)
 *
 * \li Order-statistics (CSeekCentral::ORDER_STATISTICS): the algorithm used in MEM.
 * (Adler et al, Genome Biology 2009)
 *
 * CSeekCentral can handle multiple queries at a time, but the search parameters must remain
 * the same for all queries.
 */
    class CSeekCentral {
    public:

        /*!
         * \enum SearchMode
         * \brief Search modes (see section Detailed Descriptions)
         */
        enum SearchMode {
            CV = 0, /**< Cross-validated weighting */
            EQUAL = 1, /**< Equal weighting */
            USE_WEIGHT = 2, /**< User-supplied weights */
            CV_CUSTOM = 3, /**< Cross-validated weighting,
                      but instead of using the query genes to cross-validate, use the user 
                      supplied gene-sets to validate each query partition */
            ORDER_STATISTICS = 4, /**< MEM algorithm */
            AVERAGE_Z = 5 /**< Average z-scores between query, SPELL algorithm */
        };

        /*!
         * \brief Constructor
         */
        CSeekCentral();

        /*!
         * \brief Destructor
         */
        ~CSeekCentral();

        /*!
         * \brief Initialize function
         *
         * Performs the following operations:
         * \li Read the search parameters
         * \li Read the gene mapping \c gene_map.txt
         * \li Read a list of queries
         * \li Read the dataset mapping and the search datasets
         * \li Read the CDatabaselets (ie, the gene-gene \a correlations for the
         * query genes)
         *
         * \param gene The gene mapping file name, \c gene_map.txt
         * \param quant The quant file name
         * \param dset The dataset mapping file name, \c dataset_platform.txt
         * \param search_dset The file which contains the dataset names to be used for the search
         * \param query The query file name
         * \param platform The platform directory, which contains the platform
         * \a correlation averages and standard deviations
         * \param db The CDatabaselet directory, which contains the
         * gene-centric compendium-wide \a correlations, \c *.db files
         * \param prep The Prep directory, which contains the gene \a correlation
         * average \c *.gavg, and the gene presence \c *.gpres.
         * \param gvar The gene variance directory, which contains the \c *.gvar files
         * \param sinfo The sinfo directory, which contains the \c *.sinfo files
         * \param num_db The total number of CDatabaselet files
         * \param buffer The number of query genes to store in the memory
         * \param output_dir The output directory
         * \param to_output_text If true, output the gene-ranking in textual format
         * \param bOutputWeightComponent If true, output the dataset weight components (ie the score of cross-validations)
         * \param bSimulateWeight If true, use simulated weight as dataset weight
         * \param dist_measure Distance measure, either CORRELATION or Z_SCORE
         * \param bSubtractAvg If true, subtract the average z-score on a per-gene basis
         * \param bNormPlatform If true, subtract the platform gene average, divide by platform gene standard deviation
         * \param bLogit If true, apply the logit transformation on the \a correlations
         * \param fCutOff Cutoff the \a correlation values
         * \param fPercentRequired The fraction of the query genes required to be present in a dataset
         * in order to consider the dataset for integration
         * \param bSquareZ If true, square the \a correlations
         * \param bRandom If true, shuffle the \a correlation vector
         * \param iNumRandom The number of random simulations to perform per query
         * \param rand The random number generator
         * \param useNibble Default to false
         *
         * \remark The word \a correlation refers to the z-scored, standardized Pearson.
         * \remark The parameters \c bSubtractAvg, \c bNormPlatform,
         * \c bLogit, and \c bSquareZ are options to transform the
         * \a correlation values.
         * \remark The \c bSimulateWeight option is for equal weighting or order statistics where the final gene ranking
         * is not derived from a weighted integration of datasets. In this case, if the user still wants to see
         * the contribution of each dataset, the simulated weight is computed from the distance of a dataset's coexpression ranking to the final gene ranking.
         * \remark This function is designed to be used by SeekMiner.
         */
        bool InitializeFromSeekMiner(
                const vector<CSeekDBSetting *> &vecDBSetting,
                const char *search_dset, const char *query,
                const char *output_dir,
                const utype buffer = 20, const bool to_output_text = false,
                const bool bOutputWeightComponent = false, const bool bSimulateWeight = false,
                const enum CSeekDataset::DistanceMeasure dist_measure = CSeekDataset::Z_SCORE,
                const bool bVariance = false,
                const bool bSubtractAvg = true, const bool bNormPlatform = false,
                const bool bLogit = false, const float fCutOff = -9999,
                const float fPercentQueryRequired = 0, const float fPercentGenomeRequired = 0,
                const bool bSquareZ = false, const bool bRandom = false, const int iNumRandom = 10,
                const bool bNegativeCor = false, const bool bCheckDsetSize = false,
                gsl_rng *rand = NULL, const bool useNibble = false, const int numThreads = 8);

        /*!
         * \brief Initialize function
         *
         * Load everything except the query, the search datasets, and the output directory
         *
         * \param gene The gene mapping file name, \c gene_map.txt
         * \param quant The quant file name
         * \param dset The dataset mapping file name, \c dataset_platform.txt
         * \param platform The platform directory, which contains the platform
         * \a correlation average and standard deviation
         * \param db The CDatabaselet directory, which contains the
         * gene-centric compendium-wide \a correlations, \c *.db files
         * \param prep The Prep directory, which contains the gene \a correlation
         * average \c *.gavg, and the gene presence \c *.gpres.
         * Divided by datasets.
         * \param gvar The gene variance directory, which contains the \c *.gvar files
         * \param sinfo The sinfo directory, which contains the \c *.sinfo files
         * \param num_db The total number of CDatabaselet files
         * \param buffer The number of query genes to store in the memory
         * \param to_output_text If true, output the gene-ranking in the textual format
         * \param bOutputWeightComponent If true, output the dataset weight components (ie the score of cross-validations)
         * \param bSimulateWeight If true, use simulated weight as dataset weight
         * \param dist_measure Distance measure, either CORRELATION or Z_SCORE
         * \param bSubtractAvg If true, subtract the average z-score on a per-gene basis
         * \param bNormPlatform If true, subtract the platform gene average, divide by platform gene standard deviation
         * \param bLogit If true, apply the logit transformation on the \a correlations
         * \param fCutOff Cutoff the \a correlations
         * \param fPercentRequired The fraction of the query genes required to be present in a dataset
         * \param bSquareZ If true, square the \a correlations
         * \param bRandom If true, shuffle the \a correlation vector
         * \param iNumRandom The number of random simulations to perform per query
         * \param rand The random number generator
         * \param useNibble Default to false
         *
         * \remark The word \a correlation refers to the z-scored, standardized Pearson.
         * \remark The parameters \c bSubtractAvg, \c bNormPlatform,
         * \c bLogit, and \c bSquareZ are options to transform the
         * \a correlation values.
         * \remark The \c bSimulateWeight option is for equal weighting or order statistics where the final gene ranking
         * is not derived from a weighted integration of datasets. In this case, if the user still wants to see
         * the contribution of each dataset, the simulated weight is computed from the distance of a dataset's coexpression ranking to the final gene ranking.
         * \remark This function is designed to be used by SeekMiner.
         */
        bool Initialize(
                const vector<CSeekDBSetting *> &vecDBSetting,
                const utype buffer = 20, const bool to_output_text = false,
                const bool bOutputWeightComponent = false, const bool bSimulateWeight = false,
                const enum CSeekDataset::DistanceMeasure dist_measure = CSeekDataset::Z_SCORE,
                const bool bVariance = false,
                const bool bSubtractAvg = true, const bool bNormPlatform = false,
                const bool bLogit = false, const float fCutOff = -9999,
                const float fPercentQueryRequired = 0, const float fPercentGenomeRequired = 0,
                const bool bSquareZ = false, const bool bRandom = false, const int iNumRandom = 10,
                const bool bNegativeCor = false, const bool bCheckDsetSize = false,
                gsl_rng *rand = NULL, const bool useNibble = false, const int numThreads = 8);

        /* Initialize a SeekCentral instance from config file settings
         *  throws init_error
         */
        void InitializeFromSeekConfig(const SeekSettings &settings);

        /*!
         * \brief Initialize function
         *
         * Prepares Seek to be used in a client-server environment
         *
         * \param output_dir The output directory
         * \param query The query file name
         * \param search_dset The file that contains the name of datasets to be used for the search
         * \param src The CSeekCentral instance, where some settings will be copied to here
         * \param iClient The client's socket connection
         * \param query_min_required The minimum number of query genes required to be present in a dataset
         * \param dist_measure Distance measure, either CORRELATION or Z_SCORE.
         * \param bSubtractAvg If true, subtract the average z-score on a per-gene basis
         * \param bNormPlatform If true, subtract the platform gene average, divide by platform gene standard deviation
         *
         * \remark This function is designed to be used by SeekServer.
         * \remark The parameters \c bSubtractAvg, \c bNormPlatform
         * are options to transform the \a correlation values.
         * \remark Assumes that the CDatabaselets have been read, and the \c *.gvar, \c *.sinfo files have been loaded.
         * \remark Assumes that the dataset and gene mapping files have been read.
         */
        bool InitializeQuery(const string &output_dir, const string &query,
                        const string &search_dset, CSeekCentral *src, const int iClient,
                        const float query_min_required = 0, const float genome_min_required = 0,
                        const enum CSeekDataset::DistanceMeasure = CSeekDataset::Z_SCORE,
                        const bool bSubtractGeneAvg = true, const bool bNormPlatform = false,
                        const bool bNegativeCor = false, const bool bCheckDsetSize = false);

        /*!
         * \brief Run Seek with the cross-validated dataset weighting
         *
         * \param rnd The random number generator
         * \param PART_M Query partition mode
         * \param FOLD Number of partitions to generate from the query
         * \param RATE The weighting parameter \a p
         *
         * \remark The random number generator is used for partitioning the query.
         * \remark Assumes that the CSeekCentral::Initialize() has been called.
         */
        bool CVSearch(gsl_rng *, const CSeekQuery::PartitionMode &, const utype &, const float &);

        /*!
         * \brief Run Seek with the custom dataset weighting
         *
         * \param newGoldStd The gold-standard gene-set that is used for weighting datasets
         * \param rnd The random number generator
         * \param PART_M Query partition mode
         * \param FOLD Number of partitions to generate from the query
         * \param RATE The weighting parameter \a p	 *
         * Same as CVSearch, except that the weighting is not based on the coexpression
         * of the query genes, but based on the similarity of the query genes to some custom gold
         * standard gene-set.
         *
         * \remark The random number generator is used for partitioning the query.
         * \remark Assumes that the CSeekCentral::Initialize() has been called.
         */
        bool CVCustomSearch(const vector <vector<string>> &, gsl_rng *,
                            const CSeekQuery::PartitionMode &, const utype &, const float &);

        /*!
         * \brief Run Seek with the equal dataset weighting
         * \remark Assumes that the CSeekCentral::Initialize() has been called.
         */
        bool EqualWeightSearch();

        /*!
         * \brief Run Seek with the user-given dataset weights
         * \param weights A two-dimensional array that stores the user-given weights
         *
         * \remark The two-dimensional array \c weights is \a Q by \a D :
         * where \a Q is the number of queries, \a D is the number of datasets. \c weights[i][j]
         * stores the weight of dataset \a j in query \a i.
         *
         * \remark Assumes that the CSeekCentral::Initialize() has been called.
         */
        bool WeightSearch(const vector <vector<float>> &);

        /*!
         * \brief Run Seek with the variance weighted search
         *
         * Same as CSeekCentral::WeightSearch(), except that the user-given weights are the query gene expression variances.
         *
         * \remark Assumes that the CSeekCentral::Initialize() has been called.
         */
        bool VarianceWeightSearch();

        /*!
         * \brief Run Seek with the SPELL search
         *
         * \remark Assumes that the CSeekCentral::Initialize() has been called.
         */
        bool AverageWeightSearch();

        /*!
         * \brief Run Seek with the order statistics dataset weighting algorithm
         *
         * \remark Assumes that the CSeekCentral::Initialize() has been called.
         */
        bool OrderStatistics();

        /*!
         * \brief Get the final gene-ranking for all the queries
         * \return A two-dimensional array that stores the gene-rankings
         */
        const vector <vector<AResultFloat>> &GetAllResult() const;

        /*!
         * \brief Get all the queries
         * \return A vector of queries.
         */
        const vector <CSeekQuery> &GetAllQuery() const;

        /*!
         * \brief Get the dataset weight vector for all the queries
         * \return A two-dimensional \c float array that stores the weights
         *
         * \remark The first dimension is the query. The second dimension is the dataset.
         */
        const vector <vector<float>> &GetAllWeight() const;

        /*!
         * \brief Get the gene-map ID for a given \a gene-name
         * \param strGene The \a gene-name as a \c string
         * \return The gene-map ID
         */
        utype GetGeneId(const string &strGene) const;

        /*!
         * \brief Get the \a gene-name for a given gene-map ID
         * \param geneID The gene-map ID
         * \return The \a gene-name as a \c string
         */
        string GetGeneName(const utype &geneID) const;

        /*!
         * \brief Destruct this search instance
         * \return True if successful.
         */
        bool Destruct();

        /*!
         * \brief Get the maximum genome coverage among the
         * datasets in the compendium.
         */
        int GetMaxGenomeCoverage();

        /*!
         * \brief Get the sorted gene scores in paired <geneName, score>
         */
        vector<StrDoublePair> & getGeneResult(uint32_t queryIndex) {
            return this->m_geneResults[queryIndex];
        }

        /*!
         * \brief Get the sorted dataset weights in paired <datasetName, weight>
         */
        vector<StrDoublePair> & getDatasetResult(uint32_t queryIndex) {
            return this->m_datasetResults[queryIndex];
        }

        /*!
         * \brief Get number of queries performed
         */
        uint32_t numQueries() {
            return m_vecstrAllQuery.size();
        }

        /*!
         * \brief When using RPC communication, status messages will accumulate in the log queue
         */
        void setUsingRPC(bool val, queue<string> &log) {
            m_useRPC = val;
            m_rpcLog = &log;
        }

        void convertGenesEntrezToSymbol(const vector<string> &entrez, vector<string> &symbols);

        void convertGenesSymbolToEntrez(const vector<string> &symbols, vector<string> &entrez);

        string entrezToSymbol(string &entrez);

        void setSimulateWeightFlag(bool bSimulateWeight) {
            m_bSimulateWeight = bSimulateWeight;
        }

    void PrintSettings();

    private:
        //network mode
        bool EnableNetwork(const int &);

        bool CheckDatasets(const bool &);

        /* Central search function */
        bool Common(CSeekCentral::SearchMode &, gsl_rng * = NULL,
                    const CSeekQuery::PartitionMode * = NULL,
                    const utype * = NULL, const float * = NULL,
                    const vector <vector<float>> * = NULL,
                    const vector <vector<string>> * = NULL);

        bool CheckWeight(const utype &i);

        bool CopyTopGenes(CSeekQuery &, const vector <AResultFloat> &,
                          const utype);

        bool SetQueryScoreNull(const CSeekQuery &);

        bool PrepareQuery(const vector <string> &, CSeekQuery &);

        bool CalculateRestart();

        bool PrepareOneQuery(CSeekQuery &, CSeekIntIntMap &, vector<float> &);

        bool AggregateThreads();

        bool FilterResults(const utype &);

        bool getSortedGeneScores(vector <AResultFloat> &);

        bool getSortedDatasetScores(uint32_t queryIndex, vector <AResultFloat> &);

        void setGenePairedResult(uint32_t queryIndex, vector <AResultFloat> &sortedGeneScore);

        void setDatasetPairedResult(uint32_t queryIndex, vector <AResultFloat> &sortedDatasetWeights);

        bool Write(const utype &);

        bool Display(CSeekQuery &, vector <AResultFloat> &);

        /* Gene, Dataset, and Platform Mapping*/
        vector <string> m_vecstrGenes;
        vector <string> m_vecstrDatasets;
        vector <string> m_vecstrDP;
        map <string, string> m_mapstrstrDatasetPlatform;
        map <string, utype> m_mapstrintDataset; // map from dataset name to index in dset file
        map <string, utype> m_mapstrintGene;  // map from geneName to index in gene_map file
        map <string, string> m_geneEntrezToSymbolMap;
        map <string, string> m_geneSymbolToEntrezMap;
        vector <vector<string>> m_vecstrSearchDatasets;
        vector<CSeekIntIntMap *> m_searchdsetMap;

        /* Datasets */
        vector<CSeekDataset *> m_vc;

        /* Output */
        bool m_bOutputText;

        /* If true, output random case (ie shuffle rankings per dataset)
           iNumRandom: number of repetitions (Oct 26, 2012) */
        bool m_bRandom;
        int m_iNumRandom;
        gsl_rng *m_randRandom;
        /* random dataset weight over all repetitions */
        //vector<vector<float> > m_vecRandWeight;
        /* random gene scores over all repetitions */
        //vector<vector<float> > m_vecRandScore;

        /* Gene-gene correlation matrix for all datasets
         Organized per thread */
        utype ***m_rData;

        /* Correlation discretization */
        vector<float> m_quant;

        /* Correlation transformation options */
        bool m_bSubtractGeneAvg;
        bool m_bNormPlatform;
        enum CSeekDataset::DistanceMeasure m_eDistMeasure;
        bool m_bLogit;
        bool m_bSquareZ;

        /* multi-threaded programming */
        float **m_master_rank_threads;
        float **m_sum_weight_threads;
        float **m_sum_sq_weight_threads;
        utype **m_counts_threads;
        vector <utype> *m_rank_normal_threads;
        vector <utype> *m_rank_threads;

        /* Essential search results */
        vector<float> m_master_rank;
        vector<float> m_sum_weight;
        vector<float> m_sum_sq_weight;
        vector <utype> m_counts;

        /* Holds results for all queries */
        vector <vector<float>> m_weight;
        vector <vector<AResultFloat>> m_final;
        // for each query, pairs of geneName and Score
        vector<vector<StrDoublePair>> m_geneResults;
        // for each query<vec>, pairs of datasetName and Weight
        vector<vector<StrDoublePair>> m_datasetResults;

        /* Query */
        // vector of queries, each query is a vector of genes to query on
        vector <vector<string>> m_vecstrAllQuery;
        vector <CSeekQuery> m_Query;

        /* Platform */
        SeekPlatforms m_seekPlatforms;

        //CDatabase reference
        vector<CDatabase *> m_vecDB;
        vector <vector<string>> m_vecDBDataset; //A list of dsets in each CDatabase

        size_t m_iDatasets;
        size_t m_iGenes;
        utype m_numThreads;

        utype m_maxNumDB;
        // map by index to sets of queries to be loaded together (databaselets)
        map <utype, vector<vector < string>> >
        m_mapLoadTime;
        bool DEBUG;

        bool m_bOutputWeightComponent;
        bool m_bSimulateWeight;

        string m_output_dir;
        float m_fScoreCutOff;
        float m_fPercentQueryAfterScoreCutOff;
        float m_fPercentGenomeRequired;

        /* for order statistics, a datasets-by-genes matrix */
        utype **m_rank_d;

        /* for network mode */
        int m_iClient;
        bool m_bEnableNetwork;
        bool m_useRPC;
        queue<string> *m_rpcLog;
        //bool m_bSharedDB; //if m_DB is shared between multiple CSeekCentral instances
        bool m_bNegativeCor;

        vector<CSeekDBSetting *> m_vecDBSetting; //DBSetting
        bool m_useNibble;

        float m_DEFAULT_NA;

        /* for specifying dataset size */
        bool m_bCheckDsetSize;
        int m_iNumSampleRequired;
        map <string, utype> m_mapstrintDatasetSize;
    };


}
#endif

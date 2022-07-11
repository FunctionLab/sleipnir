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
#include "seekcentral.h"
#include "seekerror.h"
#include <filesystem>
#include <cassert>
#include <boost/algorithm/string/predicate.hpp>
#include "stdafx.h"

namespace Sleipnir {

    CSeekCentral::CSeekCentral() {
        m_vecstrSearchDatasets.clear();
        m_vc.clear();
        m_quant.clear();
        m_vecstrAllQuery.clear();
        m_searchdsetMap.clear();
        m_vecDB.clear();
        m_rData = NULL;
        m_maxNumDB = 50;

        m_master_rank_threads = NULL;
        m_sum_weight_threads = NULL;
        m_sum_sq_weight_threads = NULL; //sum squared weight
        m_counts_threads = NULL;
        m_rank_normal_threads = NULL;
        m_rank_threads = NULL;
        m_rank_d = NULL;

        m_master_rank.clear();
        m_sum_weight.clear();
        m_sum_sq_weight.clear();
        m_weight.clear();
        m_counts.clear();
        m_mapLoadTime.clear();

        m_Query.clear();
        m_final.clear();

        m_iDatasets = 0;
        m_iGenes = 0;
        m_numThreads = 0;
        m_bSubtractGeneAvg = false;
        m_bNormPlatform = false;
        m_bLogit = false;
        m_eDistMeasure = CSeekDataset::Z_SCORE;
        m_bOutputWeightComponent = false;
        m_bSimulateWeight = false;
        m_bOutputText = false;
        m_bSquareZ = false;

        m_bNegativeCor = false;
        m_DEFAULT_NA = -320;

        m_useNibble = false;

        DEBUG = false;
        m_output_dir = "";

        m_iClient = -1;
        m_bEnableNetwork = false;
        m_hasPclInDatasetName = false;

        m_bCheckDsetSize = false;
        m_iNumSampleRequired = 10; //if checking for dataset size

        m_fPercentDatasetCoverage = 0.5;
    }

    CSeekCentral::~CSeekCentral() {
        // call common destructor
        this->Destruct();

        m_vecstrSearchDatasets.clear();

        m_quant.clear();

        //m_rData, m_master_rank_threads,
        //m_sum_weight_threads, m_counts_threads
        //m_rank_normal_threads, and m_rank_threads
        //should be already freed

        m_master_rank.clear();
        m_sum_weight.clear();
        m_sum_sq_weight.clear();
        m_counts.clear();
        m_weight.clear();
        m_final.clear();

        m_vecstrAllQuery.clear();
        m_Query.clear();

        m_iDatasets = 0;
        m_iGenes = 0;
        m_numThreads = 0;
        m_mapLoadTime.clear();
        m_output_dir = "";
        DEBUG = false;

        m_useNibble = false;
        m_bNegativeCor = false;

        m_bCheckDsetSize = false;
        m_iNumSampleRequired = 0;

        m_fPercentDatasetCoverage = 0;
    }

    bool CSeekCentral::CalculateRestart() {
        set <string> ss;
        m_mapLoadTime.clear();
        utype i, prev;
        prev = 0;
        for (i = 0; i < m_vecstrAllQuery.size(); i++) {
            if (m_vecstrAllQuery[i].size() > m_maxNumDB) {
                fprintf(stderr, "Cannot fit query %d in buffer\n", i);
                return false;
            }
        }

        // Here it creates m_mapLoadTime which combines together multiple query
        //  string genes whose data will be loaded into memory at the same time
        vector <vector<string>> *vv = NULL;
        m_mapLoadTime[prev] = vector < vector < string > > ();
        for (i = 0; i < m_vecstrAllQuery.size(); i++) {
            ss.insert(m_vecstrAllQuery[i].begin(), m_vecstrAllQuery[i].end());
            vv = &(m_mapLoadTime[prev]);
            vv->push_back(m_vecstrAllQuery[i]);
            if (ss.size() > m_maxNumDB) {
                vv->pop_back();
                ss.clear();
                ss.insert(m_vecstrAllQuery[i].begin(), m_vecstrAllQuery[i].end());
                prev = i;
                m_mapLoadTime[prev] = vector < vector < string > > ();
                vv = &(m_mapLoadTime[prev]);
                vv->push_back(m_vecstrAllQuery[i]);
            }
        }
        ss.clear();

        //check
        vector<char> vc;
        utype tot = 0;
        CSeekTools::InitVector(vc, m_vecstrAllQuery.size(), (char) 0);
        map < utype, vector < vector < string > > > ::const_iterator
        ci =
                m_mapLoadTime.begin();
        for (; ci != m_mapLoadTime.end(); ci++) tot += (ci->second).size();
        vc.clear();
        if (tot != m_vecstrAllQuery.size()) {
            fprintf(stderr, "ERROR, size do not match\n");
            return false;
        }
        return true;
    }


    // Per-query version of Initialize() called by SeekServer
    //for SeekServer
    //assume DB has been read (with gvar, sinfo information)
    //assume datasets and genes have been read
    //assume m_enableNetwork is on
    //* CDatabaselet collection is shared between multiple clients (m_bSharedDB)
    bool CSeekCentral::InitializeQuery(
            const string &output_dir, const string &query,
            const string &search_dset, CSeekCentral *src, const int iClient,
            const float query_min_required, const float genome_min_required,
            const enum CSeekDataset::DistanceMeasure eDistMeasure,
            const bool bSubtractGeneAvg, const bool bNormPlatform,
            const bool bNegativeCor, const bool bCheckDsetSize,
            const float percentDatasetCoverage) {

        g_CatSleipnir().debug("InitializeQuery: Request received from client");
        assert(m_missingInitParams == false);

        // Copy a pointer to the read only attributes
        roAttr = src->roAttr;

        m_output_dir = output_dir; //LATER, TO BE DELETED
        m_maxNumDB = src->m_maxNumDB;
        //m_bSharedDB = true;
        m_numThreads = src->m_numThreads;
        m_fScoreCutOff = src->m_fScoreCutOff;
        m_fPercentQueryAfterScoreCutOff = query_min_required;
        m_fPercentGenomeRequired = genome_min_required;
        m_bSquareZ = src->m_bSquareZ;
        m_bOutputText = src->m_bOutputText;
        m_bSubtractGeneAvg = bSubtractGeneAvg;
        m_bNormPlatform = bNormPlatform;
        m_bLogit = src->m_bLogit;
        m_eDistMeasure = eDistMeasure;
        m_bCheckDsetSize = bCheckDsetSize;
        m_iNumSampleRequired = src->m_iNumSampleRequired;
        m_fPercentDatasetCoverage = percentDatasetCoverage;

        //if negative correlation, then need to use a different null value
        m_bNegativeCor = bNegativeCor;
        if (m_bNegativeCor) {
            m_DEFAULT_NA = 320;
        } else {
            m_DEFAULT_NA = -320;
        }

        m_bOutputWeightComponent = src->m_bOutputWeightComponent;
        m_bSimulateWeight = src->m_bSimulateWeight;

        m_bRandom = false;
        m_iNumRandom = 1;
        m_randRandom = NULL;

        m_quant = src->m_quant;
        utype i, j;
        omp_set_num_threads(m_numThreads);

        m_iDatasets = roAttr->m_vecstrDatasets.size();
        m_iGenes = roAttr->m_vecstrGenes.size();

        if (m_bCheckDsetSize && roAttr->m_mapstrintDatasetSize.size() == 0) {
            throw config_error("Error: Check_dataset_size requested with no dataset size file configured");
        }

        //read search datasets
        vector <string> sd;
        CMeta::Tokenize(search_dset.c_str(), sd, "|", false);
        m_vecstrSearchDatasets.resize(sd.size());
        for (i = 0; i < sd.size(); i++) {
            m_vecstrSearchDatasets[i] = vector<string>();
            vector <string> vecsearchDset;
            if (sd[i] == "" or sd[i] == "NA") {
                // include all datasets as part of search
                for (j = 0; j < roAttr->m_vecstrDatasets.size(); j++) {
                    if (m_bCheckDsetSize) {
                        if (roAttr->m_mapstrintDatasetSize.count(roAttr->m_vecstrDatasets[j]) == 0 ||
                            roAttr->m_mapstrintDatasetSize.at(roAttr->m_vecstrDatasets[j]) <
                                            m_iNumSampleRequired) {
                            continue;
                        }
                    }
                    m_vecstrSearchDatasets[i].push_back(roAttr->m_vecstrDatasets[j]);
                }
            } else {
                CMeta::Tokenize(sd[i].c_str(), vecsearchDset, " ", false);
                for (j = 0; j < vecsearchDset.size(); j++) {
                    if (m_bCheckDsetSize) {
                        if (roAttr->m_mapstrintDatasetSize.count(vecsearchDset[j]) == 0 ||
                            roAttr->m_mapstrintDatasetSize.at(vecsearchDset[j]) <
                                            m_iNumSampleRequired) {
                            continue;
                        }
                    }
                    m_vecstrSearchDatasets[i].push_back(vecsearchDset[j]);
                }
            }
            m_vecstrSearchDatasets[i].resize(m_vecstrSearchDatasets[i].size());
        }
        //read queries
        vector <string> sq;
        CMeta::Tokenize(query.c_str(), sq, "|", false);
        m_vecstrAllQuery.resize(sq.size());
        for (i = 0; i < sq.size(); i++) {
            CMeta::Tokenize(sq[i].c_str(), m_vecstrAllQuery[i], " ", false);
        }

        m_searchdsetMap.resize(m_vecstrAllQuery.size());
        for (i = 0; i < m_vecstrAllQuery.size(); i++) {
            m_searchdsetMap[i] = new CSeekIntIntMap(roAttr->m_vecstrDatasets.size());
            for (j = 0; j < m_vecstrSearchDatasets[i].size(); j++) {
                if (roAttr->m_mapstrintDataset.count(m_vecstrSearchDatasets[i][j]) > 0) {
                    m_searchdsetMap[i]->Add(
                            roAttr->m_mapstrintDataset.at(m_vecstrSearchDatasets[i][j]));
                }
            }
        }

        m_useNibble = src->m_useNibble;

        m_vecDB.resize(src->m_vecDB.size());
        //commented out Jan 8, 2014
        //for(i=0; i<m_vecDB.size(); i++)
        //	m_vecDB[i] = src->m_vecDB[i];
        for (i = 0; i < m_vecDB.size(); i++) {
            bool res;
            m_vecDB[i] = NULL;
            m_vecDB[i] = new CDatabase(m_useNibble);
            res = m_vecDB[i]->Open(roAttr->m_vecDBSetting[i]->dbDir,
                                   roAttr->m_vecstrGenes, roAttr->m_vecDBDataset[i].size(),
                                   roAttr->m_vecDBSetting[i]->GetNumDB());
            if (res == false) {
                return false;
            }
        }

        CSeekTools::LoadDatabase(m_vecDB, m_iGenes, m_iDatasets,
                                 m_vc, src->m_vc, roAttr->m_seekPlatforms.getCSeekPlatforms(), roAttr->m_vecstrDatasets,
                                 roAttr->m_mapstrstrDatasetPlatform, roAttr->m_seekPlatforms.getPlatformMap());

        if (!CalculateRestart()) {
            fprintf(stderr, "Error occurred during CalculateRestart()\n");
            return false;
        }

        //fprintf(stderr, "Finished CalculateRestart()\n");

        if (!EnableNetwork(iClient)) {
            fprintf(stderr, "Error occurred during EnableNetwork()\n");
            return false;
        }

        //fprintf(stderr, "Finished EnableNetworks()\n");

        if (!CheckDatasets(true)) { //replace parameter is true
            fprintf(stderr, "Error occurred during CheckDatasets()\n");
            return false;
        }

        //fprintf(stderr, "Finished CheckDatasets()\n");
        return true;
    }

//network mode, meant to be run after Initialize()
    bool CSeekCentral::EnableNetwork(const int &iClient) {
        if (iClient > 0) {
            m_bEnableNetwork = true;
            m_iClient = iClient; //assume client connection is already open
        } else {
            m_bEnableNetwork = false;
            m_iClient = 0;
            if (m_useRPC) {
                // clear the log queue
                while (!m_rpcLog->empty()) {
                    m_rpcLog->dequeue();
                }
            }
        }
        return true;
    }

//optional step
//Checks how many datasets contain the query
//requires the queries and searchdatasets to be loaded!
    bool CSeekCentral::CheckDatasets(const bool &replace) {
        utype dd, j;
        utype l;

        stringstream ss; //search dataset (new!)
        stringstream sq; //query availability
        stringstream aq; //query (new!)

        int maxGCoverage = GetMaxGenomeCoverage();

        // iterate over queries
        for (l = 0; l < m_searchdsetMap.size(); l++) {
            utype iUserDatasets = m_searchdsetMap[l]->GetNumSet();
            const vector <utype> &allRDatasets = m_searchdsetMap[l]->GetAllReverse();
            vector<int> count; // stores dataset counts per query gene
            CSeekTools::InitVector(count, m_vecstrAllQuery[l].size(), (int) 0);
            bool isFirst = true;
            //fprintf(stderr, "iUserDatasets %d\n", iUserDatasets);

            // iterate over datasets
            for (dd = 0; dd < iUserDatasets; dd++) {
                utype i = allRDatasets[dd];
                CSeekIntIntMap *si = m_vc[i]->GetGeneMap();
                utype present = 0; // total genes in dataset

                // iterate over query genes, increment if gene is in dd
                for (j = 0, present = 0; j < m_vecstrAllQuery[l].size(); j++) {
                    if (roAttr->m_mapstrintGene.find(m_vecstrAllQuery[l][j]) ==
                        roAttr->m_mapstrintGene.end())
                        continue;
                    if (CSeekTools::IsNaN(si->GetForward(
                            roAttr->m_mapstrintGene.at(m_vecstrAllQuery[l][j]))))
                        continue;
                    count[j]++;
                    present++;
                }

                //datasets that contains all query genes (very stringent)
                //if(present==m_vecstrAllQuery[l].size()){

                int minRequired = 0;
                if (m_vecstrAllQuery[l].size() >= 5) {
                    minRequired = (int) (m_fPercentQueryAfterScoreCutOff *
                                         m_vecstrAllQuery[l].size());
                } else if (m_vecstrAllQuery[l].size() <= 2) {
                    minRequired = m_vecstrAllQuery[l].size();
                } else {
                    minRequired = 2;
                }

                int minGRequired = (int) (m_fPercentGenomeRequired *
                                          (float) maxGCoverage);

                //datasets containing some query genes (relaxed) [ 0 ]
                // dataset dd contains enough query genes
                if (present > 0 && present >= minRequired &&
                    si->GetNumSet() >= minGRequired) {
                    if (isFirst) {
                        isFirst = false;
                        ss << roAttr->m_vecstrDatasets[i];
                    } else {
                        ss << " " << roAttr->m_vecstrDatasets[i];
                    }
                }
            }

            if (isFirst) {
                string err = "Error: no dataset contains enough query genes";
                fprintf(stderr, "%s (%d)\n", err.c_str(), l);
                if (m_bEnableNetwork) {
                    CSeekNetwork::Send(m_iClient, err);
                }
            }

            if (l != m_searchdsetMap.size() - 1) {
                ss << "|";
            }

            //fprintf(stderr, "ss %s\n", ss.str().c_str());

            // iterate over query genes
            isFirst = true;
            for (j = 0; j < m_vecstrAllQuery[l].size(); j++) {
                sq << m_vecstrAllQuery[l][j] << ":" << count[j];
                if (j != m_vecstrAllQuery[l].size() - 1) {
                    sq << ";";
                }
                if (count[j] == 0) continue; // only add genes with datasets
                if (isFirst) {
                    isFirst = false;
                    aq << m_vecstrAllQuery[l][j];
                } else {
                    aq << " " << m_vecstrAllQuery[l][j];
                }
            }

            if (isFirst) {
                string err = "Error: no dataset contains any of the query genes";
                fprintf(stderr, "%s (%d)\n", err.c_str(), l);
                if (m_bEnableNetwork) {
                    CSeekNetwork::Send(m_iClient, err);
                }
            }

            if (l != m_searchdsetMap.size() - 1) {
                aq << "|";
                sq << "|";
            }

            //fprintf(stderr, "sq %s\n", sq.str().c_str());
            //fprintf(stderr, "aq %s\n", aq.str().c_str());
        }

        string refinedQuery = aq.str();
        string refinedSearchDataset = ss.str();
        string refinedGeneCount = sq.str();

        if (m_bEnableNetwork) {
            CSeekNetwork::Send(m_iClient, refinedSearchDataset);
            CSeekNetwork::Send(m_iClient, refinedGeneCount);
        } else if (m_useRPC) {
            m_rpcLog->enqueue(refinedSearchDataset);
            m_rpcLog->enqueue(refinedGeneCount);
        }

        if (replace) {
            vector <string> qq;
            utype i;
            CMeta::Tokenize(refinedQuery.c_str(), qq, "|", false);
            m_vecstrAllQuery.resize(qq.size());
            for (i = 0; i < qq.size(); i++) {
                m_vecstrAllQuery[i].clear();
                CMeta::Tokenize(qq[i].c_str(), m_vecstrAllQuery[i], " ", true);
            }

            if (m_bCheckDsetSize && roAttr->m_mapstrintDatasetSize.size() == 0) {
                throw config_error("Error: Check_dataset_size requested with no dataset size file configured");
            }

            //Change the search datasets
            vector <string> sd;
            CMeta::Tokenize(refinedSearchDataset.c_str(), sd, "|", false);

            if (sd.size() != qq.size()) {
              string err = "Error: dataset size and query size does not match";
              fprintf(stderr, "%s\n", err.c_str());
              throw query_error(err);
            }

            m_vecstrSearchDatasets.resize(sd.size());
            for (i = 0; i < sd.size(); i++) {
                m_vecstrSearchDatasets[i] = vector<string>();
                vector <string> vecsearchDset;
                CMeta::Tokenize(sd[i].c_str(), vecsearchDset, " ", false);
                for (j = 0; j < vecsearchDset.size(); j++) {
                    if (m_bCheckDsetSize) {
                        if (roAttr->m_mapstrintDatasetSize.count(vecsearchDset[j]) == 0 ||
                            roAttr->m_mapstrintDatasetSize.at(vecsearchDset[j]) <
                                            m_iNumSampleRequired) {
                            continue;
                        }
                    }
                    m_vecstrSearchDatasets[i].push_back(vecsearchDset[j]);
                }
                m_vecstrSearchDatasets[i].resize(m_vecstrSearchDatasets[i].size());
            }

            //fprintf(stderr, "replace datasets %d %d\n", sd.size(), m_vecstrSearchDatasets[0].size());
            //fprintf(stderr, "searchdsetMap size %d\n", m_searchdsetMap.size());

            for (i = 0; i < m_searchdsetMap.size(); i++) {
                delete m_searchdsetMap[i];
                m_searchdsetMap[i] = NULL;
            }
            m_searchdsetMap.clear();
            //fprintf(stderr, "cleared searchdsetMap\n");

            m_searchdsetMap.resize(m_vecstrAllQuery.size());
            for (i = 0; i < m_vecstrAllQuery.size(); i++) {
                m_searchdsetMap[i] = new CSeekIntIntMap(roAttr->m_vecstrDatasets.size());
                for (j = 0; j < m_vecstrSearchDatasets[i].size(); j++) {
                    if (roAttr->m_mapstrintDataset.count(m_vecstrSearchDatasets[i][j]) > 0) {
                        m_searchdsetMap[i]->Add(
                                roAttr->m_mapstrintDataset.at(m_vecstrSearchDatasets[i][j]));
                    }
                }
            }

        }

        return true;
    }

    void CSeekCentral::InitializeFromSeekConfig(const SeekSettings &settings) {
        bool bOutputWeightComponent = true;
        bool bSimulateWeight = false;
        bool bSubtractAvg = false;
        bool bNormPlatform = false;
        bool bLogit = false;
        bool bVariance = false;
        bool useNibble = false;

        if (settings.isNibble == true) {
            string errStr = "Nibble integration is not supported! Please use a non-nibble CDatabase";
            throw config_error(FILELINE + errStr);
        }
        // Commented out to allow PCL Server to use SeekCentral
        // Added checks that when m_bCheckDsetSize is true
        // then m_mapstrintDatasetSize isn't empty
        // if (settings.dbs[0]->dsetSizeFile == "NA")
        // {
        //     // Must be set so the query request can decide whether to use check dataset size
        //     string errStr = "Dataset size file is missing";
        //     throw config_error(FILELINE + errStr);
        // }
        if (settings.dbs[0]->gvarDir != "NA")
        {
            bVariance = true;
        }

        bool res;
        res = Initialize(settings.dbs,
                         settings.numBufferedDBs,
                         settings.outputAsText,
                         bOutputWeightComponent,
                         bSimulateWeight,
                         CSeekDataset::CORRELATION, //to be overwritten by individual search instance's setting
                         bVariance, //decide whether or not to load gvar
                         bSubtractAvg,
                         bNormPlatform, //to be overwritten by individual search instance's settings
                         bLogit, //always false
                         settings.scoreCutoff,
                         0.0, //min query fraction (to be overwrriten)
                         0.0, //min genome fraction (to be overwrriten))
                         settings.squareZ,
                         false,  // bRandom
                         1,  // iNumRandom
                         false, //negative cor (to be overwritten)
                         true, //check dataset size (to be overwritten)
                         NULL,  // gsl_rng
                         useNibble,
                         settings.numThreads,
                         0.5);  // min percent datasets
        if (res == false) {
            throw init_error(FILELINE + "Seek Initialize");
        }
    }

    // In-common version of Initialize(), to be called at startup time, no per-query info.
    // load everything except query, search datasets, output directory
    // called by SeekServer at startup time
    bool CSeekCentral::Initialize(const vector<CSeekDBSetting *> &vecDBSetting,
                                  const utype buffer, const bool to_output_text,
                                  const bool bOutputWeightComponent, const bool bSimulateWeight,
                                  const enum CSeekDataset::DistanceMeasure dist_measure, const bool bVariance,
                                  const bool bSubtractAvg, const bool bNormPlatform,
                                  const bool bLogit, const float fCutOff, const float fPercentQueryRequired,
                                  const float fPercentGenomeRequired,
                                  const bool bSquareZ, const bool bRandom, const int iNumRandom,
                                  const bool bNegativeCor, const bool bCheckDatasetSize,
                                  gsl_rng *rand, const bool useNibble, const int numThreads, float percentDatasetCoverage) {

        auto attr = make_shared<ReadOnlyAttributes>();

        m_maxNumDB = buffer;
        m_numThreads = numThreads; //changed from 8

        m_fScoreCutOff = fCutOff;
        m_fPercentQueryAfterScoreCutOff = fPercentQueryRequired;
        m_fPercentGenomeRequired = fPercentGenomeRequired;
        m_bSquareZ = bSquareZ;
        m_fPercentDatasetCoverage = percentDatasetCoverage;

        m_bNegativeCor = bNegativeCor;
        if (m_bNegativeCor) {
            m_DEFAULT_NA = 320;
        } else {
            m_DEFAULT_NA = -320;
        }

        m_bCheckDsetSize = bCheckDatasetSize;

        m_bOutputWeightComponent = bOutputWeightComponent;
        m_bSimulateWeight = bSimulateWeight;

        //random retrieval==========================
        m_bRandom = bRandom;
        m_iNumRandom = iNumRandom;
        m_randRandom = rand; //random-case only
        //==========================================

        if (!m_bRandom) {
            m_randRandom = NULL;
        }

        utype i, j;

        omp_set_num_threads(m_numThreads);

        m_bOutputText = to_output_text;
        m_bSubtractGeneAvg = bSubtractAvg;
        m_bNormPlatform = bNormPlatform;
        m_bLogit = bLogit;
        m_eDistMeasure = dist_measure;

        bool bCorrelation = false;

        if (dist_measure == CSeekDataset::CORRELATION) {
            bCorrelation = true;
            if (m_bSubtractGeneAvg || m_bNormPlatform || m_bLogit) {
                g_CatSleipnir().warn(
                        "Warning: setting subtract_avg, norm platform to false");
                m_bSubtractGeneAvg = false;
                m_bNormPlatform = false;
                m_bLogit = false;
            }
        }

        //read genes
        vector <string> vecstrGeneID;
        if (!CSeekTools::ReadListTwoColumns(vecDBSetting[0]->geneMapFile,
                                            vecstrGeneID, attr->m_vecstrGenes))
            return false;
        for (i = 0; i < attr->m_vecstrGenes.size(); i++)
            attr->m_mapstrintGene[attr->m_vecstrGenes[i]] = i;

        // read gene symbol map if provided
        if (vecDBSetting[0]->geneSymbolFile != "NA") {
            vector <string> geneEntrez;
            vector <string> geneSymbol;
            if (!CSeekTools::ReadListTwoColumns(vecDBSetting[0]->geneSymbolFile,
                                                geneEntrez, geneSymbol)) {
                throw init_error(FILELINE + "Unable to read gene symbol file: " + vecDBSetting[0]->geneSymbolFile);
            }
            // create forward and reverse map
            uint32_t numGenes = geneEntrez.size();
            for (int i=0; i<numGenes; i++) {
                attr->m_geneEntrezToSymbolMap.insert({geneEntrez[i], geneSymbol[i]});
                attr->m_geneSymbolToEntrezMap.insert({geneSymbol[i], geneEntrez[i]});
            }
            if (attr->m_geneEntrezToSymbolMap.size() != numGenes) {
                throw init_error(FILELINE + "entrezToSymbol map count mismatch: " +
                                 to_string(attr->m_geneEntrezToSymbolMap.size()) +
                                 ", " + to_string(numGenes));
            }
            if (attr->m_geneSymbolToEntrezMap.size() != numGenes) {
                uint32_t numDups = numGenes - attr->m_geneSymbolToEntrezMap.size();
                cerr << "Warning: duplicate symbols in symbolToEntrez map: " << numDups << endl;
            }
        }

        //read quant file
        if (vecDBSetting[0]->quantFile != "NA") {
            CSeekTools::ReadQuantFile(vecDBSetting[0]->quantFile, m_quant);
        } else {
            m_missingInitParams = true;
        }

        attr->m_vecstrDatasets.clear();
        attr->m_vecstrDP.clear();
        attr->m_mapstrintDatasetDB.clear();
        attr->m_mapstrstrDatasetPlatform.clear();
        attr->m_mapstrintDataset.clear();

        attr->m_seekPlatforms.clear();

        m_vecDB.resize(vecDBSetting.size());
        attr->m_vecDBDataset.resize(vecDBSetting.size());
        for (i = 0; i < vecDBSetting.size(); i++)
            m_vecDB[i] = NULL;

        attr->m_vecDBSetting.resize(vecDBSetting.size());
        for (i = 0; i < attr->m_vecDBSetting.size(); i++)
            attr->m_vecDBSetting[i] = new CSeekDBSetting(vecDBSetting[i]);

        m_useNibble = useNibble;

        for (i = 0; i < vecDBSetting.size(); i++) {
            if (dist_measure == CSeekDataset::CORRELATION &&
                vecDBSetting[i]->sinfoDir == "NA") {
                g_CatSleipnir().warn("WARNING: sinfo dir not specified!");
                m_missingInitParams = true;
            }

            m_vecDB[i] = new CDatabase(useNibble);
            //read datasets
            vector <string> vD, vDP;
            if (!CSeekTools::ReadListTwoColumns(vecDBSetting[i]->datasetFile, vD, vDP))
                return false;

            if (vD.size() > 0) {
                using boost::algorithm::ends_with;
                if (ends_with(vD[0], "pcl")) {
                    m_hasPclInDatasetName = true;
                } else {
                    m_hasPclInDatasetName = false;
                }
            }
            for (j = 0; j < vD.size(); j++) {
                attr->m_vecstrDatasets.push_back(vD[j]);
                attr->m_vecDBDataset[i].push_back(vD[j]);
                attr->m_vecstrDP.push_back(vDP[j]);
                attr->m_mapstrintDatasetDB[vD[j]] = (int) i;
            }

            if (vecDBSetting[i]->dsetSizeFile != "NA") {
                vector <string> col1, col2;
                // TODO GW - fix finicky: dataset_size file must have tabs not spaces
                if (!CSeekTools::ReadListTwoColumns(vecDBSetting[i]->dsetSizeFile, col1, col2))
                    return false;
                set <string> currentD(vD.begin(), vD.end());
                for (j = 0; j < col1.size(); j++) {
                    if (currentD.find(col1[j]) == currentD.end()) {
                        fprintf(stderr, "Specified dataset name %s does not match with dataset platform file\n",
                                col1[j].c_str());
                        return false;
                    }
                    attr->m_mapstrintDatasetSize[col1[j]] = (utype) atoi(col2[j].c_str());
                }
                for (j = 0; j < vD.size(); j++) {
                    if (attr->m_mapstrintDatasetSize.find(vD[j]) == attr->m_mapstrintDatasetSize.end()) {
                        fprintf(stderr, "There is no dataset size entry for %s\n", vD[j].c_str());
                        return false;
                    }
                }
            }

            if (vecDBSetting[i]->platDir != "NA") {
                // Load the new database platform statistics
                SeekPlatforms db_platforms;
                string platformDir = vecDBSetting[i]->platDir;
                assert(!platformDir.empty() && platformDir != "NA");
                db_platforms.loadPlatformDataFromFiles(platformDir);
                // Combine the new db platform statistics with the main db stats
                attr->m_seekPlatforms.combineWithPlatform(db_platforms);
            } else {
                m_missingInitParams = true;
            }
        }

        for (i = 0; i < attr->m_vecstrDatasets.size(); i++) {
            attr->m_mapstrstrDatasetPlatform[attr->m_vecstrDatasets[i]] = attr->m_vecstrDP[i];
            attr->m_mapstrintDataset[attr->m_vecstrDatasets[i]] = i;
        }

        m_iDatasets = attr->m_vecstrDatasets.size();
        m_iGenes = attr->m_vecstrGenes.size();

        for (i = 0; i < vecDBSetting.size(); i++) {
            if (vecDBSetting[i]->dbDir != "NA") {
                bool res;
                res = m_vecDB[i]->Open(vecDBSetting[i]->dbDir,
                                    attr->m_vecstrGenes, attr->m_vecDBDataset[i].size(),
                                    vecDBSetting[i]->GetNumDB());
                if (res == false) {
                    return false;
                }
            } else {
                m_missingInitParams = true;
            }
        }

        if (!m_missingInitParams) {
            CSeekTools::LoadDatabase(m_vecDB, m_iGenes, m_iDatasets,
                                    vecDBSetting, attr->m_vecstrDatasets, attr->m_mapstrstrDatasetPlatform,
                                    attr->m_seekPlatforms.getPlatformMap(), attr->m_seekPlatforms.getCSeekPlatforms(),
                                    m_vc, attr->m_vecDBDataset, attr->m_mapstrintDataset,
                                    bVariance, bCorrelation);
        }

        this->roAttr = attr;
        return true;
    }

    // Version of Initialize called by SeekMiner - calls in-common intialize() and then does query-specific prep
    bool CSeekCentral::InitializeFromSeekMiner(
            const vector<CSeekDBSetting *> &vecDBSetting,
            const char *search_dset, const char *query,
            const char *output_dir, const utype buffer, const bool to_output_text,
            const bool bOutputWeightComponent, const bool bSimulateWeight,
            const enum CSeekDataset::DistanceMeasure dist_measure, const bool bVariance,
            const bool bSubtractAvg, const bool bNormPlatform,
            const bool bLogit, const float fCutOff,
            const float fPercentQueryRequired, const float fPercentGenomeRequired,
            const bool bSquareZ, const bool bRandom, const int iNumRandom,
            const bool bNegativeCor, const bool bCheckDsetSize,
            gsl_rng *rand, const bool useNibble, const int numThreads,
            const float percentDatasetCoverage) {

        if (!CSeekCentral::Initialize(vecDBSetting, buffer, to_output_text,
                                      bOutputWeightComponent, bSimulateWeight, dist_measure, bVariance,
                                      bSubtractAvg, bNormPlatform, bLogit, fCutOff,
                                      fPercentQueryRequired, fPercentGenomeRequired,
                                      bSquareZ, bRandom, iNumRandom, bNegativeCor, bCheckDsetSize,
                                      rand, useNibble, numThreads, percentDatasetCoverage)) {
            return false;
        }

        utype i, j;
        omp_set_num_threads(m_numThreads);
        m_output_dir = output_dir;

        if (m_bCheckDsetSize && roAttr->m_mapstrintDatasetSize.size() == 0) {
            throw config_error("Error: Check_dataset_size requested with no dataset size file configured");
        }

        //fprintf(stderr, "Reading query...\n");
        //read queries
        if (!CSeekTools::ReadMultipleQueries(query, m_vecstrAllQuery))
            return false;

        //fprintf(stderr, "Reading search dataset...\n");
        //read search datasets
        string strSearchDset = search_dset;
        m_vecstrSearchDatasets.resize(m_vecstrAllQuery.size());
        //preparing...
        for (i = 0; i < m_vecstrAllQuery.size(); i++)
            m_vecstrSearchDatasets[i] = vector<string>();

        if (strSearchDset == "NA") {
            for (i = 0; i < m_vecstrAllQuery.size(); i++) {
                for (j = 0; j < roAttr->m_vecstrDatasets.size(); j++) {
                    if (m_bCheckDsetSize) {
                        if (roAttr->m_mapstrintDatasetSize.count(roAttr->m_vecstrDatasets[j]) == 0 ||
                            roAttr->m_mapstrintDatasetSize.at(roAttr->m_vecstrDatasets[j]) <
                                            m_iNumSampleRequired) {
                            continue;
                        }
                    }
                    m_vecstrSearchDatasets[i].push_back(roAttr->m_vecstrDatasets[j]);
                }
            }
        } else {
            vector <vector<string>> vecsearchDset;
            if (!CSeekTools::ReadMultipleQueries(search_dset, vecsearchDset)) {
                fprintf(stderr, "Error reading search datasets\n");
                return false;
            }
            if (vecsearchDset.size() != m_vecstrAllQuery.size()) {
                fprintf(stderr, "Search_dset file doesn't have enough lines. Remember 1 line / query!\n");
                return false;
            }
            for (i = 0; i < m_vecstrAllQuery.size(); i++) {
                for (j = 0; j < vecsearchDset.size(); j++) {
                    if (m_bCheckDsetSize) {
                        if (roAttr->m_mapstrintDatasetSize.count(vecsearchDset[i][j]) == 0 ||
                            roAttr->m_mapstrintDatasetSize.at(vecsearchDset[i][j]) <
                                            m_iNumSampleRequired) {
                            continue;
                        }
                    }
                    m_vecstrSearchDatasets[i].push_back(vecsearchDset[i][j]);
                }
            }
        }

        //free up unnecessary space
        for (i = 0; i < m_vecstrAllQuery.size(); i++)
            m_vecstrSearchDatasets[i].resize(m_vecstrSearchDatasets[i].size());

        //fprintf(stderr, "Making search dset map...\n");
        m_searchdsetMap.resize(m_vecstrAllQuery.size());
        for (i = 0; i < m_vecstrAllQuery.size(); i++) {
            m_searchdsetMap[i] = new CSeekIntIntMap(roAttr->m_vecstrDatasets.size());
            for (j = 0; j < m_vecstrSearchDatasets[i].size(); j++) {
                if (roAttr->m_mapstrintDataset.count(m_vecstrSearchDatasets[i][j]) > 0) {
                    m_searchdsetMap[i]->Add(
                            roAttr->m_mapstrintDataset.at(m_vecstrSearchDatasets[i][j]));
                }
            }
        }


        //fprintf(stderr, "Calculate restart...\n");

        if (!CalculateRestart()) {
            fprintf(stderr, "Error occurred during CalculateRestart()\n");
            return false;
        }
        //fprintf(stderr, "Finished CalculateRestart()\n");


        if (!CheckDatasets(true)) { //replace parameter is true
            fprintf(stderr, "Error occurred during CheckDatasets()\n");
            return false;
        }

        return true;
    }

    bool CSeekCentral::PrepareQuery(const vector <string> &vecstrQuery,
                                    CSeekQuery &query) {
        vector <utype> queryGenes;
        utype j;
        for (j = 0; j < vecstrQuery.size(); j++) {
            if (roAttr->m_mapstrintGene.find(vecstrQuery[j]) ==
                roAttr->m_mapstrintGene.end())
                continue;
            //size_t m = m_DB->GetGene(vecstrQuery[j]);
            //if(m==-1) continue;
            //queryGenes.push_back(m);
            queryGenes.push_back(roAttr->m_mapstrintGene.at(vecstrQuery[j]));
        }
        queryGenes.resize(queryGenes.size());
        query.InitializeQuery(queryGenes, m_iGenes);

        return true;
    }

    bool CSeekCentral::PrepareOneQuery(CSeekQuery &query,
                                       CSeekIntIntMap &dMap, vector<float> &weight) {

        assert(m_master_rank_threads == NULL && m_counts_threads == NULL &&
               m_sum_weight_threads == NULL && m_sum_sq_weight_threads == NULL);
        assert(m_rank_normal_threads == NULL && m_rank_threads == NULL);
        assert(m_rData == NULL);

        utype j;
        const vector <utype> &queryGenes = query.GetQuery();
        const vector <utype> &allRDatasets = dMap.GetAllReverse();
        utype iSearchDatasets = dMap.GetNumSet();
        utype iQuery = queryGenes.size();

        for (j = 0; j < iSearchDatasets; j++)
            m_vc[allRDatasets[j]]->InitializeQuery(queryGenes);

        m_rData = new utype **[m_numThreads];
        for (j = 0; j < m_numThreads; j++)
            m_rData[j] = CSeekTools::Init2DArray(m_iGenes, iQuery, (utype) 0);

        m_master_rank_threads =
                CSeekTools::Init2DArray(m_numThreads, m_iGenes, (float) 0);
        m_sum_weight_threads =
                CSeekTools::Init2DArray(m_numThreads, m_iGenes, (float) 0);
        m_sum_sq_weight_threads =
                CSeekTools::Init2DArray(m_numThreads, m_iGenes, (float) 0);
        m_counts_threads =
                CSeekTools::Init2DArray(m_numThreads, m_iGenes, (utype) 0);

        m_rank_normal_threads = new vector<utype>[m_numThreads];
        m_rank_threads = new vector<utype>[m_numThreads];

        for (j = 0; j < m_numThreads; j++) {
            m_rank_normal_threads[j].resize(m_iGenes);
            m_rank_threads[j].resize(m_iGenes);
            //CSeekTools::InitVector(m_rank_normal_threads[j], m_iGenes, (utype) 255);
            //CSeekTools::InitVector(m_rank_threads[j], m_iGenes, (utype) 255);
        }

        CSeekTools::InitVector(m_master_rank, m_iGenes, (float) 0);
        CSeekTools::InitVector(m_sum_weight, m_iGenes, (float) 0);
        CSeekTools::InitVector(m_sum_sq_weight, m_iGenes, (float) 0);
        CSeekTools::InitVector(m_counts, m_iGenes, (utype) 0);
        CSeekTools::InitVector(weight, m_iDatasets, (float) 0);

        return true;
    }

    bool CSeekCentral::AggregateThreads() {
        assert(m_master_rank_threads != NULL && m_counts_threads != NULL &&
               m_sum_sq_weight_threads != NULL && m_sum_weight_threads != NULL);
        assert(m_rank_normal_threads != NULL && m_rank_threads != NULL);

        //Aggregate into three vectors: m_master_rank, m_counts, m_sum_weight, m_sum_sq_weight
        utype j, k;
        for (j = 0; j < m_numThreads; j++) {
            for (k = 0; k < m_iGenes; k++) {
                m_master_rank[k] += m_master_rank_threads[j][k];
                m_counts[k] += m_counts_threads[j][k];
                m_sum_weight[k] += m_sum_weight_threads[j][k];
                m_sum_sq_weight[k] += m_sum_sq_weight_threads[j][k];
            }
        }

        CSeekTools::Free2DArray(m_master_rank_threads);
        CSeekTools::Free2DArray(m_counts_threads);
        CSeekTools::Free2DArray(m_sum_weight_threads);
        CSeekTools::Free2DArray(m_sum_sq_weight_threads);
        m_master_rank_threads = NULL;
        m_counts_threads = NULL;
        m_sum_weight_threads = NULL;
        m_sum_sq_weight_threads = NULL;

        for (j = 0; j < m_numThreads; j++) {
            m_rank_normal_threads[j].clear();
            m_rank_threads[j].clear();
        }

        delete[] m_rank_normal_threads;
        delete[] m_rank_threads;
        m_rank_normal_threads = NULL;
        m_rank_threads = NULL;

        return true;
    }

    bool CSeekCentral::FilterResults(const utype &iSearchDatasets) {
        utype j, k;

        uint32_t minPercentDatasets = (int)(m_fPercentDatasetCoverage * iSearchDatasets);

        bool DEBUG = false;
        if (DEBUG) {
            fprintf(stderr, "Aggregating genes\n");
            fprintf(stderr, "Searching %d datasets\n", iSearchDatasets);
	    fprintf(stderr, "Min datasets %d\n", minPercentDatasets);
        }



        for (j = 0; j < m_iGenes; j++) {
            //TO DO: make K=(int)(0.5*iSearchDatasets) a customizable parameter
            //TO DO: perhaps it is better to use K=(int)(0.5*(max of m_counts[]))??
            if (m_counts[j] < minPercentDatasets && m_counts[j] < 50) 
                m_master_rank[j] = m_DEFAULT_NA;
            else if (m_sum_weight[j] == 0)
                m_master_rank[j] = m_DEFAULT_NA;
            else {
                m_master_rank[j] =
                        (m_master_rank[j] - 320 * m_sum_weight[j]) / 100.0 / m_sum_weight[j];
                if (m_eDistMeasure == CSeekDataset::CORRELATION) {
                    m_master_rank[j] = m_master_rank[j] / 3.0;
                }
            }
            if (DEBUG) fprintf(stderr, "Gene %d %.5f\n", j, m_master_rank[j]);
        }
        return true;
    }

    bool CSeekCentral::getSortedGeneScores(vector <AResultFloat> &final) {
        if (DEBUG) fprintf(stderr, "Sorting genes\n");
        final.resize(m_iGenes);
        utype j;
        for (j = 0; j < m_iGenes; j++) {
            //fprintf(stderr, "%d %s\n", j, DB.GetGene((size_t) j).c_str());
            final[j].i = j;
            final[j].f = m_master_rank[j];
        }
        if (DEBUG) fprintf(stderr, "Begin Sorting genes\n");
        if (m_bNegativeCor) {
            sort(final.begin(), final.end(), AscendingFloat());
        } else {
            sort(final.begin(), final.end());
        }
        return true;
    }

    bool CSeekCentral::getSortedDatasetScores(uint32_t queryIndex, vector <AResultFloat> &final)  {
        final.resize(m_iDatasets);
        utype j;
        for (j = 0; j < m_iDatasets; j++) {
            final[j].i = j;
            final[j].f = m_weight[queryIndex][j];
        }
        sort(final.begin(), final.end());
        return true;
    }

    void CSeekCentral::setGenePairedResult(uint32_t queryIndex, vector <AResultFloat> &sortedGeneScore) {
        vector<StrDoublePair> &geneResult = this->m_geneResults[queryIndex];
        uint32_t numGenes = sortedGeneScore.size();
        geneResult.resize(numGenes);
        int i;
        for (i = 0; i < numGenes; i++) {
            double geneScore = sortedGeneScore[i].f;
            if (geneScore == m_DEFAULT_NA) {
                break;
            }
            uint32_t geneId = sortedGeneScore[i].i;
            string geneName = roAttr->m_vecstrGenes[geneId];
            geneResult[i].key = geneName;
            geneResult[i].val = geneScore;
        }
        geneResult.resize(i);
    }

    void CSeekCentral::setDatasetPairedResult(uint32_t queryIndex, vector <AResultFloat> &sortedDatasetWeight) {
        vector<StrDoublePair> &datasetResult = this->m_datasetResults[queryIndex];
        uint32_t numDatasets = sortedDatasetWeight.size();
        datasetResult.resize(numDatasets);
        int i;
        for (i = 0; i < numDatasets; i++) {
            double datasetWeight = sortedDatasetWeight[i].f;
            if (datasetWeight == 0) {
                break;
            }
            uint32_t datasetId = sortedDatasetWeight[i].i;
            string datasetName = roAttr->m_vecstrDatasets[datasetId];
            datasetResult[i].key = datasetName;
            datasetResult[i].val = datasetWeight;
        }
        datasetResult.resize(i);
    }


    bool CSeekCentral::Display(CSeekQuery &query, vector <AResultFloat> &final) {
        if (DEBUG) fprintf(stderr, "Results:\n");
        utype jj, ii;
        const vector<char> &cQuery = query.GetQueryPresence();
        for (ii = 0, jj = 0; jj < 500; ii++) {
            if (cQuery[final[ii].i] == 1) continue;
            //fprintf(stderr, "%s %.5f\n",
            //	m_DB->GetGene((size_t)final[ii].i).c_str(), final[ii].f);
            fprintf(stderr, "%s %.5f\n",
                    roAttr->m_vecstrGenes[(size_t) final[ii].i].c_str(), final[ii].f);
            jj++;
        }
        return true;
    }

    bool CSeekCentral::Write(const utype &i) {
        //assume m_bRandom = false
        char acBuffer[1024];
        if (m_output_dir.length() > 0) {
            sprintf(acBuffer, "%s/%d.query", m_output_dir.c_str(), i);
            CSeekTools::WriteArrayText(acBuffer, m_vecstrAllQuery[i]);

            sprintf(acBuffer, "%s/%d.dweight", m_output_dir.c_str(), i);
            CSeekTools::WriteArray(acBuffer, m_weight[i]);

            sprintf(acBuffer, "%s/%d.gscore", m_output_dir.c_str(), i);
            CSeekTools::WriteArray(acBuffer, m_master_rank);
        }

        //send data to client
        if (m_bEnableNetwork) {
            if (CSeekNetwork::Send(m_iClient, m_weight[i]) == -1) {
                fprintf(stderr, "Error sending message to client\n");
                return false;
            }
            if (CSeekNetwork::Send(m_iClient, m_master_rank) == -1) {
                fprintf(stderr, "Error sending message to client\n");
                return false;
            }
        }

        if (!m_bRandom && m_bOutputText) {
            const vector <utype> &allRDatasets =
                    m_searchdsetMap[i]->GetAllReverse();
            utype iSearchDatasets = m_searchdsetMap[i]->GetNumSet();
            vector <vector<string>> vecOutput;
            vecOutput.resize(2);
            vecOutput[0] = vector<string>();
            vecOutput[1] = vector<string>();
            sprintf(acBuffer, "%s/%d.results.txt", m_output_dir.c_str(), i);
            // Get sorted dataset weights
            vector <AResultFloat> w;
            getSortedDatasetScores(i, w);
            for (int j = 0; j < 200 && j < iSearchDatasets; j++) {
                if (w[j].f == 0) break;
                vecOutput[0].push_back(roAttr->m_vecstrDatasets[w[j].i]);
            }
            // Get sorted gene scores
            vector <AResultFloat> wd;
            getSortedGeneScores(wd);
            for (int j = 0; j < 2000 && j < wd.size(); j++) {
                if (wd[j].f == m_DEFAULT_NA) break;
                vecOutput[1].push_back(roAttr->m_vecstrGenes[wd[j].i]);
            }
            if (m_output_dir.length() > 0) {
                CSeekTools::Write2DArrayText(acBuffer, vecOutput);
            }
        }
        return true;
    }

    int CSeekCentral::GetMaxGenomeCoverage() {
        utype d;
        int max = 0;
        for (d = 0; d < roAttr->m_vecstrDatasets.size(); d++) {
            CSeekIntIntMap *mapG = m_vc[d]->GetGeneMap();
            if (mapG->GetNumSet() > max) {
                max = mapG->GetNumSet();
            }
        }
        return max;
    }

    void CSeekCentral::PrintSettings() {
        cout << "SeekCentral: " << endl;
        cout << "bRandom: " << m_bRandom << endl;
        cout << "iNumRandom: " << m_iNumRandom << endl;
        cout << "bSubtractGeneAvg: " << m_bSubtractGeneAvg << endl;
        cout << "bNormPlatform: " << m_bNormPlatform << endl;
        cout << "eDistMeasure: " << m_eDistMeasure << endl;
        cout << "bLogit: " << m_bLogit << endl;
        cout << "bSquareZ: " << m_bSquareZ << endl;
        cout << "Datasets: " << m_iDatasets << endl;
        cout << "iGenes: " << m_iGenes << endl;
        cout << "numThreads: " << m_numThreads << endl;

        cout << "maxNumDB: " << m_maxNumDB << endl;

        cout << "bOutputWeightComponent: " << m_bOutputWeightComponent << endl;
        cout << "bSimulateWeight: " << m_bSimulateWeight << endl;
        cout << "fScoreCutOff: " << m_fScoreCutOff << endl;
        cout << "fPercentQueryAfterScoreCutOff: " << m_fPercentQueryAfterScoreCutOff << endl;
        cout << "fPercentGenomeRequired: " << m_fPercentGenomeRequired << endl;
        cout << "bNegativeCor: " << m_bNegativeCor << endl;

        cout << "useNibble: " << m_useNibble << endl;

        cout << "DEFAULT_NA: " << m_DEFAULT_NA << endl;

        /* for specifying dataset size */
        cout << "bCheckDsetSize: " << m_bCheckDsetSize << endl;
        cout << "iNumSampleRequired: " << m_iNumSampleRequired << endl;
    }

    bool CSeekCentral::Common(CSeekCentral::SearchMode &sm,
                              gsl_rng *rnd, const CSeekQuery::PartitionMode *PART_M,
                              const utype *FOLD, const float *RATE,
                              const vector <vector<float>> *providedWeight,
                              const vector <vector<string>> *newGoldStd) {

        utype i, j, d, dd;
        int k; //keeps track of genes (for random case)
        utype l; //keeps track of random repetition (for random case)
        char acBuffer[1024];

        assert(m_missingInitParams == false);
        // PrintSettings();

        m_Query.resize(m_vecstrAllQuery.size());
        m_weight.resize(m_vecstrAllQuery.size());
        m_final.resize(m_vecstrAllQuery.size());
        m_geneResults.resize(m_vecstrAllQuery.size());
        m_datasetResults.resize(m_vecstrAllQuery.size());

        //random-ranking case =========================
        vector <vector<float>> vecRandWeight, vecRandScore;
        vecRandWeight.resize(m_iNumRandom);
        vecRandScore.resize(m_iNumRandom);
        for (l = 0; l < m_iNumRandom; l++) {
            CSeekTools::InitVector(vecRandWeight[l], m_iDatasets, (float) 0);
            CSeekTools::InitVector(vecRandScore[l], m_iGenes, (float) 0);
        }

        //NEED TO MAKE THIS A PARAMETER
        bool simulateWeight = m_bSimulateWeight;

        //output weight component (Mar 19)
        bool weightComponent = m_bOutputWeightComponent;

        l = 0;
        //oct 20, 2012: whether to redo current query with equal weighting
        int redoWithEqual = 0; //tri-mode: 0, 1, 2
        CSeekCentral::SearchMode current_sm;
        CSeekQuery equalWeightGold;

        //backup of scores (Feb 3)
        vector<float> backupScore;
        CSeekTools::InitVector(backupScore, m_iGenes, (float) 0);

        //fprintf(stderr, "0 %lu\n", CMeta::GetMemoryUsage());
        current_sm = sm;

        int maxGCoverage = GetMaxGenomeCoverage();

        //fprintf(stderr, "Min gene required %.2f %d %d\n", m_fPercentGenomeRequired,
        //	maxGCoverage, (int)(m_fPercentGenomeRequired*(float) maxGCoverage));

        // m_vecstrAllQuery is vector of queries
        for (i = 0; i < m_vecstrAllQuery.size(); i++) {
            //simulated weight case ======================
            /*if(simulateWeight && redoWithEqual>=1) //1 or 2
                current_sm = EQUAL;
            else //0
                current_sm = sm;*/
            //============================================

            g_CatSleipnir().info("Processing Seek subquery: %d", i);

            if (m_mapLoadTime.find(i) != m_mapLoadTime.end()) {
                if (!m_bRandom || l == 0) { //l==0: first random repetition
                    // load query genes into m_vc a vector of CSeekDatasets
                    CSeekTools::ReadDatabaselets(m_vecDB, m_iGenes, m_iDatasets,
                                                 m_mapLoadTime[i], m_vc, roAttr->m_mapstrintGene, roAttr->m_vecDBDataset,
                                                 roAttr->m_mapstrintDataset, m_iClient, m_bEnableNetwork);
                }
            }

            const vector <utype> &allRDatasets =
                    m_searchdsetMap[i]->GetAllReverse();
            utype iSearchDatasets = m_searchdsetMap[i]->GetNumSet();

            //fprintf(stderr, "1 %lu\n", CMeta::GetMemoryUsage());

            CSeekQuery &query = m_Query[i];
            vector<float> &weight = m_weight[i];
            vector <AResultFloat> &final = m_final[i];

            PrepareQuery(m_vecstrAllQuery[i], query);
            PrepareOneQuery(query, *(m_searchdsetMap[i]), weight);
            utype iQuery = query.GetQuery().size();

            //fprintf(stderr, "1b %lu\n", CMeta::GetMemoryUsage());

            //for CV_CUSTOM
            CSeekQuery customGoldStd;
            if (current_sm == CV_CUSTOM)
                PrepareQuery((*newGoldStd)[i], customGoldStd);

            if (current_sm == CV || current_sm == CV_CUSTOM)
                query.CreateCVPartitions(rnd, *PART_M, *FOLD);

            if (current_sm == ORDER_STATISTICS)
                m_rank_d = CSeekTools::Init2DArray(iSearchDatasets, m_iGenes,
                                                   (utype) 0);

            //For outputing component weights!
            vector<float> wc;
            if (weightComponent) {
                if (current_sm == CV || current_sm == CV_CUSTOM) {
                    wc.resize((int) query.GetNumFold() * (int) m_iDatasets);
                } else {
                    wc.resize((int) query.GetQuery().size() * (int) m_iDatasets);
                }
                fill(wc.begin(), wc.end(), (float) 0);
            }
            //fprintf(stderr, "2 %lu\n", CMeta::GetMemoryUsage());

#pragma omp parallel for \
        shared(customGoldStd) \
        private(dd, d, j) \
        firstprivate(i, iSearchDatasets, iQuery) \
        schedule(dynamic)
            for (dd = 0; dd < iSearchDatasets; dd++) {
                d = allRDatasets[dd];
                utype tid = omp_get_thread_num();
                if (DEBUG)
                    fprintf(stderr, "Dataset %d, %s\n",
                            d, roAttr->m_vecstrDatasets[d].c_str());

                CSeekIntIntMap *mapG = m_vc[d]->GetGeneMap();
                CSeekIntIntMap *mapQ = m_vc[d]->GetQueryMap();

                //if dataset contains less than required number of genes, skip
                if (mapG->GetNumSet() < (int) (m_fPercentGenomeRequired *
                                               (float) maxGCoverage)) { //10000
                    continue;
                }

                if (mapQ == NULL || mapQ->GetNumSet() == 0) {
                    if (DEBUG) fprintf(stderr, "This dataset is skipped\n");
                    continue;
                }

                vector <utype> this_q;
                for (j = 0; j < mapQ->GetNumSet(); j++)
                    this_q.push_back(mapQ->GetReverse(j));

                if (DEBUG)
                    fprintf(stderr, "Initializing %d\n",
                            (int) this_q.size());
                m_vc[d]->InitializeDataMatrix(m_rData[tid], m_quant, m_iGenes,
                                              iQuery, m_bSubtractGeneAvg, m_bNormPlatform, m_bLogit,
                                              m_eDistMeasure, m_fScoreCutOff, m_bRandom, m_randRandom);

                float w = -1;
                float report_w = -1; //for showing weight of dataset

                if (current_sm == CV || current_sm == CV_CUSTOM) {
                    if (DEBUG) fprintf(stderr, "Weighting dataset\n");
                    if (current_sm == CV)
                        CSeekWeighter::CVWeighting(query, *m_vc[d], *RATE,
                                                   m_fPercentQueryAfterScoreCutOff, m_bSquareZ,
                                                   false,
                                                   &m_rank_threads[tid]); //weighting always based on positive co-expression
                    else
                        CSeekWeighter::CVWeighting(query, *m_vc[d], *RATE,
                                                   m_fPercentQueryAfterScoreCutOff, m_bSquareZ,
                                                   false, &m_rank_threads[tid],
                                                   &customGoldStd); //weighting based on positive correlation

                    if ((w = m_vc[d]->GetDatasetSumWeight()) == -1) {
                        if (DEBUG) fprintf(stderr, "Bad weight\n");
                        continue;
                    }

                    if (weightComponent && current_sm == CV) {
                        utype numFold = query.GetNumFold();
                        float ww;
                        for (j = 0; j < numFold; j++) {
                            //if((ww=m_vc[d]->GetCVWeight(j))==-1) continue;
                            wc[(int) d * (int) numFold + (int) j] = m_vc[d]->GetCVWeight(j);
                        }
                    }
                } else if (current_sm == AVERAGE_Z) {
                    CSeekWeighter::AverageWeighting(query, *m_vc[d],
                                                    m_fPercentQueryAfterScoreCutOff, m_bSquareZ, w,
                                                    false); //weighting based on positive correlation
                    if (w == -1) continue;
                } else if (current_sm == EQUAL && redoWithEqual == 0) {
                    w = 1.0;
                } else if (current_sm == USE_WEIGHT) {
                    w = (*providedWeight)[i][d];
                }
                    //simulated weight case ======================
                else if (current_sm == EQUAL && redoWithEqual == 1) {
                    w = 1.0;
                    if (DEBUG) fprintf(stderr, "Before doing one gene weighting\n");
                    //calculate reported weight here!
                    CSeekWeighter::OneGeneWeighting(query, *m_vc[d], 0.95,
                                                    m_fPercentQueryAfterScoreCutOff, m_bSquareZ,
                                                    &m_rank_threads[tid], &equalWeightGold, m_bNegativeCor);
                    report_w = m_vc[d]->GetDatasetSumWeight();
                }
                //============================================

                if (DEBUG) fprintf(stderr, "Doing linear combination\n");

                const utype MIN_REQUIRED = max((utype) 1, (utype) (
                        m_fPercentQueryAfterScoreCutOff * this_q.size()));
                CSeekWeighter::LinearCombine(m_rank_normal_threads[tid], this_q,
                                             *m_vc[d], m_bNegativeCor, MIN_REQUIRED, m_bSquareZ);

                if (DEBUG)
                    fprintf(stderr,
                            "Adding contribution of dataset %d to master ranking: %.5f\n", d, w);

                utype iGeneSet = mapG->GetNumSet();
                const vector <utype> &allRGenes = mapG->GetAllReverse();
                vector<utype>::const_iterator iterR = allRGenes.begin();
                vector<utype>::const_iterator endR = allRGenes.begin() + iGeneSet;
                vector <utype> &Rank_Normal = m_rank_normal_threads[tid];
                float *Master_Rank = &m_master_rank_threads[tid][0];
                float *Sum_Weight = &m_sum_weight_threads[tid][0];
                float *Sum_Sq_Weight = &m_sum_sq_weight_threads[tid][0];
                utype *Counts = &m_counts_threads[tid][0];

                if (current_sm == ORDER_STATISTICS)
                    for (; iterR != endR; iterR++) {
                        //if(Rank_Normal[*iterR]==0) continue;
                        m_rank_d[dd][*iterR] = Rank_Normal[*iterR];
                        Counts[*iterR]++;
                    }
                else
                    for (; iterR != endR; iterR++) {
                        //if(Rank_Normal[*iterR]==0) continue;
                        Master_Rank[*iterR] += (float) Rank_Normal[*iterR] * w;
                        Sum_Weight[*iterR] += w;
                        Sum_Sq_Weight[*iterR] += w * w;
                        Counts[*iterR]++;
                    }

                //simulated weight case ======================
                if (current_sm == EQUAL && redoWithEqual == 1) {
                    weight[d] = report_w;
                }
                    //============================================
                else if ((current_sm == EQUAL || current_sm == ORDER_STATISTICS)
                         && redoWithEqual == 0) {
                    weight[d] = 0;
                } else {
                    weight[d] = w;
                }

            }
            //omp finishes
            //fprintf(stderr, "3 %lu\n", CMeta::GetMemoryUsage());
            for (j = 0; j < iSearchDatasets; j++)
                m_vc[allRDatasets[j]]->DeleteQuery();

            assert(m_rData != NULL);
            for (j = 0; j < m_numThreads; j++)
                CSeekTools::Free2DArray(m_rData[j]);
            delete[] m_rData;
            m_rData = NULL;

            AggregateThreads();

            if (current_sm != ORDER_STATISTICS) {
                FilterResults(iSearchDatasets);
            } else {
                CSeekWeighter::OrderStatisticsRankAggregation(iSearchDatasets,
                                                              m_iGenes, m_rank_d, m_counts, m_master_rank, m_numThreads,
                                                              m_bNegativeCor);
                CSeekTools::Free2DArray(m_rank_d);
                m_rank_d = NULL;
            }

            //Display(query, final);
            //fprintf(stderr, "4 %lu\n", CMeta::GetMemoryUsage());
            // SetQueryScoreNull(query);  // commented out 6/2016 by qzhu
            getSortedGeneScores(final);
            int ret; //for system calls

            if (m_bRandom) {
                /*utype z, cz;
                for(cz=0, z=0; z<m_iDatasets; z++)
                    if(m_weight[i][z]!=0)
                        cz++;
                fprintf(stderr, "Number of weighted dataset: %d\n", cz);
                */
            } else if (simulateWeight) {
                if ((current_sm == EQUAL || current_sm == ORDER_STATISTICS) && !CheckWeight(i)) {
                    g_CatSleipnir().debug("Calculate dataset ordering");
                    // ret = system("date +%s%N 1>&2");
                    if (m_bEnableNetwork) {
                        if (CSeekNetwork::Send(m_iClient, "Calculate dataset ordering") == -1) {
                            fprintf(stderr, "Error sending message to client\n");
                        }
                    } else if (m_useRPC) {
                        m_rpcLog->enqueue("Calculate dataset ordering");
                    }
                    CopyTopGenes(equalWeightGold, final, 100);
                    redoWithEqual = 1;
                    if (current_sm == ORDER_STATISTICS) {
                        current_sm = EQUAL;
                        //backup genes to old
                        copy(m_master_rank.begin(), m_master_rank.end(), backupScore.begin());
                    }
                    i--;
                    continue;
                } else if (current_sm == CV && !CheckWeight(i)) {
                    g_CatSleipnir().info("CV: Redo with equal weighting");
                    // ret = system("date +%s%N 1>&2");
                    if (m_bEnableNetwork) {
                        if (CSeekNetwork::Send(m_iClient, "Redo with equal weighting") == -1) {
                            fprintf(stderr, "Error sending message to client\n");
                        }
                    } else if (m_useRPC) {
                        m_rpcLog->enqueue("Redo with equal weighting");
                    }
                    current_sm = EQUAL;
                    i--;
                    continue;
                } else if (current_sm == EQUAL && CheckWeight(i)) {
                    redoWithEqual = 0;
                    current_sm = sm;
                    //copy genes from old to new
                    if (sm == ORDER_STATISTICS) {
                        copy(backupScore.begin(), backupScore.end(), m_master_rank.begin());
                    }
                }
            }

            // fprintf(stderr, "Done search\n");
            // ret = system("date +%s%N 1>&2");

            if (m_bEnableNetwork) {
                if (CSeekNetwork::Send(m_iClient, "Done Search") == -1) {
                    fprintf(stderr, "Error sending message to client\n");
                }
            } else if (m_useRPC) {
                m_rpcLog->enqueue("Done Search");
            }

            setGenePairedResult(i, final);
            vector <AResultFloat> finalDatasetWeights;
            getSortedDatasetScores(i, finalDatasetWeights);
            setDatasetPairedResult(i, finalDatasetWeights);

            //random-ranking case =========================
            if (m_bRandom) {
                if (m_bNegativeCor) {
                    fprintf(stderr, "Error! Random-ranking case does not support Negative Correlations!\n");
                    continue;
                }

                sort(m_master_rank.begin(), m_master_rank.end(), greater<float>());
                sort(weight.begin(), weight.end(), greater<float>());
                copy(m_master_rank.begin(), m_master_rank.end(), vecRandScore[l].begin());
                copy(weight.begin(), weight.end(), vecRandWeight[l].begin());
                l++;

                if (l == m_iNumRandom) { //last repetition
                    float confidence = 0.50;
                    int n = (int) (confidence * (float) m_iNumRandom);
                    vector<float> new_weight, new_score;
                    CSeekTools::InitVector(new_weight, m_iDatasets, (float) 0);
                    CSeekTools::InitVector(new_score, m_iGenes, (float) 0);
                    for (k = 0; k < m_iGenes; k++) {
                        vector<float> all_s;
                        all_s.resize(m_iNumRandom);
                        for (l = 0; l < m_iNumRandom; l++)
                            all_s[l] = vecRandScore[l][k];
                        std::nth_element(all_s.begin(), all_s.begin() + n, all_s.end());
                        new_score[k] = all_s[n];
                    }
                    /*for(k=0; k<m_iGenes; k++){
                        fprintf(stderr, "%d %.5f\n", k, new_score[k]);
                    }
                    getchar();*/
                    for (k = 0; k < m_iDatasets; k++) {
                        vector<float> all_w;
                        all_w.resize(m_iNumRandom);
                        for (l = 0; l < m_iNumRandom; l++)
                            all_w[l] = vecRandWeight[l][k];
                        std::nth_element(all_w.begin(), all_w.begin() + n, all_w.end());
                        new_weight[k] = all_w[n];
                    }
                    if (m_output_dir.length() > 0) {
                        sprintf(acBuffer, "%s/%d.rand.dweight", m_output_dir.c_str(), i);
                        CSeekTools::WriteArray(acBuffer, new_weight);
                        sprintf(acBuffer, "%s/%d.rand.gscore", m_output_dir.c_str(), i);
                        CSeekTools::WriteArray(acBuffer, new_score);
                    }
                    l = 0;
                } else {
                    i--;
                }
            }

            if (!m_bRandom) {
                //if m_bRandom, write at the very end when all repetitions are done
                // TODO GW - pass a tag in to write() to give the files distinct name
                Write(i);
                if (weightComponent) {
                    if (m_output_dir.length() > 0) {
                        sprintf(acBuffer, "%s/%d.dweight_comp", m_output_dir.c_str(), i);
                        CSeekTools::WriteArray(acBuffer, wc);
                    }
                    vector <vector<string>> vecParts;
                    vecParts.resize(query.GetNumFold());
                    utype kk;
                    string strParts = "";
                    for (j = 0; j < query.GetNumFold(); j++) {
                        vecParts[j] = vector<string>();
                        const vector <utype> &vu = query.GetCVQuery(j);
                        string s = "";
                        for (kk = 0; kk < vu.size(); kk++) {
                            vecParts[j].push_back(roAttr->m_vecstrGenes[vu[kk]]);
                            s += roAttr->m_vecstrGenes[vu[kk]];
                            if (kk != vu.size() - 1)
                                s += " ";
                        }
                        strParts += s;
                        if (j != query.GetNumFold() - 1)
                            strParts += "|";
                    }
                    if (m_output_dir.length() > 0) {
                        sprintf(acBuffer, "%s/%d.query_cvpart", m_output_dir.c_str(), i);
                        CSeekTools::Write2DArrayText(acBuffer, vecParts);
                    }
                    //send data to client
                    if (m_bEnableNetwork) {
                        if (CSeekNetwork::Send(m_iClient, wc) == -1) {
                            fprintf(stderr, "Error sending message to client\n");
                            return false;
                        }
                        if (CSeekNetwork::Send(m_iClient, strParts) == -1) {
                            fprintf(stderr, "Error sending message to client\n");
                            return false;
                        }
                    }
                }
            }
        }

        //fprintf(stderr, "4b %lu\n", CMeta::GetMemoryUsage());

        return true;
    }

    bool CSeekCentral::CheckWeight(const utype &i) {
        utype j = 0;
        bool valid = false;
        for (j = 0; j < m_iDatasets; j++) {
            if (m_weight[i][j] != 0) {
                valid = true;
                break;
            }
        }
        return valid;
    }

    bool CSeekCentral::CopyTopGenes(CSeekQuery &csq,
                                    const vector <AResultFloat> &src, const utype top) {
        utype i, j;
        vector <utype> topGenes;
        for (i = 0; i < top; i++) {
            if (src[i].f == m_DEFAULT_NA) continue;
            topGenes.push_back(src[i].i);
        }
        if (topGenes.size() == 0) {
            fprintf(stderr, "Error in CopyTopGenes!\n");
            return false;
        }
        csq.InitializeQuery(topGenes, m_iGenes);
        return true;
    }

    bool CSeekCentral::SetQueryScoreNull(const CSeekQuery &csq) {
        utype j;
        const vector <utype> &query = csq.GetQuery();
        for (j = 0; j < query.size(); j++) {
            m_master_rank[query[j]] = m_DEFAULT_NA;
        }
        return true;
    }

    bool CSeekCentral::EqualWeightSearch() {
        CSeekCentral::SearchMode sm = EQUAL;
        return CSeekCentral::Common(sm);
    }

    bool CSeekCentral::CVSearch(gsl_rng *rnd, const CSeekQuery::PartitionMode &PART_M,
                                const utype &FOLD, const float &RATE) {
        CSeekCentral::SearchMode sm = CV;
        return CSeekCentral::Common(sm, rnd, &PART_M, &FOLD, &RATE);
    }

/*	perform CVSearch, except that the weighting is based not on co-expression
	of query genes, but based on similarity of query genes to some custom gold
	standard gene-set */
    bool CSeekCentral::CVCustomSearch(const vector <vector<string>> &newGoldStd,
                                      gsl_rng *rnd, const CSeekQuery::PartitionMode &PART_M,
                                      const utype &FOLD, const float &RATE) {
        CSeekCentral::SearchMode sm = CV_CUSTOM;
        return CSeekCentral::Common(sm, rnd, &PART_M, &FOLD, &RATE,
                                    NULL, &newGoldStd);
    }

    bool CSeekCentral::WeightSearch(const vector <vector<float>> &weights) {
        CSeekCentral::SearchMode sm = USE_WEIGHT;
        return CSeekCentral::Common(sm, NULL, NULL, NULL, NULL, &weights);
    }

    bool CSeekCentral::OrderStatistics() {
        CSeekCentral::SearchMode sm = ORDER_STATISTICS;
        return CSeekCentral::Common(sm);
    }

    bool CSeekCentral::AverageWeightSearch() {
        CSeekCentral::SearchMode sm = AVERAGE_Z;
        return CSeekCentral::Common(sm, NULL, NULL, NULL, NULL);
    }

    bool CSeekCentral::VarianceWeightSearch() {
        vector <vector<float>> weights;
        weights.resize(m_vecstrAllQuery.size());
        utype i, j, k;
        for (i = 0; i < m_vecstrAllQuery.size(); i++) {
            weights[i] = vector<float>();
            CSeekTools::InitVector(weights[i], m_iDatasets, (float) 0);
            const vector <utype> &allRDatasets =
                    m_searchdsetMap[i]->GetAllReverse();
            utype iSearchDatasets = m_searchdsetMap[i]->GetNumSet();

            CSeekQuery q;
            PrepareQuery(m_vecstrAllQuery[i], q);
            const vector <utype> &qGenes = q.GetQuery();

            for (j = 0; j < iSearchDatasets; j++) {
                utype d = allRDatasets[j];
                for (k = 0; k < qGenes.size(); k++) {
                    float gv = m_vc[d]->GetGeneVariance(qGenes[k]);
                    if (CMeta::IsNaN(gv)) continue;
                    weights[i][d] += gv;
                }
            }
        }
        return WeightSearch(weights);
    }

    bool CSeekCentral::Destruct() {
        int i;
        for (i = 0; i < m_vc.size(); i++) {
            if (m_vc[i] == NULL) continue;
            delete m_vc[i];
            m_vc[i] = NULL;
        }
        m_vc.clear();

        for (i = 0; i < m_searchdsetMap.size(); i++) {
            if (m_searchdsetMap[i] == NULL) continue;
            delete m_searchdsetMap[i];
            m_searchdsetMap[i] = NULL;
        }
        m_searchdsetMap.clear();

        for (i = 0; i < m_vecDB.size(); i++) {
            if (m_vecDB[i] == NULL) continue;
            delete m_vecDB[i];
            m_vecDB[i] = NULL;
        }
        m_vecDB.clear();

        if (m_master_rank_threads != NULL) {
            CSeekTools::Free2DArray(m_master_rank_threads);
            m_master_rank_threads = NULL;
        }
        if (m_sum_weight_threads != NULL) {
            CSeekTools::Free2DArray(m_sum_weight_threads);
            m_sum_weight_threads = NULL;
        }
        if (m_sum_sq_weight_threads != NULL) {
            CSeekTools::Free2DArray(m_sum_sq_weight_threads);
            m_sum_sq_weight_threads = NULL;
        }
        if (m_counts_threads != NULL) {
            CSeekTools::Free2DArray(m_counts_threads);
            m_counts_threads = NULL;
        }
        if (m_rank_normal_threads != NULL) {
            for (i = 0; i < m_numThreads; i++)
                m_rank_normal_threads[i].clear();
            delete[] m_rank_normal_threads;
            m_rank_normal_threads = NULL;
        }
        if (m_rank_threads != NULL) {
            for (i = 0; i < m_numThreads; i++)
                m_rank_threads[i].clear();
            delete[] m_rank_threads;
            m_rank_threads = NULL;
        }

        //m_rank_d is for Order Statistics aggregation only! Should be null already.
        //note that m_rData already deleted after search
        if (m_rank_d != NULL){
            CSeekTools::Free2DArray(m_rank_d);
            m_rank_d = NULL;
        }

        return true;
    }

    const vector <vector<AResultFloat>> &CSeekCentral::GetAllResult() const {
        return m_final;
    }

    const vector <CSeekQuery> &CSeekCentral::GetAllQuery() const {
        return m_Query;
    }

    utype CSeekCentral::GetGeneId(const string &strGene) const {
        if (roAttr->m_mapstrintGene.find(strGene) == roAttr->m_mapstrintGene.end())
            return CSeekTools::GetNaN();
        return roAttr->m_mapstrintGene.find(strGene)->second;
    }

    string CSeekCentral::GetGeneName(const utype &geneID) const {
        return roAttr->m_vecstrGenes[(size_t) geneID];
    }

    const vector <vector<float>> &CSeekCentral::GetAllWeight() const {
        return m_weight;
    }

    void CSeekCentral::convertGenesEntrezToSymbol(const vector<string> &entrez, vector<string> &symbols) {
        uint32_t numGenes = entrez.size();
        symbols.resize(numGenes);
        for (int i=0; i<numGenes; i++) {
            try {
                if (roAttr->m_geneEntrezToSymbolMap.count(entrez[i]) > 0) {
                    symbols[i] = roAttr->m_geneEntrezToSymbolMap.at(entrez[i]);
                } else {
                    symbols[i] = entrez[i];
                }
            } catch(exception &err) {
                throw_with_nested(request_error(FILELINE + "geneEntrezToSymbolMap error for id:: " + entrez[i]));
            }
        }
    }

    void CSeekCentral::convertGenesSymbolToEntrez(const vector<string> &symbols, vector<string> &entrez) {
        uint32_t numGenes = symbols.size();
        entrez.resize(numGenes);
        for (int i=0; i<numGenes; i++) {
            try {
                entrez[i] = roAttr->m_geneSymbolToEntrezMap.at(symbols[i]);
            } catch(exception &err) {
                throw_with_nested(request_error(FILELINE + "symbol not found: " + symbols[i]));
            }
        }
    }

    string CSeekCentral::entrezToSymbol(string &entrez) {
        try {
            if (roAttr->m_geneEntrezToSymbolMap.count(entrez) > 0) {
                return roAttr->m_geneEntrezToSymbolMap.at(entrez);
            } else {
                return entrez;
            }
        } catch(exception &err) {
            throw_with_nested(request_error(FILELINE + "geneEntrezToSymbolMap error for id: " + entrez));
        }
    }
}

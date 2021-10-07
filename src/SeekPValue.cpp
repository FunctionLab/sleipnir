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
#include "stdafx.h"
#include "SeekPValue.h"
#include <filesystem>
namespace fs = std::filesystem;


void *do_pvalue_query(void *th_arg) {
    struct pvalue_thread_data *my = (struct pvalue_thread_data *) th_arg;
    /* A lock_guard will automatically set the completionFlag
     *  to true when the scope exits. This will allow other threads
     *  such as the main thread to know the execution has completed.
     */
    BoolFlag completionFlag(my->isComplete);
    lock_guard<BoolFlag> flag_lock(completionFlag);

    int new_fd = my->new_fd;
    float nan = my->nan;
    int rankBased = my->rankBased;
    int queryType = my->queryType;
    CSeekCentral *cc = my->seekCentral;
    PValueData *pvalueData = my->pvalueData;
    int numGenes = cc->m_vecstrGenes.size();

    vector <string> &dset = my->dset;
    vector<float> &dset_score = my->dset_score;
    vector<int> &dset_qsize = my->dset_qsize;

    vector <string> &queryGenes = my->query;
    vector <string> &geneEntrezIds = my->gene_entrezIds;
    vector<double> &geneScores = my->gene_scores;
    vector<int> &geneRanks = my->gene_ranks;
    vector<int> geneIds;

    size_t i, j, jj, k;
    vector<float> pval;

    if (queryType == 1) { //dataset
        for (i = 0; i < dset.size(); i++) {
            if (cc->m_mapstrintDataset.find(dset[i]) == cc->m_mapstrintDataset.end()) {
                fprintf(stderr, "Error: cannot find dataset %s\n", dset[i].c_str());
                pval.push_back(-1);
                continue;
            }
            if (dset_qsize[i] < 2) {
                pval.push_back(0.99);
                continue;
            }
            if (dset_score[i] <= 0) {
                pval.push_back(0.99);
                continue;
            }
            int pi = cc->m_mapstrintDataset[dset[i]];
            float sc = log(dset_score[i]);
            int qsize = dset_qsize[i];
            if (pvalueData->dsetScore[pi].size() == 0) {
                //This dataset is a dummy, null dataset (because it is bad or low quality
                pval.push_back(0.99);
                continue;
            }
            if (qsize - 2 >= pvalueData->dsetScore[pi].size()) {
                //Don't have modeling for this query size, most likely because query size > 100
                pval.push_back(0.99);
                continue;
            }
            struct parameter &par = pvalueData->dsetScore[pi][qsize - 2];
            vector<double> &quantile = par.quantile;
            //fprintf(stderr, "threshold %.2e, scale %.2e, shape %.2e, qsize %d\n", par.threshold, par.scale, par.shape, qsize);

            if ((double) sc <= par.threshold) {
                double min_diff = 9999;
                int min_diff_index = -1;
                /*float p_val = 0.05;
                for(k=quantile.size()-1; k>=0; k--){
                    if( (double) sc > quantile[k] ){
                        break;
                    }else{
                        p_val+=0.05;
                    }
                }
                if(p_val >= 0.99){
                    p_val = 0.99;
                }*/
                int cut_index = 0;
                for (k = 0; k <= quantile.size(); k++) {
                    if (quantile[k] >= par.threshold) {
                        cut_index = k;
                        break;
                    }
                }

                for (k = 0; k <= cut_index; k++) {
                    double diff = -1;
                    if (k == cut_index) {
                        diff = fabs(par.threshold - (double) sc);
                        if (diff < min_diff) {
                            min_diff = diff;
                            min_diff_index = k;
                        }
                        continue;
                    }
                    diff = fabs(quantile[k] - (double) sc);
                    if (diff < min_diff) {
                        min_diff = diff;
                        min_diff_index = k;
                    }
                    //fprintf(stderr, "Diff %d %.2f %.2f %d\n", k, quantile[k], diff, qsize);
                }
                if (min_diff_index == -1) {
                    //fprintf(stderr, "Error, negative one!\n");
                }
                float pv = 1.0 - (0.05 + (float) min_diff_index * 0.01);
                pval.push_back(pv);
                //pval.push_back(p_val);
                //fprintf(stderr, "original %.2e pval %.2e\n", sc, pv);
                continue;
            }

            double limit = -1.0 * par.scale / par.shape;
            limit = floor(limit * 100.0) / 100.0;

            double diff = (double) sc - par.threshold;
            double cdf = 0;

            if (par.shape == 0) {
                cdf = 1.0 - exp(-1.0 * diff / par.scale);
            } else {
                if (par.shape < 0 && diff >= limit) { //check that difference is within range
                    diff = limit; //if difference is bigger than it can accept, assume limit
                }
                cdf = 1.0 - pow(1.0 + par.shape * diff / par.scale, -1.0 / par.shape);
            }
            double pv = par.portion * (1.0 - cdf);
            //fprintf(stderr, "before original %.2e pval %.2e\n", sc, pv);
            if (isnan(pv)) {
                pv = 2.0e-5;
            }
            pval.push_back((float) pv);
            //fprintf(stderr, "original %.2e pval %.2e\n", sc, pv);
        }

    } else if (queryType == 0) { //genes
        try {
            if (rankBased == true) {
                if (geneRanks.size() == 0) {
                    // Sort geneScores and populate the geneRanks array
                    if (geneScores.size() != numGenes) {
                        // To calculate the geneRanks we need the scores for all genes
                        throw runtime_error("PValue: rank-based requested but num geneScores!=numGenes");
                    }
                    // Sort the geneScores to populate the geneRanks
                    vector <AResultFloat> sortedGenes;
                    sortedGenes.resize(geneScores.size());
                    for (i = 0; i < sortedGenes.size(); i++) {
                        sortedGenes[i].i = i;
                        sortedGenes[i].f = geneScores[i];
                    }
                    sort(sortedGenes.begin(), sortedGenes.end());
                    geneIds.resize(sortedGenes.size());
                    geneRanks.resize(sortedGenes.size());
                    for (i = 0; i < sortedGenes.size(); i++) {
                        geneIds[i] = sortedGenes[i].i;
                        geneRanks[i] = i;
                    }
                } else {
                    if (geneRanks.size() != geneEntrezIds.size()) {
                        throw runtime_error("PValue: rank-based requested but num geneRanks!=geneEntrezIds");
                    }
                    geneIds.resize(geneEntrezIds.size());
                    for (i = 0; i < geneIds.size(); i++) {
                        geneIds[i] = cc->m_mapstrintGene[geneEntrezIds[i]];
                    }
                }
                assert(geneRanks.size() == geneIds.size());
            } else { // score based prep
                if (geneEntrezIds.size() == 0) {
                    // assumption is that the geneScores are provided in the gene_map order
                    // and there are as many geneScores as numGenes
                    if (geneScores.size() != numGenes) {
                        throw runtime_error("PValue: score-based, no geneEntrezIds provided, numGeneScores!=numGenes");
                    }
                    geneIds.resize(geneScores.size());
                    for (i = 0; i < geneIds.size(); i++) {
                        geneIds[i] = i;
                    }
                } else { // geneEntrezIds.size() > 0
                    // convert entrezIds to geneIds
                    geneIds.resize(geneEntrezIds.size());
                    for (i=0; i<geneEntrezIds.size(); i++) {
                        geneIds[i] = cc->m_mapstrintGene[geneEntrezIds[i]];
                    }
                }
                assert(geneScores.size() == geneIds.size());
            }

            // TODO: Remove this, not used
            // vector<int> queryGeneID;
            // for (i = 0; i < queryGenes.size(); i++)
            //     queryGeneID.push_back(cc->m_mapstrintGene[queryGenes[i]]);

            //Query genes themselves have lowest score, to prevent
            //them from being counted in PR
            //(disabled 6/6/2016) want the query to have scores
            //for (i = 0; i < queryGeneID.size(); i++)
            //    sortedGenes[queryGeneID[i]].f = nan;


            //comparison
            // TODO - Remove this unused code
            // vector<int> geneRank;
            // geneRank.resize(numGenes);
            // for (jj = 0; jj < numGenes; jj++) {
            //     geneRank[sortedGenes[jj].i] = jj;
            // }

            CSeekTools::InitVector(pval, geneIds.size(), (float) nan);
            if (my->useGeneMapOrder == true && pval.size() != numGenes) {
                throw runtime_error("PValue: useGeneMapOrder result array insufficient size");
            }

            for (jj = 0; jj < pval.size(); jj++) {
                int geneId = geneIds[jj];
                vector<int> &rR = pvalueData->randomRank[geneId];
                vector<float> &rF = pvalueData->randomSc[geneId];
                int pvalIdx = jj;
                if (my->useGeneMapOrder == true) {
                    // Return results in geneMapOrder, put this result in geneId slot
                    pvalIdx = geneId;
                }
                int kk = 0;
                if (rankBased == false) {  // score based
                    // float gene_score = sortedGenes[jj].f;
                    float gene_score = geneScores[jj];
                    if (gene_score == nan) break;
                    if (gene_score >= 0) {
                        for (kk = 0; kk < rF.size(); kk++) {
                            if (gene_score >= rF[kk] || kk == rF.size() - 1)
                                pval[pvalIdx] = (float) kk / (float) rF.size();
                            //fprintf(stderr, "%s\t%d\t%d\t%.5e\t%.5e\n", vecstrGenes[gene].c_str(),
                            //gene_rank, kk, gene_score, randomSc[gene][kk]);
                            if (gene_score >= rF[kk])
                                break;
                        }
                    } else {
                        for (kk = rF.size() - 1; kk >= 0; kk--) {
                            if (gene_score <= rF[kk] || kk == 0)
                                pval[pvalIdx] = (float) (rF.size() - 1 - kk) / (float) rF.size();
                            //fprintf(stderr, "%s\t%d\t%d\t%.5e\t%.5e\n", vecstrGenes[gene].c_str(),
                            //gene_rank, rF.size()-1-kk, gene_score, randomSc[gene][kk]);
                            if (gene_score <= rF[kk])
                                break;
                        }
                    }
                } else if (rankBased == true) {
                    int gene_rank = geneRanks[jj];
                    if (gene_rank < numGenes / 2) {
                        for (kk = 0; kk < rR.size(); kk++) {
                            if (gene_rank <= rR[kk] || kk == rR.size() - 1)
                                pval[pvalIdx] = (float) kk / (float) rR.size();
                            //fprintf(stderr, "%s\t%d\t%d\t%.5e\t%.5e\n", vecstrGenes[gene].c_str(),
                            //gene_rank, kk, gene_score, randomSc[gene][kk]);
                            if (gene_rank <= rR[kk])
                                break;
                        }
                    } else {
                        for (kk = rR.size() - 1; kk >= 0; kk--) {
                            if (gene_rank >= rR[kk] || kk == 0)
                                pval[pvalIdx] = (float) (rR.size() - 1 - kk) / (float) rF.size();
                            //fprintf(stderr, "%s\t%d\t%d\t%.5e\t%.5e\n", vecstrGenes[gene].c_str(),
                            //gene_rank, rR.size()-1-kk, gene_score, randomSc[gene][kk]);
                            if (gene_rank >= rR[kk])
                                break;
                        }
                    }
                }
            }
        } catch (exception &err) {
            std::cerr << "ERROR: " << err.what() << endl;
            my->error = true;
            my->errorMsg = err.what();
        }
    }

    if (new_fd >= 0) {
        if (CSeekNetwork::Send(new_fd, pval) == -1) {
            fprintf(stderr, "Error sending message to client!\n");
        }
        close(new_fd);
    } else {
        my->resPvalues->insert(my->resPvalues->end(),
                                pval.begin(), pval.end());
    }

    return 0;
}


bool ReadParameter(const string& param_file, vector<struct parameter> &v) {
    ifstream ifsm;
    ifsm.open(param_file.c_str());
    if (!ifsm.is_open()) {
        fprintf(stderr, "Error opening file %s\n", param_file.c_str());
        return false;
    }
    const int lineSize = 6084;
    char acBuffer[lineSize];
    utype c_iBuffer = lineSize;
    v.clear();

    while (!ifsm.eof()) {
        ifsm.getline(acBuffer, c_iBuffer - 1);
        if (acBuffer[0] == 0) break;
        acBuffer[c_iBuffer - 1] = 0;
        vector <string> tok;
        CMeta::Tokenize(acBuffer, tok);
        struct parameter par;
        par.size = strtol(tok[0].c_str(), nullptr, 10);
        par.portion = strtod(tok[2].c_str(), nullptr);
        par.threshold = strtod(tok[3].c_str(), nullptr);
        par.scale = strtod(tok[5].c_str(), nullptr);
        par.shape = strtod(tok[6].c_str(), nullptr);
        par.quantile = vector<double>();
        //for(int k=8; k<tok.size(); k++){
        for (int k = 7; k < tok.size(); k++) {
            par.quantile.push_back(strtod(tok[k].c_str(), nullptr));
            //fprintf(stderr, "This value is %d %.2f\n", k, atof(tok[k].c_str()));
        }
        v.push_back(par);
    }
    ifsm.close();
    return true;
}

bool loadPvalueArrays(string dirname, PValueData &pvalueData) {
    fs::path scoreFile = dirname;
    fs::path rankFile = dirname;
    scoreFile /= "randomScoreFile.bin";
    rankFile /= "randomRankFile.bin";
    read2DVector(pvalueData.randomSc, scoreFile);
    read2DVector(pvalueData.randomRank, rankFile);
    return true;
}

bool initializePvalue(CSeekCentral &seekCentral, int numRandQueries, PValueData &pvalueData) {
    int numGenes = seekCentral.m_vecstrGenes.size();
    int numRandFiles = 0;
    vector <string> gscoreFiles;
    string random_directory = seekCentral.m_vecDBSetting[0]->randomDir;
    for (const auto & entry : fs::directory_iterator(random_directory)) {
        if (entry.path().extension() == ".gscore") {
            numRandFiles++;
            gscoreFiles.push_back(entry.path());
            // cout << entry.path() << endl;
        }
    }
    int ii, jj;
    char ac[256];
    int num_random = numRandFiles;

    if (numRandQueries > 0) {
        num_random = numRandQueries;
    }

    pvalueData.randomRank.resize(numGenes);
    pvalueData.randomSc.resize(numGenes);
    for (ii = 0; ii < numGenes; ii++) {
        pvalueData.randomRank[ii].resize(num_random);
        pvalueData.randomSc[ii].resize(num_random);
    }

    for (ii = 0; ii < num_random; ii++) {
        vector<float> randomScores;
        // sprintf(ac, "%s/%d.gscore", random_directory.c_str(), ii);
        cout << gscoreFiles[ii] << endl;
        CSeekTools::ReadArray(gscoreFiles[ii].c_str(), randomScores);
        assert (randomScores.size() == numGenes);
        /*vector<string> queryGenes;
        sprintf(ac, "%s/%d.query", random_directory.c_str(), ii);
        CSeekTools::ReadMultiGeneOneLine(ac, queryGenes);
        querySize.push_back(queryGenes.size());
        */
        vector <AResultFloat> sortedRandom;
        sortedRandom.resize(randomScores.size());
        for (jj = 0; jj < randomScores.size(); jj++) {
            sortedRandom[jj].i = jj;
            sortedRandom[jj].f = randomScores[jj];
        }
        sort(sortedRandom.begin(), sortedRandom.end());
        for (jj = 0; jj < randomScores.size(); jj++) {
            pvalueData.randomRank[sortedRandom[jj].i][ii] = jj;
            pvalueData.randomSc[sortedRandom[jj].i][ii] = sortedRandom[jj].f;
        }
    }

    for (jj = 0; jj < numGenes; jj++) {
        sort(pvalueData.randomRank[jj].begin(), pvalueData.randomRank[jj].end());
        sort(pvalueData.randomSc[jj].begin(), pvalueData.randomSc[jj].end(), std::greater<float>());
    }
    cout << "Done initializing PValue arrays" << endl;

    fs::path scoreFile = random_directory;
    fs::path rankFile = random_directory;
    scoreFile /= "randomScoreFile.bin";
    rankFile /= "randomRankFile.bin";

    write2DVector(pvalueData.randomSc, scoreFile);
    write2DVector(pvalueData.randomRank, rankFile);

// comment out loading dataset parameter file for now
#if 0
    string param_dir = sArgs.param_dir_arg;
    int numDatasets = seekCentral.m_vecstrDatasets.size();
    pvalueData.dsetScore.resize(numDatasets);
    for (i = 0; i < numDatasets; i++) {
        string param_file = param_dir + "/" + seekCentral.m_vecstrDatasets[i] + ".param";
        pvalueData.dsetScore[i] = vector<struct parameter>();
        if (!ReadParameter(param_file, pvalueData.dsetScore[i])) {
            fprintf(stderr, "Making this dataset null... (will always return insignificant)\n");
        }
    }
#endif

    return true;
}

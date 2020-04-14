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
#include "cmdline.h"

#include "statistics.h"
#include "datapair.h"
#include "file.h"

struct GeneStruct {
    size_t idx;
    float score;
};

struct SortGeneStruct {

    bool operator()(const GeneStruct rOne, const GeneStruct rTwo) const {
        return (rOne.score > rTwo.score);
    }
};

class PairStruct {
public:
    size_t g1;
    size_t g2;
    float score;

    PairStruct(size_t g1idx, size_t g2idx) {
        g1 = g1idx;
        g2 = g2idx;
        score = 0;
    }
};

struct SortPairStruct {

    bool operator()(const PairStruct *rOne, const PairStruct *rTwo) const {
        return (rOne->score < rTwo->score);
    }
};


bool open_genepairs(char *szFile, vector<PairStruct *> &vecPairs, CPCL &GDist) {
    ifstream istm;
    const char *pc;
    char *pcTail;
    char *acBuf;
    string strToken, strCache, strValue;
    size_t iOne, iTwo, i;
    const char c_acComment[] = "#";

    istm.open(szFile);

    acBuf = new char[CFile::GetBufferSize()];
    while (istm.peek() != EOF) {
        istm.getline(acBuf, CFile::GetBufferSize() - 1);
        strToken = CFile::OpenToken(acBuf, &pc);
        if (!strToken.length())
            break;
        if (strToken == c_acComment)
            continue;
        if (strToken != strCache) {
            strCache = strToken;

            iOne = GDist.GetGene(strToken);
            if (iOne == -1) {
                cerr << "Skipping gene pair, missing gene: " << strToken << endl;
                continue;
            }

            //iOne = MapGene( mapGenes, m_vecstrGenes, strToken );
            /*
            if(GDist != NULL){
          iOne = GDist.GetGene(strToken);
          if( iOne == -1){
            cerr << "Skipping gene pair, missing gene: " << strToken << endl;
            continue;
          }
            }else{
          iOne = atoi(strToken.c_str());
            }
            */
        }

        strToken = CFile::OpenToken(pc, &pc);
        if (!strToken.length()) {
            delete[] acBuf;
            return false;
        }

        iTwo = GDist.GetGene(strToken);
        if (iTwo == -1) {
            cerr << "Skipping gene pair, missing gene: " << strToken << endl;
            continue;
        }
        //iTwo = MapGene( mapGenes, m_vecstrGenes, strToken );
        /*
        if(GDist != NULL){
          iTwo = GDist.GetGene(strToken);
          if( iTwo == -1){
        cerr << "Skipping gene pair, missing gene: " << strToken << endl;
        continue;
          }
        }else{
          iTwo = atoi(strToken.c_str());
        }
        */

        vecPairs.push_back(new PairStruct(iOne, iTwo));
    }
    delete[] acBuf;
    return true;
}

void find_bunch(CPCL &GDist, vector<PairStruct *> &vecPairs, vector<PairStruct *> &vecFinalPairs) {
    float C, d, d1, d2, density;
    size_t i, j, k, p, Bi, Bj, x;
    vector < PairStruct * > result_pairs;
    vector < PairStruct * > unresolved_pairs;

    result_pairs.resize(vecPairs.size());
    unresolved_pairs.resize(vecPairs.size());

    density = CMeta::GetNaN();
    for (i = 0; i < GDist.GetGenes(); i++) {
        for (j = 0; j < GDist.GetGenes(); j++) {

            if (i == j || CMeta::IsNaN(d = GDist.Get(i, j)))
                continue;

            for (p = 0; p < vecPairs.size(); p++) {
                if (CMeta::IsNaN(d1 = GDist.Get(vecPairs[p]->g1, i)) ||
                    CMeta::IsNaN(d2 = GDist.Get(j, vecPairs[p]->g2))) {
                    vecPairs[p]->score = CMeta::GetNaN();
                    continue;
                }
                vecPairs[p]->score = d1 + d2;
            }

            // DEBUG!!! check if direction is correct
            sort(vecPairs.begin(), vecPairs.end(), SortPairStruct());

            C = d;
            for (p = 0; p < vecPairs.size(); p++) {
                // check nan
                // DEBUG is this right
                if (CMeta::IsNaN(vecPairs[p]->score))
                    break;

                C += vecPairs[p]->score;

                if (C / (p + 1.0) < density) {
                    density = C / (p + 1.0);

                    result_pairs.clear();
                    unresolved_pairs.clear();

                    // store to B the i,j and pairs up to p
                    for (x = 0; x < vecPairs.size(); x++) {
                        if (x <= p)
                            result_pairs.push_back(vecPairs[x]);
                        else
                            unresolved_pairs.push_back(vecPairs[x]);
                    }
                    Bi = i;
                    Bj = j;
                }
            }
        }
    }

    // now set up results
    vecFinalPairs.push_back(new PairStruct(Bi, Bj));
    for (i = 0; i < result_pairs.size(); i++) {
        vecFinalPairs.push_back(result_pairs[i]);
    }

    // set up teh resulting pairs for next execution
    vecPairs.clear();
    for (j = 0; j < unresolved_pairs.size(); j++) {
        vecPairs.push_back(unresolved_pairs[j]);
    }
}


int main(int iArgs, char **aszArgs) {
    gengetopt_args_info sArgs;
    CGenome Genome;
    CGenes Genes(Genome);
    ifstream ifsm;
    CDat Dat;
    CDat Sim;
    size_t i, j, k, p;
    bool fModified;
    float d, d1, d2;
    CPCL PCL;

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }

    CMeta Meta(sArgs.verbosity_arg, sArgs.random_arg);

    // open query set genes
    if (sArgs.query_given) {
        ifsm.open(sArgs.query_arg);
        if (!Genes.Open(ifsm)) {
            cerr << "Could not open: " << sArgs.query_arg << endl;
            return 1;
        }
        ifsm.close();
    }

    cerr << "Num genes: " << Genes.GetGenes() << endl;

    // open similarity
    if (sArgs.sim_given) {
        if (!Sim.Open(sArgs.sim_arg, sArgs.memmap_flag)) {
            cerr << "Could not open: " << sArgs.sim_arg << endl;
            return 1;
        }

        if (sArgs.normalize_flag) {
            Sim.Normalize(CDat::ENormalizeMinMax);
            cerr << "Normalize: " << sArgs.sim_arg << endl;
        }
    }

    // open tissue expression matrix
    if (sArgs.tissue_given) {
        if (!PCL.Open(sArgs.tissue_arg, sArgs.skip_arg)) {
            cerr << "Could not open: " << sArgs.tissue_arg << endl;
            return 1;
        }
    }

    vector <size_t> qidxs;
    vector <GeneStruct> geneweights;

    cerr << "check #genes: " << Sim.GetGenes() << endl;
    //cerr << "check: " << Sim.GetGene( "9261" ) << endl;
    //cerr << "check: " << Sim.GetGene( "1385" ) << endl;

    qidxs.resize(Genes.GetGenes());
    for (i = 0; i < Genes.GetGenes(); i++) {
        //qidxs[i] = Sim.GetGene( Genes.GetGene(i).GetName().c_str() );
        qidxs[i] = Sim.GetGene(Genes.GetGene(i).GetName());
        //cerr << Genes.GetGene(i).GetName() << "\t" << qidxs[i] << endl;
    }

    // tissue
    vector <size_t> qidxs_pcl;
    vector<float> genetissue;
    vector<float> tprofile;
    vector <size_t> sim2pcl_idx;
    if (sArgs.tissue_given) {
        float profilesum;

        qidxs_pcl.resize(Genes.GetGenes());
        for (i = 0; i < Genes.GetGenes(); i++) {
            qidxs_pcl[i] = PCL.GetGene(Genes.GetGene(i).GetName());
        }

        // construct query gene group tissue profile
        tprofile.resize(PCL.GetExperiments());
        for (j = 0; j < PCL.GetExperiments(); j++) {
            tprofile[j] = 0;
        }

        profilesum = 0;
        for (i = 0; i < Genes.GetGenes(); i++) {
            if (qidxs_pcl[i] == -1)
                continue;
            for (j = 0; j < PCL.GetExperiments(); j++) {
                // maybe check for nan
                if (!CMeta::IsNaN((d = PCL.Get(qidxs_pcl[i], j)))) {
                    tprofile[j] += d;
                    profilesum += d;
                }
            }
        }

        // normalize tissue profile
        for (j = 0; j < PCL.GetExperiments(); j++) {
            tprofile[j] = tprofile[j] / profilesum;
        }

        genetissue.resize(PCL.GetGenes());

        float tsum;
        float profile_sum;
        for (i = 0; i < PCL.GetGenes(); i++) {
            tsum = 0;
            profile_sum = 0;
            for (j = 0; j < PCL.GetExperiments(); j++) {
                // check NaN
                if (!CMeta::IsNaN((d = PCL.Get(i, j)))) {
                    tsum += (tprofile[j] * d);
                    profile_sum += d;
                }
            }
            genetissue[i] = tsum / profile_sum;
        }

        if (sArgs.sim_given) {
            sim2pcl_idx.resize(Sim.GetGenes());
            for (i = 0; i < Sim.GetGenes(); i++) {
                sim2pcl_idx[i] = PCL.GetGene(Sim.GetGene(i));
            }
        } else {
            // no imilarity matrix just print results
            geneweights.resize(PCL.GetGenes());
            for (i = 0; i < geneweights.size(); i++) {
                geneweights[i].score = genetissue[i];
                geneweights[i].idx = i;
            }
            // sort genes by weight
            sort(geneweights.begin(), geneweights.end(), SortGeneStruct());

            for (i = 0; i < PCL.GetGenes(); i++) {
                if (CMeta::IsNaN(geneweights[i].score))
                    continue;
                cout << PCL.GetGene(geneweights[i].idx) << "\t" << geneweights[i].score << endl;
            }

            return 0;
        }
    }

    if (sArgs.query_given) {

        cerr << "mapped: " << qidxs.size() << endl;

        geneweights.resize(Sim.GetGenes());
        for (i = 0; i < geneweights.size(); i++) {
            geneweights[i].score = 0.0;
            geneweights[i].idx = i;
        }

        for (i = 0; i < Sim.GetGenes(); i++) {
            for (j = 0; j < qidxs.size(); j++) {
                if (i == qidxs[j] || qidxs[j] == -1)
                    continue;
                if (!CMeta::IsNaN((d = Sim.Get(i, qidxs[j])))) {
                    if (sArgs.tissue_given)
                        if (sim2pcl_idx[i] != -1)
                            geneweights[i].score += (d * genetissue[sim2pcl_idx[i]]);
                        else // currently genes only in Sim matrix and not in PCL removed
                            geneweights[i].score = CMeta::GetNaN();
                    else
                        geneweights[i].score += d;
                }
            }
        }

        for (i = 0; i < Genes.GetGenes(); i++) {
            geneweights[qidxs[i]].score = CMeta::GetNaN();
        }

        // sort genes by weight
        sort(geneweights.begin(), geneweights.end(), SortGeneStruct());

        for (i = 0; i < Sim.GetGenes(); i++) {
            if (CMeta::IsNaN(geneweights[i].score))
                continue;
            cout << Sim.GetGene(geneweights[i].idx) << "\t" << geneweights[i].score << endl;
        }
        return 0;
    }
    //////////////////////////////////////////////

    CPCL GPCL;
    CPCL GDist;
    CPCL GNext;

    // open gene direct distance matrix
    if (sArgs.input_given) {
        if (!GPCL.Open(sArgs.input_arg, sArgs.skip_arg)) {
            cerr << "Could not open: " << sArgs.input_arg << endl;
            return 1;
        }
        cerr << "Open: " << sArgs.input_given << endl;

        if (sArgs.NegLog_flag) {
            // convert the values to -log(prob)
            for (i = 0; i < GPCL.GetGenes(); i++) {
                for (j = 0; j < GPCL.GetGenes(); j++) {
                    if (i == j)
                        GPCL.Set(i, j, CMeta::GetNaN());

                    if (CMeta::IsNaN(d = GPCL.Get(i, j)))
                        continue;

                    GPCL.Set(i, j, -log(d));
                }
            }
        }
    }

    GDist.Open(GPCL);
    cerr << "opend1: " << GDist.GetGenes() << endl;
    GNext.Open(GPCL);
    cerr << "opend2: " << GNext.GetGenes() << endl;

    // initialize two matrix
    for (i = 0; i < GDist.GetGenes(); i++) {
        //GDist.Set(i, CMeta::GetNaN());
        for (j = 0; j < GDist.GetGenes(); j++) {
            if (i == j)
                GNext.Set(i, j, CMeta::GetNaN());
            else
                GNext.Set(i, j, i);
        }
    }

    cerr << "WTF" << endl;

    for (k = 0; k < GDist.GetGenes(); k++) {

        if (k % 100 == 0)
            cerr << "k: " << k << endl;

        for (i = 0; i < GDist.GetGenes(); i++) {
            for (j = 0; j < GDist.GetGenes(); j++) {
                if (i == j)
                    continue;
                // update distance matrix
                // (i,k), (k,j)
                d = GDist.Get(i, j);

                if (CMeta::IsNaN(d1 = GDist.Get(i, k)) || CMeta::IsNaN(d2 = GDist.Get(k, j)))
                    continue;

                if (CMeta::IsNaN(d) || (d1 + d2) < d) {
                    GDist.Set(i, j, (d1 + d2));

                    // update path matrix
                    GNext.Set(i, j, GNext.Get(k, j));
                }
            }
        }
    }

    cerr << "done" << endl;

    if (sArgs.savematrix_flag) {
        std::ofstream ofsm;

        cerr << "Save output" << endl;
        std::stringstream sstmDist;
        std::stringstream sstmNext;

        std::string path(sArgs.input_arg);
        size_t pos = path.find_last_of("/");
        size_t lastindex = path.find_last_of(".");
        std::string rawname = path.substr(0, lastindex);

        sstmNext << "" << rawname << ".next.bin";
        cerr << "Save next matrix: " << sstmNext.str() << endl;
        ofsm.open(sstmNext.str().c_str());

        GNext.SaveBinary(ofsm);
        ofsm.close();

        if (!sArgs.savedab_flag) {
            sstmDist << "" << rawname << ".dist.bin";

            ofsm.open(sstmDist.str().c_str());
            cerr << "Save distance matrix: " << sstmDist.str() << endl;
            GDist.SaveBinary(ofsm);
            ofsm.close();
        } else {
            CDat GDistDat;
            GDistDat.Open(GDist.GetGeneNames());

            sstmDist << "" << rawname << ".dist.dab";

            // map PCL to dab
            for (i = 0; i < GDistDat.GetGenes(); i++)
                for (j = (i + 1); j < GDistDat.GetGenes(); j++) {
                    GDistDat.Set(i, j, (GDist.Get(i, j) + GDist.Get(j, i)) / 2.0);
                }

            cerr << "Save distance matrix: " << sstmDist.str() << endl;
            GDistDat.Save(sstmDist.str().c_str());
        }
    }

    //****************
    //* find steiner tree

    if (sArgs.genepairs_given) {
        // read the gene pair
        vector < PairStruct * > vecOrigPairs;
        vector < PairStruct * > vecPairs;
        vector < PairStruct * > vecFinalPairs;
        float C, density;

        if (!open_genepairs(sArgs.genepairs_arg, vecOrigPairs, GDist)) {
            cerr << "ERROR: failed to open " << sArgs.genepairs_arg << endl;
            return 1;
        }
        vecPairs.resize(vecOrigPairs.size());

        for (i = 0; i < vecOrigPairs.size(); i++) {
            vecPairs[i] = vecOrigPairs[i];
        }

        // might just make it that if all genes are covered halt
        // currently it is if all pairs are covered
        while (vecPairs.size() > 0) {
            find_bunch(GDist, vecPairs, vecFinalPairs);
        }

        // now vecFinalPairs = original required pairs + new gene pairs
        //  print output
    }


    return 0;
}

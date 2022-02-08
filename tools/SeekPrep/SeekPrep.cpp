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
#include "stdafx.h"
#include "cmdline.h"

enum ExportMode {
    DISTANCE_MATRIX, COUNT_MATRIX, WEIGHTSUM_MATRIX, GENE_PRESENCE_VECTOR
};

enum NormMode {
    RANK_NORM, Z_NORM
};

//tType can only be unsigned short or utype
//norm_mode = 0 (rank norm), 1 (subtract_z norm)
template<class tType>
bool CalculateMatrix(const NormMode norm_mode,
                     const ExportMode e, const vector <string> &dab_list,
                     const string &dab_dir, const string &outdir, const string &outdab,
                     const vector <string> &vecstrGenes, const vector<float> &w, const float rbp_p,
                     const int MAX_RANK, const float exp) { //w=weight of dsets
    //rbp_p, MAX_RANK for rank-based normalize
    //exp for value-based or subtract_z normalize
    const int RANK_NORM = 0;
    const int Z_NORM = 1;

    size_t i, j, k;

    if (e == GENE_PRESENCE_VECTOR) {
        vector <utype> gpres;
        CSeekTools::InitVector(gpres, vecstrGenes.size(), (utype) 0);
        for (i = 0; i < dab_list.size(); i++) {
            fprintf(stderr, "Reading %zu of %lu: %s\n", i, dab_list.size(),
                    dab_list[i].c_str());
            CSeekIntIntMap d1(vecstrGenes.size());
            string dabfile = dab_dir + "/" + dab_list[i];
            CSeekWriter::ReadSeekSparseMatrixHeader<tType>(dabfile.c_str(), d1);
            const vector <utype> &allR = d1.GetAllReverse();
            for (j = 0; j < d1.GetNumSet(); j++) {
                utype gi = allR[j];
                gpres[gi]++;
            }
        }
        string ofile = outdir + "/" + outdab + ".gene_count";
        CSeekTools::WriteArray(ofile.c_str(), gpres);

    } else if (e == DISTANCE_MATRIX) {
        CHalfMatrix<float> res;
        res.Initialize(vecstrGenes.size());
        for (i = 0; i < vecstrGenes.size(); i++)
            for (j = i + 1; j < vecstrGenes.size(); j++)
                res.Set(i, j, 0);

        for (i = 0; i < dab_list.size(); i++) {
            fprintf(stderr, "Reading %zu of %lu: %s\n", i, dab_list.size(),
                    dab_list[i].c_str());
            CSeekIntIntMap d1(vecstrGenes.size());
            string dabfile = dab_dir + "/" + dab_list[i];
            CSparseFlatMatrix<float> sm(0);
            if (norm_mode == RANK_NORM) {
                CSeekWriter::ReadSeekSparseMatrix<tType>(dabfile.c_str(), sm, d1,
                                                         MAX_RANK, rbp_p, vecstrGenes);
                //NEW
                //CSeekWriter::RemoveDominant<tType>(sm, d1, vecstrGenes);
                //CSeekWriter::RemoveDominant<tType>(sm, d1, vecstrGenes);
            } else if (norm_mode == Z_NORM) {
                CSeekWriter::ReadSeekSparseMatrix<tType>(dabfile.c_str(), sm, d1,
                                                         vecstrGenes, (int) (0.10 * vecstrGenes.size()), exp);
            }
            fprintf(stderr, "Summing...\n");
            CSeekWriter::SumSparseMatrix(sm, res, d1, w[i]);
            fprintf(stderr, "Finished Summing...\n");
        }

        /*for(i=0; i<vecstrGenes.size(); i++){
            for(j=i+1; j<vecstrGenes.size(); j++){
                float v = res.Get(i,j);
                if(CMeta::IsNaN(v)){
                    fprintf(stderr, "Error, %d %d is nan\n", i,j);
                }
            }
        }*/

        ofstream ofsm;
        string ofile = outdir + "/" + outdab + ".half";
        ofsm.open(ofile.c_str(), ios_base::binary);
        res.Save(ofsm, true);
        ofsm.close();

    } else if (e == COUNT_MATRIX) {
        CHalfMatrix<unsigned short> res;
        res.Initialize(vecstrGenes.size());
        for (i = 0; i < vecstrGenes.size(); i++)
            for (j = i + 1; j < vecstrGenes.size(); j++)
                res.Set(i, j, 0);

        for (i = 0; i < dab_list.size(); i++) {
            fprintf(stderr, "Reading %zu of %lu: %s\n", i, dab_list.size(),
                    dab_list[i].c_str());
            CSeekIntIntMap d1(vecstrGenes.size());
            string dabfile = dab_dir + "/" + dab_list[i];
            CSeekWriter::ReadSeekSparseMatrixHeader<tType>(dabfile.c_str(), d1);
            const vector <utype> &allR = d1.GetAllReverse();
            for (j = 0; j < d1.GetNumSet(); j++) {
                utype gi = allR[j];
                for (k = j + 1; k < d1.GetNumSet(); k++) {
                    utype gj = allR[k];
                    res.Set(gi, gj, res.Get(gi, gj) + 1);
                }
            }
        }
        ofstream ofsm;
        string ofile = outdir + "/" + outdab + ".pair_count";
        ofsm.open(ofile.c_str(), ios_base::binary);
        res.Save(ofsm, true);
        ofsm.close();

    } else if (e == WEIGHTSUM_MATRIX) {
        CHalfMatrix<float> res;
        res.Initialize(vecstrGenes.size());
        for (i = 0; i < vecstrGenes.size(); i++)
            for (j = i + 1; j < vecstrGenes.size(); j++)
                res.Set(i, j, 0);

        for (i = 0; i < dab_list.size(); i++) {
            fprintf(stderr, "Reading %zu of %lu: %s\n", i, dab_list.size(),
                    dab_list[i].c_str());
            CSeekIntIntMap d1(vecstrGenes.size());
            string dabfile = dab_dir + "/" + dab_list[i];
            CSeekWriter::ReadSeekSparseMatrixHeader<tType>(dabfile.c_str(), d1);
            const vector <utype> &allR = d1.GetAllReverse();
            for (j = 0; j < d1.GetNumSet(); j++) {
                utype gi = allR[j];
                for (k = j + 1; k < d1.GetNumSet(); k++) {
                    utype gj = allR[k];
                    res.Set(gi, gj, res.Get(gi, gj) + w[i]);
                }
            }
        }
        ofstream ofsm;
        string ofile = outdir + "/" + outdab + ".pair_weightsum";
        ofsm.open(ofile.c_str(), ios_base::binary);
        res.Save(ofsm, true);
        ofsm.close();
    } else {
        fprintf(stderr, "Unrecognized mode!\n");
        return false;
    }
    return true;
}


bool InitializeDataset(size_t &iDatasets, vector <string> &vecstrDatasets,
                       string &strPrepInputDirectory, vector<CSeekDataset *> &vc) {
    vc.clear();
    vc.resize(iDatasets);
    size_t i, j, k;

    for (i = 0; i < iDatasets; i++) {
        vc[i] = new CSeekDataset();
        string strFileStem = vecstrDatasets[i];
        string strAvgPath = strPrepInputDirectory + "/" + strFileStem + ".gavg";
        string strPresencePath = strPrepInputDirectory + "/" + strFileStem + ".gpres";
        vc[i]->ReadGeneAverage(strAvgPath);
        vc[i]->ReadGenePresence(strPresencePath);
        vc[i]->InitializeGeneMap();

        //DEBUGGING==============
        /*
        fprintf(stderr, "==Dataset %d==\n", i);
        CSeekIntIntMap *mapG = vc[i]->GetGeneMap();
        const vector<utype> &allR = mapG->GetAllReverse();
        for(j=0; j<mapG->GetNumSet(); j++){
            utype ja = mapG->GetReverse(j);
            fprintf(stderr, "%d\t%.3f\n", ja, vc[i]->GetGeneAverage(ja));
        }*/
    }
    return true;
}


bool InitializeDB(size_t &iDatasets, size_t &iGenes,
                  vector <string> &vecstrGenes, vector<CSeekDataset *> &vc, CDatabaselet &DBL) {

    utype i, j, k;
    vector<char> cQuery;
    CSeekTools::InitVector(cQuery, iGenes, (char) 0);
    vector <string> allQuery;
    map <string, size_t> mapstrGenes;
    for (i = 0; i < vecstrGenes.size(); i++) {
        mapstrGenes[vecstrGenes[i]] = i;
    }

    /* Databaselet mapping */
    map <string, size_t> dbmap;
    vector <utype> veciQuery;
    for (i = 0; i < DBL.GetGenes(); i++) {
        const string& strQuery = DBL.GetGene(i);
        dbmap[strQuery] = i;
        allQuery.push_back(strQuery);
        /* global mapping */
        cQuery[mapstrGenes[strQuery]] = 1;
        veciQuery.push_back((utype) mapstrGenes[strQuery]);
    }

    //fprintf(stderr, "Start initializing query map...\n");
    for (i = 0; i < iDatasets; i++) {
        vc[i]->InitializeQueryBlock(veciQuery);
    }

    //fprintf(stderr, "Finished initializing map\n");

    //fprintf(stderr, "Start making gene-centric...\n");
    for (i = 0; i < DBL.GetGenes(); i++) {
        vector<unsigned char> Q;
        /* expanded */
        DBL.Get(i, Q);
        utype m = mapstrGenes[DBL.GetGene(i)];
        CSeekIntIntMap *qu = nullptr;
        utype db;
        for (j = 0; j < iDatasets; j++) {
            if ((qu = vc[j]->GetDBMap()) == nullptr) continue;
            if (CSeekTools::IsNaN(db = qu->GetForward(m))) continue;
            unsigned char **r = vc[j]->GetMatrix();
            for (k = 0; k < iGenes; k++) {
                unsigned char c = Q[k * iDatasets + j];
                r[db][k] = c;
            }
        }
    }

    //fprintf(stderr, "Finished making gene-centric\n");

    return true;
}

bool OpenDBFiles(string &DBFile, vector<unsigned char *> &cc, bool &useNibble) {
    CDatabaselet *CD;
    (CD = new CDatabaselet(useNibble))->Open(DBFile);
    cc.push_back(CD->GetCharImage());
    return true;
}


bool OpenDB(string &DBFile, bool &useNibble, size_t &iDatasets,
            size_t &m_iGenes, vector<string> &vecstrGenes,
            map<utype, utype> &mapiPlatform, const vector<float> &quant,
            vector<CSeekDataset *> &vc, CFullMatrix<float> &platform_avg,
            CFullMatrix<float> &platform_stdev, CFullMatrix<uint32_t> &platform_count,
            vector<string> &vecstrQuery,
            const bool &logit)
{

    size_t i, j, jj, k;

    CDatabaselet CD(useNibble);
    CD.Open(DBFile);
    InitializeDB(iDatasets, m_iGenes, vecstrGenes, vc, CD);

    vector <string> presGenes;
    for (i = 0; i < CD.GetGenes(); i++)
        presGenes.push_back(CD.GetGene(i));

    size_t numPlatforms = platform_avg.GetRows();
    map <string, size_t> mapstriGenes;
    for (i = 0; i < vecstrGenes.size(); i++)
        mapstriGenes[vecstrGenes[i]] = i;

    //fprintf(stderr, "Start calculating platform average\n");
    for (i = 0; i < CD.GetGenes(); i++) {
        vector<float> sum, sq_sum, mean, stdev;
        vector<int> num;
        sum.resize(numPlatforms);
        sq_sum.resize(numPlatforms);
        mean.resize(numPlatforms);
        stdev.resize(numPlatforms);
        num.resize(numPlatforms);

        for (k = 0; k < numPlatforms; k++) {
            sum[k] = 0;
            sq_sum[k] = 0;
            mean[k] = 0;
            stdev[k] = 0;
            num[k] = 0;
        }

        const string& thisGene = CD.GetGene(i);
        size_t geneID = mapstriGenes[thisGene];
        vecstrQuery.push_back(thisGene);
        for (k = 0; k < iDatasets; k++) {
            CSeekIntIntMap *mapQ = vc[k]->GetDBMap();
            CSeekIntIntMap *mapG = vc[k]->GetGeneMap();
            const vector <utype> &allR = mapG->GetAllReverse();
            if (mapQ == nullptr) continue;
            unsigned char **f = vc[k]->GetMatrix();
            utype iQ = mapQ->GetForward(geneID);
            if (CSeekTools::IsNaN(iQ)) {
                continue;
            }
            utype platform_id = mapiPlatform[k];
            if (platform_id >= (utype) numPlatforms) {
                printf("Error, platforms are equal %d %d",
                       (int) platform_id, (int) numPlatforms);
                getchar();
            }

            for (jj = 0; jj < mapG->GetNumSet(); jj++) {
                j = mapG->GetReverse(jj);
                unsigned char uc = f[iQ][j];
                float v = 0;
                if (uc == 255) v = CMeta::GetNaN();
                else {
                    float vv = -1;
                    if (logit) vv = log(quant[uc]) - log((float) (1.0 - quant[uc]));
                    else vv = quant[uc];
                    v = vv - vc[k]->GetGeneAverage(j);
                    //fprintf(stderr, "%.5f\t%.5f\n", vv, v);
                    if (isnan(vv) || isinf(vv) || isnan(vc[k]->GetGeneAverage(j)) ||
                        isinf(vc[k]->GetGeneAverage(j))) {
                        fprintf(stderr, "%d D%d %.5f %.5f %d\n", (int) j, (int) k,
                                vv, vc[k]->GetGeneAverage(j), mapG->GetForward(j));
                    }
                    //v = quant[uc];
                    sum[platform_id] += v;
                    num[platform_id]++;
                    sq_sum[platform_id] += v * v;
                }
            }
        }

        for (k = 0; k < numPlatforms; k++) {
            if (num[k] == 0) continue;
            mean[k] = sum[k] / (float) num[k];
            stdev[k] = sq_sum[k] / (float) num[k] - mean[k] * mean[k];
            stdev[k] = sqrt(stdev[k]);
            // fprintf(stderr, "%s G%zu P%zu %.5f %.5f\n", thisGene.c_str(), geneID, k, mean[k], stdev[k]);
            platform_avg.Set(k, geneID, mean[k]);
            platform_stdev.Set(k, geneID, stdev[k]);
            platform_count.Set(k, geneID, num[k]);
        }
    }

    //fprintf(stderr, "Finished calculating platform average\n");
    //fprintf(stderr, "Start deleting\n");

    for (i = 0; i < iDatasets; i++)
        vc[i]->DeleteQueryBlock();
    //fprintf(stderr, "Finished deleting\n");
    return true;
}

string getPclFileStem(string fileName) {
    // match on substring '.pcl' will also match on '.pcl.bin'
    string fileStem = fileName;
    int pos;
    pos = fileStem.rfind(".pcl");
    fileStem = fileStem.substr(0, pos);
    pos = fileStem.rfind(".bin");
    fileStem = fileStem.substr(0, pos);
    return fileStem;
}


int main(int iArgs, char **aszArgs) {
    static const size_t c_iBuffer = 1024;
#ifdef WIN32
    pthread_win32_process_attach_np( );
#endif // WIN32
    gengetopt_args_info sArgs{};
    ifstream ifsm;
    istream *pistm;
    vector <string> vecstrLine, vecstrGenes, vecstrDBs, vecstrQuery;
    char acBuffer[c_iBuffer];
    utype i, j, k, l;

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }

    if (sArgs.datasetsize_flag == 1) {
        // open pcl file and ouput number of expiements
        string pclfile = sArgs.pclinput_arg;
        if (pclfile.empty() || !CMeta::IsExtension(pclfile, ".bin")) {
          cerr << "Missing pcl input or input file is not bin type!" << endl;
          return -1;
        }

        CPCL pcl;
        if (!pcl.Open(pclfile.c_str())) {
          cerr << "Error opening pcl file " << pclfile << endl;
          return -1;
        }
        string dsetName = getPclFileStem(CMeta::Basename(pclfile.c_str()));
        size_t numExp = pcl.GetExperiments();
        cout << dsetName << "\t" << numExp << endl;
        return 0;
    }

    /* reading gene-mapping file */
    if (sArgs.input_arg) {
        ifsm.open(sArgs.input_arg);
        pistm = &ifsm;
    } else
        pistm = &cin;

    map <string, size_t> mapstriGenes;
    while (!pistm->eof()) {
        pistm->getline(acBuffer, c_iBuffer - 1);
        acBuffer[c_iBuffer - 1] = 0;
        vecstrLine.clear();
        CMeta::Tokenize(acBuffer, vecstrLine);
        if (vecstrLine.size() < 2) {
            //cerr << "Ignoring line: " << acBuffer << endl;
            continue;
        }
        if (!(i = strtol(vecstrLine[0].c_str(), nullptr, 10))) {
            cerr << "Illegal gene ID: " << vecstrLine[0] << " for "
                 << vecstrLine[1] << endl;
            return 1;
        }
        i--;
        if (vecstrGenes.size() <= i)
            vecstrGenes.resize(i + 1);
        vecstrGenes[i] = vecstrLine[1];
        mapstriGenes[vecstrGenes[i]] = i;
    }

    if (sArgs.input_arg) ifsm.close();

    if ((sArgs.dab_flag == 1 && sArgs.norm_flag == 1) ||
        sArgs.dabset_flag == 1) {
        if (sArgs.default_type_arg != 0 && sArgs.default_type_arg != 1) {
            fprintf(stderr, "Error, invalid default_type argument!\n");
            return -1;
        }
    }

    /* PCL mode */
    if (sArgs.pclbin_flag == 1) {

        if (sArgs.sinfo_flag == 1) {
            string fileName = CMeta::Basename(sArgs.pclinput_arg);
            string fileStem = getPclFileStem(fileName);

            char outFile[125];
            sprintf(outFile, "%s/%s.sinfo", sArgs.dir_out_arg, fileStem.c_str());

            string pclfile = sArgs.pclinput_arg;
            if (!CMeta::IsExtension(pclfile, ".bin")) {
                fprintf(stderr, "Input file is not bin type!\n");
                return 1;
            }

            CPCL pcl;
            if (!pcl.Open(pclfile.c_str())) {
                fprintf(stderr, "Error opening file\n");
            }
            CMeasurePearNorm pn;
            CDat Dat;
            int numG = pcl.GetGeneNames().size();

            vector<int> veciGenes;
            veciGenes.resize(numG);
            for (j = 0; j < veciGenes.size(); ++j)
                veciGenes[j] = j;
            int iOne, iTwo;
            Dat.Open(pcl.GetGeneNames());

            for (k = 0; k < Dat.GetGenes(); ++k)
                for (j = (k + 1); j < Dat.GetGenes(); ++j)
                    Dat.Set(k, j, CMeta::GetNaN());

            omp_set_num_threads(8);

#pragma omp parallel for \
            shared(pcl, Dat, veciGenes, pn) \
            private(k, iOne, j, iTwo) \
            firstprivate(numG) \
            schedule(dynamic)
            for (k = 0; k < numG; ++k) {
                if ((iOne = veciGenes[k]) == -1)
                    continue;
                float *adOne = &pcl.Get(iOne)[0];
                for (j = (k + 1); j < numG; ++j) {
                    if ((iTwo = veciGenes[j]) != -1) {
                        auto x = (float) pn.Measure(adOne,
                                                    pcl.GetExperiments(), &pcl.Get(iTwo)[0],
                                                    pcl.GetExperiments());
                        Dat.Set(k, j, x);
                        //fprintf(stderr, "%d %d %.5f\n", k, j, x);
                    }
                }
            }

            double gmean, gstdev;
            size_t iN;

            Dat.AveStd(gmean, gstdev, iN);
            fprintf(stderr, "%.5f %.5f\n", (float) gmean, (float) gstdev);
            vector<float> vv;
            vv.resize(2);
            vv[0] = (float) gmean;
            vv[1] = (float) gstdev;
            CSeekTools::WriteArray(outFile, vv);

        }

            //if calculating gene variance per dataset
        else if (sArgs.gexpvarmean_flag == 1) {
            string fileName = CMeta::Basename(sArgs.pclinput_arg);
            string fileStem = getPclFileStem(fileName);
            char outFile[125];

            string pclfile = sArgs.pclinput_arg;
            if (!CMeta::IsExtension(pclfile, ".bin")) {
                fprintf(stderr, "Input file is not bin type!\n");
                return 1;
            }

            CPCL pcl;
            if (!pcl.Open(pclfile.c_str())) {
                fprintf(stderr, "Error opening file\n");
                return 1;
            }

            vector<float> var;
            vector<float> avg;

            CSeekTools::InitVector(var, vecstrGenes.size(), (float) CMeta::GetNaN());
            CSeekTools::InitVector(avg, vecstrGenes.size(), (float) CMeta::GetNaN());

            int totNumExperiments = pcl.GetExperiments();
            if (totNumExperiments <= 2) {
                fprintf(stderr, "This dataset is skipped because it has <=2 columns\n");
                fprintf(stderr, "An empty vector will be returned\n");
            } else {
                for (j = 0; j < vecstrGenes.size(); j++) {
                    utype gene = pcl.GetGene(vecstrGenes[j]);
                    if (CSeekTools::IsNaN(gene)) continue; //gene does not exist in the dataset
                    float *val = pcl.Get(gene);
                    vector<float> rowVal;
                    for (k = 0; k < totNumExperiments; k++)
                        rowVal.push_back(val[k]);
                    float mean = 0;
                    float variance = 0;
                    for (k = 0; k < rowVal.size(); k++)
                        mean += rowVal[k];
                    mean /= rowVal.size();
                    for (k = 0; k < rowVal.size(); k++)
                        variance += (rowVal[k] - mean) * (rowVal[k] - mean);
                    variance /= rowVal.size();
                    var[j] = variance;
                    avg[j] = mean;
                    //fprintf(stderr, "%.5f %.5f\n", mean, variance);
                }
                //fprintf(stderr, "done\n");
            }
            //DEBUG
            fprintf(stderr, "UCSC genes\tVariance\tAverage\n");
            for (j = 0; j < vecstrGenes.size(); j++) {
                utype gene = pcl.GetGene(vecstrGenes[j]);
                if (CSeekTools::IsNaN(gene)) continue;
                fprintf(stderr, "%s\t%0.4f\t%.4f\n", vecstrGenes[j].c_str(),
                        var[j], avg[j]);
            }
            //fprintf(stderr, "G\n");
            sprintf(outFile, "%s/%s.gexpvar", sArgs.dir_out_arg,
                    fileStem.c_str());
            CSeekTools::WriteArray(outFile, var);
            sprintf(outFile, "%s/%s.gexpmean", sArgs.dir_out_arg,
                    fileStem.c_str());
            CSeekTools::WriteArray(outFile, avg);
        }

    }


    /* DB mode */
    if (sArgs.db_flag == 1) {
        if (sArgs.combineplat_flag == 1) {
            // combine the platform statistics (avg, stdev, count) from in_dir1 and in_dir2 and write to dir_out
            assert(sArgs.in_dir1_arg != NULL);
            assert(sArgs.in_dir2_arg != NULL);
            assert(sArgs.dir_out_arg != NULL);
            combinePlatformStatistics(sArgs.in_dir1_arg, sArgs.in_dir2_arg, sArgs.dir_out_arg);
        } else if (sArgs.gplat_flag == 1) {
            if (!sArgs.quant_arg)
            {
                fprintf(stderr, "Must give quant file\n");
                return -1;
            }

            vector<float> quant;
            string strQuantFile = sArgs.quant_arg;
            //fprintf(stderr, "%s\n", strQuantFile.c_str());
            CSeekTools::ReadQuantFile(strQuantFile, quant);

            //fprintf(stderr, "quant %d %.5f\n", 170, quant[170]);
            //getchar();


            bool logit = false;
            if (sArgs.logit_flag == 1) logit = true;

            map <string, utype> mapstriPlatform;
            map <utype, string> mapistrPlatform;
            map <utype, utype> mapiPlatform;
            vector <string> vecstrPlatforms;

            vector <string> vecstrDatasets;

            ifsm.open(sArgs.dset_arg);
            i = 0;
            while (!pistm->eof()) {
                pistm->getline(acBuffer, c_iBuffer - 1);
                if (acBuffer[0] == 0) break;
                acBuffer[c_iBuffer - 1] = 0;
                vecstrLine.clear();
                CMeta::Tokenize(acBuffer, vecstrLine);
                /* read dataset name */
                vecstrDatasets.push_back(vecstrLine[0]);
                /* just read the platform information */
                string pl = vecstrLine[1];
                map<string, utype>::const_iterator iter;
                iter = mapstriPlatform.find(pl);
                if (iter == mapstriPlatform.end()) {
                    utype s = mapstriPlatform.size();
                    mapstriPlatform[pl] = s;
                    mapistrPlatform[s] = pl;
                    vecstrPlatforms.push_back(pl);
                }
                utype platform_id = mapstriPlatform[pl];
                mapiPlatform[i] = platform_id;
                i++;
            }
            ifsm.close();

            vector <string> dblist;
            ifsm.open(sArgs.dblist_arg);
            i = 0;
            while (!pistm->eof()) {
                pistm->getline(acBuffer, c_iBuffer - 1);
                if (acBuffer[0] == 0) break;
                acBuffer[c_iBuffer - 1] = 0;
                dblist.emplace_back(acBuffer);
            }
            dblist.resize(dblist.size());
            ifsm.close();

            bool useNibble = false;
            if (sArgs.useNibble_flag == 1) useNibble = true;

            size_t numPlatforms = mapstriPlatform.size();
            size_t iDatasets = vecstrDatasets.size();
            size_t m_iGenes = vecstrGenes.size();

            // Create a SeekPlatforms object to hold all the platform data
            SeekPlatforms platforms(numPlatforms, m_iGenes);
            platforms.setPlatformNames(vecstrPlatforms);
            platforms.setPlatformNameMap(mapstriPlatform);

            //printf("Size: %d %d\n", numPlatforms, m_iGenes); getchar();

            int numThreads = omp_get_max_threads();
            if (numThreads > 16) {
                numThreads = 16;
            }
            omp_set_num_threads(numThreads);

            /*if(iDatasets<numThreads){
                numThreads = iDatasets;
                omp_set_num_threads(numThreads);
            }*/

            if (dblist.size() < numThreads) {
                numThreads = dblist.size();
                omp_set_num_threads(numThreads);
            }
            cout << "Calc Platform Stats: Num threads: " << numThreads << endl;

            time_t startTime, endTime;
            startTime = time(nullptr);

            string strPrepInputDirectory = sArgs.dir_prep_in_arg;
            auto *vc = new vector<CSeekDataset *>[numThreads];
            auto *platform_avg_threads =
                    new CFullMatrix<float>[numThreads];
            auto *platform_stdev_threads =
                    new CFullMatrix<float>[numThreads];
            CFullMatrix<uint32_t> *platform_count_threads =
                new CFullMatrix<uint32_t>[numThreads];

            // Initialize per-thread data
#pragma omp parallel for \
            shared(iDatasets, vecstrDatasets, strPrepInputDirectory, vc, numThreads, \
            platform_avg_threads, platform_stdev_threads, platform_count_threads, \
            numPlatforms, m_iGenes) \
            private(i) schedule(dynamic)
            for (i = 0; i < numThreads; i++) {
                InitializeDataset(iDatasets, vecstrDatasets,
                                  strPrepInputDirectory, vc[i]);
                platform_avg_threads[i].Initialize(numPlatforms, m_iGenes);
                platform_stdev_threads[i].Initialize(numPlatforms, m_iGenes);
                platform_count_threads[i].Initialize(numPlatforms, m_iGenes);
                for (int j = 0; j < numPlatforms; j++) {
                    for (int k = 0; k < m_iGenes; k++) {
                        platform_avg_threads[i].Set(j, k, CMeta::GetNaN());
                        platform_stdev_threads[i].Set(j, k, CMeta::GetNaN());
                        platform_count_threads[i].Set(j, k, 0);
                        //platform_avg_threads[i].Set(j, k, (float) 0);
                        //platform_stdev_threads[i].Set(j, k, (float) 1.0);
                    }
                }
            }

            endTime = time(nullptr);
            printf("init time %ld sec\n", endTime - startTime);
            startTime = time(nullptr);

            //printf("Dataset initialized"); getchar();
            vector <string> localVectrQuery;

            // Open the db files and accumulate statistics per thread for those file
#pragma omp parallel for \
			shared(vc, dblist, iDatasets, m_iGenes, vecstrGenes, mapiPlatform, quant, \
			platform_avg_threads, platform_stdev_threads, localVectrQuery, logit) \
			private(i) firstprivate(useNibble) schedule(dynamic)
            for (i = 0; i < dblist.size(); i++) {
                int tid = omp_get_thread_num();
                string DBFile = dblist[i];
                fprintf(stderr, "opening db file %s\n", DBFile.c_str());
                OpenDB(DBFile, useNibble, iDatasets, m_iGenes,
                       vecstrGenes, mapiPlatform, quant, vc[tid],
                       platform_avg_threads[tid], platform_stdev_threads[tid],
                       platform_count_threads[tid],
                       localVectrQuery, logit);
                fprintf(stderr, "finished opening db file %s\n",
                        DBFile.c_str());
            }

            endTime = time(nullptr);
            printf("thread calc time %ld sec\n", endTime - startTime);
            startTime = time(nullptr);

            // combine the per-thread avg and stdev stats together
#pragma omp parallel for \
            shared(platform_avg_threads, platform_stdev_threads, platform_count_threads, \
            platforms, numPlatforms, m_iGenes) \
            private(j) schedule(dynamic)
            for (j = 0; j < numPlatforms; j++) {
                int tid = omp_get_thread_num();
                // printf("tid: %d, j: %d\n", tid, j);
                for (int k = 0; k < m_iGenes; k++) {
                    uint32_t runningCount = 0;
                    double runningAvg = 0;
                    double runningVar = 0;
                    for (int i = 0; i < numThreads; i++) {
                        double avg = platform_avg_threads[i].Get(j, k);
                        double stdev = platform_stdev_threads[i].Get(j, k);
                        uint32_t cnt = platform_count_threads[i].Get(j, k);
                        if (avg == CMeta::GetNaN() || stdev == CMeta::GetNaN() || cnt == 0) {
                            continue;
                        }
                        double var = pow(stdev, 2);
                        double newAvg = (runningAvg * runningCount + avg * cnt) / (runningCount + cnt);
                        double diffRunningAvg = runningAvg - newAvg;
                        double diffTAvg = avg - newAvg;
                        double newVar = (runningCount * (runningVar + pow(diffRunningAvg, 2)) +
                                         cnt * (var + pow(diffTAvg, 2))) / (runningCount + cnt);
                        runningCount += cnt;
                        runningAvg = newAvg;
                        runningVar = newVar;
                    }
                    if (runningCount > 0) {
                        platforms.platformAvgMatrix.Set(j, k, runningAvg);
                        platforms.platformStdevMatrix.Set(j, k, sqrt(runningVar));
                        platforms.platformCountMatrix.Set(j, k, runningCount);
                    }
                }
            }

            endTime = time(nullptr);
            printf("thread combine time %ld sec\n", endTime - startTime);

            /*
            for(i=0; i<numPlatforms; i++){
                printf("Platform %s\n", mapistrPlatform[i].c_str());
                for(j=0; j<localVectrQuery.size(); j++){
                    size_t iGene = mapstriGenes[localVectrQuery[j]];
                    printf("Gene %s %.5f %.5f\n", localVectrQuery[j].c_str(),
                        platform_avg.Get(i, iGene), platform_stdev.Get(i,iGene));
                }
            }*/

            platforms.savePlatformDataToFiles(sArgs.dir_out_arg);
        }
    } else if (sArgs.dab_flag == 1) {

        string norm_mode = sArgs.norm_mode_arg;
        if (sArgs.norm_flag == 1 && norm_mode == "rank") {
            if (sArgs.default_type_arg == -1) {
                fprintf(stderr, "Please supply parameter --default_type\n");
                return 1;
            }
            if (sArgs.max_rank_arg == -1) {
                fprintf(stderr, "Please supply parameter --max_rank\n");
                return 1;
            }

            CDataPair Dat;
            char outFile[1024];
            fprintf(stderr, "Opening file...\n");
            //if(!Dat.Open(sArgs.dabinput_arg, false, 2, false, false, true)){
            if (!Dat.Open(sArgs.dabinput_arg, false, false, 2, false, false)) {
                cerr << "error opening file" << endl;
                return 1;
            }
            fprintf(stderr, "Finished opening file\n");
            string fileName = CMeta::Basename(sArgs.dabinput_arg);
            string fileStem = CMeta::Deextension(fileName);
            sprintf(outFile, "%s/%s.2.dab", sArgs.dir_out_arg,
                    fileStem.c_str());
            int max_rank = sArgs.max_rank_arg;
            fprintf(stderr, "Using max_rank: %d\n", max_rank);
            //cutoff, expTransform, divideNorm, subtractNorm
            //CSeekWriter::NormalizeDAB(Dat, vecstrGenes, true, false, true, false);
            //CSeekWriter::RankNormalizeDAB(Dat, vecstrGenes, max_rank, rbp_p);
            //Dat.Save(outFile);
            if (sArgs.default_type_arg == 0) {
                vector <map<utype, unsigned short>> umat;
                CSeekWriter::GetSparseRankMatrix<utype>(Dat, umat, max_rank,
                                                        vecstrGenes);
                CSeekWriter::WriteSparseMatrix<utype>(Dat, umat,
                                                      vecstrGenes, outFile);
            } else if (sArgs.default_type_arg == 1) {
                vector <map<unsigned short, unsigned short>> umat;
                CSeekWriter::GetSparseRankMatrix<unsigned short>(Dat, umat, max_rank,
                                                                 vecstrGenes);
                CSeekWriter::WriteSparseMatrix<unsigned short>(Dat, umat,
                                                               vecstrGenes, outFile);
            } else {
                fprintf(stderr, "Invalid default type!\n");
                return -1;
            }

            //TEST========================================
            /*CSeekIntIntMap d1(vecstrGenes.size());
            CSparseFlatMatrix<float> sm(0);
            float rbp_p = sArgs.rbp_p_arg;
            CSeekWriter::ReadSeekSparseMatrix<unsigned short>(outFile, sm, d1,
            max_rank, rbp_p, vecstrGenes);*/
            //============================================

            /*fprintf(stderr, "Begin\n");
            vector<unsigned short> l;
            vector<vector<float> > mat;
            CSeekWriter::ReadSparseMatrixAsArray(l, outFile);
            fprintf(stderr, "Begin 2\n");
            CSeekWriter::ReadSparseMatrix(l, mat, 0.99, vecstrGenes);*/
        }

        if (sArgs.view_flag == 1) {
            //fprintf(stderr, "Operation not implemented yet!\n");
            CDataPair Dat;
            char outFile[125];
            if (!Dat.Open(sArgs.dabinput_arg, false, false)) {
                cerr << "error opening file" << endl;
                return 1;
            }
            vector<unsigned int> veciGenes;
            veciGenes.resize(vecstrGenes.size());
            for (i = 0; i < vecstrGenes.size(); i++)
                veciGenes[i] = (unsigned int) Dat.GetGeneIndex(vecstrGenes[i]);

            unsigned int s, t;
            vector<float> d;
            for (i = 0; i < vecstrGenes.size(); i++) {
                if ((s = veciGenes[i]) == (unsigned int) -1) continue;
                for (j = 0; j < vecstrGenes.size(); j++) {
                    if ((t = veciGenes[j]) == (unsigned int) -1) continue;
                    d.push_back(Dat.Get(s, t));
                }
            }
            sort(d.begin(), d.end());
            float ff[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
            for (i = 0; i < 9; i++) {
                auto ind = (size_t)((float) d.size() * ff[i]);
                fprintf(stderr, "%.3f %.3e\n", ff[i], d[ind]);
            }
            return 1;
        }

        if (sArgs.norm_flag == 1 && norm_mode == "subtract_z") {
            if (sArgs.cutoff_value_arg == -1.0) {
                fprintf(stderr, "Please supply parameter --cutoff_value\n");
                return 1;
            }
            if (sArgs.default_type_arg == -1) {
                fprintf(stderr, "Please supply parameter --default_type\n");
                return 1;
            }

            CDataPair Dat;
            char outFile[125];
            if (!Dat.Open(sArgs.dabinput_arg, false, false)) {
                cerr << "error opening file" << endl;
                return 1;
            }
            fprintf(stderr, "Finished opening file\n");
            string fileName = CMeta::Basename(sArgs.dabinput_arg);
            string fileStem = CMeta::Deextension(fileName);
            sprintf(outFile, "%s/%s.2.dab", sArgs.dir_out_arg, fileStem.c_str());
            //expTransform, divideNorm, subtractNorm
            CSeekWriter::NormalizeDAB(Dat, vecstrGenes, false, false, true);
            float cutoff = sArgs.cutoff_value_arg;
            if (sArgs.default_type_arg == 0) {
                //unsigned int
                vector <map<utype, unsigned short>> umat;
                CSeekWriter::ConvertToSparseMatrix<utype>(Dat, umat, vecstrGenes,
                                                          cutoff);
                CSeekWriter::WriteSparseMatrix<utype>(Dat, umat, vecstrGenes,
                                                      outFile);
            } else if (sArgs.default_type_arg == 1) {
                //unsigned short
                vector <map<unsigned short, unsigned short>> umat;
                CSeekWriter::ConvertToSparseMatrix<unsigned short>(Dat, umat, vecstrGenes,
                                                                   cutoff);
                CSeekWriter::WriteSparseMatrix<unsigned short>(Dat, umat, vecstrGenes,
                                                               outFile);
            } else {
                fprintf(stderr, "Error, unsupported type --default_type_arg\n");
                return 1;
            }
        }

        if (sArgs.norm_flag == 1 && norm_mode == "topological_overlap") {
            omp_set_num_threads(8);
            CDataPair Dat;
            char outFile[125];
            if (!Dat.Open(sArgs.dabinput_arg, false, false)) {
                cerr << "error opening file" << endl;
                return 1;
            }
            fprintf(stderr, "Finished opening file\n");
            string fileName = CMeta::Basename(sArgs.dabinput_arg);
            string fileStem = CMeta::Deextension(fileName);
            sprintf(outFile, "%s/%s.dab", sArgs.dir_out_arg, fileStem.c_str());
            CSeekWriter::TopologicalOverlap(Dat, vecstrGenes);
            Dat.Save(outFile);
        }

        if (sArgs.gavg_flag == 1) {
            bool logit = false;
            if (sArgs.logit_flag == 1) logit = true;

            CDataPair Dat;
            char outFile[125];
            if (!Dat.Open(sArgs.dabinput_arg, false, false)) {
                cerr << "error opening file" << endl;
                return 1;
            }
            vector<float> vecGeneAvg;
            string fileName = CMeta::Basename(sArgs.dabinput_arg);
            string fileStem = CMeta::Deextension(fileName);
            sprintf(outFile, "%s/%s.gavg", sArgs.dir_out_arg,
                    fileStem.c_str());
            CSeekWriter::GetGeneAverage(Dat, vecstrGenes, vecGeneAvg, logit,
                                        sArgs.top_avg_percent_arg);

            //DEBUGGING
            //for(i=0; i<vecGeneAvg.size(); i++){
            //	fprintf(stderr, "%s\t%.3f\n", vecstrGenes[i].c_str(), vecGeneAvg[i]);
            //}

            CSeekTools::WriteArray(outFile, vecGeneAvg);
        } else if (sArgs.gpres_flag == 1) {
            CDataPair Dat;
            char outFile[125];
            if (!Dat.Open(sArgs.dabinput_arg, false, false)) {
                cerr << "error opening file" << endl;
                return 1;
            }
            vector<char> vecGenePresence;
            string fileName = CMeta::Basename(sArgs.dabinput_arg);
            string fileStem = CMeta::Deextension(fileName);
            sprintf(outFile, "%s/%s.gpres", sArgs.dir_out_arg,
                    fileStem.c_str());
            CSeekWriter::GetGenePresence(Dat, vecstrGenes, vecGenePresence);

            //DEBUGGING
            //for(i=0; i<vecGenePresence.size(); i++){
            //	fprintf(stderr, "%s\t%d\n", vecstrGenes[i].c_str(), vecGenePresence[i]);
            //}

            CSeekTools::WriteArray(outFile, vecGenePresence);
        }

    }

    if (sArgs.dabset_flag == 1) {
        NormMode n;
        string norm_mode = sArgs.norm_mode_arg;
        if (norm_mode == "NA") {
            fprintf(stderr, "Error, please supply --norm_mode\n");
            return 1;
        }
        if (norm_mode == "subtract_z") {
            n = Z_NORM;
            if (sArgs.exp_arg == -1) {
                fprintf(stderr, "Error, please supply --exp\n");
                return 1;
            }
        } else if (norm_mode == "rank") {
            n = RANK_NORM;
            if (sArgs.rbp_p_arg == -1 || sArgs.max_rank_arg == -1) {
                fprintf(stderr, "Error, Need both --rbp_p and --max_rank\n");
                return 1;
            }
        }

        vector <string> dab_list;
        int numGenes = vecstrGenes.size();
        string dab_dir = sArgs.dab_dir_arg;
        string outdir = sArgs.dir_out_arg;
        string outdab = sArgs.out_dab_arg;
        if (dab_dir == "NA" || outdab == "NA") {
            fprintf(stderr, "Required args: dab_dir, out_dab");
            return -1;
        }
        string dw = sArgs.dataset_w_arg;

        CSeekTools::ReadListOneColumn(sArgs.dablist_arg, dab_list);

        vector<float> dweight;
        if (dw == "NA") {
            CSeekTools::InitVector(dweight, dab_list.size(), (float) 1.0);
        } else {
            CSeekTools::ReadArray(dw.c_str(), dweight);
            if (dweight.size() != dab_list.size()) {
                fprintf(stderr, "Incorrect Size!\n");
                return -1;
            }
        }

        float rbp_p = sArgs.rbp_p_arg;
        int max_rank = sArgs.max_rank_arg;
        float exp = sArgs.exp_arg;
        fprintf(stderr, "Using rbp_p: %.3f, max_rank: %d, exp: %.3f\n", rbp_p, max_rank, exp);

        if (sArgs.default_type_arg == 0) {
            CalculateMatrix<utype>(n, GENE_PRESENCE_VECTOR, dab_list, dab_dir, outdir,
                                   outdab, vecstrGenes, dweight, rbp_p, max_rank, exp);
            CalculateMatrix<utype>(n, COUNT_MATRIX, dab_list, dab_dir, outdir, outdab,
                                   vecstrGenes, dweight, rbp_p, max_rank, exp);
            //CalculateMatrix<utype>(n, WEIGHTSUM_MATRIX, dab_list, dab_dir, outdir, outdab,
            //vecstrGenes, dweight, rbp_p, max_rank, exp);
            CalculateMatrix<utype>(n, DISTANCE_MATRIX, dab_list, dab_dir, outdir, outdab,
                                   vecstrGenes, dweight, rbp_p, max_rank, exp);
        } else if (sArgs.default_type_arg == 1) {
            CalculateMatrix<unsigned short>(n, GENE_PRESENCE_VECTOR, dab_list, dab_dir, outdir,
                                            outdab, vecstrGenes, dweight, rbp_p, max_rank, exp);
            CalculateMatrix<unsigned short>(n, COUNT_MATRIX, dab_list, dab_dir, outdir, outdab,
                                            vecstrGenes, dweight, rbp_p, max_rank, exp);
            //CalculateMatrix<unsigned short>(n, WEIGHTSUM_MATRIX, dab_list, dab_dir, outdir, outdab,
            //vecstrGenes, dweight, rbp_p, max_rank);
            CalculateMatrix<unsigned short>(n, DISTANCE_MATRIX, dab_list, dab_dir, outdir, outdab,
                                            vecstrGenes, dweight, rbp_p, max_rank, exp);
        } else {
            fprintf(stderr, "Invalid default type!\n");
            return -1;
        }
    }

    if (sArgs.combined_dab_flag == 1) {
        utype numGenes = vecstrGenes.size();
        string dab_base = sArgs.dab_basename_arg;
        string dab_dir = sArgs.dab_dir2_arg;
        string outdir = sArgs.dir_out_arg;

        if (dab_base == "NA" || dab_dir == "NA") {
            fprintf(stderr, "Required args: dab_dir2, dab_basename");
            return -1;
        }

        string file1 = dab_dir + "/" + dab_base + ".half";
        string file2 = dab_dir + "/" + dab_base + ".pair_count";
        string file3 = dab_dir + "/" + dab_base + ".gene_count";
        string file4 = dab_dir + "/" + dab_base + ".pair_weightsum"; //optional

        vector <utype> gene_count;
        CSeekTools::ReadArray(file3.c_str(), gene_count);

        //normalize by counts
        float max_count = (float) *max_element(gene_count.begin(), gene_count.end());
        float min_count_required = max_count * 0.50;
        CSeekIntIntMap d1(vecstrGenes.size());
        vector <string> validGenes;
        int numValidGenes = 0;
        for (i = 0; i < vecstrGenes.size(); i++) {
            if (gene_count[i] < min_count_required) continue;
            numValidGenes++;
            validGenes.push_back(vecstrGenes[i]);
            d1.Add(i);
        }

        fprintf(stderr, "%.2f %.2f %d\n", max_count, min_count_required,
                numValidGenes);

        uint32_t dim;
        ifstream istm1, istm2, istm3;
        istm1.open(file1.c_str(), ios_base::binary);
        istm1.read((char *) (&dim), sizeof(dim));
        istm2.open(file2.c_str(), ios_base::binary);
        istm2.read((char *) (&dim), sizeof(dim));
        istm3.open(file4.c_str(), ios_base::binary); //test if file3 exists
        bool enableWeight = false;
        if (!istm3.fail()) {
            istm3.read((char *) (&dim), sizeof(dim)); //read dimension if exists
            enableWeight = true;
        }

        /*CHalfMatrix<float> res;
        res.Initialize(dim);
        if(dim!=numGenes){
            fprintf(stderr, "Error in number of dimensions\n");
            return -1;
        }*/

        fprintf(stderr, "Begin Reading\n");
        if (enableWeight) {
            fprintf(stderr, "Enable weight: True\n");
        } else {
            fprintf(stderr, "Enable weight: False\n");
        }
        CDat CD;
        CD.Open(validGenes);
        for (i = 0; i < validGenes.size(); i++)
            for (j = i + 1; j < validGenes.size(); j++)
                CD.Set(i, j, CMeta::GetNaN());

        auto *adScores = new float[dim - 1];
        auto *uCounts = new unsigned short[dim - 1];
        auto *fSum = new float[dim - 1];

        for (i = 0; (i + 1) < dim; i++) {
            if (i % 2000 == 0) {
                fprintf(stderr, "Finished %d\n", i);
            }
            istm1.read((char *) adScores, sizeof(*adScores) * (dim - i - 1));
            istm2.read((char *) uCounts, sizeof(*uCounts) * (dim - i - 1));
            if (enableWeight)
                istm3.read((char *) fSum, sizeof(*fSum) * (dim - i - 1));

            utype gi = d1.GetForward(i);
            if (CSeekTools::IsNaN(gi))
                continue;
            for (j = 0; j < dim - i - 1; j++) {
                //reverse of HalfIndex function from halfmatrixi.h
                //to get the original coordinate for j:
                //original j = new j + 1 + i
                utype gj = d1.GetForward(j + i + 1);
                if (CSeekTools::IsNaN(gj))
                    continue;
                if ((float) uCounts[j] < min_count_required)
                    continue;
                if (enableWeight)
                    adScores[j] /= fSum[j];
                else
                    adScores[j] /= (float) uCounts[j];
                CD.Set(gi, gj, adScores[j]);
            }
            //res.Set(i, adScores);
        }
        delete[] adScores;
        delete[] uCounts;
        delete[] fSum;
        istm1.close();
        istm2.close();
        if (enableWeight)
            istm3.close();
        string file_out = outdir + "/" + dab_base + ".dab";
        CD.Save(file_out.c_str());

    }

#ifdef WIN32
    pthread_win32_process_detach_np( );
#endif // WIN32
    return 0;
}

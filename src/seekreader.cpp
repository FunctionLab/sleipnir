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
#include "seekreader.h"
#include "seekerror.h"
#include <regex>

namespace Sleipnir {

    string CSeekTools::ConvertInt(const int &number) {
        stringstream ss;//create a stringstream
        ss << number;//add number to the stream
        return ss.str();//return a string with the contents of the stream
    }

    bool CSeekTools::IsNaN(const utype &v) {
        if (v == MAX_UTYPE) return true;
        return false;
    }

    utype CSeekTools::GetNaN() {
        return MAX_UTYPE;
    }

    bool CSeekTools::ReadDatabaselets(const vector<CDatabase *> &DB,
                                      const size_t &iGenes, const size_t &iDatasets,
                                      const vector <vector<string>> &vecstrAllQuery,
                                      vector<CSeekDataset *> &vc, const map <string, utype> &mapstriGenes,
                                      const vector <vector<string>> &dbDatasets,
                                      const map <string, utype> &mapstriDatasets,
            //network mode (data sent to client)
                                      const int &iClient, const bool &bNetwork) {

        //requires LoadDatabase to be called beforehand
        size_t i, j, k;
        vector<char> cAllQuery;

        // cAllQuery indicates which genes are part of the query by gene map index
        // It has one byte per gene in the collection if that byte is 0 that gene
        // is not present in the query, if 1 it is present.
        CSeekTools::InitVector(cAllQuery, iGenes, (char) 0);

        // loop through all queries
        for (i = 0; i < vecstrAllQuery.size(); i++) {
            // in each query loop through all genes
            for (j = 0; j < vecstrAllQuery[i].size(); j++) {
                // is the gene name contained in the gene_map file
                if (mapstriGenes.find(vecstrAllQuery[i][j]) == mapstriGenes.end()) continue;
                // get the gene_map index value associated with the gene name
                utype k = mapstriGenes.find(vecstrAllQuery[i][j])->second;
                // indicate this gene (by gene index) is present
                cAllQuery[k] = 1;
            }
        }

        // create a list of just the query genes from the cAllQuery presence map
        vector <utype> allQ;
        for (i = 0; i < cAllQuery.size(); i++) if (cAllQuery[i] == 1) allQ.push_back(i);
        allQ.resize(allQ.size());

        //for now
        // reset all seek datasets that were previously loaded for other queries
        for (i = 0; i < iDatasets; i++) {
            if (vc[i]->GetDBMap() != NULL) {
                vc[i]->DeleteQueryBlock();
            }
        }

        int ret; //system call return

        fprintf(stderr, "Initializing query map\n");
        //ret = system("date +%s%N 1>&2");
        /*if(bNetwork && CSeekNetwork::Send(iClient, "Initializing query map")==-1){
            fprintf(stderr, "Error sending client message\n");
            return false;
        }*/

#pragma omp parallel for \
    shared(allQ) private(i) schedule(dynamic)
        for (i = 0; i < iDatasets; i++) {
            // initialize each dataset with the presence map of query genes
            vc[i]->InitializeQueryBlock(allQ);
        }

        fprintf(stderr, "Done initializing query map\n");
        //ret = system("date +%s%N 1>&2");
        /*if(bNetwork && CSeekNetwork::Send(iClient,
            "Done initializing query map")==-1){
            fprintf(stderr, "Error sending client message\n");
            return false;
        }*/

        fprintf(stderr, "Reading %lu query genes' correlations\n",
                allQ.size());
        ret = system("date +%s%N 1>&2");
        if (bNetwork && CSeekNetwork::Send(iClient, "Reading " +
                                                    CSeekTools::ConvertInt(allQ.size()) +
                                                    " query genes' correlations (estimated time " +
                                                    CSeekTools::ConvertInt(allQ.size() * 2) +
                                                    "s; stay on this page)") == -1) {
            fprintf(stderr, "Error sending client message\n");
            return false;
        }

        size_t m, d;

        for (i = 0; i < allQ.size(); i++) {
            m = allQ[i];
            for (d = 0; d < DB.size(); d++) { //number of CDatabase collections
                vector<unsigned char> Qi;
                // reads into Qi this m gene's data (pair values for all datasets)
                //  in same layout as on disk
                if (!DB[d]->GetGene(m, Qi)) {
                    cerr << "Gene does not exist" << endl;
                    continue;
                }
                utype db;
                CSeekIntIntMap *qu = NULL;
                unsigned char **r = NULL;
                vector <utype> vecDatasetID;
                // loop through each dataset in this DB collection
                for (j = 0; j < dbDatasets[d].size(); j++) {
                    // get the index number of this dataset name add it to vecDatasetId
                    utype qq = mapstriDatasets.find(dbDatasets[d][j])->second;
                    vecDatasetID.push_back(qq);
                }
#pragma omp parallel for \
            shared(Qi) private(j, k) \
            firstprivate(m, qu, r, db) schedule(dynamic)
                // loop through each dataset index
                for (j = 0; j < vecDatasetID.size(); j++) {
                    // get the gene presence map for the dataset into qu
                    if ((qu = vc[vecDatasetID[j]]->GetDBMap()) == NULL) continue;
                    // get the index position of the query gene m in the query (db = qu->GetForward(m))
                    if (CSeekTools::IsNaN(db = (qu->GetForward(m)))) continue;
                    // fill the qgene-gene correlation values
                    for (r = vc[vecDatasetID[j]]->GetMatrix(), k = 0; k < iGenes; k++)
                        // data order is for query gene (db), all correlations with other genes for this dataset (j)
                        //   followed by all gene correlations for the next dataset, etc.
                        // GRW - change needed for different seek_db layout
                        r[db][k] = Qi[k * vecDatasetID.size() + j];
                }
                Qi.clear();
            }
        }

        fprintf(stderr, "Finished reading query genes' correlations\n");
        ret = system("date +%s%N 1>&2");
        if (bNetwork && CSeekNetwork::Send(iClient,
                                           "Finished reading. Now doing search (estimated time " +
                                           CSeekTools::ConvertInt(allQ.size()) +
                                           "s; stay on this page)") == -1) {
            fprintf(stderr, "Error sending client message\n");
            return false;
        }

        return true;
    }

    bool CSeekTools::ReadQuantFile(const string &strFile, vector<float> &quant,
                                   const int lineSize) {
        return CSeekTools::ReadQuantFile(strFile.c_str(), quant, lineSize);
    }

    bool CSeekTools::ReadQuantFile(const char *file, vector<float> &quant, const int lineSize) {
        ifstream ifsm;
        ifsm.open(file);
        char acBuffer[lineSize];
        utype c_iBuffer = lineSize;
        vector <string> vecstrLine;

        ifsm.getline(acBuffer, c_iBuffer - 1);
        //fprintf(stderr, "%s\n", acBuffer);
        CMeta::Tokenize(acBuffer, vecstrLine, " ", false);
        quant.clear();
        utype i;
        for (i = 0; i < vecstrLine.size(); i++) {
            quant.push_back(atof(vecstrLine[i].c_str()));
            //fprintf(stderr, "%.5f\n", atof(vecstrLine[i].c_str()));
        }
        quant.resize(quant.size());
        ifsm.close();
        return true;
    }

    // Called by the per-query version of CSeekCentral::Initialize()
    bool CSeekTools::LoadDatabase(const vector<CDatabase *> &DB,
                                  const size_t &iGenes, const size_t &iDatasets,
                                  vector<CSeekDataset *> &vc, const vector<CSeekDataset *> &vc_src,
                                  vector <CSeekPlatform> &vp,
                                  const vector <string> &vecstrDatasets,
                                  const map <string, string> &mapstrstrDatasetPlatform,
                                  const map <string, utype> &mapstriPlatform) {

        size_t i, j, k;

        vc.clear();
        vc.resize(iDatasets);

        int ret; //system call returns

        fprintf(stderr, "Initializing gene map\n");
        //ret = system("date +%s%N 1>&2");
#pragma omp parallel for \
    private(i) schedule(dynamic)
        for (i = 0; i < iDatasets; i++) {
            vc[i] = new CSeekDataset();
            vc[i]->Copy(vc_src[i]);
            string strFileStem = vecstrDatasets[i];
            string strPlatform =
                    mapstrstrDatasetPlatform.find(strFileStem)->second;
            utype platform_id = mapstriPlatform.find(strPlatform)->second;
            vc[i]->SetPlatform(vp[platform_id]);
        }

        fprintf(stderr, "Done initializing gene map\n");
        //ret = system("date +%s%N 1>&2");
        return true;
    }

    // Called by in-common (startup) version of CSeekCentral::Initialize()
    bool CSeekTools::LoadDatabase(const vector<CDatabase *> &DB,
                                  const size_t &iGenes, const size_t &iDatasets,
                                  const vector<CSeekDBSetting *> &DBSetting,
                                  const vector <string> &vecstrDatasets,
                                  const map <string, string> &mapstrstrDatasetPlatform,
                                  const map <string, utype> &mapstriPlatform, vector <CSeekPlatform> &vp,
                                  vector<CSeekDataset *> &vc, const vector <vector<string>> &dbDataset,
                                  const map <string, utype> &mapstriDataset,
                                  const bool bVariance, const bool bCorrelation) {

        size_t i, j, k;
        vc.clear();
        vc.resize(iDatasets);

        if (bCorrelation) {
            for (i = 0; i < DB.size(); i++) {
                if (DBSetting[i]->sinfoDir == "NA") {
                    fprintf(stderr, "sinfo parameter must be given.\n");
                    return false;
                }
            }
        }

        if (bVariance) {
            for (i = 0; i < DB.size(); i++) {
                if (DBSetting[i]->gvarDir == "NA") {
                    fprintf(stderr, "gene variance parameter must be given.\n");
                    return false;
                }
            }
        }

        int ret; //system call return

        fprintf(stderr, "Start reading average and presence files\n");
        ret = system("date +%s%N 1>&2");
        for (i = 0; i < DB.size(); i++) {
            const vector <string> &dset = dbDataset[i];
            string strPrepInputDirectory = DBSetting[i]->prepDir;
            string strGvarInputDirectory = DBSetting[i]->gvarDir;
            string strSinfoInputDirectory = DBSetting[i]->sinfoDir;

            for (j = 0; j < dset.size(); j++) {
                utype d = mapstriDataset.find(dset[j])->second;
                vc[d] = new CSeekDataset();
                string strFileStem = dset[j];
                string strAvgPath = strPrepInputDirectory + "/" +
                                    strFileStem + ".gavg";
                string strPresencePath = strPrepInputDirectory + "/" +
                                         strFileStem + ".gpres";
                vc[d]->ReadGeneAverage(strAvgPath);
                vc[d]->ReadGenePresence(strPresencePath);
                if (bVariance) {
                    string strVariancePath = strGvarInputDirectory + "/" +
                                             strFileStem + ".gexpvar";
                    vc[d]->ReadGeneVariance(strVariancePath);
                }
                if (bCorrelation) {
                    string strSinfoPath = strSinfoInputDirectory + "/" +
                                          strFileStem + ".sinfo";
                    vc[d]->ReadDatasetAverageStdev(strSinfoPath);
                }
                string strPlatform =
                        mapstrstrDatasetPlatform.find(strFileStem)->second;
                utype platform_id = mapstriPlatform.find(strPlatform)->second;
                vc[d]->SetPlatform(vp[platform_id]);
            }
        }

        fprintf(stderr, "Done reading average and presence files\n");
        ret = system("date +%s%N 1>&2");

        fprintf(stderr, "Initializing gene map\n");
        ret = system("date +%s%N 1>&2");
#pragma omp parallel for \
    shared(vc, iDatasets) private(i) schedule(dynamic)
        for (i = 0; i < iDatasets; i++) vc[i]->InitializeGeneMap();

        fprintf(stderr, "Done initializing gene map\n");
        ret = system("date +%s%N 1>&2");
        return true;
    }

    bool CSeekTools::ReadListTwoColumns(const string &strFile,
                                        vector <string> &vecstrList1, vector <string> &vecstrList2,
                                        const int lineSize) {
        return CSeekTools::ReadListTwoColumns(strFile.c_str(),
                                              vecstrList1, vecstrList2, lineSize);
    }

    bool CSeekTools::ReadListTwoColumns(const char *file,
                                        vector <string> &vecstrList1, vector <string> &vecstrList2,
                                        const int lineSize) {
        ifstream ifsm;
        ifsm.open(file);
        if (!ifsm.is_open()) {
            fprintf(stderr, "Error opening file %s\n", file);
            return false;
        }
        char acBuffer[lineSize];
        utype c_iBuffer = lineSize;
        vecstrList1.clear();
        vecstrList2.clear();
        regex ws_re("\\s+");  // match on one or more space type values [tab, space]

        while (!ifsm.eof()) {
            ifsm.getline(acBuffer, c_iBuffer - 1);
            if (acBuffer[0] == 0) break;
            acBuffer[c_iBuffer - 1] = 0;
            string line(acBuffer);
            vector<string> tok {
                sregex_token_iterator(line.begin(), line.end(), ws_re, -1), {}
            };
            if (tok.size() != 2) {
                throw init_error(FILELINE + "Expecting two values per line in file: " + file + ", line: " + to_string(vecstrList1.size()));
            }
            vecstrList1.push_back(tok[0]);
            vecstrList2.push_back(tok[1]);
        }
        vecstrList1.resize(vecstrList1.size());
        vecstrList2.resize(vecstrList2.size());
        ifsm.close();
        return true;
    }

    bool CSeekTools::ReadListOneColumn(const string &strFile,
                                       vector <string> &vecstrList, CSeekStrIntMap &mapstriList, const int lineSize) {
        return CSeekTools::ReadListOneColumn(strFile.c_str(),
                                             vecstrList, mapstriList, lineSize);
    }


    bool CSeekTools::ReadListOneColumn(const char *file,
                                       vector <string> &vecstrList, CSeekStrIntMap &mapstriList, const int lineSize) {
        ifstream ifsm;
        ifsm.open(file);
        if (!ifsm.is_open()) {
            fprintf(stderr, "Error opening file %s\n", file);
            return false;
        }

        char acBuffer[lineSize];
        utype c_iBuffer = lineSize;
        vecstrList.clear();

        int i = 0;
        while (!ifsm.eof()) {
            ifsm.getline(acBuffer, c_iBuffer - 1);
            if (acBuffer[0] == 0) break;
            acBuffer[c_iBuffer - 1] = 0;
            string line = acBuffer;
            vecstrList.push_back(line);
            mapstriList.Set(line, i);
            i++;
        }
        vecstrList.resize(vecstrList.size());
        ifsm.close();
        return true;
    }

    bool CSeekTools::ReadMultipleQueries(const string &strFile,
                                         vector <vector<string>> &qList, const int lineSize) {
        return CSeekTools::ReadMultipleQueries(strFile.c_str(), qList, lineSize);
    }

    bool CSeekTools::ReadMultipleQueries(const char *file,
                                         vector <vector<string>> &qList, const int lineSize) {
        qList.clear();
        FILE *infile;
        if ((infile = fopen(file, "r")) == NULL) {
            fprintf(stderr, "Error opening file %s\n", file);
            return false;
        }

        char *acBuffer;
        int MAX_CHAR_PER_LINE = lineSize;
        int lineLen = MAX_CHAR_PER_LINE;
        acBuffer = (char *) malloc(lineLen);
        // This is looping through the file to determine if lineLen is
        //   ever exceeded and increasing/realloc lineLen if needed
        while (fgets(acBuffer, lineLen, infile) != NULL) {
            while (strlen(acBuffer) == lineLen - 1) {
                int len = strlen(acBuffer);
                fseek(infile, -len, SEEK_CUR);
                lineLen += MAX_CHAR_PER_LINE;
                acBuffer = (char *) realloc(acBuffer, lineLen);
                char *ret = fgets(acBuffer, lineLen, infile);
            }
        }
        rewind(infile);

        // Use regex to separate tokens on any number or type of spaces
        regex ws_re("\\s+");
        while (fgets(acBuffer, lineLen, infile) != NULL) {
            char *p = strtok(acBuffer, "\n");
            if (p == NULL) continue;  // skip empty lines
            string line(p);
            vector<string> tok {
                sregex_token_iterator(line.begin(), line.end(), ws_re, -1), {}
            };
            qList.push_back(tok);
        }
        qList.resize(qList.size());
        free(acBuffer);

        fclose(infile);
        return true;
    }

    bool CSeekTools::ReadMultipleNotQueries(const char *file,
                                            vector<vector<vector<string>>> &qList,
                                            const int lineSize)
    {
        qList.clear();
        FILE *infile;
        if ((infile = fopen(file, "r")) == NULL)
        {
            fprintf(stderr,
                    "Error opening file %s\n", file);
            return false;
        }

        char *acBuffer;
        int MAX_CHAR_PER_LINE = lineSize;
        int lineLen = MAX_CHAR_PER_LINE;
        acBuffer = (char *)malloc(lineLen);
        while (fgets(acBuffer, lineLen, infile) != NULL) {
            while (strlen(acBuffer) == lineLen - 1) {
                int len = strlen(acBuffer);
                fseek(infile,
                      -len, SEEK_CUR);
                lineLen +=
                    MAX_CHAR_PER_LINE;
                acBuffer = (char *)realloc(acBuffer, lineLen);
                char *ret = fgets(acBuffer, lineLen, infile);
            }
        }
        rewind(infile);

        while (fgets(acBuffer, lineLen, infile) != NULL) {
            char *p = strtok(acBuffer, "\n");
            if (p == NULL) continue;  // skip empty lines
            vector<vector<string>> aQ;
            vector<string> tok;
            CMeta::Tokenize(p, tok, "|");
            aQ.resize(tok.size());
            int i;
            for ( i = 0; i < tok. size(); i++) {
                vector<string> tmp;
                CMeta::Tokenize(tok[i].c_str(), tmp, " ");
                aQ[i] = tmp;
            }
            qList.push_back(aQ);
        }
        qList.resize(qList.size());
        free(acBuffer);
        fclose(infile);
        return true;
    }

bool CSeekTools::ReadMultiGeneOneLine(const string &strFile,
                                      vector <string> &list, const int lineSize) {
    return CSeekTools::ReadMultiGeneOneLine(strFile.c_str(), list, lineSize);
}

bool CSeekTools::ReadMultiGeneOneLine(const char *file,
                                      vector <string> &list, const int lineSize) {
    list.clear();
    ifstream ifsm;
    ifsm.open(file);
    if (!ifsm.is_open()) {
        fprintf(stderr, "Error opening file %s\n", file);
        return false;
    }

    char *acBuffer = (char *) malloc(lineSize);
    int c_iBuffer = lineSize;
    int i = 0;

    ifsm.getline(acBuffer, c_iBuffer - 1);
    acBuffer[c_iBuffer - 1] = 0;
    vector <string> tok;
    CMeta::Tokenize(acBuffer, tok, " ");
    for (i = 0; i < tok.size(); i++) {
        list.push_back(tok[i]);
    }

    list.resize(list.size());
    ifsm.close();
    free(acBuffer);
    return true;
}

bool CSeekTools::ReadListOneColumn(const string &strFile,
                                   vector <string> &vecstrList, const int lineSize) {
    return CSeekTools::ReadListOneColumn(strFile.c_str(), vecstrList, lineSize);
}

bool CSeekTools::ReadListOneColumn(const char *file,
                                   vector <string> &vecstrList, const int lineSize) {
    ifstream ifsm;
    ifsm.open(file);

    if (!ifsm.is_open()) {
        fprintf(stderr, "Error opening file %s\n", file);
        return false;
    }
    char acBuffer[lineSize];
    utype c_iBuffer = lineSize;
    vecstrList.clear();
    int i = 0;
    while (!ifsm.eof()) {
        ifsm.getline(acBuffer, c_iBuffer - 1);
        if (acBuffer[0] == 0) break;
        acBuffer[c_iBuffer - 1] = 0;
        string line = acBuffer;
        vecstrList.push_back(line);
    }
    vecstrList.resize(vecstrList.size());
    ifsm.close();
    return true;
}

}

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
#include <filesystem>
#include <unistd.h>
#include <boost/algorithm/string/predicate.hpp>
#include "seekerror.h"
#include "PclQuery.h"

using namespace Sleipnir;
using namespace std;


string stripPclExtensions(string pclDsetName) {
    // Remove .pcl extensions
    size_t pos;
    pos = pclDsetName.find(".pcl");
    if (pos != string::npos) {
        pclDsetName = pclDsetName.substr(0, pos);
    }
    return pclDsetName;
}


string stripBinExtensions(string pclDsetName) {
    // Remove .bin extensions
    size_t pos;
    pos = pclDsetName.find(".bin");
    if (pos != string::npos) {
        pclDsetName = pclDsetName.substr(0, pos);
    }
    return pclDsetName;
}


void *do_pcl_query(void *th_arg) {
    struct pcl_thread_data *my = (struct pcl_thread_data *) th_arg;
    vector <string> &datasetNames = my->datasetNames;
    vector <string> &geneName = my->geneNames;
    vector <string> &queryName = my->queryGeneNames;
    CSeekCentral *cc = my->seekCentral;
    bool outputNormalized = my->outputNormalized;
    bool outputCoexpression = my->outputCoexpression;
    bool outputQueryCoexpression = my->outputQueryCoexpression;
    bool outputExpression = my->outputExpression;
    bool outputQueryExpression = my->outputQueryExpression;
    int new_fd = my->new_fd;
    float RBP_P = my->rbp_p;

    /* A lock_guard will automatically set the completionFlag
     *  to true when the scope exits. This will allow other threads
     *  such as the main thread to know the execution has completed.
     */
    BoolFlag completionFlag(my->isComplete);
    lock_guard<BoolFlag> flag_lock(completionFlag);

    fprintf(stderr, "start processing...\n");

    string pcl_input_dir = cc->m_vecDBSetting[0]->pclDir;

    // GW New Caching Algo
    int numDatasets = datasetNames.size();
    // A vector to hold the datasets we need for the query
    vector<PclPtrS> dsetPcls(datasetNames.size(), nullptr);
    // Check for datasets already cached before adding new ones which might evict ones we need already loaded
    // for (auto & dataset : datasetNames) {
    for (int idx = 0; idx < numDatasets; idx++) {
        // strip any trailing .pcl or .pcl.bin for cache lookup
        // Note: don't change the datasetNames[] themselves as those will be needed for opening
        //  other types of files and there can be many types of extension.
        string dsetCacheName = stripPclExtensions(datasetNames[idx]);
        dsetCacheName = stripBinExtensions(dsetCacheName);
        if (my->pclCache->get(dsetCacheName, dsetPcls[idx]) == false) {
            dsetPcls[idx] = nullptr;
        }
    }
    for (int idx = 0; idx < numDatasets; idx++) {
        string dsetCacheName = stripPclExtensions(datasetNames[idx]);
        dsetCacheName = stripBinExtensions(dsetCacheName);
        if (dsetPcls[idx] == nullptr) {
            // check cache again in case it was added in the mean time
            if (my->pclCache->get(dsetCacheName, dsetPcls[idx]) == true) {
                continue;
            }
            // open the new pcl in the cache
            bool res;
            // TODO - Note, opening CPCP with /pcl instead of /pclbin gives
            //  different expression results - bug in PCLServer algorithm?
            // Check if pclbin file name is the dsetname.bin
            string pcl_path = pcl_input_dir + "/" + datasetNames[idx] + ".bin";
            if (!filesystem::exists(pcl_path)) {
                // if datasetName doesn't have .pcl in it, try adding that
                if (datasetNames[idx].find(".pcl") == string::npos) {
                    pcl_path = pcl_input_dir + "/" + datasetNames[idx] + ".pcl.bin";
                    if (!filesystem::exists(pcl_path)) {
                        cerr << "Failed to find PCL file (.bin or .pcl.bin) " << pcl_path << endl;
                        throw argument_error("Failed to find PCL (.bin or .pcl.bin) " + pcl_path);
                    }
                }
            }
            PclPtrS dsetPcl = make_shared<CPCL>();
            res = dsetPcl->Open(pcl_path.c_str());
            if (res == false) {
                cerr << "Failed to open CPCL " << pcl_path << endl;
                throw argument_error("Failed to open CPCL " + pcl_path);
            }
            my->pclCache->set(dsetCacheName, dsetPcl);
            dsetPcls[idx] = move(dsetPcl);
        }
    }

    int genes = geneName.size();
    int queries = queryName.size(); //(EXTRA)
    int datasets = datasetNames.size();

    vector<float> vecG, vecQ, vecCoexpression, vecqCoexpression;
    vector<float> sizeD;

    sizeD.resize(datasets);

    vector <vector<float>> d_vecG, d_vecQ, d_vecCoexpression, d_vecqCoexpression;
    d_vecG.resize(datasets);
    d_vecQ.resize(datasets);
    d_vecCoexpression.resize(datasets);
    d_vecqCoexpression.resize(datasets);

    //Qian added
    //vector<CSeekDataset *> vd; //(EXTRA)
    //vd.resize(vc.size()); //(EXTRA)
    //CFullMatrix<float> *vCoexpression = new CFullMatrix<float>(); //(EXTRA)
    //vCoexpression->Initialize(genes, datasets); //(EXTRA)

    fprintf(stderr, "Reading data...\n");

    float NaN = -9999;

    map<string, utype> & platformMap = cc->m_seekPlatforms.getPlatformMap();
    vector<CSeekPlatform> & seekPlatforms = cc->m_seekPlatforms.getCSeekPlatforms();

    size_t i;

    bool hasError = false;
    string anError;

#pragma omp parallel for \
    private(i) \
    firstprivate(genes, queries, datasets, outputCoexpression, outputQueryCoexpression, outputNormalized, outputExpression, \
    outputQueryExpression, NaN) \
    shared(sizeD, datasetNames, queryName, geneName, d_vecG, d_vecQ, d_vecCoexpression, d_vecqCoexpression, \
    cc, platformMap, seekPlatforms, hasError, anError) \
    schedule(dynamic)
    for (i = 0; i < datasets; i++) {
        // CPCL *pp = vc[i];
        PclPtrS pp = dsetPcls[i];
        size_t j, k;
        int ps = pp->GetExperiments();
        int gs = pp->GetExperiments();

        CFullMatrix<float> *fq = NULL;

        CFullMatrix<float> *fq_RBP = NULL;

        CFullMatrix<float> *ff = NULL;
        vector<float> vCoexpression;
        vCoexpression.resize(genes);
        vector<float> vqCoexpression;
        vqCoexpression.resize(queryName.size());

        CSeekDataset *vd = NULL;
        CSeekPlatform *pl = NULL;
        sizeD[i] = (float) ps;
        d_vecG[i] = vector<float>();
        d_vecQ[i] = vector<float>();
        d_vecCoexpression[i] = vector<float>();
        d_vecqCoexpression[i] = vector<float>();

        set <string> absent;

        if (outputQueryExpression || outputQueryCoexpression) {
            fq = new CFullMatrix<float>();
            fq->Initialize(queryName.size(), ps);
        }

        if (outputCoexpression || outputQueryCoexpression) {
            vd = new CSeekDataset();
            // string strFileStem = datasetNames[i].substr(0, datasetNames[i].find(".bin")); //for human-SEEK
            // string strFileStem = datasetNames[i]; //for model-organism-SEEK
            string dsetName = datasetNames[i];
            dsetName = stripBinExtensions(dsetName);

            auto dbid_iter = cc->roAttr->m_mapstrintDatasetDB.find(dsetName);
            if (dbid_iter == cc->roAttr->m_mapstrintDatasetDB.end()) {
                // Try stripping .pcl extension
                string tmpName = stripPclExtensions(dsetName);
                dbid_iter = cc->roAttr->m_mapstrintDatasetDB.find(tmpName);
                if (dbid_iter == cc->roAttr->m_mapstrintDatasetDB.end()) {
                    hasError = true;
                    anError = "DbID not found in m_mapstrintDatasetDB: " + dsetName;
                    continue;
                }
            }
            int dbID = dbid_iter->second;

            string strAvgPath =
                    cc->m_vecDBSetting[dbID]->prepDir + "/" + dsetName + ".gavg"; //avg and prep path share same directory
            string strPresencePath = cc->m_vecDBSetting[dbID]->prepDir + "/" + dsetName + ".gpres";
            string strSinfoPath = cc->m_vecDBSetting[dbID]->sinfoDir + "/" + dsetName + ".sinfo";

            if (!filesystem::exists(strAvgPath) ||  !filesystem::exists(strPresencePath) || 
                !filesystem::exists(strSinfoPath)) {
                // Try stripping .pcl extension
                string tmpName = stripPclExtensions(dsetName);
                strAvgPath = cc->m_vecDBSetting[dbID]->prepDir + "/" + tmpName + ".gavg";
                strPresencePath = cc->m_vecDBSetting[dbID]->prepDir + "/" + tmpName + ".gpres";
                strSinfoPath = cc->m_vecDBSetting[dbID]->sinfoDir + "/" + tmpName + ".sinfo";
                if (!filesystem::exists(strAvgPath) ||  !filesystem::exists(strPresencePath) || 
                    !filesystem::exists(strSinfoPath)) {
                    hasError = true;
                    anError = "Missing one of files: " + strAvgPath + ", " +
                              strPresencePath + ", " + strSinfoPath;
                    continue;
                }
            }
            vd->ReadGeneAverage(strAvgPath);
            vd->ReadGenePresence(strPresencePath);
            vd->ReadDatasetAverageStdev(strSinfoPath);
            vd->InitializeGeneMap();

            auto dp_iter = cc->roAttr->m_mapstrstrDatasetPlatform.find(dsetName);
            if (dp_iter == cc->roAttr->m_mapstrstrDatasetPlatform.end()) {
                // Try stripping .pcl extension
                string tmpName = stripPclExtensions(dsetName);
                dp_iter = cc->roAttr->m_mapstrstrDatasetPlatform.find(tmpName);
                if (dp_iter == cc->roAttr->m_mapstrstrDatasetPlatform.end()) {
                    hasError = true;
                    anError = "Dataset not found in datasetPlatformMap: " + dsetName;
                    continue;
                }
            }
            string strPlatform = dp_iter->second;

            auto p_iter = platformMap.find(strPlatform);
            if (p_iter == platformMap.end()) {
                hasError = true;
                anError = "Platform not found in PlatformMap: " + strPlatform;
                continue;
            }
            utype platform_id = p_iter->second;

            vd->SetPlatform(seekPlatforms[platform_id]);
            pl = &vd->GetPlatform();

            //NEW 3/5/2013====================================
            if (outputQueryCoexpression) {
                const vector <string> gNames = pp->GetGeneNames();
                fq_RBP = new CFullMatrix<float>();
                fq_RBP->Initialize(gNames.size(), ps);
                //fprintf(stderr, "%d %d X1\n", i, gNames.size());
                //#pragma omp parallel for \
				shared(pp, fq_RBP) private(k, j) firstprivate(gs) \
				schedule(dynamic)
                /*for(k=0; k<gNames.size(); k++){
                    int g = pp->GetGene(gNames[k]);
                    float *vv = pp->Get(g);
                    for(j=0; j<gs; j++)
                        fq_RBP->Set(g, j, vv[j]);
                    float mean = 0;
                    for(j=0; j<gs; j++)
                        mean+=fq_RBP->Get(g, j);
                    mean/=(float) (gs);
                    float stdev = 0;
                    for(j=0; j<gs; j++)
                        stdev+=(fq_RBP->Get(g, j) - mean) * (fq_RBP->Get(g, j) - mean);
                    stdev/=(float) (gs);
                    if(stdev<=0){
                        absent.insert(gNames[k]);
                    }
                    stdev = sqrt(stdev);
                    if(isnan(stdev) || isinf(stdev)){
                        fprintf(stderr, "Error:, standard deviation is zero\n");
                    }
                    for(j=0; j<gs; j++){
                        float t1 = fq_RBP->Get(g, j);
                        fq_RBP->Set(g, j, (t1 - mean) / stdev);

                    }
                }*/
                //fprintf(stderr, "%d X2\n", i);
            }
            //====================================================
        }

        if (outputQueryExpression) {
            for (k = 0; k < queryName.size(); k++) {
                int g = pp->GetGene(queryName[k]);
                if (g == -1) { //does not contain the gene in the dataset
                    for (j = 0; j < gs; j++) {
                        fq->Set(k, j, NaN);
                        d_vecQ[i].push_back(fq->Get(k, j));
                    }
                    continue;
                }
                float *vv = pp->Get(g);
                for (j = 0; j < gs; j++)
                    fq->Set(k, j, vv[j]);
                if (!outputNormalized) {
                    for (j = 0; j < gs; j++)
                        d_vecQ[i].push_back(fq->Get(k, j));
                }

                //normalize
                float mean = 0;
                for (j = 0; j < gs; j++)
                    mean += fq->Get(k, j);
                mean /= (float) (gs);
                float stdev = 0;
                for (j = 0; j < gs; j++)
                    stdev += (fq->Get(k, j) - mean) * (fq->Get(k, j) - mean);
                stdev /= (float) (gs);
                stdev = sqrt(stdev);
                for (j = 0; j < gs; j++) {
                    float t1 = fq->Get(k, j);
                    fq->Set(k, j, (t1 - mean) / stdev);
                }

                if (outputNormalized) {
                    for (j = 0; j < gs; j++)
                        d_vecQ[i].push_back(fq->Get(k, j));
                }
            }
        }

        fprintf(stderr, "allocating space %lu %d...\n", geneName.size(),
                ps);
        ff = new CFullMatrix<float>();
        ff->Initialize(genes, ps);
        fprintf(stderr, "done allocating space.\n");

        if (outputExpression) {
            for (k = 0; k < geneName.size(); k++) {
                int g = pp->GetGene(geneName[k]);
                if (g == -1) {
                    for (j = 0; j < gs; j++) {
                        ff->Set(k, j, NaN);
                        d_vecG[i].push_back(ff->Get(k, j));
                    }
                    continue;
                }
                float *vv = pp->Get(g);
                for (j = 0; j < gs; j++)
                    ff->Set(k, j, vv[j]);

                if (!outputNormalized) {
                    for (j = 0; j < gs; j++) {
                        d_vecG[i].push_back(ff->Get(k, j));
                    }
                }

                //normalize
                float mean = 0;
                for (j = 0; j < gs; j++)
                    mean += ff->Get(k, j);
                mean /= (float) (gs);
                float stdev = 0;
                for (j = 0; j < gs; j++)
                    stdev += (ff->Get(k, j) - mean) * (ff->Get(k, j) - mean);
                stdev /= (float) (gs);
                stdev = sqrt(stdev);
                for (j = 0; j < gs; j++) {
                    float t1 = ff->Get(k, j);
                    ff->Set(k, j, (t1 - mean) / stdev);
                }

                if (outputNormalized) {
                    for (j = 0; j < gs; j++) {
                        d_vecG[i].push_back(ff->Get(k, j));
                    }
                }
            }
        }

        if (outputCoexpression) {
            int kk;
            CMeasurePearNorm pn;
            for (k = 0; k < geneName.size(); k++) {
                int g = pp->GetGene(geneName[k]);
                if (g == -1) {
                    vCoexpression[k] = NaN;
                    d_vecCoexpression[i].push_back(vCoexpression[k]);
                    continue;
                }
                float avgP = 0;
                int totalQueryPresent = queryName.size();
                for (kk = 0; kk < queryName.size(); kk++) {
                    int gg = pp->GetGene(queryName[kk]);
                    if (gg == -1) {
                        totalQueryPresent--;
                        continue;
                    }

                    float *x1 = NULL;
                    float *x2 = NULL;
                    if (g < gg) {
                        x1 = pp->Get(g);
                        x2 = pp->Get(gg);
                    } else {
                        x1 = pp->Get(gg);
                        x2 = pp->Get(g);
                    }

                    float p = (float) pn.Measure(x1, ps, x2, ps,
                                                 IMeasure::EMapCenter, NULL, NULL);

                    p = (p - vd->GetDatasetAverage()) / vd->GetDatasetStdev();

                    int gID = 0; // previous default case when []operator was used for m_mapstrintGene
                    if (cc->roAttr->m_mapstrintGene.count(geneName[k]) > 0) {
                        gID = cc->roAttr->m_mapstrintGene.at(geneName[k]);
                    }
                    int qID = 0;  // previous default case when []operator was used for m_mapstrintGene
                    if (cc->roAttr->m_mapstrintGene.count(queryName[kk]) > 0) {
                        qID = cc->roAttr->m_mapstrintGene.at(queryName[kk]);
                    }

                    int qb = CMeta::Quantize(p, cc->m_quant);
                    p = cc->m_quant[qb];

                    //fprintf(stderr, "Correlation %d %d %.5f %.5f %.5f %.5f\n", qID, gID, p, vd->GetGeneAverage(gID),
                    //	pl->GetPlatformAvg(qID), pl->GetPlatformStdev(qID));

                    //subtract hubbiness
                    p = (p - vd->GetGeneAverage(gID) - pl->GetPlatformAvg(qID)) / pl->GetPlatformStdev(qID);
                    //p = p - vd->GetGeneAverage(gID);
                    p = max((float) min(p, (float) 3.2), (float) -3.2);

                    avgP += p;
                }
                if (totalQueryPresent == 0)
                    avgP = NaN;
                else
                    avgP /= (float) (totalQueryPresent);

                vCoexpression[k] = avgP;
                d_vecCoexpression[i].push_back(avgP);
            }
        }

        if (outputQueryCoexpression) {
            int kk = 0;

            const vector <string> gNames = pp->GetGeneNames();
            vector<char> qMap;
            CSeekTools::InitVector(qMap, cc->roAttr->m_vecstrGenes.size(), (char) 0);

            int totQuery = 0;
            for (k = 0; k < queryName.size(); k++) {
                int g = pp->GetGene(queryName[k]);
                if (g == -1) continue;
                int mg = 0;  // previous default case when []operator was used for m_mapstrintGene
                if (cc->roAttr->m_mapstrintGene.count(queryName[k]) > 0) {
                    mg = cc->roAttr->m_mapstrintGene.at(queryName[k]);
                }
                qMap[mg] = 1;
                totQuery++;
            }

            CMeasurePearNorm pn;

            for (k = 0; k < queryName.size(); k++) {
                int g = pp->GetGene(queryName[k]);
                if (g == -1) {
                    vqCoexpression[k] = NaN;
                    d_vecqCoexpression[i].push_back(vqCoexpression[k]);
                    continue;
                }
                vector <AResult> ar;
                ar.resize(gNames.size());

                for (kk = 0; kk < gNames.size(); kk++) {
                    int gg = pp->GetGene(gNames[kk]);
                    int id = 0; // previous default case when []operator was used for m_mapstrintGene
                    if (cc->roAttr->m_mapstrintGene.count(gNames[kk]) > 0) {
                        id = (utype) cc->roAttr->m_mapstrintGene.at(gNames[kk]);
                    }
                    ar[kk].i = (utype) id;
                    if (g == gg) {
                        ar[kk].f = 0;
                        continue;
                    }

                    float *x1 = NULL;
                    float *x2 = NULL;
                    if (g < gg) {
                        x1 = pp->Get(g);
                        x2 = pp->Get(gg);
                    } else {
                        x1 = pp->Get(gg);
                        x2 = pp->Get(g);
                    }

                    float p = (float) pn.Measure(x1, ps, x2, ps,
                                                 IMeasure::EMapCenter, NULL, NULL);

                    float px = p;
                    if (!(p < 5.0 && p > -5.0)) {
                        ar[kk].f = 0;
                        continue;
                    }

                    //get z-score (dataset wide)
                    p = (p - vd->GetDatasetAverage()) / vd->GetDatasetStdev();
                    int gID = 0; // previous default case when []operator was used for m_mapstrintGene
                    if (cc->roAttr->m_mapstrintGene.count(gNames[kk]) > 0) {
                        gID = cc->roAttr->m_mapstrintGene.at(gNames[kk]);
                    }
                    int qID = 0;  // previous default case
                    if (cc->roAttr->m_mapstrintGene.count(queryName[k]) > 0) {
                        qID = cc->roAttr->m_mapstrintGene.at(queryName[k]);
                    }

                    int qb = CMeta::Quantize(p, cc->m_quant);
                    p = cc->m_quant[qb];

                    //fprintf(stderr, "Correlation %d %d %.5f %.5f %.5f %.5f\n", qID, gID, p, vd->GetGeneAverage(gID),
                    //	pl->GetPlatformAvg(qID), pl->GetPlatformStdev(qID));

                    //subtract hubbiness
                    p = (p - vd->GetGeneAverage(gID) - pl->GetPlatformAvg(qID)) / pl->GetPlatformStdev(qID);
                    //p = p - vd->GetGeneAverage(gID);
                    p = max((float) min(p, (float) 3.2), (float) -3.2);

                    //fprintf(stderr, "%d error, infinite or nan! %.3f, %.3f, %.3f\n", i,
                    //	vd->GetDatasetAverage(), vd->GetDatasetStdev(),
                    //	vd->GetGeneAverage(mapstrintGene[gNames[kk]]));
                    //fprintf(stderr, "Correlation %d %d %.5f\n", qID, gID, p);
                    ar[kk].f = (utype)(p * 100.0 + 320);

                }

                int TOP = 1000;
                if (ar.size() < TOP) {
                    TOP = ar.size();
                }
                //fprintf(stderr, "%d H2\n", i);
                nth_element(ar.begin(), ar.begin() + TOP, ar.end());
                sort(ar.begin(), ar.begin() + TOP);
                //fprintf(stderr, "%d H3\n", i);

                float rbp = 0;
                //utype jj = 0;
                for (kk = 0; kk < TOP; kk++) {
                    if (qMap[ar[kk].i] == 0) continue;
                    if (ar[kk].f == 0) break;
                    rbp += pow(RBP_P, kk);
                    // fprintf(stderr, "Sorted %zu %d %d %.5f\n", i, kk, ar[kk].i, (ar[kk].f - 320) / 100.0f);
                    //jj++;
                    //fprintf(stderr, "Good %d %d\n", i, kk);
                    //jj++;
                }
                rbp *= (1.0 - RBP_P);

                rbp = rbp / totQuery * 1000;
                fprintf(stderr, "%zu %.3e\n", i, rbp);
                vqCoexpression[k] = rbp;
                d_vecqCoexpression[i].push_back(rbp);
            }


            /*for(k=0; k<queryName.size(); k++){
                int g = pp->GetGene(queryName[k]);
                if(g==-1){
                    vqCoexpression[k] = NaN;
                    vecqCoexpression.push_back(vqCoexpression[k]);
                    continue;
                }
                float avgP = 0;
                int totalQueryPresent = queryName.size() - 1;

                for(kk=0; kk<queryName.size(); kk++){
                    if(kk==k) continue;
                    int gg = pp->GetGene(queryName[kk]);
                    if(gg==-1){
                        totalQueryPresent--;
                        continue;
                    }

                    float *x1 = NULL;
                    float *x2 = NULL;
                    if(g<gg){
                        x1 = pp->Get(g);
                        x2 = pp->Get(gg);
                    }else{
                        x1 = pp->Get(gg);
                        x2 = pp->Get(g);
                    }

                    float p = (float) pn.Measure(x1, ps, x2, ps,
                        IMeasure::EMapCenter, NULL, NULL);
                    float px = p;

                    if(!(p<5.0 && p>-5.0)){
                        continue;
                    }

                    //get z-score (dataset wide)
                    p = (p - vd->GetDatasetAverage()) / vd->GetDatasetStdev();
                    int gID = mapstrintGene[queryName[kk]];
                    int qID = mapstrintGene[queryName[k]];

                    int qb = CMeta::Quantize(p, quant);
                    p = quant[qb];

                    //fprintf(stderr, "Correlation %d %d %.5f %.5f %.5f %.5f\n", qID, gID, p, vd->GetGeneAverage(gID),
                    //	pl->GetPlatformAvg(qID), pl->GetPlatformStdev(qID));

                    //subtract hubbiness
                    p = (p - vd->GetGeneAverage(gID) - pl->GetPlatformAvg(qID)) / pl->GetPlatformStdev(qID) ;
                    //p = p - vd->GetGeneAverage(gID);
                    p = max((float) min(p, (float) 3.2), (float)-3.2);

                    //float p = 0;
                    //for(j=2; j<gs; j++)
                    //	p+= fq->Get(k, j-2)*fq->Get(kk, j-2);
                    //p/=(float)(gs-2);
                    //p = 0.5 * log((1.0+p)/(1.0-p));
                    //p = max((float) min(p, (float) 3.2), (float)-3.2);
                    //get z-score (dataset wide)
                    //p = (p - vd->GetDatasetAverage()) / vd->GetDatasetStdev();
                    //subtract hubbiness
                    //p = p - vd->GetGeneAverage(mapstrintGene[queryName[kk]]);


                    avgP+=p;
                }
                if(totalQueryPresent==0)
                    avgP = NaN;
                else
                    avgP/=(float)(totalQueryPresent);

                vqCoexpression[k] = avgP;
                d_vecqCoexpression[i].push_back(avgP);
            }*/

        }

        if (outputCoexpression || outputQueryCoexpression) {
            delete vd;
            if (outputQueryCoexpression)
                delete fq_RBP;
        }
        if (outputQueryExpression || outputQueryCoexpression) {
            delete fq;
        }

        delete ff;
    }

    if (hasError == true) {
        throw named_error(anError);
    }

    for (i = 0; i < datasets; i++) {
        if (new_fd >= 0) {
            vecG.insert(vecG.end(), d_vecG[i].begin(), d_vecG[i].end());
            vecQ.insert(vecQ.end(), d_vecQ[i].begin(), d_vecQ[i].end());
            vecCoexpression.insert(vecCoexpression.end(),
                                d_vecCoexpression[i].begin(), d_vecCoexpression[i].end());
            vecqCoexpression.insert(vecqCoexpression.end(),
                                    d_vecqCoexpression[i].begin(), d_vecqCoexpression[i].end());
        } else {
            my->resGeneExpression->insert(my->resGeneExpression->end(), 
                                          d_vecG[i].begin(),
                                          d_vecG[i].end());
            my->resGeneCoexpression->insert(my->resGeneCoexpression->end(),
                                            d_vecCoexpression[i].begin(),
                                            d_vecCoexpression[i].end());
            my->resQueryExpression->insert(my->resQueryExpression->end(),
                                           d_vecQ[i].begin(), d_vecQ[i].end());
            my->resQueryCoexpression->insert(my->resQueryCoexpression->end(),
                                             d_vecqCoexpression[i].begin(),
                                             d_vecqCoexpression[i].end());
            if (i == 0) {
                // Copy over the dataset sizes once
                my->resDatasetSizes->insert(my->resDatasetSizes->end(),
                                            sizeD.begin(),
                                            sizeD.end());
                // *my->resDatasetSizes = sizeD;  // type conversion error
            }
        }
    }

    if (new_fd >= 0) {
        if (CSeekNetwork::Send(new_fd, sizeD) != 0) {
            fprintf(stderr, "Error sending messages\n");
        }

        if (outputExpression) {
            if (CSeekNetwork::Send(new_fd, vecG) != 0) {
                fprintf(stderr, "Error sending messages\n");
            }
        }

        if (outputQueryExpression) {
            if (CSeekNetwork::Send(new_fd, vecQ) != 0) {
                fprintf(stderr, "Error sending messages\n");
            }
        }

        if (outputCoexpression) {
            if (CSeekNetwork::Send(new_fd, vecCoexpression) != 0) {
                fprintf(stderr, "Error sending messages\n");
            }
        }

        if (outputQueryCoexpression) {
            if (CSeekNetwork::Send(new_fd, vecqCoexpression) != 0) {
                fprintf(stderr, "Error sending messages\n");
            }
        }
        close(new_fd);
    }

    return 0;
}

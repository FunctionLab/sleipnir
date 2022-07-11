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

float get_mean(vector<float> &f) {
    int i;
    float sum = 0;
    for (i = 0; i < f.size(); i++) {
        sum += f[i];
    }
    return sum / (float) (f.size());
}

float get_stdev(vector<float> &f, float mean) {
    int i;
    float sum_dev = 0;
    float diff = 0;
    for (i = 0; i < f.size(); i++) {
        diff = f[i] - mean;
        sum_dev += diff * diff;
    }
    float dev = sqrt(sum_dev / (float) f.size());
    return dev;
}


int main(int iArgs, char **aszArgs) {
    static const size_t c_iBuffer = 1024;
#ifdef WIN32
    pthread_win32_process_attach_np( );
#endif // WIN32
    gengetopt_args_info sArgs;
    ifstream ifsm;
    istream *pistm;
    vector <string> vecstrLine, vecstrGenes, vecstrDatasets, vecstrUserDatasets;
    char acBuffer[c_iBuffer];
    size_t i;

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }

    vector <string> vecstrGeneID;
    map <string, utype> mapstrintGene;
    if (!CSeekTools::ReadListTwoColumns(sArgs.input_arg, vecstrGeneID, vecstrGenes))
        return -1;

    for (i = 0; i < vecstrGenes.size(); i++)
        mapstrintGene[vecstrGenes[i]] = i;

    fprintf(stderr, "Finished reading gene map\n");

    if (sArgs.read_bin_vec_flag == 1) {
        vector<float> vv;
        CSeekTools::ReadArray(sArgs.input_vec_arg, vv);
        for (i = 0; i < vv.size(); i++) {
            fprintf(stderr, "%.5f\n", vv[i]);
        }
        return 0;
    }

    if (sArgs.add_gscore_flag == 1) {
        vector <string> vecstrGScore;
        if (!CSeekTools::ReadListOneColumn(sArgs.gscore_list_arg, vecstrGScore))
            return -1;
        string gscoreDir = sArgs.gscore_dir_arg;
        vector<float> totalScore;
        totalScore.resize(vecstrGeneID.size());
        for (int j = 0; j < totalScore.size(); j++) {
            totalScore[j] = 0;
        }
        for (int j = 0; j < vecstrGScore.size(); j++) {
            string path = gscoreDir + "/" + vecstrGScore[j];
            vector<float> v1;
            CSeekTools::ReadArray(path.c_str(), v1);
            for (int k = 0; k < v1.size(); k++) {
                if (v1[k] < -320) {
                    totalScore[k] = -320;
                } else {
                    totalScore[k] += v1[k];
                }
            }
        }
        for (int j = 0; j < totalScore.size(); j++) {
            if (totalScore[j] == -320) {
                continue;
            }
            totalScore[j] /= (float) totalScore.size();
        }
        CSeekTools::WriteArray(sArgs.gscore_output2_arg, totalScore);
        return 0;
    }


    if (sArgs.increase_gscore_flag == 1) {
        vector<float> v1;
        CSeekTools::ReadArray(sArgs.gscore_file_arg, v1);

        vector<float> v2;
        CSeekTools::ReadArray(sArgs.gscore_file_2_arg, v2);

        int num_genes_1 = 0;
        int num_genes_2 = 0;
        for (int j = 0; j < v1.size(); j++) {
            if (v1[j] < -320) continue;
            num_genes_1++;
        }
        for (int j = 0; j < v2.size(); j++) {
            if (v2[j] < -320) continue;
            num_genes_2++;
        }

        vector <CPair<float>> cp1, cp2;
        cp1.resize(v1.size());
        cp2.resize(v2.size());
        for (int j = 0; j < v1.size(); j++) {
            cp1[j].i = (utype) j;
            cp1[j].v = v1[j];
            cp2[j].i = (utype) j;
            cp2[j].v = v2[j];
        }
        sort(cp1.begin(), cp1.end(), CDescendingValue<float>());
        sort(cp2.begin(), cp2.end(), CDescendingValue<float>());

        for (int j = 0; j < v1.size(); j++) {
            if (cp1[j].v < -320) {
            } else {
                v1[cp1[j].i] = (float) (num_genes_1 - j) / num_genes_1;
            }
            if (cp2[j].v < -320) {
            } else {
                v2[cp2[j].i] = (float) (num_genes_2 - j) / num_genes_2;
            }
        }

        vector<float> sum_v;
        sum_v.resize(v1.size());
        for (int j = 0; j < v1.size(); j++) {
            if (v1[j] < -320 || v2[j] < -320) {
                sum_v[j] = -320;
            } else {
                sum_v[j] = v1[j] + v2[j];
            }
        }

        CSeekTools::WriteArray(sArgs.gscore_output_arg, sum_v);
        return 0;
    }

    if (sArgs.weight2_flag == 1) {
        vector <string> vecstrDataset;
        if (!CSeekTools::ReadListOneColumn(sArgs.dweight_map_arg, vecstrDataset))
            return -1;
        vector <vector<float>> vec_score, orig_score;
        utype i, j;
        int num_query = sArgs.dweight_num_arg; //random query
        orig_score.resize(num_query);

        vec_score.resize(vecstrDataset.size());
        for (j = 0; j < vec_score.size(); j++)
            vec_score[j].resize(num_query);

        char x[256];
        for (i = 0; i < num_query; i++) { //i is query id
            vector<float> v;
            sprintf(x, "%s/%d.dweight", sArgs.dweight_dir_arg, i);
            CSeekTools::ReadArray(x, v);
            orig_score[i] = v;
            for (j = 0; j < vec_score.size(); j++) //j is dataset id
                vec_score[j][i] = v[j];
        }

        for (j = 0; j < vec_score.size(); j++)
            sort(vec_score[j].begin(), vec_score[j].end());

        int test_num_query = sArgs.dweight_test_num_arg;
        for (i = 0; i < test_num_query; i++) { //i is query id
            fprintf(stderr, "Query %d\n", i);
            vector<float> v;
            sprintf(x, "%s/%d.dweight", sArgs.dweight_test_dir_arg, i);
            CSeekTools::ReadArray(x, v);
            for (j = 0; j < v.size(); j++) {
                utype k = 0;
                for (k = 0; k < num_query; k++) {
                    if (v[j] < vec_score[j][k])
                        break;
                }
                fprintf(stderr, "%s\t%.3e\t%d\n", vecstrDataset[j].c_str(), v[j], k);
            }
        }
        fprintf(stderr, "Done!\n");
        return 0;
    }

    if (sArgs.combine_pcl_flag == 1) {
        vector <string> vecstrPCL;
        if (!CSeekTools::ReadListOneColumn(sArgs.pcl_list_arg, vecstrPCL))
            return -1;
        vector < CPCL * > vc;
        vc.resize(vecstrPCL.size());
        utype i, j, k;

        float **new_m;
        int totExperiments = 0;
        int totDatasets = vecstrPCL.size();
        vector<int> geneFreq;
        CSeekTools::InitVector(geneFreq, vecstrGenes.size(), (int) 0);

        for (i = 0; i < vecstrPCL.size(); i++) {
            fprintf(stderr, "Reading %s...\n", vecstrPCL[i].c_str());
            vc[i] = new CPCL();
            CPCL *cp = vc[i];
            cp->Open(vecstrPCL[i].c_str());
            totExperiments += cp->GetExperiments();
            for (j = 0; j < vecstrGenes.size(); j++) {
                int g = cp->GetGene(vecstrGenes[j]);
                if (g == -1) continue;
                geneFreq[j]++;
            }
        }

        vector <string> pr; //all present genes
        int totGenes = 0;
        for (j = 0; j < vecstrGenes.size(); j++) {
            if (geneFreq[j] == totDatasets) {
                totGenes++;
                pr.push_back(vecstrGenes[j]);
            }
        }

        //option 1 (Normalize within dataset)
        vector < vector < vector < float > > > mm;
        vector <vector<int>> good_experiment_id;

        mm.resize(vc.size());
        good_experiment_id.resize(vc.size());
        for (i = 0; i < vecstrPCL.size(); i++) {
            CPCL *cp = vc[i];
            mm[i].resize(totGenes);

            float **new_m = new float *[totGenes];
            for (j = 0; j < totGenes; j++)
                new_m[j] = new float[cp->GetExperiments()];

            vector<float> exp_avg;
            vector<float> exp_stdev;
            CSeekTools::InitVector(exp_avg, cp->GetExperiments(), (float) 0);
            CSeekTools::InitVector(exp_stdev, cp->GetExperiments(), (float) 0);

            for (j = 0; j < pr.size(); j++) {
                int g = cp->GetGene(pr[j]);
                float *vv = cp->Get(g);
                for (k = 0; k < cp->GetExperiments(); k++) {
                    new_m[j][k] = vv[k];
                }
            }

            /*
            set<int> badExperiments;
            for(k=0; k<cp->GetExperiments(); k++){
                vector<float> vf;
                for(j=0; j<pr.size(); j++){
                    vf.push_back(new_m[j][k]);
                }
                exp_avg[k] = get_mean(vf);
                exp_stdev[k] = get_stdev(vf, exp_avg[k]);
                //float fl = (new_m[j][k] - exp_avg[k]) / exp_stdev[k];
                //if(isinf(exp_avg[k]) || isnan(exp_avg[k])){
                //	badExperiments.insert(k);
                //	continue;
                //}
                //if(isinf(exp_stdev[k]) || isnan(exp_stdev[k])){
                //	badExperiments.insert(k);
                //	continue;
                //}
            }*/


            good_experiment_id[i] = vector<int>();
            for (k = 0; k < cp->GetExperiments(); k++) {
                //if(badExperiments.find(k)==badExperiments.end()){
                good_experiment_id[i].push_back(k);
                //}
            }

            for (j = 0; j < totGenes; j++) {
                vector<float> gv;
                for (k = 0; k < cp->GetExperiments(); k++) {
                    //if(badExperiments.find(k)==badExperiments.end()){
                    gv.push_back(new_m[j][k]);
                    //}
                }
                float mean = get_mean(gv);
                float stdev = get_stdev(gv, mean);
                mm[i][j] = vector<float>();

                for (k = 0; k < cp->GetExperiments(); k++) {
                    //if(badExperiments.find(k)==badExperiments.end()){
                    //z-scoring
                    //float fl = new_m[j][k];
                    float fl = (new_m[j][k] - mean) / stdev; //row normalization
                    if (isnan(fl) || isinf(fl)) {
                        fprintf(stderr, "A %.2f\t%.2f\t%.2f\t%.2f\tD%d\tG%d\tE%d/%zu\n", fl, mean,
                                stdev, new_m[j][k], i, j, k, cp->GetExperiments());
                    }
                    mm[i][j].push_back(fl);
                    //original expression values
                    //mm[i][j].push_back(new_m[j][k]);
                    //}
                }
            }

            //Column normalization
            for (k = 0; k < cp->GetExperiments(); k++) {
                vector<float> cv;
                //do not use badExperiments
                for (j = 0; j < totGenes; j++) {
                    cv.push_back(mm[i][j][k]);
                }
                float mean = get_mean(cv);
                float stdev = get_stdev(cv, mean);
                //column normalization
                for (j = 0; j < totGenes; j++) {
                    float fl = (mm[i][j][k] - mean) / stdev;
                    if (isnan(fl) || isinf(fl)) {
                        fprintf(stderr, "B %.2f\t%.2f\t%.2f\t%.2f\tD%d\tG%d\tE%d/%zu\n", fl, mean,
                                stdev, mm[i][j][k], i, j, k, cp->GetExperiments());
                    }
                    mm[i][j][k] = fl;
                }
            }

            for (j = 0; j < totGenes; j++)
                delete[] new_m[j];
            delete[] new_m;
        }

        FILE *pFile;
        pFile = fopen(sArgs.output_pcl_arg, "w");

        for (j = 0; j < totGenes; j++) {
            vector<float> vv;
            for (i = 0; i < vecstrPCL.size(); i++) {
                for (k = 0; k < mm[i][j].size(); k++) {
                    vv.push_back(mm[i][j][k]);
                }
            }
            for (i = 0; i < vv.size(); i++) {
                //int vi = 0;
                //binarization
                //if(vv[i]>=1.0 || vv[i] <= -1.0)
                //	vi = 1;

                //show integer z-scores
                //int vi = (int) vv[i];
                //fprintf(pFile, "%d", vi);

                //show original values (float)
                float vi = vv[i];
                fprintf(pFile, "%.2f", vi);

                if (i == vv.size() - 1) {
                    fprintf(pFile, "\n");
                } else {
                    fprintf(pFile, "\t");
                }
            }
        }
        fclose(pFile);

        string mms = sArgs.output_pcl_arg;
        mms += "_map";
        pFile = fopen(mms.c_str(), "w");
        for (i = 0; i < vecstrPCL.size(); i++) {
            for (k = 0; k < good_experiment_id[i].size(); k++) {
                fprintf(pFile, "Dataset %d\tExperiment %d\n", i, good_experiment_id[i][k]);
            }
        }
        for (j = 0; j < totGenes; j++) {
            fprintf(pFile, "Gene %s\n", pr[j].c_str());
        }
        fclose(pFile);

        //fprintf(stderr, "Number of datasets: %d. Number of genes with full coverage: %d.\n",
        //totDatasets, totGenes);

        /*
        new_m = new float*[totGenes];
        for(i=0; i<totGenes; i++){
            new_m[i] = new float[totExperiments];
        }

        int kk = 0;
        for(i=0; i<vecstrPCL.size(); i++){
            CPCL *cp = vc[i];
            for(j=0; j<pr.size(); j++){
                int g = cp->GetGene(pr[j]);
                float *vv = cp->Get(g);
                for(k=0; k<cp->GetExperiments(); k++){
                    new_m[j][kk+k] = vv[k];
                }
            }
            kk+=cp->GetExperiments();
        }

        vector<float> allV;
        allV.resize(totGenes);
        set<int> badExperiments;
        vector<float> averages;
        averages.resize(totExperiments);
        float mean_of_mean = 0;

        for(i=0; i<totExperiments; i++){
            for(j=0; j<totGenes; j++){
                allV[j] = new_m[j][i];
            }
            float mean = get_mean(allV);
            float stdev = get_stdev(allV, mean);
            if(isinf(mean) || isnan(mean)){
                badExperiments.insert(i);
                continue;
            }
            averages[i] = mean;
            mean_of_mean += mean;
            //fprintf(stderr, "%d\t%.2f\t%.2f\n", i, mean, stdev);
        }

        mean_of_mean /= (float) (totExperiments - badExperiments.size());
        //fprintf(stderr, "Mean of mean %.2f\n", mean_of_mean);

        //adjustment by mean
        vector<int> goodExperiments;
        for(i=0; i<totExperiments; i++){
            if(badExperiments.find(i)==badExperiments.end()){
                float adj = averages[i] - mean_of_mean;
                for(j=0; j<totGenes; j++){
                    new_m[j][i] -= adj;
                }
                goodExperiments.push_back(i);
            }
        }
        fprintf(stderr, "Number of experiments: %d\n", goodExperiments.size());
        */

        //Option 2 (normalize across all datasets)
        /*
        vector<float> mean_gene;
        vector<float> stdev_gene;
        CSeekTools::InitVector(mean_gene, pr.size(), (float) 0); //present genes
        CSeekTools::InitVector(stdev_gene, pr.size(), (float) 0); //present genes
        //calculate mean
        for(i=0; i<goodExperiments.size(); i++){
            int ii = goodExperiments[i];
            for(j=0; j<totGenes; j++){
                mean_gene[j] += new_m[j][ii];
            }
        }
        //calculate mean and standard deviation
        for(j=0; j<totGenes; j++)
            mean_gene[j] /= (float) goodExperiments.size();
        for(i=0; i<goodExperiments.size(); i++){
            int ii = goodExperiments[i];
            for(j=0; j<totGenes; j++){
                float diff = new_m[j][ii] - mean_gene[j];
                stdev_gene[j] += diff * diff;
            }
        }
        for(j=0; j<totGenes; j++)
            stdev_gene[j] = sqrt(stdev_gene[j] / (float) goodExperiments.size());
        for(j=0; j<totGenes; j++){
            for(i=0; i<goodExperiments.size(); i++){
                float v = (new_m[j][goodExperiments[i]] - mean_gene[j]) / stdev_gene[j];
                int vi = 0;
                if(v>=1.0)
                    vi = 1;
                fprintf(stdout, "%d", vi);
                if(i==goodExperiments.size()-1){
                    fprintf(stdout, "\n");
                }else{
                    fprintf(stdout, "\t");
                }
            }
        }
        */

        //original values
        /*
        for(j=0; j<totGenes; j++){
            for(i=0; i<goodExperiments.size(); i++){
                fprintf(stdout, "%d", (int) new_m[j][goodExperiments[i]]);
                if(i==goodExperiments.size()-1){
                    fprintf(stdout, "\n");
                }else{
                    fprintf(stdout, "\t");
                }
            }
        }

        for(i=0; i<totGenes; i++){
            delete new_m[i];
        }
        delete new_m;
        */
        for (i = 0; i < vecstrPCL.size(); i++)
            delete vc[i];
        vc.clear();

        return 0;
    }

    if (sArgs.convert_aracne_flag == 1) {
        int lineLen = 1024;
        char *acBuffer = (char *) malloc(lineLen);
        FILE *infile;
        if ((infile = fopen(sArgs.aracne_file_arg, "r")) == NULL) {
            fprintf(stderr, "Error no such file %s\n", sArgs.aracne_file_arg);
            return -1;
        }
        while (fgets(acBuffer, lineLen, infile) != NULL) {
            while (strlen(acBuffer) == lineLen - 1) {
                int len = strlen(acBuffer);
                fseek(infile, -len, SEEK_CUR);
                lineLen += 1024;
                acBuffer = (char *) realloc(acBuffer, lineLen);
                char *ret = fgets(acBuffer, lineLen, infile);
            }
        }
        fclose(infile);
        fprintf(stderr, "Max line length: %d\n", lineLen);
        free(acBuffer);

        acBuffer = (char *) malloc(lineLen);
        ifstream ia;
        ia.open(sArgs.aracne_file_arg);
        if (!ia.is_open()) {
            fprintf(stderr, "Error opening file %s\n", sArgs.aracne_file_arg);
            return -1;
        }
        set <string> allGenes;
        size_t ci = 0;
        while (!ia.eof()) {
            ia.getline(acBuffer, lineLen - 1);
            if (acBuffer[0] == 0) break;
            acBuffer[lineLen - 1] = 0;
            if (acBuffer[0] == '>') continue;
            vector <string> tok;
            CMeta::Tokenize(acBuffer, tok);
            allGenes.insert(tok[0]);
            size_t si;
            for (si = 1; si < tok.size(); si += 2) {
                allGenes.insert(tok[si]);
            }
            if (ci % 100 == 0) {
                fprintf(stderr, "Gene Number %zu\n", ci);
            }
            ci++;
        }
        ia.close();

        vector <string> vecGenes;
        copy(allGenes.begin(), allGenes.end(), back_inserter(vecGenes));
        map <string, size_t> mapiGenes;
        for (i = 0; i < vecGenes.size(); i++)
            mapiGenes[vecGenes[i]] = i;

        CDat Dat;
        Dat.Open(vecGenes);
        ci = 0;
        ia.open(sArgs.aracne_file_arg);
        while (!ia.eof()) {
            ia.getline(acBuffer, lineLen - 1);
            if (acBuffer[0] == 0) break;
            acBuffer[lineLen - 1] = 0;
            if (acBuffer[0] == '>') continue;
            vector <string> tok;
            CMeta::Tokenize(acBuffer, tok);
            size_t tG = mapiGenes[tok[0]];
            size_t si;
            for (si = 1; si < tok.size(); si += 2) {
                size_t oG = mapiGenes[tok[si]];
                float fG = atof(tok[si + 1].c_str());
                Dat.Set(tG, oG, fG);
            }
            if (ci % 100 == 0) {
                fprintf(stderr, "Gene Number %zu\n", ci);
            }
            ci++;
        }
        ia.close();
        free(acBuffer);
        Dat.Save(sArgs.output_dab_file_arg);
        return 0;
    }

    if (sArgs.weight_flag == 1) {
        vector <string> vecstrDataset;
        if (!CSeekTools::ReadListOneColumn(sArgs.dweight_map_arg, vecstrDataset))
            return -1;

        vector <vector<float>> vec_score;
        vector <vector<float>> orig_score;
        utype i, j;
        int num_query = sArgs.dweight_num_arg; //random query
        orig_score.resize(num_query);

        vec_score.resize(vecstrDataset.size());
        for (j = 0; j < vec_score.size(); j++)
            vec_score[j].resize(num_query);

        char x[256];
        for (i = 0; i < num_query; i++) { //i is query id
            vector<float> v;
            sprintf(x, "%s/%d.dweight", sArgs.dweight_dir_arg, i);
            CSeekTools::ReadArray(x, v);
            orig_score[i] = v;
            for (j = 0; j < vec_score.size(); j++) //j is dataset id
                vec_score[j][i] = v[j];
        }

        vector<float> score_cutoff;
        score_cutoff.resize(vec_score.size());
        for (j = 0; j < vec_score.size(); j++) {
            sort(vec_score[j].begin(), vec_score[j].end());
            score_cutoff[j] = vec_score[j][899];
            //fprintf(stderr, "Dataset %d: %.4e\n", j, vec_score[j][899]);
        }

        sprintf(x, "/tmp/dataset.cutoff");
        CSeekTools::WriteArray(x, score_cutoff);

        vector<int> numGoodDataset;
        CSeekTools::InitVector(numGoodDataset, num_query, (int) 0);
        for (i = 0; i < num_query; i++) //i is query id
            for (j = 0; j < vecstrDataset.size(); j++)
                if (orig_score[i][j] > vec_score[j][899])
                    numGoodDataset[i]++;

        sort(numGoodDataset.begin(), numGoodDataset.end());
        fprintf(stderr, "10 percentile %d\n", numGoodDataset[99]);
        fprintf(stderr, "90 percentile %d\n", numGoodDataset[899]);

        int test_num_query = sArgs.dweight_test_num_arg;
        for (i = 0; i < test_num_query; i++) { //i is query id
            vector<float> v;
            sprintf(x, "%s/%d.dweight", sArgs.dweight_test_dir_arg, i);
            CSeekTools::ReadArray(x, v);
            int numGood = 0;
            for (j = 0; j < v.size(); j++)
                if (v[j] > vec_score[j][899])
                    numGood++;
            fprintf(stderr, "Query %d: ", i);
            if (numGood > numGoodDataset[899])
                fprintf(stderr, "Unique (upper)\n");
            else if (numGood < numGoodDataset[99])
                fprintf(stderr, "Unique (lower)\n");
            else
                fprintf(stderr, "Not unique\n");
        }
        fprintf(stderr, "Done!\n");
        return 0;
    }


    if (sArgs.convert_dab_flag == 1) {
        string dab_file_name = sArgs.dab_file_arg;
        if (dab_file_name.find(".2.dab") != string::npos) {
            CSeekIntIntMap d1(vecstrGenes.size());
            CSparseFlatMatrix<float> sm(0);
            CSeekWriter::ReadSeekSparseMatrix<unsigned short>(sArgs.dab_file_arg,
                                                              sm, d1, 1000, 0.99, vecstrGenes);
            utype ii = 0;
            const vector <utype> &allRGenes = d1.GetAllReverse();
            FILE *pFile;
            pFile = fopen(sArgs.output_matrix_arg, "w");
            fprintf(pFile, "Genes\t");
            for (i = 0; i < d1.GetNumSet(); i++) {
                string g = vecstrGenes[allRGenes[i]];
                fprintf(pFile, "%s", g.c_str());
                if (i < d1.GetNumSet() - 1)
                    fprintf(pFile, "\t");
            }
            fprintf(pFile, "\n");

            unsigned int j, t;
            for (i = 0; i < d1.GetNumSet(); i++) {
                string g = vecstrGenes[allRGenes[i]];
                fprintf(pFile, "%s\t", g.c_str());
                vector<float> vf;
                vf.resize(vecstrGenes.size());
                for (j = 0; j < vecstrGenes.size(); j++)
                    vf[j] = 0;
                utype ii = allRGenes[i];
                vector < CPair < float > > ::iterator
                row_it;
                float rv;
                for (row_it = sm.RowBegin(ii); row_it != sm.RowEnd(ii); row_it++) {
                    j = (size_t) row_it->i;
                    rv = row_it->v;
                    vf[j] = rv;
                }

                for (j = 0; j < d1.GetNumSet(); j++) {
                    t = allRGenes[j];
                    fprintf(pFile, "%.3e", vf[t]);
                    if (j < d1.GetNumSet() - 1)
                        fprintf(pFile, "\t");
                }
                fprintf(pFile, "\n");
            }
            fclose(pFile);
            return 0;
        }

        CDataPair Dat;
        fprintf(stderr, "Opening file...\n");
        if (!Dat.Open(sArgs.dab_file_arg, false, false, 2, false, false)) {
            cerr << "error opening file" << endl;
            return -1;
        }

        vector<unsigned int> veciGenes;
        veciGenes.resize(vecstrGenes.size());
        for (i = 0; i < vecstrGenes.size(); i++)
            veciGenes[i] = (unsigned int) Dat.GetGeneIndex(vecstrGenes[i]);

        unsigned int s, t, j, ss, tt;
        float d;
        CSeekIntIntMap m1(vecstrGenes.size());
        for (i = 0; i < vecstrGenes.size(); i++) {
            if ((s = veciGenes[i]) == (unsigned int) -1) continue;
            m1.Add(i);
        }

        set <utype> badIndex;
        const vector <utype> &allRGenes1 = m1.GetAllReverse();
        for (i = 0; i < m1.GetNumSet(); i++) {
            s = veciGenes[allRGenes1[i]];
            int countNaN = 0;
            for (j = 0; j < m1.GetNumSet(); j++) {
                t = veciGenes[allRGenes1[j]];
                if (s == t) continue;
                if (CMeta::IsNaN(d = Dat.Get(s, t))) {
                    countNaN++;
                }
            }
            if (countNaN == m1.GetNumSet() - 1) {
                badIndex.insert(allRGenes1[i]);
            }
        }

        //fprintf(stderr, "Bad indices: %d\n", badIndex.size());
        //getchar();

        CSeekIntIntMap m(vecstrGenes.size());
        for (i = 0; i < m1.GetNumSet(); i++) {
            utype s = allRGenes1[i];
            if (badIndex.find(s) == badIndex.end()) {
                m.Add(s);
            } else {
                fprintf(stderr, "Gene %s is deleted\n", vecstrGenes[s].c_str());
            }
        }

        const vector <utype> &allRGenes = m.GetAllReverse();
        FILE *pFile;
        pFile = fopen(sArgs.output_matrix_arg, "w");
        fprintf(pFile, "Genes\t");
        for (i = 0; i < m.GetNumSet(); i++) {
            string g = vecstrGenes[allRGenes[i]];
            fprintf(pFile, "%s", g.c_str());
            if (i < m.GetNumSet() - 1) {
                fprintf(pFile, "\t");
            }
        }
        fprintf(pFile, "\n");

        for (i = 0; i < m.GetNumSet(); i++) {
            s = veciGenes[allRGenes[i]];
            string g = vecstrGenes[allRGenes[i]];
            fprintf(pFile, "%s\t", g.c_str());
            for (j = 0; j < m.GetNumSet(); j++) {
                t = veciGenes[allRGenes[j]];
                if (CMeta::IsNaN(d = Dat.Get(s, t)))
                    d = 0;
                fprintf(pFile, "%.3e", d);
                if (j < m.GetNumSet() - 1) {
                    fprintf(pFile, "\t");
                }
            }
            fprintf(pFile, "\n");
        }
        fclose(pFile);
        return 0;
    }


    if (sArgs.limit_hub_flag == 1) {
        float per = 0.30;
        CDataPair Dat;
        char outFile[1024];
        fprintf(stderr, "Opening file...\n");
        if (!Dat.Open(sArgs.dabinput_arg, false, false, 2, false, false)) {
            cerr << "error opening file" << endl;
            return -1;
        }
        vector<unsigned int> veciGenes;
        veciGenes.resize(vecstrGenes.size());
        for (i = 0; i < vecstrGenes.size(); i++)
            veciGenes[i] = (unsigned int) Dat.GetGeneIndex(vecstrGenes[i]);

        unsigned int s, t, j, ss, tt;
        float d;
        CSeekIntIntMap m(vecstrGenes.size());
        for (i = 0; i < vecstrGenes.size(); i++) {
            if ((s = veciGenes[i]) == (unsigned int) -1) continue;
            m.Add(i);
        }
        vector<float> vecSum;
        CSeekTools::InitVector(vecSum, vecstrGenes.size(), (float) 0);
        const vector <utype> &allRGenes = m.GetAllReverse();
        for (i = 0; i < m.GetNumSet(); i++) {
            s = veciGenes[allRGenes[i]];
            for (j = i + 1; j < m.GetNumSet(); j++) {
                t = veciGenes[allRGenes[j]];
                if (CMeta::IsNaN(d = Dat.Get(s, t))) continue;
                vecSum[allRGenes[i]] += pow(abs(d), 9);
                vecSum[allRGenes[j]] += pow(abs(d), 9);
            }
        }
        vector<float> backupSum;
        for (i = 0; i < vecSum.size(); i++)
            backupSum.push_back(vecSum[i]);
        sort(vecSum.begin(), vecSum.end(), greater<float>());
        int index = (int) (per * (float) m.GetNumSet());
        float percentile = vecSum[index];
        vector <string> limitedGenes;
        for (i = 0; i < backupSum.size(); i++) {
            if (backupSum[i] > percentile) {
                limitedGenes.push_back(vecstrGenes[i]);
            }
        }
        fprintf(stderr, "%zu / %d genes to be written!\n", limitedGenes.size(), m.GetNumSet());
        CDat NewDat;
        NewDat.Open(limitedGenes);
        for (i = 0; i < limitedGenes.size(); i++) {
            s = (unsigned int) Dat.GetGeneIndex(limitedGenes[i]);
            ss = NewDat.GetGeneIndex(limitedGenes[i]);
            for (j = i + 1; j < limitedGenes.size(); j++) {
                t = (unsigned int) Dat.GetGeneIndex(limitedGenes[j]);
                tt = NewDat.GetGeneIndex(limitedGenes[j]);
                d = Dat.Get(s, t);
                NewDat.Set(ss, tt, d);
            }
        }
        NewDat.Save(sArgs.hub_dab_output_arg);
        return 0;
    }

    if (sArgs.comp_ranking_flag == 1) {
        utype i, j;
        int num_query = sArgs.gscore_num1_arg; //random query
        char x[256];
        for (i = 0; i < num_query; i++) { //i is query id
            vector<float> v1, v2;
            sprintf(x, "%s/%d.gscore", sArgs.gscore_dir1_arg, i);
            CSeekTools::ReadArray(x, v1);
            sprintf(x, "%s/%d.gscore", sArgs.gscore_dir2_arg, i);
            CSeekTools::ReadArray(x, v2);
            vector <CPair<float>> cp1, cp2;
            cp1.resize(v1.size());
            cp2.resize(v2.size());
            for (j = 0; j < v1.size(); j++) {
                cp1[j].i = (utype) j;
                cp1[j].v = v1[j];
                cp2[j].i = (utype) j;
                cp2[j].v = v2[j];
            }
            sort(cp1.begin(), cp1.end(), CDescendingValue<float>());
            sort(cp2.begin(), cp2.end(), CDescendingValue<float>());
            vector<char> presence;
            CSeekTools::InitVector(presence, v1.size(), (char) 0);
            for (j = 0; j < 500; j++) {
                presence[cp1[j].i]++;
                presence[cp2[j].i]++;
            }
            int count = 0;
            for (j = 0; j < v1.size(); j++) {
                if (presence[j] == 2) {
                    count++;
                }
            }
            fprintf(stderr, "Query %d %d\n", i, count);
        }
        return 0;
    }

    if (sArgs.dataset_flag == 1) {
        string db = sArgs.db_arg;
        string dset_list = sArgs.dset_list_arg;
        string dir_in = sArgs.dir_in_arg;
        string dir_prep = sArgs.dir_prep_in_arg;
        if (db == "NA" || dset_list == "NA" || dir_in == "NA" ||
            dir_prep == "NA") {
            fprintf(stderr, "Requires: -x, -X, -d -p\n");
            return -1;
        }

        vector <string> vecstrDP, vecstrUserDP;
        //dataset-platform mapping (required)
        if (!CSeekTools::ReadListTwoColumns(sArgs.db_arg, vecstrDatasets, vecstrDP))
            return -1;

        //dataset filter
        if (!CSeekTools::ReadListTwoColumns(sArgs.dset_list_arg, vecstrUserDatasets, vecstrUserDP))
            return -1;

        map <string, utype> mapstrintDataset;
        map <string, string> mapstrstrDatasetPlatform;
        utype i;
        for (i = 0; i < vecstrDatasets.size(); i++) {
            mapstrstrDatasetPlatform[vecstrDatasets[i]] = vecstrDP[i];
            mapstrintDataset[vecstrDatasets[i]] = i;
        }

        CSeekIntIntMap *mapUserDatasets = new CSeekIntIntMap(vecstrDatasets.size());
        for (i = 0; i < vecstrUserDatasets.size(); i++) {
            mapUserDatasets->Add(mapstrintDataset[vecstrUserDatasets[i]]);
        }

        //fprintf(stderr, "Finished reading dataset\n");

        size_t iDatasets = vecstrDatasets.size();

        SeekPlatforms platforms;
        platforms.loadPlatformDataFromFiles(sArgs.platform_dir_arg);
        const vector<CSeekPlatform> &vp = platforms.getCSeekPlatforms();
        const map<string, utype> &mapstriPlatform = platforms.getPlatformMap();
        //fprintf(stderr, "Finished reading platform\n");

        vector < CSeekDataset * > vc;
        vc.resize(iDatasets);
        string strPrepInputDirectory = sArgs.dir_prep_in_arg;
        string strGvarInputDirectory = sArgs.dir_gvar_in_arg;
        string strSinfoInputDirectory = sArgs.dir_sinfo_in_arg;

        for (i = 0; i < iDatasets; i++) {
            vc[i] = new CSeekDataset();
            string strFileStem = vecstrDatasets[i];
            string strAvgPath = strPrepInputDirectory + "/" +
                                strFileStem + ".gavg";
            string strPresencePath = strPrepInputDirectory + "/" +
                                     strFileStem + ".gpres";

            if (strGvarInputDirectory != "NA") {
                string strGvarPath = strGvarInputDirectory + "/" +
                                     strFileStem + ".gexpvar";
                vc[i]->ReadGeneVariance(strGvarPath);
            }
            if (strSinfoInputDirectory != "NA") {
                string strSinfoPath = strSinfoInputDirectory + "/" +
                                      strFileStem + ".sinfo";
                vc[i]->ReadDatasetAverageStdev(strSinfoPath);
            }

            vc[i]->ReadGeneAverage(strAvgPath);
            vc[i]->ReadGenePresence(strPresencePath);
            string strPlatform =
                    mapstrstrDatasetPlatform.find(strFileStem)->second;
            utype platform_id = mapstriPlatform.find(strPlatform)->second;
            vc[i]->SetPlatform(vp[platform_id]);
        }

        //fprintf(stderr, "Finished reading prep\n");

        for (i = 0; i < iDatasets; i++) vc[i]->InitializeGeneMap();

        for (i = 0; i < vecstrGenes.size(); i++) {
            utype ii = mapstrintGene[vecstrGenes[i]];
            fprintf(stderr, "Gene %s ", vecstrGenes[i].c_str());
            utype j = 0;
            vector<float> va;
            for (j = 0; j < iDatasets; j++) {
                CSeekIntIntMap *gm = vc[j]->GetGeneMap();
                utype ij = gm->GetForward(ii);
                if (CSeekTools::IsNaN(ij)) {
                    continue;
                }
                float a = vc[j]->GetGeneAverage(ii);
                va.push_back(a);
                //fprintf(stderr, "%.2f ", a);
            }
            sort(va.begin(), va.end());
            utype g = 0;
            for (g = 0; g < 20; g++) {
                utype ik = (utype)((float) 0.05 * (float) (g + 1) * (float) va.size() - 1.0);
                fprintf(stderr, "%.2f ", va[ik]);
            }
            fprintf(stderr, "\n");
        }

        fprintf(stderr, "Done\n");
        return 1;

        vector <vector<string>> vecstrQueries;
        string multiQuery = sArgs.multi_query_arg;

        if (multiQuery == "NA") {
            vecstrQueries.resize(1);
            if (!CSeekTools::ReadMultiGeneOneLine(sArgs.single_query_arg, vecstrQueries[0]))
                return -1;
        } else {
            if (!CSeekTools::ReadMultipleQueries(multiQuery, vecstrQueries))
                return -1;
        }

        bool toOutput = false;
        string output_file = sArgs.output_file_arg;
        FILE *out = NULL;

        if (output_file == "NA") {
            toOutput = false;
        } else {
            toOutput = true;
            out = fopen(output_file.c_str(), "w");
        }

        utype ii = 0;
        float cutoff = sArgs.gvar_cutoff_arg;
        bool toCutoff = false;
        if (cutoff < 0) {
            toCutoff = false;
        } else {
            toCutoff = true;
        }

        if (sArgs.order_stat_single_gene_query_flag == 1) {
            if (strGvarInputDirectory == "NA") {
                fprintf(stderr, "Order statistics mode, but need to provide gvar!\n");
                return -1;
            }
            if (cutoff < 0) {
                fprintf(stderr, "Need to provide positive Gvar cutoff\n");
                return -1;
            }
        }

        for (ii = 0; ii < vecstrQueries.size(); ii++) {

            vector <string> &vecstrQuery = vecstrQueries[ii];

            //fprintf(stderr, "Finished reading query\n");
            bool isFirst = true;
            vector<int> count;
            utype j;
            CSeekTools::InitVector(count, vecstrQuery.size(), (int) 0);

            //only analyze and report on user datasets
            const vector <utype> &allRDatasets = mapUserDatasets->GetAllReverse();
            utype iUserDatasets = mapUserDatasets->GetNumSet();
            utype dd;
            for (dd = 0; dd < iUserDatasets; dd++) {
                i = allRDatasets[dd];
                //fprintf(stdout, "%d %d %s\n", i, dd, vecstrDatasets[i].c_str());
                CSeekIntIntMap *si = vc[i]->GetGeneMap();

                if (toCutoff && sArgs.order_stat_single_gene_query_flag == 1) {
                    if (mapstrintGene.find(vecstrQuery[0]) == mapstrintGene.end())
                        continue;
                    if (CSeekTools::IsNaN(si->GetForward(
                            mapstrintGene[vecstrQuery[0]])))
                        continue;
                    float gene_var = vc[i]->GetGeneVariance(mapstrintGene[vecstrQuery[0]]);
                    if (gene_var < cutoff) continue;
                    if (isFirst) {
                        isFirst = false;
                        if (toOutput) {
                            fprintf(out, "%s", vecstrDatasets[i].c_str());
                        } else {
                            fprintf(stdout, "%s", vecstrDatasets[i].c_str());
                        }
                    } else {
                        if (toOutput) {
                            fprintf(out, " %s", vecstrDatasets[i].c_str());
                        } else {
                            fprintf(stdout, " %s", vecstrDatasets[i].c_str());
                        }
                    }

                } else {
                    utype present = 0;
                    for (j = 0, present = 0; j < vecstrQuery.size(); j++) {
                        if (mapstrintGene.find(vecstrQuery[j]) == mapstrintGene.end())
                            continue;
                        if (CSeekTools::IsNaN(si->GetForward(
                                mapstrintGene[vecstrQuery[j]])))
                            continue;
                        count[j]++;
                        present++;
                    }
                    if (present == vecstrQuery.size()) {
                        if (isFirst) {
                            isFirst = false;
                            if (toOutput) {
                                fprintf(out, "%s", vecstrDatasets[i].c_str());
                            } else {
                                fprintf(stdout, "%s", vecstrDatasets[i].c_str());
                            }
                        } else {
                            if (toOutput) {
                                fprintf(out, " %s", vecstrDatasets[i].c_str());
                            } else {
                                fprintf(stdout, " %s", vecstrDatasets[i].c_str());
                            }
                        }
                        //fprintf(stderr, "%s\t%s\t%d\t%d\n",
                        //	vecstrDatasets[i].c_str(), vecstrDP[i].c_str(),
                        //	present, si->GetNumSet());
                    }
                }
            }

            for (j = 0; j < vecstrQuery.size(); j++)
                fprintf(stderr, "Gene %s: %d\n", vecstrQuery[j].c_str(), count[j]);

            if (toOutput) {
                fprintf(out, "\n");
            } else {
                fprintf(stdout, "\n");
            }

        }

        if (toOutput) {
            fclose(out);
        }
        return 0;
    }

    if (sArgs.dataset2_flag == 1) {
        //string dset_list = sArgs.dset_list_arg;
        vector <string> vecstrDP;
        if (!CSeekTools::ReadListTwoColumns(sArgs.db_arg, vecstrDatasets, vecstrDP))
            return -1;
        size_t iDatasets = vecstrDatasets.size();
        string strSinfoInputDirectory = sArgs.dir_sinfo_in_arg;
        if (strSinfoInputDirectory == "NA") {
            fprintf(stderr, "Error: requires sinfo input directory\n");
            return -1;
        }
        for (i = 0; i < iDatasets; i++) {
            CSeekDataset *vc = new CSeekDataset();
            string strFileStem = vecstrDatasets[i];
            string strSinfoPath = strSinfoInputDirectory + "/" +
                                  strFileStem + ".sinfo";
            vc->ReadDatasetAverageStdev(strSinfoPath);
            fprintf(stderr, "%s\n", strFileStem.c_str());
            fprintf(stderr, "%.2f\t%.2f\n", vc->GetDatasetAverage(), vc->GetDatasetStdev());
            delete vc;
        }
        return 0;
    }

    if (sArgs.databaselet_flag == 1) {
        string db = sArgs.db_arg;
        string dset_list = sArgs.dset_list_arg;
        string dir_in = sArgs.dir_in_arg;
        string dir_prep = sArgs.dir_prep_in_arg;
        string single_query = sArgs.single_query_arg;

        if (db == "NA" || dset_list == "NA" || dir_in == "NA" ||
            dir_prep == "NA" || single_query == "NA") {
            fprintf(stderr, "Requires: -x, -X, -d -p -q\n");
            return -1;
        }

        bool useNibble = false;
        if (sArgs.is_nibble_flag == 1) useNibble = true;
        CDatabase DB(useNibble);
        vector <string> vecstrDP;
        vector <string> vecstrQuery;
        if (!CSeekTools::ReadListTwoColumns(sArgs.db_arg, vecstrDatasets, vecstrDP))
            return -1;
        if (!CSeekTools::ReadMultiGeneOneLine(sArgs.single_query_arg, vecstrQuery))
            return -1;

        string strInputDirectory = sArgs.dir_in_arg;
        DB.Open(strInputDirectory);

        size_t iDatasets = DB.GetDatasets();
        size_t iGenes = DB.GetGenes();

        size_t j, k;
        vector < CSeekDataset * > vc;
        vc.clear();
        vc.resize(iDatasets);
        for (i = 0; i < iDatasets; i++) {
            vc[i] = new CSeekDataset();
            string strPrepInputDirectory = sArgs.dir_prep_in_arg;
            string strFileStem = vecstrDatasets[i];
            string strAvgPath = strPrepInputDirectory + "/" +
                                strFileStem + ".gavg";
            string strPresencePath = strPrepInputDirectory + "/" +
                                     strFileStem + ".gpres";
            vc[i]->ReadGeneAverage(strAvgPath);
            vc[i]->ReadGenePresence(strPresencePath);
        }

        vector<char> cQuery;
        vector <utype> allQ;
        CSeekTools::InitVector(cQuery, iGenes, (char) 0);

        for (i = 0; i < vecstrQuery.size(); i++) {
            k = DB.GetGene(vecstrQuery[i]);
            if (k == -1) continue;
            cQuery[k] = 1;
        }

        for (k = 0; k < iGenes; k++)
            if (cQuery[k] == 1)
                allQ.push_back(k);

        for (i = 0; i < iDatasets; i++) {
            vc[i]->InitializeGeneMap();
            vc[i]->InitializeQueryBlock(allQ);
        }

        vector<unsigned char> *Q =
                new vector<unsigned char>[vecstrQuery.size()];

        for (i = 0; i < vecstrQuery.size(); i++)
            if (!DB.GetGene(vecstrQuery[i], Q[i]))
                cerr << "Gene does not exist" << endl;

        //printf("Before"); getchar();
        for (i = 0; i < vecstrQuery.size(); i++) {
            if (DB.GetGene(vecstrQuery[i]) == -1) continue;
            size_t m = DB.GetGene(vecstrQuery[i]);
            size_t l = 0;
            for (j = 0; j < iDatasets; j++) {
                CSeekIntIntMap *qu = vc[j]->GetDBMap();
                if (qu == NULL) continue;
                unsigned char **r = vc[j]->GetMatrix();
                utype query = qu->GetForward(m);
                if (CSeekTools::IsNaN(query)) continue;
                for (k = 0; k < iGenes; k++) {
                    unsigned char c = Q[i][k * iDatasets + j];
                    r[query][k] = c;
                }
            }
        }

        for (i = 0; i < iDatasets; i++) {
            printf("Dataset %ld\n", i);
            unsigned char **r = vc[i]->GetMatrix();
            CSeekIntIntMap *mapG = vc[i]->GetGeneMap();
            CSeekIntIntMap *mapDB = vc[i]->GetDBMap();
            if (mapDB == NULL) continue;
            for (j = 0; j < mapDB->GetNumSet(); j++) {
                if (vecstrDatasets[i] == "GSE19470.GPL5175.pcl") {
                    printf("Row %ld\n", j);
                    for (k = 0; k < mapG->GetNumSet(); k++) {
                        printf("%d ", r[j][mapG->GetReverse(k)]);
                    }
                    printf("\n");
                }
            }
        }

    }

#ifdef WIN32
    pthread_win32_process_detach_np( );
#endif // WIN32
    return 0;
}

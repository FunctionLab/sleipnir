#include "stdafx.h"
#include "cmdline.h"

bool transfer(CDat &Dat, vector <vector<float>> &mat, CSeekIntIntMap *geneMap,
              vector <string> &vecstrGenes) {

    int i, j;
    vector <utype> veciGenes;
    veciGenes.clear();
    veciGenes.resize(vecstrGenes.size());
    for (i = 0; i < vecstrGenes.size(); ++i)
        veciGenes[i] = Dat.GetGene(vecstrGenes[i]);

    mat.resize(vecstrGenes.size());
    for (i = 0; i < vecstrGenes.size(); i++)
        CSeekTools::InitVector(mat[i], vecstrGenes.size(), CMeta::GetNaN());

    for (i = 0; i < vecstrGenes.size(); i++) {
        utype s = veciGenes[i];
        if (CSeekTools::IsNaN(s)) continue;
        geneMap->Add(i);
        float *v = Dat.GetFullRow(s);
        for (j = 0; j < vecstrGenes.size(); j++) {
            utype t = veciGenes[j];
            if (CSeekTools::IsNaN(t)) continue;
            if (CMeta::IsNaN(v[t])) continue;
            mat[i][j] = Dat.Get(s, t);
        }
        free(v);
    }
    return true;
}

//Add scores from two vectors and integrate them using alpha
bool integrate(vector<float> &d1, vector<float> &d2,
               vector<float> &dest, CSeekIntIntMap *geneMap, int k, float alpha) {
    utype i;
    CSeekTools::InitVector(dest, geneMap->GetSize(), (float) CMeta::GetNaN());
    const vector <utype> &allGenes = geneMap->GetAllReverse();
    for (i = 0; i < geneMap->GetNumSet(); i++) {
        utype gi = allGenes[i];
        dest[gi] = (1.0 - alpha) * d1[gi] + alpha / k * d2[gi];
    }
    return true;
}

//initialize score vector
bool init_score(vector<float> &dest, CSeekIntIntMap *geneMap) {
    CSeekTools::InitVector(dest, geneMap->GetSize(), (float) CMeta::GetNaN());
    utype i;
    const vector <utype> &allGenes = geneMap->GetAllReverse();
    for (i = 0; i < geneMap->GetNumSet(); i++)
        dest[allGenes[i]] = 0;
    return true;
}

bool add_score(vector<float> &src, vector<float> &dest, CSeekIntIntMap *geneMap) {
    const vector <utype> &allGenes = geneMap->GetAllReverse();
    utype i, j;
    for (i = 0; i < geneMap->GetNumSet(); i++) {
        utype gi = allGenes[i];
        dest[gi] += src[gi];
    }
    return true;
}

bool search_one_dab(vector<float> &gene_score,
                    CDat &mat, const size_t numGenes,
                    CSeekIntIntMap &d1,
                    vector <utype> &indexConvReverse,
                    vector<float> &q_weight,
                    vector <vector<float>> &nq_weight
) {
    //q_weight is query presence - 1.0 if gene at the index i is a query gene

    vector<float> gene_count;
    CSeekTools::InitVector(gene_score, numGenes, (float) CMeta::GetNaN());
    CSeekTools::InitVector(gene_count, numGenes, (float) 0);
    utype qqi;
    utype kk, k;

    const vector <utype> &allGenes = d1.GetAllReverse();
    for (kk = 0; kk < d1.GetNumSet(); kk++) {
        utype k = allGenes[kk];
        utype gi = indexConvReverse[k];
        gene_score[gi] = 0;
    }

    for (qqi = 0; qqi < d1.GetNumSet(); qqi++) {
        utype qi = allGenes[qqi];
        utype qq = indexConvReverse[qi];
        if (q_weight[qq] == 0) //not a query gene
            continue;
        //now a query gene
        float *vc = mat.GetFullRow(qi);
        for (kk = 0; kk < d1.GetNumSet(); kk++) {
            utype k = allGenes[kk];
            utype gi = indexConvReverse[k];
            float fl = vc[k];
            gene_score[gi] += fl;
            gene_count[gi] += 1.0;
        }
        delete[] vc;
    }

    for (kk = 0; kk < d1.GetNumSet(); kk++) {
        utype k = allGenes[kk];
        utype gi = indexConvReverse[k];
        gene_score[gi] /= gene_count[gi];
    }

    utype ind;
    for (ind = 0; ind < nq_weight.size(); ind++) {
        vector<float> ngene_count, ngene_score;
        CSeekTools::InitVector(ngene_score, numGenes, (float) CMeta::GetNaN());
        CSeekTools::InitVector(ngene_count, numGenes, (float) 0);

        for (kk = 0; kk < d1.GetNumSet(); kk++)
            ngene_score[(utype) indexConvReverse[(utype) allGenes[kk]]] = 0;

        for (qqi = 0; qqi < d1.GetNumSet(); qqi++) {
            utype qi = allGenes[qqi];
            utype qq = indexConvReverse[qi];
            if (nq_weight[ind][qq] == 0) //not a NOT query gene
                continue;
            //now a NOT query gene
            float *vc = mat.GetFullRow(qi);
            for (kk = 0; kk < d1.GetNumSet(); kk++) {
                utype k = allGenes[kk];
                utype gi = indexConvReverse[k];
                float fl = vc[k];
                ngene_score[gi] += fl;
                ngene_count[gi] += 1.0;
            }
            delete[] vc;
        }
        for (kk = 0; kk < d1.GetNumSet(); kk++) {
            utype gi = indexConvReverse[(utype) allGenes[kk]];
            ngene_score[gi] /= ngene_count[gi];
            gene_score[gi] -= ngene_score[gi]; //subtract the correlation to NOT genes
        }
    }

    return true;
}

bool get_score(vector<float> &gene_score,
               CSparseFlatMatrix<float> &mat,
               CSeekIntIntMap *geneMap, vector<float> &q_weight) {
    vector<float> gene_count;
    int numGenes = geneMap->GetSize();
    CSeekTools::InitVector(gene_score, numGenes, (float) CMeta::GetNaN());
    CSeekTools::InitVector(gene_count, numGenes, (float) 0);
    int qi = 0;
    const vector <utype> &allGenes = geneMap->GetAllReverse();
    utype kk, k;

    for (kk = 0; kk < geneMap->GetNumSet(); kk++) {
        utype gi = allGenes[kk];
        gene_score[gi] = 0;
    }
    for (qi = 0; qi < geneMap->GetNumSet(); qi++) {
        utype qq = allGenes[qi];
        if (q_weight[qq] == 0)
            continue;
        const vector <CPair<float>> &vc = mat.GetRow(qq);
        for (kk = 0; kk < vc.size(); kk++) {
            float fl = vc[kk].v;
            utype gi = vc[kk].i;
            gene_score[gi] += fl * q_weight[qq];
        }
        for (kk = 0; kk < geneMap->GetNumSet(); kk++) {
            utype gi = allGenes[kk];
            gene_count[gi] += q_weight[qq];
        }
        /*for(kk=0; kk<geneMap->GetNumSet(); kk++){
            utype gi = allGenes[kk];
            map<utype,float>::iterator it;
            float fl = 0;
            if((it=mat[qq].find(gi))==mat[qq].end()){
                fl = 0;
            }else{
                fl = it->second;
            }
            //float fl = mat[qq][gi];
            gene_score[gi] += fl * q_weight[qq];
            gene_count[gi] += q_weight[qq];
        }*/
    }

    for (k = 0; k < geneMap->GetNumSet(); k++) {
        utype gi = allGenes[k];
        gene_score[gi] /= gene_count[gi];
        if (gene_count[gi] == 0) {
            gene_score[gi] = 0;
        }
    }
    return true;
}

bool get_score(vector<float> &gene_score, CFullMatrix<float> &mat,
               CSeekIntIntMap *geneMap, vector<float> &q_weight, vector<float> &geneAverage) {
    bool enableGeneAverage = true;
    if (geneAverage.size() == 0)
        enableGeneAverage = false;

    vector<float> gene_count;
    int numGenes = geneMap->GetSize();
    CSeekTools::InitVector(gene_score, numGenes, (float) CMeta::GetNaN());
    CSeekTools::InitVector(gene_count, numGenes, (float) 0);
    int qi = 0;
    const vector <utype> &allGenes = geneMap->GetAllReverse();
    utype kk, k;
    for (kk = 0; kk < geneMap->GetNumSet(); kk++) {
        utype gi = allGenes[kk];
        gene_score[gi] = 0;
    }
    for (qi = 0; qi < geneMap->GetNumSet(); qi++) {
        utype qq = allGenes[qi];
        if (q_weight[qq] == 0)
            continue;
        for (kk = 0; kk < geneMap->GetNumSet(); kk++) {
            utype gi = allGenes[kk];
            float fl = mat.Get(qq, gi);
            fl = max(min(fl, 3.20f), -3.20f);

            if (enableGeneAverage) {
                float fl_avg = geneAverage[gi];
                if (fl_avg == CMeta::GetNaN()) {
                    fprintf(stderr, "Warning: average value for %d NaN\n", gi);
                }
                fl = fl - fl_avg;
            }

            gene_score[gi] += fl * q_weight[qq];
        }
        for (kk = 0; kk < geneMap->GetNumSet(); kk++) {
            utype gi = allGenes[kk];
            gene_count[gi] += q_weight[qq];
        }
    }
    for (k = 0; k < geneMap->GetNumSet(); k++) {
        utype gi = allGenes[k];
        gene_score[gi] /= gene_count[gi];
        if (gene_count[gi] == 0)
            gene_score[gi] = 0;
    }
    return true;
}

/*bool iterate(vector<float> &src, vector<float> &dest, 
	vector<vector<float> > &tr, CSeekIntIntMap *geneMap,
	float alpha, int numIterations){
	
	const vector<utype> &allGenes = geneMap->GetAllReverse();
	dest.resize(src.size());
	vector<float> backup;
	CSeekTools::InitVector(backup, src.size(), CMeta::GetNaN());
	int i, j, k;
	for(i=0; i<dest.size(); i++)
		dest[i] = src[i];	
		
	for(i=0; i<numIterations; i++){
		for(j=0; j<geneMap->GetNumSet(); j++){
			utype gi = allGenes[j];
			backup[gi] = dest[gi];
		}
		get_score(dest, tr, geneMap, backup);
		for(j=0; j<geneMap->GetNumSet(); j++){
			utype gi = allGenes[j];
			dest[gi] = (1.0 - alpha) * src[j] + alpha * dest[gi];
		}
	}
	return true;
}*/

//is_query is like a gold-standard gene list
bool weight_fast(vector<char> &is_query, vector<float> &d1,
                 CSeekIntIntMap *geneMap, float &w) {
    utype i;
    w = 0;
    const vector <utype> &allGenes = geneMap->GetAllReverse();
    for (i = 0; i < geneMap->GetNumSet(); i++) {
        utype gi = allGenes[i];
        if (is_query[gi] == 1) {
            //if(CMeta::IsNaN(d1[gi])) continue;
            w += d1[gi];
        }
    }
    return true;
}

bool weight(vector<char> &is_query, vector<float> &d1,
            CSeekIntIntMap *geneMap, float rbp_p, float &w) {
    vector <AResultFloat> ar;
    ar.resize(geneMap->GetNumSet());
    utype i;
    const vector <utype> &allGenes = geneMap->GetAllReverse();
    for (i = 0; i < geneMap->GetNumSet(); i++) {
        utype gi = allGenes[i];
        ar[i].i = gi;
        ar[i].f = d1[gi];
    }
    int MAX = 1000;
    nth_element(ar.begin(), ar.begin() + MAX, ar.end());
    sort(ar.begin(), ar.begin() + MAX);
    w = 0;
    for (i = 0; i < MAX; i++)
        if (is_query[ar[i].i] == 1)
            w += pow(rbp_p, i);
    w *= (1.0 - rbp_p);
    return true;
}

bool cv_weight_LOO(vector <utype> &query, CSparseFlatMatrix<float> &mat,
                   CSeekIntIntMap *geneMap, float rbp_p, float &tot_w) {
    //leave one out cross-validation
    utype i, j;
    int numGenes = geneMap->GetSize();
    tot_w = 0;

    for (i = 0; i < query.size(); i++) {
        vector<char> is_query;
        vector<float> q_weight;
        float w = 0;
        CSeekTools::InitVector(is_query, numGenes, (char) 0);
        CSeekTools::InitVector(q_weight, numGenes, (float) 0);
        for (j = 0; j < query.size(); j++) {
            if (i == j) continue;
            q_weight[query[j]] = 1.0;
        }
        is_query[query[i]] = 1;
        vector<float> gene_score;
        get_score(gene_score, mat, geneMap, q_weight);
        for (j = 0; j < query.size(); j++) {
            if (i == j) continue;
            gene_score[query[j]] = 0;
        }
        //weight_fast(is_query, gene_score, geneMap, w);
        weight(is_query, gene_score, geneMap, rbp_p, w);
        tot_w += w;
    }
    tot_w /= query.size();
    return true;
}

bool cv_weight(vector <utype> &query, CSparseFlatMatrix<float> &mat,
               CSeekIntIntMap *geneMap, float rbp_p, float &tot_w) {
    //leave one in cross-validation
    utype i, j;
    int numGenes = geneMap->GetSize();
    tot_w = 0;

    for (i = 0; i < query.size(); i++) {
        vector<char> is_query;
        vector<float> q_weight;
        CSeekTools::InitVector(is_query, numGenes, (char) 0);
        CSeekTools::InitVector(q_weight, numGenes, (float) 0);
        q_weight[query[i]] = 1.0;
        float w = 0;
        for (j = 0; j < query.size(); j++) {
            if (i == j) continue;
            is_query[query[j]] = 1;
        }
        vector<float> gene_score;
        get_score(gene_score, mat, geneMap, q_weight);
        gene_score[query[i]] = 0;
        weight(is_query, gene_score, geneMap, rbp_p, w);
        //fprintf(stderr, "Q%d %.3f\n", i, w);
        //weight_fast(is_query, gene_score, geneMap, w);
        tot_w += w;
    }
    tot_w /= query.size();
    return true;
}

//accepts fullmatrix, leave one in cross-validation
bool cv_weight(vector <utype> &query, CFullMatrix<float> &mat,
               CSeekIntIntMap *geneMap, float rbp_p, float &tot_w, vector<float> &geneAverage) {
    utype i, j;
    int numGenes = geneMap->GetSize();
    tot_w = 0;
    for (i = 0; i < query.size(); i++) {
        vector<char> is_query;
        vector<float> q_weight;
        CSeekTools::InitVector(is_query, numGenes, (char) 0);
        CSeekTools::InitVector(q_weight, numGenes, (float) 0);
        q_weight[query[i]] = 1.0;
        float w = 0;
        for (j = 0; j < query.size(); j++) {
            if (i == j) continue;
            is_query[query[j]] = 1;
        }
        vector<float> gene_score;
        get_score(gene_score, mat, geneMap, q_weight, geneAverage);
        gene_score[query[i]] = 0;
        weight(is_query, gene_score, geneMap, rbp_p, w);
        tot_w += w;
    }
    tot_w /= query.size();
    return true;
}

bool spell_weight(vector <utype> &query, CFullMatrix<float> &mat,
                  CSeekIntIntMap *geneMap, float &tot_w) {
    utype i, j;
    int numGenes = geneMap->GetSize();
    tot_w = 0;
    int count = 0;
    for (i = 0; i < query.size(); i++) {
        utype ui = geneMap->GetForward(query[i]);
        if (ui == (utype) - 1) {
            continue;
        }
        float w = 0;
        for (j = 0; j < query.size(); j++) {
            if (i == j) continue;
            utype uj = geneMap->GetForward(query[j]);
            if (uj == (utype) - 1) {
                continue;
            }
            w = mat.Get(ui, uj);
            if (w < 0) {
                w = 0;
            }
            if (isnan(w) || isinf(w)) {
                fprintf(stderr, "Warning: nan or inf detected for Mat(%d, %d)\n", ui, uj);
                continue;
            }
            tot_w += w;
            count++;
        }
    }
    tot_w /= (float) count;
    return true;
}

int main(int iArgs, char **aszArgs) {
    static const size_t c_iBuffer = 1024;
    char acBuffer[c_iBuffer];

    gengetopt_args_info sArgs;
    ifstream ifsm;
    istream *pistm;
    vector <string> vecstrLine, vecstrGenes, vecstrDBs, vecstrQuery;
    utype i, qi, j, k, l, kk;

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }

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
        if (!(i = atoi(vecstrLine[0].c_str()))) {
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

    fprintf(stderr, "Finished reading gene map\n");

    string dab_dir = sArgs.dab_dir_arg;
    string output_dir = sArgs.dir_out_arg;
    int numGenes = vecstrGenes.size();
    vector <vector<string>> vecstrAllQuery;
    fprintf(stderr, "Reading queries\n");
    if (!CSeekTools::ReadMultipleQueries(sArgs.query_arg, vecstrAllQuery))
        return -1;
    fprintf(stderr, "Finished reading queries\n");

    vector <vector<utype>> qu;
    qu.resize(vecstrAllQuery.size());
    for (i = 0; i < vecstrAllQuery.size(); i++) {
        qu[i] = vector<utype>();
        for (j = 0; j < vecstrAllQuery[i].size(); j++)
            qu[i].push_back(mapstriGenes[vecstrAllQuery[i][j]]);
    }

    vector < vector < vector < string > > > vecstrAllNotQuery;
    vecstrAllNotQuery.resize(vecstrAllQuery.size());
    string not_query = sArgs.not_query_arg;
    if (not_query != "NA") {
        fprintf(stderr, "Reading NOT queries\n");
        if (!CSeekTools::ReadMultipleNotQueries(sArgs.not_query_arg,
                                                vecstrAllNotQuery))
            return -1;
    }

    if (sArgs.visualize_flag == 1) {
        string dab_base = sArgs.dab_basename_arg;
        string file1 = dab_dir + "/" + dab_base + ".dab";
        float cutoff_par = sArgs.cutoff_arg;
        string genome = sArgs.genome_arg;
        vector <string> s1, s2;
        CSeekTools::ReadListTwoColumns(genome.c_str(), s1, s2);

        CGenome g;
        g.Open(s1);
        for (i = 0; i < s2.size(); i++) {
            CGene &g1 = g.GetGene(g.FindGene(s1[i]));
            g.AddSynonym(g1, s2[i]);
        }

        CDat CD;
        CD.Open(file1.c_str(), false, 2, false, false, false);

        CSeekIntIntMap d1(CD.GetGenes());
        vector <utype> indexConvReverse;
        CSeekTools::InitVector(indexConvReverse, vecstrGenes.size(), (utype) - 1);
        for (i = 0; i < CD.GetGenes(); i++) {
            map<string, size_t>::iterator it = mapstriGenes.find(CD.GetGene(i));
            if (it == mapstriGenes.end()) continue;
            indexConvReverse[i] = it->second;
            d1.Add(i);
        }

        //Visualize
        for (j = 0; j < vecstrAllQuery.size(); j++) {
            vector <string> vec_s;
            vector <utype> vec_g;
            for (k = 0; k < vecstrAllQuery[j].size(); k++) {
                size_t ind = CD.GetGeneIndex(vecstrAllQuery[j][k]);
                if (ind == (size_t) - 1) continue;
                vec_g.push_back(CD.GetGeneIndex(vecstrAllQuery[j][k]));
                vec_s.push_back(vecstrAllQuery[j][k]);
            }
            CDat V;
            V.Open(vec_s);
            for (k = 0; k < vec_s.size(); k++) {
                for (l = k + 1; l < vec_s.size(); l++) {
                    V.Set(k, l, CD.Get(vec_g[k], vec_g[l]));
                }
            }
            fprintf(stderr, "Query %d\n", j);

            if (sArgs.print_distr_flag == 1) {
                float min = 9999;
                float max = -1;
                for (k = 0; k < vec_s.size(); k++) {
                    for (l = k + 1; l < vec_s.size(); l++) {
                        float v = CD.Get(vec_g[k], vec_g[l]);
                        if (CMeta::IsNaN(v)) continue;
                        if (v > max) max = v;
                        if (v < min) min = v;
                    }
                }
                float init = min;
                float step = (max - min) / 100.0;
                float upper = max;
                float cutoff = init;
                while (cutoff < upper) {
                    int count = 0;
                    for (k = 0; k < vec_s.size(); k++) {
                        for (l = k + 1; l < vec_s.size(); l++) {
                            if (!CMeta::IsNaN(V.Get(k, l)) && V.Get(k, l) >= cutoff) {
                                count++;
                            }
                        }
                    }
                    fprintf(stderr, "%.5f\t%d\n", cutoff, count);
                    cutoff += step;
                }
            }

            sprintf(acBuffer, "%s/%d.dot", output_dir.c_str(), j);
            ofstream ot(acBuffer);
            V.SaveDOT(ot, cutoff_par, &g, false, true, NULL, NULL);
        }

    }

    if (sArgs.combined_flag == 1) {
        string dab_base = sArgs.dab_basename_arg;
        string file1 = dab_dir + "/" + dab_base + ".dab";

        vector <vector<float>> q_weight;
        vector < vector < vector < float > > > nq_weight;

        q_weight.resize(vecstrAllQuery.size());
        for (i = 0; i < vecstrAllQuery.size(); i++) {
            CSeekTools::InitVector(q_weight[i], numGenes, (float) 0);
            for (j = 0; j < vecstrAllQuery[i].size(); j++)
                q_weight[i][mapstriGenes[vecstrAllQuery[i][j]]] = 1;
        }

        nq_weight.resize(vecstrAllNotQuery.size());
        if (not_query != "NA") {
            for (i = 0; i < vecstrAllNotQuery.size(); i++) {
                nq_weight[i].resize(vecstrAllNotQuery[i].size());
                for (j = 0; j < vecstrAllNotQuery[i].size(); j++) {
                    CSeekTools::InitVector(nq_weight[i][j], numGenes, (float) 0);
                    for (k = 0; k < vecstrAllNotQuery[i][j].size(); k++)
                        nq_weight[i][j][mapstriGenes[vecstrAllNotQuery[i][j][k]]] = 1;
                }
            }
        }

        fprintf(stderr, "Start reading dab\n");
        CDat CD;
        CD.Open(file1.c_str(), false, 2, false, false, false);
        fprintf(stderr, "Finished reading dab\n");

        CSeekIntIntMap d1(CD.GetGenes());
        vector <utype> indexConvReverse;
        CSeekTools::InitVector(indexConvReverse, vecstrGenes.size(), (utype) - 1);
        for (i = 0; i < CD.GetGenes(); i++) {
            map<string, size_t>::iterator it = mapstriGenes.find(CD.GetGene(i));
            if (it == mapstriGenes.end()) continue;
            indexConvReverse[i] = it->second;
            d1.Add(i);
        }

        vector <vector<float>> final_score;
        final_score.resize(vecstrAllQuery.size());
        for (j = 0; j < vecstrAllQuery.size(); j++)
            CSeekTools::InitVector(final_score[j], vecstrGenes.size(),
                                   (float) CMeta::GetNaN());

        const vector <utype> &allGenes = d1.GetAllReverse();
        for (j = 0; j < vecstrAllQuery.size(); j++) {
            vector<float> tmp_score;
            search_one_dab(tmp_score, CD, vecstrGenes.size(), d1, indexConvReverse,
                           q_weight[j], nq_weight[j]);
            for (kk = 0; kk < d1.GetNumSet(); kk++) {
                utype k = allGenes[kk];
                utype gi = indexConvReverse[k];
                final_score[j][gi] = tmp_score[gi];
            }
        }
        fprintf(stderr, "Writing output\n");

        for (j = 0; j < vecstrAllQuery.size(); j++) {
            for (k = 0; k < final_score[j].size(); k++) {
                if (CMeta::IsNaN(final_score[j][k])) {
                    final_score[j][k] = -320;
                    continue;
                }
            }
            sprintf(acBuffer, "%s/%d.query", output_dir.c_str(), j);
            CSeekTools::WriteArrayText(acBuffer, vecstrAllQuery[j]);
            sprintf(acBuffer, "%s/%d.gscore", output_dir.c_str(), j);
            CSeekTools::WriteArray(acBuffer, final_score[j]);
        }

        fprintf(stderr, "Finished with search\n");
        if (sArgs.print_distr_flag == 0 && sArgs.generate_dot_flag == 0) {
            return 0;
        }

        float cutoff_par = sArgs.cutoff_arg;
        string genome = sArgs.genome_arg;
        vector <string> s1, s2;
        CSeekTools::ReadListTwoColumns(genome.c_str(), s1, s2);

        CGenome g;
        g.Open(s1);
        for (i = 0; i < s2.size(); i++) {
            CGene &g1 = g.GetGene(g.FindGene(s1[i]));
            g.AddSynonym(g1, s2[i]);
        }

        //Visualize
        for (j = 0; j < vecstrAllQuery.size(); j++) {
            //fprintf(stderr, "Query %d\n", j);
            vector <AResultFloat> ar;
            ar.resize(final_score[j].size());
            for (k = 0; k < final_score[j].size(); k++) {
                ar[k].i = k;
                ar[k].f = final_score[j][k];
            }
            sort(ar.begin(), ar.end());
            vector <utype> vec_g;
            vector <string> vec_s;
            vector<float> node_w;
            int FIRST = sArgs.top_genes_arg;
            for (k = 0; k < FIRST; k++) {
                if (ar[k].f == -320) break;
                vec_g.push_back(CD.GetGeneIndex(vecstrGenes[ar[k].i]));
                vec_s.push_back(vecstrGenes[ar[k].i]);
                node_w.push_back(0);
            }
            //if(vec_g.size()!=FIRST) continue;
            for (k = 0; k < vecstrAllQuery[j].size(); k++) {
                size_t ind = CD.GetGeneIndex(vecstrAllQuery[j][k]);
                if (ind == (size_t) - 1) continue;
                vec_g.push_back(ind);
                vec_s.push_back(vecstrAllQuery[j][k]);
                node_w.push_back(1.0);
            }

            CDat V;
            V.Open(vec_s);
            for (k = 0; k < vec_s.size(); k++) {
                for (l = k + 1; l < vec_s.size(); l++) {
                    float v = CD.Get(vec_g[k], vec_g[l]);
                    V.Set(k, l, v);
                }
            }

            if (sArgs.print_distr_flag == 1) {
                float min = 9999;
                float max = -1;
                for (k = 0; k < vec_s.size(); k++) {
                    for (l = k + 1; l < vec_s.size(); l++) {
                        float v = CD.Get(vec_g[k], vec_g[l]);
                        if (CMeta::IsNaN(v)) continue;
                        if (v > max) max = v;
                        if (v < min) min = v;
                    }
                }

                float init = min;
                float step = (max - min) / 100.0;
                float upper = max;

                float cutoff = init;
                while (cutoff < upper) {
                    int count = 0;
                    for (k = 0; k < vec_s.size(); k++) {
                        for (l = k + 1; l < vec_s.size(); l++) {
                            if (!CMeta::IsNaN(V.Get(k, l)) && V.Get(k, l) >= cutoff) {
                                count++;
                            }
                        }
                    }
                    fprintf(stderr, "%.5f\t%d\n", cutoff, count);
                    cutoff += step;
                }
            }

            if (sArgs.generate_dot_flag == 1) {
                sprintf(acBuffer, "%s/%d.dot", output_dir.c_str(), j);
                ofstream ot(acBuffer);
                V.SaveDOT(ot, cutoff_par, &g, false, true, &node_w, NULL);
                //V.SaveDOT(ot, 0.0001, NULL, true, false, NULL, NULL);
            }
        }

    }

    if (sArgs.test_flag == 1) {
        string dab_base = sArgs.dab_basename_arg;
        string file1 = dab_dir + "/" + dab_base + ".half";

        vector <vector<float>> q_weight;
        q_weight.resize(vecstrAllQuery.size());
        for (i = 0; i < vecstrAllQuery.size(); i++) {
            CSeekTools::InitVector(q_weight[i], numGenes, (float) 0);
            for (j = 0; j < vecstrAllQuery[i].size(); j++)
                q_weight[i][mapstriGenes[vecstrAllQuery[i][j]]] = 1;
        }

        CHalfMatrix<float> res;
        res.Initialize(vecstrGenes.size());
        ifstream istm1;
        uint32_t dim;
        istm1.open(file1.c_str(), ios_base::binary);
        istm1.read((char *) (&dim), sizeof(dim));
        float *adScores = new float[dim - 1];
        for (i = 0; (i + 1) < dim; i++) {
            istm1.read((char *) adScores, sizeof(*adScores) * (dim - i - 1));
            res.Set(i, adScores);
        }
        delete[] adScores;
        istm1.close();

        string s1 = "uc003vor.2";
        vector <string> r1;
        r1.push_back("uc003vos.2");
        r1.push_back("uc003vop.1");
        r1.push_back("uc011jwz.1");
        r1.push_back("uc002rgw.1");

        for (i = 0; i < r1.size(); i++) {
            size_t i1 = mapstriGenes[s1];
            size_t j1 = mapstriGenes[r1[i]];
            fprintf(stderr, "%s %s %.5f\n", s1.c_str(), r1[i].c_str(), res.Get(i1, j1));
        }
    }

    if (sArgs.testcount_flag == 1) {
        string dab_base = sArgs.dab_basename_arg;
        string file1 = dab_dir + "/" + dab_base + ".pair_count";

        vector <vector<float>> q_weight;
        q_weight.resize(vecstrAllQuery.size());
        for (i = 0; i < vecstrAllQuery.size(); i++) {
            CSeekTools::InitVector(q_weight[i], numGenes, (float) 0);
            for (j = 0; j < vecstrAllQuery[i].size(); j++)
                q_weight[i][mapstriGenes[vecstrAllQuery[i][j]]] = 1;
        }

        CHalfMatrix<unsigned short> res;
        res.Initialize(vecstrGenes.size());
        ifstream istm1;
        uint32_t dim;
        istm1.open(file1.c_str(), ios_base::binary);
        istm1.read((char *) (&dim), sizeof(dim));
        unsigned short *adScores = new unsigned short[dim - 1];
        for (i = 0; (i + 1) < dim; i++) {
            istm1.read((char *) adScores, sizeof(*adScores) * (dim - i - 1));
            res.Set(i, adScores);
        }
        delete[] adScores;
        istm1.close();

        string s1 = "uc003vor.2";
        vector <string> r1;
        r1.push_back("uc003vos.2");
        r1.push_back("uc003vop.1");
        r1.push_back("uc011jwz.1");
        r1.push_back("uc002rgw.1");

        for (i = 0; i < r1.size(); i++) {
            size_t i1 = mapstriGenes[s1];
            size_t j1 = mapstriGenes[r1[i]];
            fprintf(stderr, "%s %s %d\n", s1.c_str(), r1[i].c_str(), res.Get(i1, j1));
        }
    }

    if (sArgs.testcombined_flag == 1) {
        string dab_base = sArgs.dab_basename_arg;
        string file3 = dab_dir + "/" + dab_base + ".gene_count";
        string file2 = dab_dir + "/" + dab_base + ".pair_count";
        string file1 = dab_dir + "/" + dab_base + ".half";

        vector <utype> gene_count;
        CSeekTools::ReadArray(file3.c_str(), gene_count);

        float max_count = *max_element(gene_count.begin(), gene_count.end());
        float min_count_required = max_count * 0.50;
        CSeekIntIntMap d1(vecstrGenes.size());
        vector <string> validGenes;
        for (i = 0; i < vecstrGenes.size(); i++) {
            if (gene_count[i] < min_count_required) continue;
            validGenes.push_back(vecstrGenes[i]);
            d1.Add(i);
        }

        CDat CD;
        CD.Open(validGenes);
        for (i = 0; i < validGenes.size(); i++) {
            for (j = i + 1; j < validGenes.size(); j++) {
                CD.Set(i, j, CMeta::GetNaN());
            }
        }
        //CHalfMatrix<float> res;
        //res.Initialize(vecstrGenes.size());
        ifstream istm1, istm2;
        uint32_t dim;
        istm1.open(file1.c_str(), ios_base::binary);
        istm1.read((char *) (&dim), sizeof(dim));
        istm2.open(file2.c_str(), ios_base::binary);
        istm2.read((char *) (&dim), sizeof(dim));
        float *adScores = new float[dim - 1];
        unsigned short *adScores2 = new unsigned short[dim - 1];
        for (i = 0; (i + 1) < dim; i++) {
            istm1.read((char *) adScores, sizeof(*adScores) * (dim - i - 1));
            istm2.read((char *) adScores2, sizeof(*adScores2) * (dim - i - 1));
            if (i % 2000 == 0)
                fprintf(stderr, "%d Finished\n", i);
            utype gi = d1.GetForward(i);
            if (CSeekTools::IsNaN(gi)) continue;

            for (j = 0; j < dim - i - 1; j++) {
                utype gj = d1.GetForward(j + i + 1);
                if (CSeekTools::IsNaN(gj) ||
                    (float) adScores2[j] < min_count_required)
                    continue;
                adScores[j] = adScores[j] / (float) adScores2[j];
                CD.Set(gi, gj, adScores[j]);
            }

            //res.Set(i, adScores);
        }
        delete[] adScores;
        delete[] adScores2;
        istm1.close();
        istm2.close();

        string s1 = "uc003vor.2";
        vector <string> r1;
        r1.push_back("uc003vos.2");
        r1.push_back("uc003vop.1");
        r1.push_back("uc011jwz.1");
        r1.push_back("uc002rgw.1");

        for (i = 0; i < r1.size(); i++) {
            //size_t i1 = mapstriGenes[s1];
            //size_t j1 = mapstriGenes[r1[i]];
            size_t i1 = CD.GetGeneIndex(s1);
            size_t j1 = CD.GetGeneIndex(r1[i]);
            fprintf(stderr, "%s %s %.5f\n", s1.c_str(), r1[i].c_str(), CD.Get(i1, j1));
        }
    }

    //Traditional DAB mode (.dab file)
    if (sArgs.tdab_flag == 1) {
        string search_mode = sArgs.tsearch_mode_arg;
        if (search_mode == "NA") {
            fprintf(stderr, "Please specify a search mode!\n");
            return 1;
        }
        vector <string> dab_list;
        CSeekTools::ReadListOneColumn(sArgs.tdab_list_arg, dab_list);

        fprintf(stderr, "Finished reading dablist\n");
        //preparing query
        vector <vector<float>> q_weight;
        q_weight.resize(vecstrAllQuery.size());
        for (i = 0; i < vecstrAllQuery.size(); i++) {
            CSeekTools::InitVector(q_weight[i], numGenes, (float) 0);
            for (j = 0; j < vecstrAllQuery[i].size(); j++) {
                map<string, size_t>::iterator it = mapstriGenes.find(vecstrAllQuery[i][j]);
                if (it != mapstriGenes.end()) {
                    q_weight[i][it->second] = 1;
                }
            }
        }

        //preparing query2
        vector <vector<unsigned int>> qq;
        qq.resize(vecstrAllQuery.size());
        for (i = 0; i < vecstrAllQuery.size(); i++) {
            qq[i] = vector<unsigned int>();
            for (j = 0; j < vecstrAllQuery[i].size(); j++) {
                map<string, size_t>::iterator it = mapstriGenes.find(vecstrAllQuery[i][j]);
                if (it != mapstriGenes.end()) {
                    qq[i].push_back(it->second);
                }
            }
        }

        //selected datasets for each query
        vector <vector<char>> selectedDataset;
        selectedDataset.resize(vecstrAllQuery.size());
        for (i = 0; i < vecstrAllQuery.size(); i++)
            CSeekTools::InitVector(selectedDataset[i], dab_list.size(), (char) 0);

        fprintf(stderr, "Reading DAB\n");
        vector <vector<float>> final_score, count, dweight;
        vector <vector<int>> freq;
        final_score.resize(vecstrAllQuery.size());
        count.resize(vecstrAllQuery.size());
        freq.resize(vecstrAllQuery.size());
        dweight.resize(vecstrAllQuery.size());
        for (j = 0; j < vecstrAllQuery.size(); j++) {
            CSeekTools::InitVector(final_score[j], vecstrGenes.size(), (float) CMeta::GetNaN());
            CSeekTools::InitVector(count[j], vecstrGenes.size(), (float) 0);
            CSeekTools::InitVector(freq[j], vecstrGenes.size(), (int) 0);
            CSeekTools::InitVector(dweight[j], dab_list.size(), (float) 0);
        }

        float threshold_g = sArgs.threshold_g_arg;
        float threshold_q = sArgs.threshold_q_arg;
        bool DEBUG = !!sArgs.debug_flag;

        size_t l;
        for (l = 0; l < dab_list.size(); l++) {
            CDat Dat;
            fprintf(stderr, "Reading %zu: %s\n", l, dab_list[l].c_str());
            string dabfile = dab_dir + "/" + dab_list[l];
            Dat.Open(dabfile.c_str(), false, 2, false, false, false);

            vector<float> geneAverage;
            size_t dab_found = dab_list[l].find(".dab");
            string gavg_dir = sArgs.gavg_dir_arg;
            if (gavg_dir == "NA") {
                fprintf(stderr, "INFO: Gene average is disabled.\n");
            } else {
                if (dab_found != string::npos) {
                    string gavgfile = gavg_dir + "/" + dab_list[l].substr(0, dab_found) + ".gavg";
                    if (!CSeekTools::ReadArray(gavgfile.c_str(), geneAverage)) {
                        fprintf(stderr, "Error encountered.\n");
                        return 1;
                    }
                } else {
                    fprintf(stderr, "Error: file name must contain .dab suffix\n");
                    return 1;
                }
            }


            fprintf(stderr, "Finished reading DAB\n");
            vector<unsigned int> veciGenes;
            veciGenes.resize(vecstrGenes.size());
            unsigned int ki;
            for (ki = 0; ki < vecstrGenes.size(); ki++)
                veciGenes[ki] = (unsigned int) Dat.GetGeneIndex(vecstrGenes[ki]);
            unsigned int s, t;
            float d;
            CSeekIntIntMap m(vecstrGenes.size());
            for (i = 0; i < vecstrGenes.size(); i++) {
                if ((s = veciGenes[i]) == (unsigned int) -1) continue;
                m.Add(i);
            }

            //Copy to matrix sm
            CFullMatrix<float> sm;
            sm.Initialize(vecstrGenes.size(), vecstrGenes.size());
            const vector <utype> &allRGenes = m.GetAllReverse();
            for (i = 0; i < m.GetNumSet(); i++) {
                unsigned int si = allRGenes[i];
                s = veciGenes[si];
                for (j = i + 1; j < m.GetNumSet(); j++) {
                    unsigned int tj = allRGenes[j];
                    t = veciGenes[allRGenes[j]];
                    if (CMeta::IsNaN(d = Dat.Get(s, t))) {
                        sm.Set(si, tj, 0);
                        sm.Set(tj, si, 0);
                    } else {
                        sm.Set(si, tj, d);
                        sm.Set(tj, si, d);
                    }
                }
                sm.Set(si, si, 0);
            }

            fprintf(stderr, "Finished copying matrix\n");

            for (j = 0; j < vecstrAllQuery.size(); j++) {
                int numPresent = 0;
                for (k = 0; k < qq[j].size(); k++) {
                    if (m.GetForward(qq[j][k]) == (unsigned int) -1) continue;
                    numPresent++;
                }
                if (search_mode == "eq" && numPresent == 0) {
                    continue;
                } else if (search_mode == "cv_loi" && numPresent <= 1) {
                    continue;
                } else if (search_mode == "spell" && numPresent <= 1) {
                    continue;
                }
                int minRequired = 1;
                if (search_mode == "cv_loi" || search_mode == "spell")
                    minRequired = 2;
                int numThreshold = (int) (threshold_q * qq[j].size());
                numThreshold = max(minRequired, numThreshold);
                if (numPresent > numThreshold)
                    selectedDataset[j][l] = 1;
            }

            for (j = 0; j < vecstrAllQuery.size(); j++) {
                //not enough query genes present
                if (selectedDataset[j][l] == 0) continue;

                float dw = 1.0;
                if (search_mode == "eq") {
                    dw = 1.0;
                    //}else{ //cv_loi, rbp_p = 0.99
                    //	cv_weight(qu[j], sm, &m, 0.99, dw);
                    //}
                } else if (search_mode == "cv_loi") { //cv_loi, rbp_p = 0.99
                    //fprintf(stderr, "Weighting query %d\n", (int) j);
                    cv_weight(qu[j], sm, &m, 0.99, dw, geneAverage); //returns weight to dw
                } else if (search_mode == "spell") {
                    spell_weight(qu[j], sm, &m, dw);
                }
                dweight[j][l] = dw;
                vector<float> tmp_score;
                get_score(tmp_score, sm, &m, q_weight[j], geneAverage);
                for (k = 0; k < m.GetNumSet(); k++) {
                    utype gi = allRGenes[k];
                    if (CMeta::IsNaN(final_score[j][gi]))
                        final_score[j][gi] = 0;
                    final_score[j][gi] += tmp_score[gi] * dw;
                    count[j][gi] += dw;
                    freq[j][gi]++;
                }
            }

        }

        for (j = 0; j < vecstrAllQuery.size(); j++) {

            int countSelected = 0;
            for (k = 0; k < selectedDataset[j].size(); k++)
                countSelected += selectedDataset[j][k];

            //int minRequired = (int) ((float) dab_list.size() * 0.50);
            int minRequired = (int) ((float) countSelected * threshold_g);

            if (DEBUG) {
                int nG = 0;
                for (k = 0; k < final_score[j].size(); k++) {
                    if (freq[j][k] >= minRequired) {
                        nG++;
                    }
                }
                fprintf(stderr, "Query %d numSelectedDataset %d numGenesIntegrated %d\n", j, countSelected, nG);
            }

            for (k = 0; k < final_score[j].size(); k++) {
                if (freq[j][k] < minRequired) {
                    final_score[j][k] = -320;
                    continue;
                }
                if (CMeta::IsNaN(final_score[j][k])) {
                    final_score[j][k] = -320;
                    continue;
                }
                final_score[j][k] /= count[j][k];
            }
            sprintf(acBuffer, "%s/%d.query", output_dir.c_str(), j);
            CSeekTools::WriteArrayText(acBuffer, vecstrAllQuery[j]);
            sprintf(acBuffer, "%s/%d.gscore", output_dir.c_str(), j);
            CSeekTools::WriteArray(acBuffer, final_score[j]);
            sprintf(acBuffer, "%s/%d.dweight", output_dir.c_str(), j);
            CSeekTools::WriteArray(acBuffer, dweight[j]);
        }
        return 0;
    }


    //DAB mode (sparse DAB: .2.dab)
    if (sArgs.dab_flag == 1) {
        float threshold_g = sArgs.threshold_g_arg;
        float threshold_q = sArgs.threshold_q_arg;
        bool DEBUG = !!sArgs.debug_flag;

        string search_mode = sArgs.tsearch_mode_arg;
        if (search_mode == "NA") {
            fprintf(stderr, "Please specify a search mode!\n");
            return 1;
        } else if (search_mode == "spell") {
            fprintf(stderr, "Note that spell is not supported for this mode!\n");
            return 1;
        }

        string norm_mode = sArgs.norm_mode_arg;
        if (sArgs.default_type_arg != 0 && sArgs.default_type_arg != 1) {
            fprintf(stderr, "Error, invalid type!\n");
            return -1;
        }

        if (norm_mode == "NA") {
            fprintf(stderr, "Error, please specify a norm mode!\n");
            return -1;
        }

        float exp = sArgs.exp_arg;
        if (norm_mode == "rank") {
            if (sArgs.max_rank_arg == -1) {
                fprintf(stderr, "Error, please supply the max rank flag.\n");
                return -1;
            }
            if (sArgs.rbp_p_arg == -1) {
                fprintf(stderr, "Error, please supply the rbp_p flag.\n");
                return -1;
            }
        } else if (norm_mode == "subtract_z") {
            if (exp == -1) {
                fprintf(stderr, "Invalid exponent!\n");
                return -1;
            }
            if (sArgs.rbp_p_arg == -1) {
                fprintf(stderr, "Error, please supply the rbp_p flag.\n");
                return -1;
            }
        }

        vector<float> score_cutoff;
        bool bDatasetCutoff = false;
        string dset_cutoff_file = sArgs.dset_cutoff_file_arg;
        if (dset_cutoff_file != "NA") {
            CSeekTools::ReadArray(dset_cutoff_file.c_str(), score_cutoff);
            bDatasetCutoff = true;
            fprintf(stderr, "Filtered mode is on!\n");
        }

        float rbp_p = sArgs.rbp_p_arg;
        int max_rank = sArgs.max_rank_arg;
        int num_iter = sArgs.num_iter_arg;
        vector <string> dab_list;
        CSeekTools::ReadListOneColumn(sArgs.dab_list_arg, dab_list);
        vector < CSeekIntIntMap * > dm;
        dm.resize(dab_list.size());

        //selected datasets for each query
        vector <vector<char>> selectedDataset;
        selectedDataset.resize(vecstrAllQuery.size());
        for (i = 0; i < vecstrAllQuery.size(); i++)
            CSeekTools::InitVector(selectedDataset[i], dab_list.size(), (char) 0);

        //MODE 1 - Normal search:
        //preparing query
        vector <vector<float>> q_weight;
        q_weight.resize(vecstrAllQuery.size());
        for (i = 0; i < vecstrAllQuery.size(); i++) {
            CSeekTools::InitVector(q_weight[i], numGenes, (float) 0);
            for (j = 0; j < vecstrAllQuery[i].size(); j++)
                q_weight[i][mapstriGenes[vecstrAllQuery[i][j]]] = 1;
        }

        //preparing query2
        vector <vector<unsigned int>> qq;
        qq.resize(vecstrAllQuery.size());
        for (i = 0; i < vecstrAllQuery.size(); i++) {
            qq[i] = vector<unsigned int>();
            for (j = 0; j < vecstrAllQuery[i].size(); j++)
                qq[i].push_back(mapstriGenes[vecstrAllQuery[i][j]]);
        }

        fprintf(stderr, "Reading sparse DAB\n");
        vector <vector<float>> final_score, count;
        vector <vector<float>> dweight;
        vector <vector<int>> freq;
        final_score.resize(vecstrAllQuery.size());
        count.resize(vecstrAllQuery.size());
        freq.resize(vecstrAllQuery.size());
        dweight.resize(vecstrAllQuery.size());
        for (j = 0; j < vecstrAllQuery.size(); j++) {
            CSeekTools::InitVector(final_score[j], vecstrGenes.size(), (float) CMeta::GetNaN());
            CSeekTools::InitVector(count[j], vecstrGenes.size(), (float) 0);
            CSeekTools::InitVector(freq[j], vecstrGenes.size(), (int) 0);
            CSeekTools::InitVector(dweight[j], dab_list.size(), (float) 0);
        }

        if (bDatasetCutoff) fprintf(stderr, "Dataset cutoff is on!\n");
        else fprintf(stderr, "Dataset cutoff is off!\n");


        for (i = 0; i < dab_list.size(); i++) {
            fprintf(stderr, "Reading %d: %s\n", i, dab_list[i].c_str());
            CSeekIntIntMap d1(vecstrGenes.size());
            string dabfile = dab_dir + "/" + dab_list[i];
            CSparseFlatMatrix<float> sm(0);

            if (norm_mode == "rank") {
                if (sArgs.default_type_arg == 0) //utype
                    CSeekWriter::ReadSeekSparseMatrix<utype>(dabfile.c_str(), sm, d1,
                                                             max_rank, rbp_p, vecstrGenes);
                else
                    CSeekWriter::ReadSeekSparseMatrix<unsigned short>(dabfile.c_str(),
                                                                      sm, d1, max_rank, rbp_p, vecstrGenes);
            } else if (norm_mode == "subtract_z") {
                if (sArgs.default_type_arg == 0) //utype
                    CSeekWriter::ReadSeekSparseMatrix<utype>(dabfile.c_str(), sm, d1,
                                                             vecstrGenes, 1000, exp);
                else
                    CSeekWriter::ReadSeekSparseMatrix<unsigned short>(dabfile.c_str(), sm, d1,
                                                                      vecstrGenes, 1000, exp);
            }
            const vector <utype> &allGenes = d1.GetAllReverse();

            for (j = 0; j < vecstrAllQuery.size(); j++) {
                int numPresent = 0;
                for (k = 0; k < qq[j].size(); k++) {
                    if (d1.GetForward(qq[j][k]) == (unsigned int) -1) continue;
                    numPresent++;
                }
                if (search_mode == "eq" && numPresent == 0) {
                    continue;
                } else if (search_mode == "cv_loi" && numPresent <= 1) {
                    continue;
                } else if (search_mode == "spell" && numPresent <= 1) {
                    continue;
                }
                int numThreshold = (int) (threshold_q * qq[j].size());
                if (numPresent > numThreshold)
                    selectedDataset[j][i] = 1;
            }

#pragma omp parallel for \
            shared(qu, sm, d1, dweight, final_score, count, freq, score_cutoff) \
            private(j, k) firstprivate(bDatasetCutoff, i) schedule(dynamic)
            for (j = 0; j < vecstrAllQuery.size(); j++) {
                //not enough query genes present
                if (selectedDataset[j][i] == 0) continue;

                float dw = 1.0;
                if (search_mode == "eq") {
                    dw = 1.0;
                } else if (search_mode == "cv_loi") {
                    //cv_weight_LOO(qu[j], sm, &d1, rbp_p, dw);
                    cv_weight(qu[j], sm, &d1, rbp_p, dw);
                } else if (search_mode == "spell") {
                    //spell_weight(qu[j], sm, &d1);
                }

                if (bDatasetCutoff) {
                    if (score_cutoff[i] > dw)
                        dw = 0;
                    //fprintf(stderr, "%.3e %.3e\n", score_cutoff[i], dw);
                }
                //fprintf(stderr, "%.3e\n", dw);
                dweight[j][i] = dw;
                //if(dw==0)
                //	continue;
                vector<float> tmp_score;
                get_score(tmp_score, sm, &d1, q_weight[j]);
                for (k = 0; k < d1.GetNumSet(); k++) {
                    utype gi = allGenes[k];
                    if (CMeta::IsNaN(final_score[j][gi]))
                        final_score[j][gi] = 0;
                    final_score[j][gi] += tmp_score[gi] * dw;
                    count[j][gi] += dw;
                    freq[j][gi]++;
                }
            }
        }

        /*
        ofstream ofsm;
        ofsm.open("/memex/qzhu/p4/concatenate_tumor_network.half", ios_base::binary);
        res.Save(ofsm, true);
        */

        for (j = 0; j < vecstrAllQuery.size(); j++) {
            int countSelected = 0;
            for (k = 0; k < selectedDataset[j].size(); k++)
                countSelected += selectedDataset[j][k];

            int minRequired = (int) ((float) countSelected * threshold_g);

            if (DEBUG) {
                int nG = 0;
                for (k = 0; k < final_score[j].size(); k++) {
                    if (freq[j][k] >= minRequired) {
                        nG++;
                    }
                }
                fprintf(stderr, "Query %d numSelectedDataset %d numGenesIntegrated %d\n", j, countSelected, nG);
            }

            for (k = 0; k < final_score[j].size(); k++) {
                if (freq[j][k] < minRequired) {
                    final_score[j][k] = -320;
                    continue;
                }
                if (CMeta::IsNaN(final_score[j][k])) {
                    final_score[j][k] = -320;
                    continue;
                }
                final_score[j][k] /= count[j][k];
            }

            sprintf(acBuffer, "%s/%d.query", output_dir.c_str(), j);
            CSeekTools::WriteArrayText(acBuffer, vecstrAllQuery[j]);
            sprintf(acBuffer, "%s/%d.gscore", output_dir.c_str(), j);
            CSeekTools::WriteArray(acBuffer, final_score[j]);
            sprintf(acBuffer, "%s/%d.dweight", output_dir.c_str(), j);
            CSeekTools::WriteArray(acBuffer, dweight[j]);
        }

        //MODE 2
        /*vector<vector<utype> > qu;
        qu.resize(vecstrAllQuery.size());
        for(i=0; i<vecstrAllQuery.size(); i++){
            qu[i] = vector<utype>();
            for(j=0; j<vecstrAllQuery[i].size(); j++)
                qu[i].push_back(mapstriGenes[vecstrAllQuery[i][j]]);
        }
        vector<vector<float> > q_weight;
        q_weight.resize(vecstrAllQuery.size());
        for(i=0; i<vecstrAllQuery.size(); i++){
            CSeekTools::InitVector(q_weight[i], numGenes, (float) 0);
            for(j=0; j<vecstrAllQuery[i].size(); j++)
                q_weight[i][mapstriGenes[vecstrAllQuery[i][j]]] = 1;
        }

        float alpha1 = 0.1; //inter - propagation
        float alpha2 = 0.1; //intra - propagation
        //float alpha = 0;
        bool label_propagate = false;

        for(i=0; i<dab_list.size(); i++){
            CSeekIntIntMap d1(vecstrGenes.size());
            vector<vector<float> > sm1;
            vector<utype> l1;
            string dabfile1 = dab_dir + "/" + dab_list[i];
            CSeekWriter::ReadSparseMatrixAsArray(l1, dabfile1.c_str());
            CSeekWriter::ReadSparseMatrix(l1, sm1, d1, 1000, 0.99, vecstrGenes);

            cerr << "1 " << dabfile1 << endl;
            vector<vector<float> > final_score;
            final_score.resize(vecstrAllQuery.size());

            //vector<vector<float> > tmp_score;
            //tmp_score.resize(vecstrAllQuery.size());

            vector<vector<float> > tmp_score_2;
            tmp_score_2.resize(vecstrAllQuery.size());

            //this score
            //component 1
            for(j=0; j<vecstrAllQuery.size(); j++){
                //get_score(tmp_score[j], sm1, &d1, q_weight[j]);
                init_score(tmp_score_2[j], &d1);
                init_score(final_score[j], &d1);
            }

            for(j=0; label_propagate && j<dab_list.size(); j++){
                if(i==j) continue;
                CSeekIntIntMap d2(vecstrGenes.size());
                vector<vector<float> > sm2;
                vector<utype> l2;
                string dabfile2 = dab_dir + "/" + dab_list[j];
                cerr << "2 " << dabfile2 << endl;
                CSeekWriter::ReadSparseMatrixAsArray(l2, dabfile2.c_str());
                CSeekWriter::ReadSparseMatrix(l2, sm2, d2, 1000, 0.99, vecstrGenes);

                //similarity matrix
                vector<vector<float> > sim;
                CSeekWriter::ProductNorm(sm1, sm2, d1, d2, sim);

                utype kj, kk;
                for(k=0; k<vecstrAllQuery.size(); k++){
                    //component 2
                    cerr << k << " of " << vecstrAllQuery.size() << endl;
                    vector<float> intra;
                    vector<float> inter;
                    get_score(intra, sm2, &d2, q_weight[k]);
                    get_score(inter, sim, &d2, intra);
                    //get_score(inter, sim, &d2, q_weight[k]);
                    add_score(inter, tmp_score_2[k], &d2);
                }
            }

            for(j=0; j<vecstrAllQuery.size(); j++){
                if(label_propagate){
                    vector<float> tmp3;
                    integrate(q_weight[j], tmp_score_2[j], tmp3, &d1, dab_list.size()-1, alpha1);
                    iterate(tmp3, final_score[j], sm1, &d1, alpha2, 1);
                }else{
                    //integrate(q_weight[j], tmp_score_2[j], tmp3, &d1, dab_list.size()-1, alpha1);
                    iterate(q_weight[j], final_score[j], sm1, &d1, alpha2, 1);
                }

                for(k=0; k<final_score[j].size(); k++){
                    if(CMeta::IsNaN(final_score[j][k])){
                        final_score[j][k] = -320;
                        continue;
                    }
                }

                for(k=0; k<qu[j].size(); k++){
                    final_score[j][qu[j][k]] = -320;
                }

                sprintf(acBuffer, "%s/%s/%d.query", output_dir.c_str(), dab_list[i].c_str(), j);
                CSeekTools::WriteArrayText(acBuffer, vecstrAllQuery[j]);
                sprintf(acBuffer, "%s/%s/%d.gscore", output_dir.c_str(), dab_list[i].c_str(), j);
                CSeekTools::WriteArray(acBuffer, final_score[j]);
            }
        }*/

    }

}

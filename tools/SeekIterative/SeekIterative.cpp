#include "stdafx.h"
#include "cmdline.h"

bool transfer(CDat &Dat, vector<vector<float> > &mat, CSeekIntIntMap *geneMap,
	vector<string> &vecstrGenes){

	int i, j;
	vector<utype> veciGenes;
	veciGenes.clear();
	veciGenes.resize(vecstrGenes.size());
	for( i = 0; i < vecstrGenes.size( ); ++i )
		veciGenes[ i ] = Dat.GetGene( vecstrGenes[i] );
	
	mat.resize(vecstrGenes.size());
	for(i=0; i<vecstrGenes.size(); i++)
		CSeekTools::InitVector(mat[i], vecstrGenes.size(), CMeta::GetNaN());

	for(i=0; i<vecstrGenes.size(); i++){
		utype s = veciGenes[i];
		if(CSeekTools::IsNaN(s)) continue;
		geneMap->Add(i);
		float *v = Dat.GetFullRow(s);
		for(j=0; j<vecstrGenes.size(); j++){
			utype t = veciGenes[j];
			if(CSeekTools::IsNaN(t)) continue;
			if(CMeta::IsNaN(v[t])) continue;
			mat[i][j] = Dat.Get(s, t);
		}
		free(v);
	}
	return true;
}

//Add scores from two vectors and integrate them using alpha
bool integrate(vector<float> &d1, vector<float> &d2, 
	vector<float> &dest, CSeekIntIntMap *geneMap, int k, float alpha){
	utype i;
	CSeekTools::InitVector(dest, geneMap->GetSize(), (float) CMeta::GetNaN());
	const vector<utype> &allGenes = geneMap->GetAllReverse();
	for(i=0; i<geneMap->GetNumSet(); i++){
		utype gi = allGenes[i];
		dest[gi] = (1.0 - alpha) * d1[gi] + alpha / k * d2[gi];
	}
	return true;
}

//initialize score vector
bool init_score(vector<float> &dest, CSeekIntIntMap *geneMap){
	CSeekTools::InitVector(dest, geneMap->GetSize(), (float) CMeta::GetNaN());
	utype i;	
	const vector<utype> &allGenes = geneMap->GetAllReverse();
	for(i=0; i<geneMap->GetNumSet(); i++)
		dest[allGenes[i]] = 0;
	return true;
}

bool add_score(vector<float> &src, vector<float> &dest, CSeekIntIntMap *geneMap){
	const vector<utype> &allGenes = geneMap->GetAllReverse();
	utype i, j;
	for(i=0; i<geneMap->GetNumSet(); i++){
		utype gi = allGenes[i];
		dest[gi] += src[gi];
	}
	return true;
}

bool search_one_dab(vector<float> &gene_score, 
	CDat &mat, const size_t numGenes,
	CSeekIntIntMap &d1,
	vector<utype> &indexConvReverse,
	vector<float> &q_weight){ 
	//q_weight is query presence - 1.0 if gene at the index i is a query gene

	vector<float> gene_count;
	CSeekTools::InitVector(gene_score, numGenes, (float)CMeta::GetNaN());
	CSeekTools::InitVector(gene_count, numGenes, (float)0);
	utype qqi;
	utype kk, k;

	const vector<utype> &allGenes = d1.GetAllReverse();	
	for(kk=0; kk<d1.GetNumSet(); kk++){
		utype k = allGenes[kk];
		utype gi = indexConvReverse[k];
		gene_score[gi] = 0;
	}

	for(qqi=0; qqi<d1.GetNumSet(); qqi++){
		utype qi = allGenes[qqi];
		utype qq = indexConvReverse[qi];
		if(q_weight[qq]==0) //not a query gene
			continue;
		//now a query gene
		float *vc = mat.GetFullRow(qi);
		for(kk=0; kk<d1.GetNumSet(); kk++){
			utype k = allGenes[kk];
			utype gi = indexConvReverse[k];
			float fl = vc[k];
			gene_score[gi] += fl;
			gene_count[gi] += 1.0;
		}
		delete[] vc;
	}

	for(kk=0; kk<d1.GetNumSet(); kk++){
		utype k = allGenes[kk];
		utype gi = indexConvReverse[k];
		gene_score[gi] /= gene_count[gi];
	}

	return true;
}

bool get_score(vector<float> &gene_score, 
	//vector<vector<float> > &mat, 
	CSparseFlatMatrix<float> &mat,
	CSeekIntIntMap *geneMap, vector<float> &q_weight){
	vector<float> gene_count;
	int numGenes = geneMap->GetSize();
	CSeekTools::InitVector(gene_score, numGenes, (float)CMeta::GetNaN());
	CSeekTools::InitVector(gene_count, numGenes, (float)0);
	int qi=0;
	const vector<utype> &allGenes = geneMap->GetAllReverse();
	utype kk, k;
	
	for(kk=0; kk<geneMap->GetNumSet(); kk++){
		utype gi = allGenes[kk];
		gene_score[gi] = 0;
	}
	for(qi=0; qi<geneMap->GetNumSet(); qi++){
		utype qq = allGenes[qi];
		if(q_weight[qq]==0) 
			continue;
		const vector<CPair<float> > &vc = mat.GetRow(qq);
		for(kk=0; kk<vc.size(); kk++){
			float fl = vc[kk].v;
			utype gi = vc[kk].i;
			gene_score[gi] += fl * q_weight[qq];
		}
		for(kk=0; kk<geneMap->GetNumSet(); kk++){
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

	for(k=0; k<geneMap->GetNumSet(); k++){
		utype gi = allGenes[k];
		gene_score[gi] /= gene_count[gi];
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

bool weight(vector<char> &is_query, vector<float> &d1,
	CSeekIntIntMap *geneMap, float rbp_p, float &w){
	vector<AResultFloat> ar;
	ar.resize(geneMap->GetSize());
	utype i;
	for(i=0; i<ar.size(); i++)
		ar[i].i = i;
	
	const vector<utype> &allGenes = geneMap->GetAllReverse();
	for(i=0; i<geneMap->GetNumSet(); i++){
		utype gi = allGenes[i];
		ar[gi].f = d1[gi];
	}
	
	int MAX = 5000;
	nth_element(ar.begin(), ar.begin()+MAX, ar.end());
	sort(ar.begin(), ar.begin()+MAX);
	
	w = 0;
	for(i=0; i<MAX; i++)
		if(is_query[ar[i].i]==1)
			w += pow(rbp_p, i);
	w *= (1.0 - rbp_p);	
	return true;
}

bool cv_weight(vector<utype> &query, 
	//vector<vector<float> > &mat,
	CSparseFlatMatrix<float> &mat,
	CSeekIntIntMap *geneMap, float rbp_p, float &tot_w){
	//leave one in cross-validation
	utype i, j;
	int numGenes = geneMap->GetSize();
	tot_w = 0;
	
	for(i=0; i<query.size(); i++){
		vector<char> is_query;
		vector<float> q_weight;
		CSeekTools::InitVector(is_query, numGenes, (char) 0);
		CSeekTools::InitVector(q_weight, numGenes, (float) 0);
		q_weight[query[i]] = 1.0;
		float w = 0;
		
		for(j=0; j<query.size(); j++){
			if(i==j) continue;
			is_query[query[j]] = 1;
		}
		vector<float> gene_score;
		get_score(gene_score, mat, geneMap, q_weight);
		gene_score[query[i]] = 0;
		weight(is_query, gene_score, geneMap, rbp_p, w);
		tot_w += w;
	}
	tot_w /= query.size();
	
	return true;	
}

int main(int iArgs, char **aszArgs){
	static const size_t c_iBuffer   = 1024;
	char acBuffer[c_iBuffer];
	
	gengetopt_args_info sArgs;
	ifstream ifsm;
	istream *pistm;
	vector<string> vecstrLine, vecstrGenes, vecstrDBs, vecstrQuery;
	utype i, qi, j, k, l, kk;
	
	if(cmdline_parser(iArgs, aszArgs, &sArgs)){
		cmdline_parser_print_help();
		return 1;
	}

	if( sArgs.input_arg ) {
		ifsm.open( sArgs.input_arg );
		pistm = &ifsm; 
	}else
		pistm = &cin;
	
	map<string, size_t> mapstriGenes;
	while( !pistm->eof( ) ) {
		pistm->getline( acBuffer, c_iBuffer - 1 );
		acBuffer[ c_iBuffer - 1 ] = 0;
		vecstrLine.clear( );
		CMeta::Tokenize( acBuffer, vecstrLine );
		if( vecstrLine.size( ) < 2 ) {
			//cerr << "Ignoring line: " << acBuffer << endl;
			continue;
		}
		if( !( i = atoi( vecstrLine[ 0 ].c_str( ) ) ) ) {
			cerr << "Illegal gene ID: " << vecstrLine[ 0 ] << " for "
				<< vecstrLine[ 1 ] << endl;
			return 1;
		}
		i--;
		if( vecstrGenes.size( ) <= i )
			vecstrGenes.resize( i + 1 );
		vecstrGenes[ i ] = vecstrLine[ 1 ];
		mapstriGenes[vecstrGenes[i]] = i;
	}

	fprintf(stderr, "Finished reading gene map\n");

	string dab_dir = sArgs.dab_dir_arg;
	string output_dir = sArgs.dir_out_arg;
	int numGenes = vecstrGenes.size();
	vector< vector<string> > vecstrAllQuery;
		
	fprintf(stderr, "Reading queries\n");
	if(!CSeekTools::ReadMultipleQueries(sArgs.query_arg, vecstrAllQuery))
		return -1;

	vector<vector<utype> > qu;
	qu.resize(vecstrAllQuery.size());
	for(i=0; i<vecstrAllQuery.size(); i++){
		qu[i] = vector<utype>();
		for(j=0; j<vecstrAllQuery[i].size(); j++)
			qu[i].push_back(mapstriGenes[vecstrAllQuery[i][j]]);
	}

	if(sArgs.combined_flag==1){
		string dab_base = sArgs.dab_basename_arg;
		string file1 = dab_dir + "/" + dab_base + ".dab";

		vector<vector<float> > q_weight;
		q_weight.resize(vecstrAllQuery.size());
		for(i=0; i<vecstrAllQuery.size(); i++){
			CSeekTools::InitVector(q_weight[i], numGenes, (float) 0);
			for(j=0; j<vecstrAllQuery[i].size(); j++)
				q_weight[i][mapstriGenes[vecstrAllQuery[i][j]]] = 1;
		}

		CDat CD;
		CD.Open(file1.c_str(), false, 2, false, false, false);

		CSeekIntIntMap d1(CD.GetGenes());
		vector<utype> indexConvReverse;
		CSeekTools::InitVector(indexConvReverse, vecstrGenes.size(), (utype) -1);	
		for(i=0; i<CD.GetGenes(); i++){
			map<string,size_t>::iterator it = mapstriGenes.find(CD.GetGene(i));
			if(it==mapstriGenes.end()) continue;
			indexConvReverse[i] = it->second;
			d1.Add(i);
		}

		vector<vector<float> > final_score;
		final_score.resize(vecstrAllQuery.size());
		for(j=0; j<vecstrAllQuery.size(); j++)
			CSeekTools::InitVector(final_score[j], vecstrGenes.size(), 
			(float)CMeta::GetNaN());

		const vector<utype> &allGenes = d1.GetAllReverse();			
		for(j=0; j<vecstrAllQuery.size(); j++){
			vector<float> tmp_score;
			search_one_dab(tmp_score, CD, vecstrGenes.size(), d1, indexConvReverse,
			q_weight[j]);
			for(kk=0; kk<d1.GetNumSet(); kk++){
				utype k = allGenes[kk];
				utype gi = indexConvReverse[k];
				final_score[j][gi] = tmp_score[gi];
			}
		}

		for(j=0; j<vecstrAllQuery.size(); j++){
			for(k=0; k<final_score[j].size(); k++){
				if(CMeta::IsNaN(final_score[j][k])){
					final_score[j][k] = -320;
					continue;
				}
			}
			sprintf(acBuffer, "%s/%d.query", output_dir.c_str(), j);
			CSeekTools::WriteArrayText(acBuffer, vecstrAllQuery[j]);
			sprintf(acBuffer, "%s/%d.gscore", output_dir.c_str(), j);
			CSeekTools::WriteArray(acBuffer, final_score[j]);
		}

		//Visualize
		for(j=0; j<vecstrAllQuery.size(); j++){
			vector<AResultFloat> ar;
			ar.resize(final_score[j].size());
			for(k=0; k<final_score[j].size(); k++){
				ar[k].i = k;
				ar[k].f = final_score[j][k];
			}
			sort(ar.begin(), ar.end());
			vector<utype> vec_g;
			vector<string> vec_s;
			int FIRST = 100;
			for(k=0; k<FIRST; k++){
				if(ar[k].f==-320) break;
				vec_g.push_back(CD.GetGeneIndex(vecstrGenes[ar[k].i]));
				vec_s.push_back(vecstrGenes[ar[k].i]);
			}
			if(vec_g.size()!=FIRST) continue;
			CDat V;
			V.Open(vec_s);
			for(k=0; k<vec_s.size(); k++){
				for(l=k+1; l<vec_s.size(); l++){
					V.Set(k, l, CD.Get(vec_g[k], vec_g[l]));
				}
			}

			sprintf(acBuffer, "%s/%d.dot", output_dir.c_str(), j);			
			ofstream ot(acBuffer);
			V.SaveDOT(ot, 0.0001, NULL, true, false, NULL, NULL);
		}
		

	}

	if(sArgs.test_flag==1){
		string dab_base = sArgs.dab_basename_arg;
		string file1 = dab_dir + "/" + dab_base + ".half";

		vector<vector<float> > q_weight;
		q_weight.resize(vecstrAllQuery.size());
		for(i=0; i<vecstrAllQuery.size(); i++){
			CSeekTools::InitVector(q_weight[i], numGenes, (float) 0);
			for(j=0; j<vecstrAllQuery[i].size(); j++)
				q_weight[i][mapstriGenes[vecstrAllQuery[i][j]]] = 1;
		}

		CHalfMatrix<float> res;
		res.Initialize(vecstrGenes.size());
		ifstream istm1;
		uint32_t dim;
		istm1.open(file1.c_str(), ios_base::binary);
		istm1.read((char*)(&dim), sizeof(dim));
		float *adScores = new float[dim-1];
		for(i=0; (i+1)<dim; i++){
			istm1.read((char*)adScores, sizeof(*adScores) * (dim - i - 1));
			res.Set(i, adScores);
		}
		delete[] adScores;
		istm1.close();
	
		string s1 = "uc003vor.2";
		vector<string> r1;
		r1.push_back("uc003vos.2");
		r1.push_back("uc003vop.1");
		r1.push_back("uc011jwz.1");
		r1.push_back("uc002rgw.1");

		for(i=0; i<r1.size(); i++){
			size_t i1 = mapstriGenes[s1];
			size_t j1 = mapstriGenes[r1[i]];
			fprintf(stderr, "%s %s %.5f\n", s1.c_str(), r1[i].c_str(), res.Get(i1, j1));
		}
	}

	if(sArgs.testcount_flag==1){
		string dab_base = sArgs.dab_basename_arg;
		string file1 = dab_dir + "/" + dab_base + ".pair_count";

		vector<vector<float> > q_weight;
		q_weight.resize(vecstrAllQuery.size());
		for(i=0; i<vecstrAllQuery.size(); i++){
			CSeekTools::InitVector(q_weight[i], numGenes, (float) 0);
			for(j=0; j<vecstrAllQuery[i].size(); j++)
				q_weight[i][mapstriGenes[vecstrAllQuery[i][j]]] = 1;
		}

		CHalfMatrix<unsigned short> res;
		res.Initialize(vecstrGenes.size());
		ifstream istm1;
		uint32_t dim;
		istm1.open(file1.c_str(), ios_base::binary);
		istm1.read((char*)(&dim), sizeof(dim));
		unsigned short *adScores = new unsigned short[dim-1];
		for(i=0; (i+1)<dim; i++){
			istm1.read((char*)adScores, sizeof(*adScores) * (dim - i - 1));
			res.Set(i, adScores);
		}
		delete[] adScores;
		istm1.close();
	
		string s1 = "uc003vor.2";
		vector<string> r1;
		r1.push_back("uc003vos.2");
		r1.push_back("uc003vop.1");
		r1.push_back("uc011jwz.1");
		r1.push_back("uc002rgw.1");

		for(i=0; i<r1.size(); i++){
			size_t i1 = mapstriGenes[s1];
			size_t j1 = mapstriGenes[r1[i]];
			fprintf(stderr, "%s %s %d\n", s1.c_str(), r1[i].c_str(), res.Get(i1, j1));
		}
	}

	if(sArgs.testcombined_flag==1){
		string dab_base = sArgs.dab_basename_arg;
		string file3 = dab_dir + "/" + dab_base + ".gene_count";
		string file2 = dab_dir + "/" + dab_base + ".pair_count";
		string file1 = dab_dir + "/" + dab_base + ".half";

		vector<utype> gene_count;
		CSeekTools::ReadArray(file3.c_str(), gene_count);

		float max_count = *max_element(gene_count.begin(), gene_count.end());
		float min_count_required = max_count*0.50;
		CSeekIntIntMap d1(vecstrGenes.size());
		vector<string> validGenes;
		for(i=0; i<vecstrGenes.size(); i++){
			if(gene_count[i]<min_count_required) continue;
			validGenes.push_back(vecstrGenes[i]);
			d1.Add(i);
		}
		
		CDat CD;
		CD.Open(validGenes);
		for(i=0; i<validGenes.size(); i++){
			for(j=i+1; j<validGenes.size(); j++){
				CD.Set(i,j,CMeta::GetNaN());
			}
		}
		//CHalfMatrix<float> res;
		//res.Initialize(vecstrGenes.size());
		ifstream istm1, istm2;
		uint32_t dim;
		istm1.open(file1.c_str(), ios_base::binary);
		istm1.read((char*)(&dim), sizeof(dim));
		istm2.open(file2.c_str(), ios_base::binary);
		istm2.read((char*)(&dim), sizeof(dim));
		float *adScores = new float[dim-1];
		unsigned short *adScores2 = new unsigned short[dim-1];
		for(i=0; (i+1)<dim; i++){
			istm1.read((char*)adScores, sizeof(*adScores) * (dim - i - 1));
			istm2.read((char*)adScores2, sizeof(*adScores2) * (dim - i - 1));
			if(i%2000==0)
				fprintf(stderr, "%d Finished\n", i);
			utype gi = d1.GetForward(i);
			if(CSeekTools::IsNaN(gi)) continue;

			for(j=0; j<dim-i-1; j++){
				utype gj = d1.GetForward(j+i+1);
				if(CSeekTools::IsNaN(gj) || 
					(float) adScores2[j]<min_count_required) continue;
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
		vector<string> r1;
		r1.push_back("uc003vos.2");
		r1.push_back("uc003vop.1");
		r1.push_back("uc011jwz.1");
		r1.push_back("uc002rgw.1");

		for(i=0; i<r1.size(); i++){
			//size_t i1 = mapstriGenes[s1];
			//size_t j1 = mapstriGenes[r1[i]];
			size_t i1 = CD.GetGeneIndex(s1);
			size_t j1 = CD.GetGeneIndex(r1[i]);
			fprintf(stderr, "%s %s %.5f\n", s1.c_str(), r1[i].c_str(), CD.Get(i1, j1));
		}
	}
	//DAB mode
	if(sArgs.dab_flag==1){
		if(sArgs.default_type_arg!=0 && sArgs.default_type_arg!=1){
			fprintf(stderr, "Error, invalid type!\n");
			return -1;
		}

		if(sArgs.max_rank_arg==-1){
			fprintf(stderr, "Error, please supply the max rank flag.\n");
			return -1;
		}

		if(sArgs.rbp_p_arg==-1){
			fprintf(stderr, "Error, please supply the rbp_p flag.\n");
			return -1;
		}

		float rbp_p = sArgs.rbp_p_arg;
		int max_rank = sArgs.max_rank_arg;

		int num_iter = sArgs.num_iter_arg;
		vector<string> dab_list;
		CSeekTools::ReadListOneColumn(sArgs.dab_list_arg, dab_list);
		vector<CSeekIntIntMap*> dm;
		dm.resize(dab_list.size());
		
		//MODE 1 - Normal search:
		vector<vector<float> > q_weight;
		q_weight.resize(vecstrAllQuery.size());
		for(i=0; i<vecstrAllQuery.size(); i++){
			CSeekTools::InitVector(q_weight[i], numGenes, (float) 0);
			for(j=0; j<vecstrAllQuery[i].size(); j++)
				q_weight[i][mapstriGenes[vecstrAllQuery[i][j]]] = 1;
		}
		
		
		fprintf(stderr, "Reading sparse DAB\n");
		vector<vector<float> > final_score, count;
		vector<vector<int> > freq;
		final_score.resize(vecstrAllQuery.size());
		count.resize(vecstrAllQuery.size());
		freq.resize(vecstrAllQuery.size());
		for(j=0; j<vecstrAllQuery.size(); j++){
			CSeekTools::InitVector(final_score[j], vecstrGenes.size(), (float)CMeta::GetNaN());
			CSeekTools::InitVector(count[j], vecstrGenes.size(), (float) 0);
			CSeekTools::InitVector(freq[j], vecstrGenes.size(), (int) 0);
		}
	
		//CSparseFlatHalfMatrix<float> res(0);
		/*CHalfMatrix<float> res;
		res.Initialize(vecstrGenes.size());
		for(i=0; i<vecstrGenes.size(); i++)
			for(j=i+1; j<vecstrGenes.size(); j++)
				res.Set(i, j, 0);
		*/
		for(i=0; i<dab_list.size(); i++){
			fprintf(stderr, "Reading %d: %s\n", i, dab_list[i].c_str());
			//vector<vector<float> > sm;
			//vector<map<utype,float> > sm;
			CSeekIntIntMap d1(vecstrGenes.size());
			string dabfile = dab_dir + "/" + dab_list[i];
			//CDat Dat;
			//Dat.Open(dabfile.c_str(), false, 2, false, false, false);
			//CSeekWriter::NormalizeDAB(Dat, vecstrGenes, false, false, false, true);
			//transfer(Dat, sm, &d1, vecstrGenes);
			
			CSparseFlatMatrix<float> sm (0);
			//CSeekWriter::ReadSparseMatrix(dabfile.c_str(), sm, d1, 3000, 0.997, vecstrGenes);

			if(sArgs.default_type_arg==0) //utype
				CSeekWriter::ReadSeekSparseMatrix<utype>(dabfile.c_str(), sm, d1, max_rank, rbp_p, vecstrGenes);
			else
				CSeekWriter::ReadSeekSparseMatrix<unsigned short>(dabfile.c_str(), sm, d1, max_rank, rbp_p, vecstrGenes);
			/*if(i==0){
				fprintf(stderr, "Summing...\n");
				res.Copy(sm);
				fprintf(stderr, "Finished Summing...\n");
			}else{*/
			//	fprintf(stderr, "Summing...\n");
			//	CSeekWriter::SumSparseMatrix(sm, res, d1, 1.0);
			//	fprintf(stderr, "Finished Summing...\n");
			/*}
			*/
			//fprintf(stderr, "Finished!\n");	
			
			const vector<utype> &allGenes = d1.GetAllReverse();			
			for(j=0; j<vecstrAllQuery.size(); j++){
				float dw = 1.0;
				//cv_weight(qu[j], sm, &d1, 0.99, dw);
				//fprintf(stderr, "%.3e\n", dw);
				vector<float> tmp_score;
				get_score(tmp_score, sm, &d1, q_weight[j]);
				for(k=0; k<d1.GetNumSet(); k++){
					utype gi = allGenes[k];
					if(CMeta::IsNaN(final_score[j][gi]))
						final_score[j][gi] = 0;
					final_score[j][gi] += tmp_score[gi] * dw;
					count[j][gi]+=dw;
					freq[j][gi]++;
				}	
			}
		}

		/*
		ofstream ofsm;
		ofsm.open("/memex/qzhu/p4/concatenate_tumor_network.half", ios_base::binary);
		res.Save(ofsm, true);
		*/
		
		int minRequired = (int) ((float) dab_list.size() * 0.50);
		for(j=0; j<vecstrAllQuery.size(); j++){
			for(k=0; k<final_score[j].size(); k++){
				if(freq[j][k]<minRequired){
					final_score[j][k] = -320;
					continue;
				}
				if(CMeta::IsNaN(final_score[j][k])){
					final_score[j][k] = -320;
					continue;
				}
				final_score[j][k] /= count[j][k];
			}
			sprintf(acBuffer, "%s/%d.query", output_dir.c_str(), j);
			CSeekTools::WriteArrayText(acBuffer, vecstrAllQuery[j]);
			sprintf(acBuffer, "%s/%d.gscore", output_dir.c_str(), j);
			CSeekTools::WriteArray(acBuffer, final_score[j]);
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

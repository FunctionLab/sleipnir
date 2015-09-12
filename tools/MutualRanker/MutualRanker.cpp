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

bool weight(vector<char> &is_query, vector<float> &d1,
CSeekIntIntMap *geneMap, float rbp_p, float &w){
	vector<AResultFloat> ar;
	ar.resize(geneMap->GetNumSet());
	utype i;
	const vector<utype> &allGenes = geneMap->GetAllReverse();
	for(i=0; i<geneMap->GetNumSet(); i++){
		utype gi = allGenes[i];
		ar[i].i = gi;
		ar[i].f = d1[gi];
	}
	int MAX = 1000;
	nth_element(ar.begin(), ar.begin()+MAX, ar.end(), AscendingFloat());
	sort(ar.begin(), ar.begin()+MAX, AscendingFloat());
	w = 0;
	for(i=0; i<MAX; i++)
		if(is_query[ar[i].i]==1)
			w += pow(rbp_p, i);
	w *= (1.0 - rbp_p);	
	return true;
}
bool get_score(vector<float> &gene_score, CFullMatrix<float> &mat,
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
		if(q_weight[qq]==0) //if this query gene does not exist
			continue;
		for(kk=0; kk<geneMap->GetNumSet(); kk++){
			utype gi = allGenes[kk];
			float fl = mat.Get(qq, gi);
			/*if(fl==99999){
				fprintf(stderr, "Bad %d %d \n", qi, kk);
			}*/
			gene_score[gi] += log(fl) * q_weight[qq];
			//gene_score[gi] += fl * q_weight[qq];
		}
		for(kk=0; kk<geneMap->GetNumSet(); kk++){
			utype gi = allGenes[kk];
			gene_count[gi] += q_weight[qq];
		}
	}
	for(k=0; k<geneMap->GetNumSet(); k++){
		utype gi = allGenes[k];
		gene_score[gi] /= gene_count[gi];
		if(gene_count[gi] == 0){
			gene_score[gi] = CMeta::GetNaN();
		}else{
			//fprintf(stderr, "Good %d %.2f\n", k, gene_score[gi]);
		}
	}
	return true;
}
//accepts fullmatrix, leave one in cross-validation
bool cv_weight(vector<utype> &query, CFullMatrix<float> &mat,
CSeekIntIntMap *geneMap, float rbp_p, float &tot_w){
	utype i, j;
	int numGenes = geneMap->GetSize();
	tot_w = 0;
	const vector<utype> &allRGenes = geneMap->GetAllReverse();
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
		for(j=0; j<geneMap->GetNumSet(); j++){
			utype gi = allRGenes[j];
			if(CMeta::IsNaN(gene_score[gi])){
				gene_score[gi] = 99999;
			}
		}
		gene_score[query[i]] = 99999;
		weight(is_query, gene_score, geneMap, rbp_p, w);
		tot_w += w;
	}
	tot_w /= query.size();
	return true;	
}

int main( int iArgs, char** aszArgs ) {
	static const size_t	c_iBuffer	= 1024;
#ifdef WIN32
	pthread_win32_process_attach_np( );
#endif // WIN32
	gengetopt_args_info	sArgs;
	ifstream			ifsm;
	istream*			pistm;
	vector<string>		vecstrLine, vecstrGenes, vecstrDatasets, vecstrUserDatasets;
	char				acBuffer[ c_iBuffer ];
	size_t				i, j, k;

	omp_set_num_threads(8);

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }

	vector<string> vecstrGeneID;
	map<string, size_t> mapstriGenes;
	if(!CSeekTools::ReadListTwoColumns(sArgs.input_arg, vecstrGeneID, vecstrGenes))
		return false;

	for(i=0; i<vecstrGenes.size(); i++)
		mapstriGenes[vecstrGenes[i]] = i;

	fprintf(stderr, "Finished reading gene map\n");

	bool convertRank = !!sArgs.convert_rank_flag;

	if(convertRank){
		string dab_file_name = sArgs.input_file_arg;
		CDataPair Dat;
		fprintf(stderr, "Opening file...\n");
		if(!Dat.Open(sArgs.input_file_arg, false, false, 2, false, false)){
			cerr << "error opening file" << endl;
			return 1;
		}

		vector<string> geneNames = Dat.GetGeneNames();

		sort(geneNames.begin(), geneNames.end());
		sort(vecstrGenes.begin(), vecstrGenes.end());
	
		vector<string> intersect(vecstrGenes.size() + geneNames.size());
		vector<string>::iterator it;
		it=set_intersection(geneNames.begin(), geneNames.end(), 
			vecstrGenes.begin(), vecstrGenes.end(), intersect.begin());
		intersect.resize(it - intersect.begin());
	
		int nG = intersect.size();
		float **mat = new float*[nG];
		mat[0] = new float[nG * nG];
		for(i=1; i<nG; i++){
			mat[i] = mat[i-1] + nG;
		}
		for(i=0; i<nG; i++){
			for(j=0; j<nG; j++){
				mat[i][j] = -999.0;
			}
		}

		unsigned int s, t;

		#pragma omp parallel for \
		shared(intersect, Dat, mat) \
		private(i, j, s, t) schedule(dynamic)
		for(i=0; i<intersect.size(); i++){	 //intersect genes
			if(i%100==0){
				fprintf(stderr, "Current progress %d of %d\n", i, intersect.size());
			}
			s = (unsigned int) Dat.GetGeneIndex(intersect[i]);
			//if(s==(unsigned int)-1) continue;
			vector<struct AResultFloat> af;
			af.resize(nG);
			for(j=0; j<intersect.size(); j++){
				t = (unsigned int) Dat.GetGeneIndex(intersect[j]);
				af[j].i = j;
				af[j].f = -999.0f;
				if(t==(unsigned int)-1) continue;
				af[j].f = Dat.Get(s, t);
			}
			sort(af.begin(), af.end());
			for(j=0; j<intersect.size(); j++){
				mat[i][af[j].i] = j + 1;
			}
		}

		CDat NewDat;
		NewDat.Open(intersect);
		#pragma omp parallel for \
		shared(intersect, NewDat, mat) \
		private(i, j, s, t) schedule(dynamic)
		for(i=0; i<intersect.size(); i++){	 //intersect genes
			s = (unsigned int) NewDat.GetGeneIndex(intersect[i]);
			for(j=i+1; j<intersect.size(); j++){
				t = (unsigned int) NewDat.GetGeneIndex(intersect[j]);
				float new_f = sqrt(mat[i][j] * mat[j][i]);
				NewDat.Set(s, t, new_f);
			}
		}
		NewDat.Save(sArgs.output_rank_file_arg);
		return 0;
	}


	string output_dir = sArgs.dir_out_arg;
	float threshold_g = sArgs.threshold_g_arg;
	float threshold_q = sArgs.threshold_q_arg;

	bool mergeComponent = !!sArgs.merge_component_flag;
	if(mergeComponent){
		vector<string> alist;
		CSeekTools::ReadListOneColumn(sArgs.merge_gscore_list_arg, alist);
		for(j=0; j<alist.size(); j++){
			vector<float> fs; //final score
			vector<float> ct; //count
			vector<int> fq; //freq
			vector<int> ct_selected; //count_selected
			vector<string> aQuery;

			int last_index = alist[j].find_last_of(".");
			string base = alist[j].substr(0, last_index);

			sprintf(acBuffer, "%s.query", base.c_str());
			CSeekTools::ReadMultiGeneOneLine(acBuffer, aQuery, 1024);
			sprintf(acBuffer, "%s.gscore_comp", base.c_str());
			CSeekTools::ReadArray(acBuffer, fs);
			sprintf(acBuffer, "%s.count_comp", base.c_str());
			CSeekTools::ReadArray(acBuffer, ct);
			sprintf(acBuffer, "%s.freq_comp", base.c_str());
			CSeekTools::ReadArray(acBuffer, fq);
			sprintf(acBuffer, "%s.count_selected", base.c_str());
			CSeekTools::ReadArray(acBuffer, ct_selected);

			int count_selected = ct_selected[0];
			int minRequired = (int) ((float) count_selected * threshold_g);
			minRequired = max(1, minRequired);

			for(k=0; k<fs.size(); k++){
				if(fq[k]<minRequired){
					fs[k] = 99999;
					continue;
				}
				if(CMeta::IsNaN(fs[k])){
					fs[k] = 99999;
					continue;
				}
				fs[k] /= ct[k];
			}
			vector<struct AResultFloat> af;
			af.resize(vecstrGenes.size());
			for(k=0; k<vecstrGenes.size(); k++){
				af[k].i = k;
				af[k].f = fs[k];
			}
			sort(af.begin(), af.end(), AscendingFloat()); //sort mutual ranks in ascending order
			int num_good = 0;
			for(k=0; k<vecstrGenes.size(); k++){
				if(af[k].f==99999) break;
			}
			num_good = k;
			for(k=0; k<num_good; k++){
				fs[af[k].i] = af[num_good - k - 1].f;
			}
			for(k=num_good; k<vecstrGenes.size(); k++){
				fs[af[k].i] = -320;
			}
			sprintf(acBuffer, "%s/%d.query", output_dir.c_str(), j);
			CSeekTools::WriteArrayText(acBuffer, aQuery);
			sprintf(acBuffer, "%s/%d.gscore", output_dir.c_str(), j);
			CSeekTools::WriteArray(acBuffer, fs);
		}
		return 0;
	}


	string search_mode = sArgs.search_mode_arg;
	//Doing search
	int numGenes = vecstrGenes.size();

	string dab_dir = sArgs.dab_dir_arg;
	vector<string> dab_list;
	CSeekTools::ReadListOneColumn(sArgs.tdab_list_arg, dab_list);

	vector< vector<string> > vecstrAllQuery;
	fprintf(stderr, "Reading queries\n");
	if(!CSeekTools::ReadMultipleQueries(sArgs.query_arg, vecstrAllQuery))
		return -1;
	fprintf(stderr, "Finished reading queries\n");

	//query in index
	vector<vector<utype> > qu;
	qu.resize(vecstrAllQuery.size());
	for(i=0; i<vecstrAllQuery.size(); i++){
		qu[i] = vector<utype>();
		for(j=0; j<vecstrAllQuery[i].size(); j++)
			qu[i].push_back(mapstriGenes[vecstrAllQuery[i][j]]);
	}

	//preparing query
	vector<vector<float> > q_weight;
	q_weight.resize(vecstrAllQuery.size());
	for(i=0; i<vecstrAllQuery.size(); i++){
		CSeekTools::InitVector(q_weight[i], numGenes, (float) 0);
		for(j=0; j<vecstrAllQuery[i].size(); j++){
			map<string, size_t>::iterator it = mapstriGenes.find(vecstrAllQuery[i][j]);
			if(it!=mapstriGenes.end()){
				q_weight[i][it->second] = 1;
			}
		}
	}

	//preparing query2
	vector<vector<unsigned int> > qq;
	qq.resize(vecstrAllQuery.size());
	for(i=0; i<vecstrAllQuery.size(); i++){
		qq[i] = vector<unsigned int>();
		for(j=0; j<vecstrAllQuery[i].size(); j++){
			map<string, size_t>::iterator it = mapstriGenes.find(vecstrAllQuery[i][j]);
			if(it!=mapstriGenes.end()){
				qq[i].push_back(it->second);
			}
		}
	}

	//selected datasets for each query
	vector<vector<char> > selectedDataset;
	selectedDataset.resize(vecstrAllQuery.size());
	for(i=0; i<vecstrAllQuery.size(); i++)
		CSeekTools::InitVector(selectedDataset[i], dab_list.size(), (char)0);

	vector<vector<float> > final_score, count;
	vector<vector<int> > freq;
	int nq = vecstrAllQuery.size();
	final_score.resize(nq);
	count.resize(nq);
	freq.resize(nq);
	for(j=0; j<vecstrAllQuery.size(); j++){
		CSeekTools::InitVector(final_score[j], vecstrGenes.size(), (float)CMeta::GetNaN());
		CSeekTools::InitVector(count[j], vecstrGenes.size(), (float) 0);
		CSeekTools::InitVector(freq[j], vecstrGenes.size(), (int) 0);
	}

	size_t l;
	for(l=0; l<dab_list.size(); l++){
		CDat Dat;
		fprintf(stderr, "Reading %d: %s\n", l, dab_list[l].c_str());
		string dabfile = dab_dir + "/" + dab_list[l];
		Dat.Open(dabfile.c_str(), false, 2, false, false, false);

		vector<unsigned int> veciGenes;
		veciGenes.resize(vecstrGenes.size());
		unsigned int ki;
		for(ki=0; ki<vecstrGenes.size(); ki++)
			veciGenes[ki] = (unsigned int) Dat.GetGeneIndex(vecstrGenes[ki]);
		unsigned int s,t;
		float d;
		CSeekIntIntMap m(vecstrGenes.size());
		for(i=0; i<vecstrGenes.size(); i++){
			if((s=veciGenes[i])==(unsigned int)-1) continue;
			m.Add(i);
		}

		//Copy to matrix sm
		CFullMatrix<float> sm;
		sm.Initialize(vecstrGenes.size(), vecstrGenes.size());
		const vector<utype> &allRGenes = m.GetAllReverse();
		for(i=0; i<m.GetNumSet(); i++){
			unsigned int si = allRGenes[i];
			s = veciGenes[si];
			for(j=i+1; j<m.GetNumSet(); j++){
				unsigned int tj = allRGenes[j];
				t = veciGenes[allRGenes[j]];
				if(CMeta::IsNaN(d = Dat.Get(s,t))){
					sm.Set(si, tj, 99999);
					sm.Set(tj, si, 99999);
				}else{
					sm.Set(si, tj, d);
					sm.Set(tj, si, d);
				}
			}
			sm.Set(si, si, 0);
		}

		fprintf(stderr, "Finished copying matrix\n");

		for(j=0; j<vecstrAllQuery.size(); j++){
			int numPresent = 0;
			for(k=0; k<qq[j].size(); k++){
				if(m.GetForward(qq[j][k])==(unsigned int)-1) continue;
				numPresent++;
			}
			if(search_mode=="cv_loi" && numPresent<=1){
				continue;
			}else if(search_mode=="eq" && numPresent==0){
				continue;
			}
			int minRequired = 1;
			if(search_mode=="cv_loi") 
				minRequired = 2;
			int numThreshold = (int) (threshold_q * qq[j].size());
			numThreshold = max(minRequired, numThreshold);
			if(numPresent>=numThreshold)
				selectedDataset[j][l] = 1;
		}
	
		#pragma omp parallel for \
		shared(qu, sm, final_score, count, freq, selectedDataset, q_weight, vecstrAllQuery, m) \
		private(j, k) firstprivate(l) schedule(dynamic)
		for(j=0; j<vecstrAllQuery.size(); j++){ //each query
			//not enough query genes present
			if(selectedDataset[j][l]==0) continue;
			float dw = 1.0;
			if(search_mode=="eq"){
				dw = 1.0;
			}else if(search_mode=="cv_loi"){ //cv_loi, rbp_p = 0.99
				cv_weight(qu[j], sm, &m, 0.99, dw); //returns weight to dw
			}
			vector<float> tmp_score;
			get_score(tmp_score, sm, &m, q_weight[j]);
			for(k=0; k<m.GetNumSet(); k++){ //each gene in the result
				utype gi = allRGenes[k];
				if(CMeta::IsNaN(tmp_score[gi])){
					//skip
					//should not come here!
				}else{
					if(CMeta::IsNaN(final_score[j][gi])){
						final_score[j][gi] = 0;
					}
					final_score[j][gi] += tmp_score[gi] * dw; //summed over all dataset
					count[j][gi]+=dw; //summed over all dataset
					freq[j][gi]++; //summed over all dataset
				}
			}
		}
	}

	bool outputComponent = !!sArgs.output_component_flag;
	if(outputComponent){
		for(j=0; j<vecstrAllQuery.size(); j++){
			int countSelected = 0;
			for(k=0; k<selectedDataset[j].size(); k++) //k is dataset iterator
				countSelected+=selectedDataset[j][k];
			vector<int> vecCountSelected;
			vecCountSelected.push_back(countSelected);
			sprintf(acBuffer, "%s/%d.query", output_dir.c_str(), j);
			CSeekTools::WriteArrayText(acBuffer, vecstrAllQuery[j]);
			sprintf(acBuffer, "%s/%d.gscore_comp", output_dir.c_str(), j);
			CSeekTools::WriteArray(acBuffer, final_score[j]);
			sprintf(acBuffer, "%s/%d.count_comp", output_dir.c_str(), j);
			CSeekTools::WriteArray(acBuffer, count[j]);
			sprintf(acBuffer, "%s/%d.freq_comp", output_dir.c_str(), j);
			CSeekTools::WriteArray(acBuffer, freq[j]);
			sprintf(acBuffer, "%s/%d.count_selected", output_dir.c_str(), j);
			CSeekTools::WriteArrayText(acBuffer, vecCountSelected);
		}
		return 0;
	}


	for(j=0; j<vecstrAllQuery.size(); j++){
		int countSelected = 0;
		for(k=0; k<selectedDataset[j].size(); k++) //k is dataset iterator
			countSelected+=selectedDataset[j][k];
		int minRequired = (int) ((float) countSelected * threshold_g);
		minRequired = max(1, minRequired);
		for(k=0; k<final_score[j].size(); k++){
			if(freq[j][k]<minRequired){
				final_score[j][k] = 99999;
				continue;
			}
			if(CMeta::IsNaN(final_score[j][k])){
				final_score[j][k] = 99999;
				continue;
			}
			final_score[j][k] /= count[j][k];
		}
		vector<struct AResultFloat> af;
		af.resize(vecstrGenes.size());
		for(k=0; k<vecstrGenes.size(); k++){
			af[k].i = k;
			af[k].f = final_score[j][k];
		}
		sort(af.begin(), af.end(), AscendingFloat()); //sort mutual ranks in ascending order
		int num_good = 0;
		for(k=0; k<vecstrGenes.size(); k++){
			if(af[k].f==99999) break;
		}
		num_good = k;
		for(k=0; k<num_good; k++){
			final_score[j][af[k].i] = af[num_good - k - 1].f;
			//fprintf(stderr, "%d %.2f\n", k, af[num_good - k - 1].f);
		}
		for(k=num_good; k<vecstrGenes.size(); k++){
			final_score[j][af[k].i] = -320;
		}
		sprintf(acBuffer, "%s/%d.query", output_dir.c_str(), j);
		CSeekTools::WriteArrayText(acBuffer, vecstrAllQuery[j]);
		sprintf(acBuffer, "%s/%d.gscore", output_dir.c_str(), j);
		CSeekTools::WriteArray(acBuffer, final_score[j]);
	}

#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; 
}

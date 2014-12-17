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
#include "seekwriter.h"
#include "seekreader.h"

namespace Sleipnir {

//mat a symmetric matrix
bool CSeekWriter::ReadSparseMatrix(const char *fileName,
	vector<map<utype,float> > &mat, CSeekIntIntMap &m, 
	const int maxRank, const float rbp_p,
	const vector<string> &vecstrGenes){

	FILE *f = fopen(fileName, "rb");
	if(f==NULL){
		cerr << "File not found" << endl;
		return false;
	}

	utype numGenes = 0;
	utype numPresent = 0;
	utype i, j;
	int ret;
	mat.clear();

	mat.resize(vecstrGenes.size());
	for(i=0; i<vecstrGenes.size(); i++)
		mat[i] = map<utype,float>();

	//m need to be initialized to size vecstrGenes.size() first!
	ret = fread((char*) (&numPresent), 1, sizeof(numPresent), f);
	for(j=0; j<numPresent; j++){
		utype val;
		ret = fread((char*)(&val), 1, sizeof(val), f);
		m.Add(val);
	}

	ret = fread((char*) (&numGenes), 1, sizeof(numGenes), f);

	vector<float> rbp_score;
	rbp_score.resize(maxRank);
	for(i=0; i<maxRank; i++)
		rbp_score[i] = (1.0 - rbp_p) * pow(rbp_p, (float)i);

	for(i=0; i<numGenes; i++){
		utype id, id2;
		unsigned short numEntries;
		unsigned short val;
		ret = fread((char*)(&id), 1, sizeof(id), f);
		ret = fread((char*)(&numEntries), 1, sizeof(numEntries), f);
		for(j=0; j<numEntries; j++){
			ret = fread((char*)(&id2),1,sizeof(id2),f);
			ret = fread((char*)(&val),1,sizeof(val),f);
			utype first = id;
			utype second = id2;
			if(first>=second){
				first = id2;
				second = id;
			}
			mat[first][second] = rbp_score[val];
		}
	}
	fclose(f);

	utype ii, jj;
	const vector<utype> &allRGenes = m.GetAllReverse();
	fprintf(stderr, "Begin calculating row sum\n");
	vector<float> vecSum;
	CSeekTools::InitVector(vecSum, vecstrGenes.size(), (float) 0);
	for(ii=0; ii<m.GetNumSet(); ii++){
		i = allRGenes[ii];
		map<utype,float>::iterator it;
		for(it=mat[i].begin(); it!=mat[i].end(); it++){
			j = it->first;
			float second = it->second;
			vecSum[i] += second;
			vecSum[j] += second;
		}
	}

	vector<float> vecSqrtSum;
	CSeekTools::InitVector(vecSqrtSum, vecstrGenes.size(), (float) 0);

	for(ii=0; ii<m.GetNumSet(); ii++){
		i = allRGenes[ii];
		if(vecSum[i]==0) continue;
		vecSqrtSum[i] = sqrtf(vecSum[i]);
	}

	fprintf(stderr, "Begin normalization using row sum\n");
	for(ii=0; ii<m.GetNumSet(); ii++){
		i = allRGenes[ii];
		map<utype,float>::iterator it;
		for(it=mat[i].begin(); it!=mat[i].end(); it++){
			j = it->first;
			if(vecSqrtSum[i]==0 || vecSqrtSum[j]==0) continue;
			it->second = it->second / vecSqrtSum[i] / vecSqrtSum[j];
			//symmetric matrix
			mat[j][i] = it->second;
		}
	}
	return true;
}

//add this matrix with weight w
bool CSeekWriter::SumSparseMatrix(CSparseFlatMatrix<float> &mat1,
	CSparseFlatHalfMatrix<float> &res, const CSeekIntIntMap &mi, const float w){
	utype i, j, ii, jj;
	//suppose res is already initialized
	const vector<utype> &allR = mi.GetAllReverse();
	for(ii=0; ii<mi.GetNumSet(); ii++){
		i = allR[ii];
		vector<CPair<float> >::iterator row_it;
		vector<CPair<float> > create;
		for(row_it = mat1.RowBegin(i); row_it!=mat1.RowEnd(i); row_it++){
			utype j = row_it->i;
			float rv = row_it->v;
			if(j<=i) continue; //only add if j is greater than i
			CPair<float> *pp = res.GetElement(i,j);
			if(pp==NULL){		
				CPair<float> cp;
				cp.i = j;
				cp.v = rv;
				create.push_back(cp);
				continue;
			}
			pp->v += rv * w;
		}
		for(jj=0; jj<create.size(); jj++)
			res.Add(i, create[jj].i, create[jj].v * w);
		if(create.size()>0)
			res.SortRow(i);
	}
	return true;
}

bool CSeekWriter::SumSparseMatrix(CSparseFlatMatrix<float> &mat1,
	CHalfMatrix<float> &res, const CSeekIntIntMap &mi, const float w){
	utype i, ii;
	//suppose res is already initialized
	const vector<utype> &allR = mi.GetAllReverse();
	for(ii=0; ii<mi.GetNumSet(); ii++){
		i = allR[ii];
		vector<CPair<float> >::iterator row_it;
		for(row_it = mat1.RowBegin(i); row_it!=mat1.RowEnd(i); row_it++){
			utype j = row_it->i;
			float rv = row_it->v;
			if(j<=i) continue; //only add if j is greater than i
			res.Set(i, j, res.Get(i, j) + rv * w);
		}
	}
	return true;
}
//Calculate the similarity of two distance matrices
//by simply taking product of two matrix for corresponding entries
bool CSeekWriter::ProductNorm(const vector<map<utype,float> > &mat1,
	const vector<map<utype,float> > &mat2, const CSeekIntIntMap &m1, 
	const CSeekIntIntMap &m2, vector<map<utype,float> > &re){

	utype ii, jj;
	utype i, j;

	re.resize(mat1.size());
	for(i=0; i<mat1.size(); i++)
		re[i] = map<utype,float>();

	const vector<utype> &allRGenes1 = m1.GetAllReverse();
	CSeekIntIntMap mi(mat1.size());
	for(ii=0; ii<m1.GetNumSet(); ii++){
		i = allRGenes1[ii];
		if(CSeekTools::IsNaN(m2.GetForward(i))) continue;
		mi.Add(i);
	}

	const vector<utype> &allR = mi.GetAllReverse();
	fprintf(stderr, "Begin calculating row sum\n");
	vector<float> vecSum;
	CSeekTools::InitVector(vecSum, mat1.size(), (float) 0);
	for(ii=0; ii<mi.GetNumSet(); ii++){
		i = allR[ii];
		map<utype,float>::const_iterator it;
		for(it=mat1[i].begin(); it!=mat1[i].end(); it++){
			j = it->first;
			float f1 = it->second;
			map<utype,float>::const_iterator it2;
			if((it2 = mat2[i].find(j))==mat2[i].end()) continue;
			float f2 = it2->second;
			re[i][j] = sqrtf(f1*f2);
			vecSum[i] += re[i][j];
			vecSum[j] += re[i][j];
		}
	}

	vector<float> vecSqrtSum;
	CSeekTools::InitVector(vecSqrtSum, mat1.size(), (float)0);
	for(ii=0; ii<mi.GetNumSet(); ii++){
		i = allR[ii];
		if(vecSum[i]==0) continue;
		vecSqrtSum[i] = sqrtf(vecSum[i]);
	}

	fprintf(stderr, "Begin normalization using row sum\n");
	for(ii=0; ii<mi.GetNumSet(); ii++){
		i = allR[ii];
		map<utype,float>::iterator it;
		for(it=re[i].begin(); it!=re[i].end(); it++){
			j = it->first;
			if(vecSqrtSum[i]==0 || vecSqrtSum[j]==0) continue;
			it->second = it->second / vecSqrtSum[i] / vecSqrtSum[j];
		}
	}
	return true;
}

bool CSeekWriter::TopologicalOverlap(CDataPair &Dat,
const vector<string> &vecstrGenes){
	size_t i, j;
	vector<unsigned int> veciGenes;
	veciGenes.resize(vecstrGenes.size());
	for(i=0; i<vecstrGenes.size(); i++)
		veciGenes[i] = (unsigned int) Dat.GetGeneIndex(vecstrGenes[i]);

	unsigned int s,t;
	float d;
	CSeekIntIntMap m(vecstrGenes.size());
	for(i=0; i<vecstrGenes.size(); i++){
		if((s=veciGenes[i])==(unsigned int)-1) continue;
		m.Add(i);
	}

	size_t trueSize = m.GetNumSet();
	vector<float> inner(trueSize);
	vector<vector<float> > fs(trueSize, inner);
	
	//float* fs = new float[trueSize*trueSize];
	const vector<utype> &allRGenes = m.GetAllReverse();

	for(i=0; i<m.GetNumSet(); i++){	
		s = veciGenes[allRGenes[i]];
		for(j=i+1; j<m.GetNumSet(); j++){
			t = veciGenes[allRGenes[j]];
			if(CMeta::IsNaN(d = Dat.Get(s,t))){
				fs[i][j] = 0;
				fs[j][i] = 0;
				fprintf(stderr, "Warning i, j is NaN, set to 0!\n", i, j);
			}else{
				fs[i][j] = d;
				fs[j][i] = d;
			}
		}
		fs[i][i] = 0;
	}

	//first step: transform z-scores back to correlation 
	//(exp(z*2) = (1+r)/(1-r) -> exp(2z) - exp(2z)*r = 1+r -> exp(2z) - 1 = exp(2z)*r + r = (exp(2z) + 1)*r->
	//(exp(2z) - 1) / (exp(2z) + 1) = r
	//second step: transform by abs, then take it to the exponent 9
	for(i=0; i<m.GetNumSet(); i++){	
		for(j=i+1; j<m.GetNumSet(); j++){
			if((d = fs[i][j])==0) continue;
			d = (expf(2.0*d) - 1.0) / (expf(2.0*d) + 1.0);
			d = pow(abs(d), 9);
			fs[i][j] = d;
			fs[j][i] = d;
			//fprintf(stderr, "%.3e\n", d);
		}
	}

	fprintf(stderr, "Finished step 1: tranform z-score back to pearson\n");
	vector<float> vecSum;
	CSeekTools::InitVector(vecSum, trueSize, CMeta::GetNaN());
	for(i=0; i<m.GetNumSet(); i++)
		vecSum[i] = 0;

	for(i=0; i<m.GetNumSet(); i++){	
		for(j=i+1; j<m.GetNumSet(); j++){
			vecSum[i] += fs[i][j];
			vecSum[j] += fs[i][j];
		}
	}

	//duplicate of fs
	vector<vector<float> > fs2(trueSize, inner);
	for(i=0; i<m.GetNumSet(); i++){	
		for(j=i+1; j<m.GetNumSet(); j++){
			fs2[i][j] = fs[i][j];
			fs2[j][i] = fs[i][j];
		}
		fs2[i][i] = 0;
	}


	//temporary storage matrix
	CHalfMatrix<float> res;
	res.Initialize(trueSize);
	for(i=0; i<m.GetNumSet(); i++){	
		for(j=i+1; j<m.GetNumSet(); j++){
			res.Set(i, j, 0);
		}
	}

	//result of multiplication
	fprintf(stderr, "Begin!\n");
	vector<vector<float> > fs_result(trueSize, inner);
	CStrassen::strassen(fs, fs2, fs_result, trueSize);
	fprintf(stderr, "Done!\n");

	size_t k;
	unsigned int u;
	//size_t ii = 0;
	for(i=0; i<m.GetNumSet(); i++){	
		for(j=i+1; j<m.GetNumSet(); j++){
			float tsum = fs_result[i][j];
			/*float *pi = &fs[i*trueSize];
			float *pj = &fs[j*trueSize];
			for(k=0; k<m.GetNumSet(); k++){
				tsum += pi[k] * pj[k];
			}*/
			tsum -= fs[i][i] * fs[j][i];
			tsum -= fs[i][j] * fs[j][j];
			float tmin = 0;
			if(vecSum[i] < vecSum[j])
				tmin = vecSum[i];
			else
				tmin = vecSum[j];
			float to = (tsum + d) / (tmin + 1.0 - d);
			res.Set(i, j, (float) to); //temporary matrix to store the results
			//fprintf(stderr, "%.3e\n", to);	
		}
		if(i%100==0){
			fprintf(stderr, "Doing topological overlap calculation, current %d\n", i);
		}
	}

	fprintf(stderr, "Finished step 2: topological overlap calculation\n");

	for(i=0; i<m.GetNumSet(); i++){	
		s = veciGenes[allRGenes[i]];
		for(j=i+1; j<m.GetNumSet(); j++){
			t = veciGenes[allRGenes[j]];
			Dat.Set(s, t, res.Get(i, j));
		}
		Dat.Set(s, s, 1.0);
	}

}


bool CSeekWriter::NormalizeDAB(CDataPair &Dat,
const vector<string> &vecstrGenes, 
//bool cutoff, float cutoff_val,
bool expTransform, bool divideNorm, bool subtractNorm){
	//default cutoff_val is 0

	size_t i, j;
	vector<unsigned int> veciGenes;
	veciGenes.clear();
	veciGenes.resize(vecstrGenes.size());
	for(i=0; i<vecstrGenes.size(); i++)
		veciGenes[i] = (unsigned int) Dat.GetGeneIndex(vecstrGenes[i]);

	vector<float> vecSum;
	vector<int> vecNum;
	CSeekTools::InitVector(vecSum, vecstrGenes.size(), CMeta::GetNaN());
	CSeekTools::InitVector(vecNum, vecstrGenes.size(), (int)-9999);

	unsigned int s,t;
	for(i=0; i<vecstrGenes.size(); i++){
		if((s=veciGenes[i])==(unsigned int)-1) continue;
		vecSum[i] = 0;
		vecNum[i] = 0;
	}

	if(divideNorm && subtractNorm){
		fprintf(stderr, "Error: both divideNorm and subtractNorm are true\n");
		return false;
	}else if(!divideNorm && !subtractNorm){
		fprintf(stderr, "Error: both divideNorm and subtractNorm are false\n");
		return false;
	}

	float d = -1;
	float r = -1;
	for(i=0; i<vecstrGenes.size(); i++){
		if((s=veciGenes[i])==(unsigned int)-1) continue;
		for(j=i+1; j<vecstrGenes.size(); j++){
			if((t=veciGenes[j])==(unsigned int)-1) continue;
			if(CMeta::IsNaN(d = Dat.Get(s,t))) continue;
			/*if(cutoff){
				if(d>cutoff_val){
					if(expTransform)
						r = expf(-1.0*d*d/2.0);
					else
						r = d;
					vecSum[i] += r;
					vecSum[j] += r;
					vecNum[i]++;
					vecNum[j]++;
				}
			}
			else{*/
				//fprintf(stderr, "Warning: Negative Z-Scores");
				if(expTransform)
					r = expf(-1.0*d*d/2.0);
				else
					r = d;
				vecSum[i] += r;
				vecSum[j] += r;
				vecNum[i]++;
				vecNum[j]++;
			//}	
		}
	}

	for(i=0; i<vecstrGenes.size(); i++){
		if((s=veciGenes[i])==(unsigned int)-1) continue;
		for(j=i+1; j<vecstrGenes.size(); j++){
			if((t=veciGenes[j])==(unsigned int)-1) continue;
			if(CMeta::IsNaN(d = Dat.Get(s,t))) continue;
			/*if(cutoff){
				if(d>cutoff_val){
					if(expTransform){
						if(divideNorm)
							r=expf(-1.0*d*d/2.0)/sqrtf(vecSum[i])/sqrtf(vecSum[j]);
						else if(subtractNorm)
							r=expf(-1.0*d*d/2.0)-vecSum[i]/vecNum[i]-vecSum[j]/vecNum[j];
					}else{
						if(divideNorm)
							r=d/sqrtf(vecSum[i])/sqrtf(vecSum[j]);
						else if(subtractNorm)
							r=d-vecSum[i]/vecNum[i]-vecSum[j]/vecNum[j];
					}
				}else{
					r=0; //default value
				}
				Dat.Set(s, t, r);
			}
			else{*/
				if(expTransform){
					if(divideNorm)
						r=expf(-1.0*d*d/2.0)/sqrtf(vecSum[i])/sqrtf(vecSum[j]);
					else if(subtractNorm)
						r=expf(-1.0*d*d/2.0)-vecSum[i]/vecNum[i]-vecSum[j]/vecNum[j];
				}else{
					if(divideNorm){
						//DANGEROUS
						if(vecSum[i]<=0){
							fprintf(stderr, "Warning, divide sqrt(z), when z<=0\n");
							r=0; //default value
						}else
							r=d/sqrtf(vecSum[i])/sqrtf(vecSum[j]);
					}else if(subtractNorm){
						r=d-0.5*(vecSum[i]/vecNum[i]+vecSum[j]/vecNum[j]);
					}
				}
				Dat.Set(s, t, r);
			//}
		}
	}

	//Plot a distribution
	/*vector<unsigned long> bins;
	bins.resize(55);
	float upper = 5.0; //assume z scores
	float lower = -5.0;
	float bin_size = (upper - lower) / 50;
	for(i=0; i<55; i++)
		bins[i] = 0;
	for(i=0; i<Dat.GetGenes(); i++){
		for(j=i+1; j<Dat.GetGenes(); j++){
			d = Dat.Get(i,j);
			if(CMeta::IsNaN(d)) continue;
			int b = (int) ((d - lower) / bin_size);
			if(b<0){
				bins[0]++;
				continue;
			}
			if(b>=55){
				bins[54]++;
				continue;
			}
			bins[b]++;
		}
	}
	fprintf(stderr, 
	"Distances: bin size: %.5f, num of bins: %d, min bin val: %.5f, max bin val: %.5f\n",
	bin_size, 55, lower, upper);
	for(i=0; i<55; i++){
		fprintf(stderr, "%lu\t%lu\n", i, bins[i]);
	}
	*/
	return true;
}

bool CSeekWriter::GetGeneAverage(CDataPair &Dat,
	const vector<string> &vecstrGenes,
	vector<float> &vecResult, bool logit, float top_percent){

	/* assume datapair is already opened */
	utype i, j;
	vector<utype> veciGenes;
	veciGenes.clear();
	veciGenes.resize(vecstrGenes.size());
	for( i = 0; i < vecstrGenes.size( ); ++i )
		veciGenes[ i ] = Dat.GetGene( vecstrGenes[i] );

	CSeekTools::InitVector(vecResult, vecstrGenes.size(), CMeta::GetNaN());
	for(i=0; i<vecstrGenes.size(); i++){
		utype s = veciGenes[i];
		if(CSeekTools::IsNaN(s)) continue;
		float *v = Dat.GetFullRow(s);
		float sum = 0;
		utype num = 0;
		vector<float> all;
		for(j=0; j<vecstrGenes.size(); j++){
			utype t = veciGenes[j];
			if(CSeekTools::IsNaN(t)) continue;
			if(CMeta::IsNaN(v[t])) continue;
			if(logit){
				//sum+=log(v[t]) - log((float) (1.0-v[t]));
				all.push_back(log(v[t]) - log((float) (1.0-v[t])));
			}else{
				//sum+=v[t];
				all.push_back(v[t]);
			}
			//num++;
		}
		sort(all.begin(), all.end());
		int top_start = (int) (((float)1.0 - top_percent)*(float)all.size());

		if(top_start<0){
			top_start = 0;
		}
		for(j=top_start; j<all.size(); j++){
			sum+=all[j];
			num++;
		}
		vecResult[i] = sum / (float) num;
		//fprintf(stderr, "%.2f\n", vecResult[i]);
		free(v);
	}
	return true;
}

bool CSeekWriter::GetGenePresence(CDataPair &Dat,
	const vector<string> &vecstrGenes,
	vector<char> &vecResult){
	/* assume datapair is already opened */
	utype i, j;
	vector<utype> veciGenes;
	veciGenes.clear();
	veciGenes.resize(vecstrGenes.size());
	for( i = 0; i < vecstrGenes.size( ); ++i )
		veciGenes[ i ] = Dat.GetGene( vecstrGenes[i] );

	CSeekTools::InitVector(vecResult, vecstrGenes.size(), (char) 0);

	for(i=0; i<vecstrGenes.size(); i++){
		if(CSeekTools::IsNaN(veciGenes[i])) continue;
		vecResult[i]=1;
	}
	return true;
}

bool CSeekWriter::GetDatasetSinfo(CDataPair &Dat,
	float &mean, float &stdev){
	utype i, j;
	mean = CMeta::GetNaN();
	stdev = CMeta::GetNaN();

	utype iGenes = Dat.GetGenes();

	unsigned int num = 0;
	float sum = 0;

	for(i=0; i<iGenes; i++){
		utype s = i;
		float *v = Dat.GetFullRow(s);
		for(j=0; j<iGenes; j++){
			if(CMeta::IsNaN(v[j])) continue;
			sum+=v[j];
			num++;
		}
	}

	if(num==0) return true;

	mean = sum / (float) num;
	float diff = 0;
	for(i=0; i<iGenes; i++){
		utype s = i;
		float *v = Dat.GetFullRow(s);
		for(j=0; j<iGenes; j++){
			if(CMeta::IsNaN(v[j])) continue;
			diff += (v[j] - mean) * (v[j] - mean);
		}
	}
	diff /= (float) num;
	stdev = sqrt(diff);

	return true;
}

}

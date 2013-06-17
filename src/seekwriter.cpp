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

bool CSeekWriter::ReadSparseMatrix(const char *fileName,
	vector<vector<float> > &mat, CSeekIntIntMap &m, 
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
		CSeekTools::InitVector(mat[i], vecstrGenes.size(), (float) 0);

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
		rbp_score[i] = (1.0 - rbp_p) * pow(rbp_p, i);

	fprintf(stderr, "Begin assigning rbp scores\n");

	for(i=0; i<numGenes; i++){
		utype id, id2;
		unsigned short numEntries;
		unsigned short val;
		ret = fread((char*)(&id), 1, sizeof(id), f);
		ret = fread((char*)(&numEntries), 1, sizeof(numEntries), f);
		for(j=0; j<numEntries; j++){
			ret = fread((char*)(&id2),1,sizeof(id2),f);
			ret = fread((char*)(&val),1,sizeof(val),f);
			mat[id][id2] = rbp_score[val];
		}
	}
	fclose(f);

	fprintf(stderr, "Filling zero rbp score\n");
	/*for(i=0; i<vecstrGenes.size(); i++){
		if(CSeekTools::IsNaN(m.GetForward(i))){
			for(j=0; j<vecstrGenes.size(); j++){
				mat[i][j] = CMeta::GetNaN();
			}
			continue;
		}
		for(j=0; j<vecstrGenes.size(); j++){
			if(CSeekTools::IsNaN(m.GetForward(j))){
				mat[i][j] = CMeta::GetNaN();
			}
		}
	}*/

	/*
	for(ii=0; ii<m.GetNumSet(); ii++){
		i = allRGenes[ii];
		for(jj=ii+1; jj<m.GetNumSet(); jj++){
			j = allRGenes[jj];
			if(CMeta::IsNaN(mat[i][j])){
				mat[i][j] = 0;
				mat[j][i] = 0;
			}
		}
	}*/

	utype ii, jj;
	const vector<utype> &allRGenes = m.GetAllReverse();
	fprintf(stderr, "Begin calculating row sum\n");
	vector<float> vecSum;
	CSeekTools::InitVector(vecSum, vecstrGenes.size(), (float) 0);
	for(ii=0; ii<m.GetNumSet(); ii++){
		i = allRGenes[ii];
		for(jj=ii+1; jj<m.GetNumSet(); jj++){
			j = allRGenes[jj];
			//if(CMeta::IsNaN(mat[i][j])) continue;
			//if(mat[i][j] < 0 || mat[i][j]>1){
			//	fprintf(stderr, "Should not happen, error!\n");
			//}
			if(mat[i][j]==0) continue;
			vecSum[i] += mat[i][j];
			vecSum[j] += mat[i][j];
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
		for(jj=ii+1; jj<m.GetNumSet(); jj++){
			j = allRGenes[jj];
			if(mat[i][j]==0 || vecSqrtSum[i]==0 || vecSqrtSum[j]==0) continue;
			mat[i][j] = mat[i][j] / vecSqrtSum[i] / vecSqrtSum[j];
			mat[j][i] = mat[i][j];
		}
	}
	return true;
}

bool CSeekWriter::ProductNorm(const vector<vector<float> > &mat1,
	const vector<vector<float> > &mat2, const CSeekIntIntMap &m1, 
	const CSeekIntIntMap &m2, vector<vector<float> > &re){

	utype ii, jj;
	utype i, j;

	re.resize(mat1.size());
	for(i=0; i<mat1.size(); i++)
		CSeekTools::InitVector(re[i], mat1.size(), (float)0);

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
		for(jj=ii+1; jj<mi.GetNumSet(); jj++){
			j = allR[jj];
			if(mat1[i][j]==0 || mat2[i][j]==0) continue;
			re[i][j] = sqrtf(mat1[i][j] * mat2[i][j]);
			re[j][i] = re[i][j];
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

	utype numNonZero = 0;
	fprintf(stderr, "Begin normalization using row sum\n");
	vector<float> vf;
	for(ii=0; ii<mi.GetNumSet(); ii++){
		i = allR[ii];
		for(jj=ii+1; jj<mi.GetNumSet(); jj++){
			j = allR[jj];
			if(mat1[i][j]==0 || mat2[i][j]==0) continue;
			if(vecSqrtSum[i]==0 || vecSqrtSum[j]==0) continue;
			re[i][j] = re[i][j] / vecSqrtSum[i] / vecSqrtSum[j];
			re[j][i] = re[i][j];
			//vf.push_back(rx);
			numNonZero++;
			//fprintf(stderr, "%.3e\n", re[i][j]);
		}
	}
	//sort(vf.begin(), vf.end(), greater<float>());
	//int xi;
	//for(xi=0; xi<vf.size(); xi++){
		//if(isinf(vf[xi]) || isnan(vf[xi])){
	//		fprintf(stderr, "%.3e\n", vf[xi]);
		//}
	//}
	//fprintf(stderr, "Non Zero: %d\n", numNonZero);
	return true;
}

bool CSeekWriter::WriteSparseMatrix(vector<vector<unsigned short> > &umat,
	int maxRank, const vector<string> &vecstrGenes, const char *fileName){

	FILE *f = fopen(fileName, "wb");
	if(f==NULL){
		cerr << "File not found!" << endl;
		return false;
	}
	utype numGenes = 0;
	utype i, j;

	CSeekIntIntMap mm(vecstrGenes.size());
	for(i=0; i<vecstrGenes.size(); i++){
		for(j=0; j<vecstrGenes.size(); j++)
			if(!CSeekTools::IsNaN(umat[i][j])) break;
		if(j!=vecstrGenes.size()){
			mm.Add(i);
		}
	}

	utype numPresent = mm.GetNumSet();
	//1 utype
	fwrite((char*) (&numPresent), 1, sizeof(numPresent), f);
	const vector<utype> &allR = mm.GetAllReverse();
	//numPresent utype
	for(i=0; i<numPresent; i++)
		fwrite((char*) (&allR[i]), 1, sizeof(allR[i]), f);

	for(i=0; i<vecstrGenes.size(); i++){
		for(j=i+1; j<vecstrGenes.size(); j++)
			if(!CSeekTools::IsNaN(umat[i][j]) && umat[i][j]!=maxRank) 
				break;
		if(j==vecstrGenes.size()) 
			continue;
		numGenes++;
	}

	//1 utype
	fwrite((char*) (&numGenes), 1, sizeof(numGenes), f);

	for(i=0; i<vecstrGenes.size(); i++){
		unsigned short numEntries = 0; //should be 1000
		for(j=i+1; j<vecstrGenes.size(); j++){
			if(CSeekTools::IsNaN(umat[i][j]) || umat[i][j]==maxRank)
				continue;
			numEntries++;
		}
		if(numEntries==0) 
			continue;
		//1 utype
		fwrite((char*) (&i), 1, sizeof(i), f);
		//1 unsigned short
		fwrite((char*) (&numEntries), 1, sizeof(numEntries), f);
		for(j=i+1; j<vecstrGenes.size(); j++){
			if(CSeekTools::IsNaN(umat[i][j]) || umat[i][j]==maxRank)
				continue;
			//1 utype
			fwrite((char*) (&j), 1, sizeof(j), f);
			//1 unsigned short
			fwrite((char*) (&umat[i][j]), 1, sizeof(umat[i][j]), f);
		}
	}

	fclose(f);
	return true;
}

bool CSeekWriter::GetSparseRankMatrix(CDat &Dat,
	vector<vector<unsigned short> > &umat, const unsigned short nullValue,
	int maxRank, //1000
	const vector<string> &vecstrGenes){

	utype i, j;
	vector<utype> veciGenes;
	veciGenes.clear();
	veciGenes.resize(vecstrGenes.size());
	for( i = 0; i < vecstrGenes.size( ); ++i )
		veciGenes[ i ] = Dat.GetGene( vecstrGenes[i] );
	umat.resize(vecstrGenes.size());

	for(i=0; i<vecstrGenes.size(); i++){
		CSeekTools::InitVector(umat[i], vecstrGenes.size(), nullValue);
		utype s = veciGenes[i];
		if(CSeekTools::IsNaN(s)) continue;

		float *v = Dat.GetFullRow(s);
		vector<AResultFloat> vv;
		vv.resize(vecstrGenes.size());
		for(j=0; j<vecstrGenes.size(); j++){
			utype t = veciGenes[j];
			vv[j].i = j;
			if(CSeekTools::IsNaN(t) || CMeta::IsNaN(v[t])){
				vv[j].f = -9999;
				continue;
			}
			vv[j].f = v[t];
		}
		nth_element(vv.begin(), vv.begin()+maxRank, vv.end());
		sort(vv.begin(), vv.begin()+maxRank);
		for(j=0; j<vecstrGenes.size(); j++){
			if(j<maxRank){
				umat[i][vv[j].i] = j;
			}else if(vv[j].f!=-9999){
				umat[i][vv[j].i] = maxRank;
			}
		}
		free(v);
	}
	for(i=0; i<vecstrGenes.size(); i++){
		utype s = veciGenes[i];
		if(CSeekTools::IsNaN(s)) continue;
		for(j=i+1; j<vecstrGenes.size(); j++){
			utype t = veciGenes[j];
			if(CSeekTools::IsNaN(t)) continue;
			umat[i][j] = std::min(umat[i][j], umat[j][i]);
			umat[j][i] = umat[i][j];
		}
	}
	return true;
}

bool CSeekWriter::RankNormalizeDAB(CDat &Dat,
	const vector<string> &vecstrGenes, int max_rank, float rbp_p){

	utype i, j;
	vector<utype> veciGenes;
	veciGenes.clear();
	veciGenes.resize(vecstrGenes.size());
	for( i = 0; i < vecstrGenes.size( ); ++i )
		veciGenes[ i ] = Dat.GetGene( vecstrGenes[i] );

	vector<float> vecSum;
	vector<int> vecNum;
	CSeekTools::InitVector(vecSum, vecstrGenes.size(), CMeta::GetNaN());
	CSeekTools::InitVector(vecNum, vecstrGenes.size(), (int)-9999);

	vector<vector<float> > mat;
	mat.resize(vecstrGenes.size());
	int max = max_rank;
	//float rbp_p = 0.99;

	bool expTransform = true;
	for(i=0; i<vecstrGenes.size(); i++){
		CSeekTools::InitVector(mat[i], vecstrGenes.size(), CMeta::GetNaN());

		utype s = veciGenes[i];
		if(CSeekTools::IsNaN(s)) continue;
		float *v = Dat.GetFullRow(s);
		vector<AResultFloat> vv;
		vv.resize(vecstrGenes.size());
		int numV = 0;		

		for(j=0; j<vecstrGenes.size(); j++){
			utype t = veciGenes[j];
			vv[j].i = j;
			if(CSeekTools::IsNaN(t) || CMeta::IsNaN(v[t])){
				vv[j].f = -9999;
				continue;
			}
			vv[j].f = v[t];
			numV++;
		}

		if(expTransform){
			nth_element(vv.begin(), vv.begin()+max, vv.end());
			sort(vv.begin(), vv.begin()+max);
			for(j=0; j<vecstrGenes.size(); j++){
				if(j<max){
					float rank = (1.0 - rbp_p) * pow(rbp_p, j);
					mat[i][vv[j].i] = rank;
				}else if(vv[j].f!=-9999){
					mat[i][vv[j].i] = 0;
				}
			}
		}else{
			sort(vv.begin(), vv.end());
			for(j=0; j<vecstrGenes.size(); j++){
				if(vv[j].f!=-9999){
					mat[i][vv[j].i] = numV - j;
				}
			}
		}

		free(v);
	}

	for(i=0; i<vecstrGenes.size(); i++){
		utype s = veciGenes[i];
		if(CSeekTools::IsNaN(s)) continue;
		for(j=i+1; j<vecstrGenes.size(); j++){
			utype t = veciGenes[j];
			if(CSeekTools::IsNaN(t)) continue;
			if(CMeta::IsNaN(mat[i][j]) || CMeta::IsNaN(mat[j][i])){
				fprintf(stderr, "%.3e %.3e\n", mat[i][j], mat[j][i]);
			}
			mat[i][j] = std::max(mat[i][j], mat[j][i]);
			mat[j][i] = mat[i][j];
		}
	}

	for(i=0; i<vecstrGenes.size(); i++){
		utype s = veciGenes[i];
		if(CSeekTools::IsNaN(s)) continue;
		vecSum[i] = 0;
		vecNum[i] = 0;
		for(j=0; j<vecstrGenes.size(); j++){
			utype t = veciGenes[j];
			if(CSeekTools::IsNaN(t)) continue;
			if(CMeta::IsNaN(mat[i][j])) continue;
			vecSum[i] += mat[i][j];
			vecNum[i]++;
		}
		//fprintf(stderr, "%.3e\n", vecSum[i]);
	}

	for(i=0; i<vecstrGenes.size(); i++){
		utype s = veciGenes[i];
		if(CSeekTools::IsNaN(s)) continue;
		for(j=0; j<vecstrGenes.size(); j++){
			utype t = veciGenes[j];
			if(CSeekTools::IsNaN(t)) continue;
			//fprintf(stderr, "%.3e %.3e\n", vecSum[i], vecSum[j]);
			float r = mat[i][j] / sqrtf(vecSum[i]) / sqrtf(vecSum[j]);
			Dat.Set(s, t, r);
		}
	}
	
	return true;
}



bool CSeekWriter::NormalizeDAB(CDat &Dat,
	const vector<string> &vecstrGenes,
	bool cutoff, bool expTransform, bool divideNorm, bool subtractNorm){

	utype i, j;
	vector<utype> veciGenes;
	veciGenes.clear();
	veciGenes.resize(vecstrGenes.size());
	for( i = 0; i < vecstrGenes.size( ); ++i )
		veciGenes[ i ] = Dat.GetGene( vecstrGenes[i] );

	vector<float> vecSum;
	vector<int> vecNum;
	CSeekTools::InitVector(vecSum, vecstrGenes.size(), CMeta::GetNaN());
	CSeekTools::InitVector(vecNum, vecstrGenes.size(), (int)-9999);

	for(i=0; i<vecstrGenes.size(); i++){
		utype s = veciGenes[i];
		if(CSeekTools::IsNaN(s)) continue;
		float *v = Dat.GetFullRow(s);
		float sum = 0;
		int num = 0;
		vector<float> all;
		for(j=0; j<vecstrGenes.size(); j++){
			utype t = veciGenes[j];
			if(CSeekTools::IsNaN(t)) continue;
			if(CMeta::IsNaN(v[t])) continue;
			if(cutoff){
				if(v[t]>0){
					if(expTransform)
						all.push_back(expf(-1.0*v[t]*v[t]/2.0));
					else
						all.push_back(v[t]);
				}
			}
			else{
				//fprintf(stderr, "Warning: Negative Z-Scores");
				if(expTransform)
					all.push_back(expf(-1.0*v[t]*v[t]/2.0));
				else
					all.push_back(v[t]);
			}	
		}

		for(j=0; j<all.size(); j++){
			sum+=all[j];
			num++;
		}
		vecSum[i] = sum;
		vecNum[i] = num;
		free(v);
	}

	for(i=0; i<vecstrGenes.size(); i++){
		utype s = veciGenes[i];
		if(CSeekTools::IsNaN(s)) continue;
		float *v = Dat.GetFullRow(s);

		for(j=0; j<vecstrGenes.size(); j++){
			utype t = veciGenes[j];
			if(CSeekTools::IsNaN(t)) continue;
			if(CMeta::IsNaN(v[t])) continue;
			if(cutoff){
				if(v[t]>0){
					if(expTransform){
						if(divideNorm){
							float r = expf(-1.0*v[t]*v[t]/2.0) / sqrtf(vecSum[i]) / sqrtf(vecSum[j]);
							Dat.Set(s, t, r);
						}else if(subtractNorm){
							float r = expf(-1.0*v[t]*v[t]/2.0) - vecSum[i] / vecNum[i] - vecSum[j] / vecNum[j];
							Dat.Set(s, t, r);
						}
					}else{
						if(divideNorm){
							float r = v[t] / sqrtf(vecSum[i]) / sqrtf(vecSum[j]);
							Dat.Set(s, t, r);
						}else if(subtractNorm){
							float r = v[t] - vecSum[i] / vecNum[i] - vecSum[j] / vecNum[j];
							Dat.Set(s, t, r);
						}
					}
				}else{
					Dat.Set(s, t, 0);
				}
			}
			else{
				if(expTransform){
					if(divideNorm){
						float r = expf(-1.0*v[t]*v[t]/2.0) / sqrtf(vecSum[i]) / sqrtf(vecSum[j]);
						Dat.Set(s, t, r);
					}else if(subtractNorm){
						float r = expf(-1.0*v[t]*v[t]/2.0) - vecSum[i] / vecNum[i] - vecSum[j] / vecNum[j];
						Dat.Set(s, t, r);
					}
				}else{
					if(divideNorm){
						float r = 0;
						//DANGEROUS
						if(vecSum[i]<=0){
							fprintf(stderr, "Warning, Dangerous, divide sqrt(z), where z could be negative\n");
							r = 0;
						}else{
							r = v[t] / sqrtf(vecSum[i]) / sqrtf(vecSum[j]);
						}
						Dat.Set(s, t, r);
					}else if(subtractNorm){
						float r = v[t] - vecSum[i] / vecNum[i] - vecSum[j] / vecNum[j];
						Dat.Set(s, t, r);
					}

				}
			}
		}
		free(v);
	}

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

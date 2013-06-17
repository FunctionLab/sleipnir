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

bool CSeekWriter::WriteSparseMatrix(CDataPair &Dat,
	vector< map<utype,unsigned short> > &umat, 
	int maxRank, const vector<string> &vecstrGenes, const char *fileName){

	FILE *f = fopen(fileName, "wb");
	if(f==NULL){
		cerr << "File not found!" << endl;
		return false;
	}
	utype numGenes = 0;
	utype i, j;

	vector<utype> veciGenes;
	veciGenes.clear();
	veciGenes.resize(vecstrGenes.size());
	for( i = 0; i < vecstrGenes.size( ); ++i )
		veciGenes[ i ] = Dat.GetGeneIndex( vecstrGenes[i] );

	CSeekIntIntMap mm(vecstrGenes.size());
	for(i=0; i<vecstrGenes.size(); i++)
		if(!CSeekTools::IsNaN(veciGenes[i]))
			mm.Add(i);

	utype numPresent = mm.GetNumSet();
	//1 utype
	fwrite((char*) (&numPresent), 1, sizeof(numPresent), f);
	const vector<utype> &allR = mm.GetAllReverse();
	//numPresent utype
	for(i=0; i<numPresent; i++)
		fwrite((char*) (&allR[i]), 1, sizeof(allR[i]), f);

	for(i=0; i<vecstrGenes.size(); i++){
		if(umat[i].size()==0) continue;
		numGenes++;
	}

	//1 utype
	fwrite((char*) (&numGenes), 1, sizeof(numGenes), f);

	for(i=0; i<vecstrGenes.size(); i++){
		unsigned short numEntries = umat[i].size(); //should be 1000
		if(numEntries==0) 
			continue;
		//1 utype
		fwrite((char*) (&i), 1, sizeof(i), f);
		//1 unsigned short
		fwrite((char*) (&numEntries), 1, sizeof(numEntries), f);
		map<utype,unsigned short>::iterator it;
		for(it=umat[i].begin(); it!=umat[i].end(); it++){
			utype first = it->first;
			unsigned short second = it->second;
			//1 utype
			fwrite((char*) (&first), 1, sizeof(first), f);
			//1 unsigned short
			fwrite((char*) (&second), 1, sizeof(second), f);
		}
	}
	fclose(f);
	return true;
}

bool CSeekWriter::GetSparseRankMatrix(CDataPair &Dat,
	vector< map<utype,unsigned short> > &umat, 
	int maxRank, //1000
	const vector<string> &vecstrGenes){

	utype i, j;
	vector<utype> veciGenes;
	veciGenes.clear();
	veciGenes.resize(vecstrGenes.size());
	for( i = 0; i < vecstrGenes.size( ); ++i )
		veciGenes[ i ] = Dat.GetGeneIndex( vecstrGenes[i] );
	umat.resize(vecstrGenes.size());
	for(i=0; i<vecstrGenes.size(); i++){
		umat[i] = map<utype, unsigned short>();
	}

	fprintf(stderr, "Start reading DAB...\n");
	for(i=0; i<vecstrGenes.size(); i++){
		utype s = veciGenes[i];
		if(CSeekTools::IsNaN(s)) continue;
		if(i%1000==0)
			fprintf(stderr, "Start reading gene %d...\n", i);

		//float *v = Dat.GetRowSeek(s);
		float *v = Dat.GetFullRow(s);
		vector<AResultFloat> vv;
		vv.resize(vecstrGenes.size());

		for(j=0; j<vecstrGenes.size(); j++){
			utype t = veciGenes[j];
			vv[j].i = j;
			vv[j].f = -9999;
			if(CSeekTools::IsNaN(t)) continue;
			float d = v[t];
			if(CMeta::IsNaN(d)) continue;
			vv[j].f = d;
		}

		nth_element(vv.begin(), vv.begin()+maxRank, vv.end());
		sort(vv.begin(), vv.begin()+maxRank);

		for(j=0; j<vecstrGenes.size(); j++){
			if(j<maxRank){
				utype first = i;
				utype second = vv[j].i;
				if(i >= vv[j].i){
					first = vv[j].i;
					second = i;
				}
				map<utype,unsigned short>::iterator it;
				if((it=umat[first].find(second))==umat[first].end())
					umat[first][second] = (unsigned short) j;
				else
					umat[first][second] = std::min(it->second, (unsigned short) j);
			}
		}
		free(v);
	}
	fprintf(stderr, "Finished reading DAB\n");
	return true;
}

bool CSeekWriter::NormalizeDAB(CDataPair &Dat,
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
		float sum = 0;
		int num = 0;
		vector<float> all;
		for(j=0; j<vecstrGenes.size(); j++){
			utype t = veciGenes[j];
			float d = Dat.Get(s,t);
			if(CSeekTools::IsNaN(t)) continue;
			if(CMeta::IsNaN(d)) continue;
			if(cutoff){
				if(d>0){
					if(expTransform)
						all.push_back(expf(-1.0*d*d/2.0));
					else
						all.push_back(d);
				}
			}
			else{
				//fprintf(stderr, "Warning: Negative Z-Scores");
				if(expTransform)
					all.push_back(expf(-1.0*d*d/2.0));
				else
					all.push_back(d);
			}	
		}

		for(j=0; j<all.size(); j++){
			sum+=all[j];
			num++;
		}
		vecSum[i] = sum;
		vecNum[i] = num;
	}

	for(i=0; i<vecstrGenes.size(); i++){
		utype s = veciGenes[i];
		if(CSeekTools::IsNaN(s)) continue;
		float *v = Dat.GetFullRow(s);

		for(j=0; j<vecstrGenes.size(); j++){
			utype t = veciGenes[j];
			float d = v[t];
			if(CSeekTools::IsNaN(t)) continue;
			if(CMeta::IsNaN(d)) continue;
			if(cutoff){
				if(d>0){
					if(expTransform){
						if(divideNorm){
							float r = expf(-1.0*d*d/2.0) / sqrtf(vecSum[i]) / sqrtf(vecSum[j]);
							Dat.Set(s, t, r);
						}else if(subtractNorm){
							float r = expf(-1.0*d*d/2.0) - vecSum[i] / vecNum[i] - vecSum[j] / vecNum[j];
							Dat.Set(s, t, r);
						}
					}else{
						if(divideNorm){
							float r = d / sqrtf(vecSum[i]) / sqrtf(vecSum[j]);
							Dat.Set(s, t, r);
						}else if(subtractNorm){
							float r = d - vecSum[i] / vecNum[i] - vecSum[j] / vecNum[j];
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
						float r = expf(-1.0*d*d/2.0) / sqrtf(vecSum[i]) / sqrtf(vecSum[j]);
						Dat.Set(s, t, r);
					}else if(subtractNorm){
						float r = expf(-1.0*d*d/2.0) - vecSum[i] / vecNum[i] - vecSum[j] / vecNum[j];
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
							r = d / sqrtf(vecSum[i]) / sqrtf(vecSum[j]);
						}
						Dat.Set(s, t, r);
					}else if(subtractNorm){
						float r = d - vecSum[i] / vecNum[i] - vecSum[j] / vecNum[j];
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

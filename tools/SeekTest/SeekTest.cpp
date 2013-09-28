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
#include <iomanip>


float** LoadGenes(const vector<string> struserGenes, 
	const vector<utype> &veciGenes, const vector<string> &vecstrGenes, 
	const vector<utype> &veciallGenes, CSeekIntIntMap *gmap, 
	CSeekDataset *vcd, float **vall){

	float** v2 = CSeekTools::Init2DArray(struserGenes.size(), 
		vecstrGenes.size(), CMeta::GetNaN());
	size_t i, j;

	for(i=0; i<struserGenes.size(); i++){
		utype s = veciGenes[i]; //s is user gene ID
		if(CSeekTools::IsNaN(s)) continue;

		for(j=0; j<vecstrGenes.size(); j++){
			utype t = veciallGenes[j];
			if(CSeekTools::IsNaN(t)) continue;
			if(CMeta::IsNaN(vall[s][t])) continue;
			if(CSeekTools::IsNaN(gmap->GetForward(t))) continue;
			v2[i][t] = vall[s][t];
		}
	}
	return v2;
}

float** ReadUserGenes(const vector<string> &struserGenes, 
	const vector<utype> &veciGenes, const vector<string> &vecstrGenes, 
	const vector<utype> &veciallGenes, CSeekIntIntMap *gmap, 
	CSeekDataset *vcd, CDataPair &Dat){

	float** v2 = CSeekTools::Init2DArray(struserGenes.size(), 
		vecstrGenes.size(), CMeta::GetNaN());
	size_t i, j;

	for(i=0; i<struserGenes.size(); i++){
		utype s = veciGenes[i]; //s is user gene ID
		if(CSeekTools::IsNaN(s)) continue;
		float *v = Dat.GetFullRow(s);

		for(j=0; j<vecstrGenes.size(); j++){
			utype t = veciallGenes[j];
			if(CSeekTools::IsNaN(t)) continue;
			if(CMeta::IsNaN(v[t])) continue;
			if(CSeekTools::IsNaN(gmap->GetForward(t))) continue;
			//Disable gene average subtraction!
			//v[t] = v[t] - vcd->GetGeneAverage(t);
			v2[i][t] = v[t];
			//fprintf(stderr, "%.2f\n", v2[i][t]);
		}
		free(v);
	}	
	return v2;
}

float** ReadAllGenes(const vector<string> &vecstrGenes, 
	const vector<utype> &veciallGenes, CSeekIntIntMap *gmap, 
	CSeekDataset *vcd, CDataPair &Dat){

	float** v2 = CSeekTools::Init2DArray(vecstrGenes.size(), 
		vecstrGenes.size(), CMeta::GetNaN());
	size_t i, j;
	size_t ss = Dat.GetGenes();

	map<string, utype> mapstrintGenes;
	for(i=0; i<vecstrGenes.size(); i++)
		mapstrintGenes[vecstrGenes[i]] = i;


	//fprintf(stderr, "thread limit is %d\n", omp_get_max_threads());

	//getchar();
	#pragma omp parallel for \
	private(i, j) \
	shared(Dat, gmap, mapstrintGenes, vcd, v2) \
	firstprivate(ss) \
	schedule(dynamic)
	for(i=0; i<ss; i++){
		utype tid = omp_get_thread_num();
		string s = Dat.GetGene(i);
		utype ii = mapstrintGenes[s];
		if(CSeekTools::IsNaN(gmap->GetForward(ii))) continue;
		float *v = Dat.GetFullRow(i);

		for(j=0; j<ss; j++){
			string t = Dat.GetGene(j);
			utype jj = mapstrintGenes[t];
			if(CSeekTools::IsNaN(gmap->GetForward(jj))) continue;
			v[j] = v[j] - vcd->GetGeneAverage(jj);
			v2[ii][jj] = v[j];
			//fprintf(stderr, "%.2f\n", v2[s][t]);
		}
		//fprintf(stderr, "Done row %d", i);
		//getchar();
		free(v);
		//fprintf(stderr, "TID %d Gene %d is done.\n", tid, i);
		//fprintf(stderr, "Gene %d / %d is done.\n", i, ss);
		//fprintf(stderr, "Done freeing row %d", i); getchar();
	}	
	return v2;
}

bool GetRandomGenes(gsl_rng *r, int size, vector<string> &struserGenes, 
	vector<utype> &veciGenes, const vector<string> &vecstrGenes,
	const vector<utype> &veciallGenes, CSeekIntIntMap *gmap){

	struserGenes.clear();
	veciGenes.clear();
	struserGenes.resize(size);
	veciGenes.resize(size);

	utype i;
	int totSize = 0;

	for(i=0; i<veciallGenes.size(); i++){
		utype t = veciallGenes[i];
		if(CSeekTools::IsNaN(t)) continue;
		if(CSeekTools::IsNaN(gmap->GetForward(t))) continue;
		totSize++;
	}

	int *gsample = (int*)malloc(size*sizeof(int));		
	int *gs = (int*)malloc(totSize*sizeof(int));
	utype ti = 0;

	for(i=0; i<veciallGenes.size(); i++){
		utype t = veciallGenes[i];
		if(CSeekTools::IsNaN(t)) continue;
		if(CSeekTools::IsNaN(gmap->GetForward(t))) continue;
		gs[ti] = i;
		ti++;
	}

	gsl_ran_choose(r, gsample, size, gs, totSize, sizeof(int));

	for(i=0; i<size; i++){
		veciGenes[i] = veciallGenes[gsample[i]];
		struserGenes[i] = vecstrGenes[gsample[i]];
	}

	free(gsample);
	free(gs);
	return true;
}

bool GetMeanStdev(const vector<float> &f, float &mean, float &stdev){
	size_t i, j;
	mean = 0;
	stdev = 0;
	for(i=0; i<f.size(); i++){
		mean+=f[i];
	}
	mean/=f.size();
	for(i=0; i<f.size(); i++){
		stdev+=(f[i] - mean) * (f[i] - mean);
	}
	stdev/=f.size();
	stdev=sqrt(stdev);
	return true;
}

bool CalculateWelch(const float &mean1, const float &stdev1, const int &size1, 
	const float &mean2, const float &stdev2, const int &size2, 
	double &pvalt, double &t){

	double v1 = (double) stdev1 * (double) stdev1;
	double v2 = (double) stdev2 * (double) stdev2;

	int n1 = size1;
	int n2 = size2;
	t = (mean1 - mean2) /sqrt(v1/n1 + v2/n2);
	double x1 = v1/n1;
	double x2 = v2/n2;
	double i1 = x1 + x2;
	double df = i1*i1 / (x1*x1/(n1-1) + x2*x2/(n2-1));
	if(t>0){
		//fprintf(stderr, "t:%.2f df:%.2f n1:%d n2:%d mean1:%.2f mean2:%.2f v1:%.2f v2:%.2f\n", 
		//	t, df, n1, n2, mean1, mean2, v1, v2);
		double pval1 = gsl_cdf_tdist_Q(t, df);
		double pval2 = gsl_cdf_tdist_P(-1.0*t, df);
		pvalt = pval1 + pval2;

	}else{
		//fprintf(stderr, "t:%.2f df:%.2f n1:%d n2:%d mean1:%.2f mean2:%.2f v1:%.2f v2:%.2f\n", 
		//	t, df, n1, n2, mean1, mean2, v1, v2);
		double pval1 = gsl_cdf_tdist_Q(-1.0*t, df);
		double pval2 = gsl_cdf_tdist_P(t, df);
		pvalt = pval1 + pval2;
	}

	return true;
}

vector<string> Do_Pair_Proportion(){


}


vector<string> Do_T_Test(
	gsl_rng *rnd,
	float **vall, 
	const vector<utype> &veciGenes, //process
	const vector<string> &vecstrProcess, //process

	const vector<utype> &veciallGenes, //background
	const vector<string> &vecstrGenes, //background
	CSeekIntIntMap *gmap
){

	unsigned int seed = gsl_rng_get(rnd);
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);

	size_t i, j, k;

	vector< vector<utype> > vi;
	vector< vector<string> > si;
	vi.resize(100);
	si.resize(100);
	int size = veciGenes.size();
	vector<string> outstr;

	for(i=0; i<100; i++){
		vi[i] = vector<utype>();
		si[i] = vector<string>();
		GetRandomGenes(r, size, si[i], vi[i], vecstrGenes, veciallGenes, gmap);
	}

	float user_mean, user_stdev;
	vector<float> rand_mean, rand_stdev;
	rand_mean.resize(100);
	rand_stdev.resize(100);

	for(i=0; i<100; i++){
		vector<float> vv; //preparation
		for(j=0; j<veciGenes.size(); j++){
			for(k=0; k<veciGenes.size(); k++){
				if(CMeta::IsNaN(vall[vi[i][j]][vi[i][k]])) continue;
				vv.push_back(vall[vi[i][j]][vi[i][k]]);
			}
		}
		GetMeanStdev(vv, rand_mean[i], rand_stdev[i]);
		//fprintf(stderr, "Random %.2f %.2f\n", rand_mean[i], rand_stdev[i]);
	}

	vector<float> vf; //preparation
	for(j=0; j<veciGenes.size(); j++){
		for(k=0; k<veciGenes.size(); k++){
			if(CMeta::IsNaN(vall[veciGenes[j]][veciGenes[k]])) continue;
			vf.push_back(vall[veciGenes[j]][veciGenes[k]]);
		}
	}
	GetMeanStdev(vf, user_mean, user_stdev);

	ostringstream ost;
	ost.setf(ios::fixed);
	ost << "User " << setprecision (2) << user_mean << " " << setprecision (2) << user_stdev;
	outstr.push_back(ost.str());
	
	for(i=0; i<100; i++){
		double pvalt, t;
		CalculateWelch(rand_mean[i], rand_stdev[i], veciGenes.size(), 
			user_mean, user_stdev, veciGenes.size(), pvalt, t);
		//printf("%.2f %.2f %.3E\n", (float) user_mean - rand_mean[i], t, pvalt);
		ostringstream ost2;
		ost2.setf(ios::fixed);
		ost2 << setprecision(2) << (float) user_mean - rand_mean[i] << " " << setprecision(2) << t << " " << setprecision(2) << scientific << pvalt;
		outstr.push_back(ost2.str());
	}
	return outstr;
}

bool Get_Mann_Whitney_U_Statistics(const vector<float> &v1, const vector<float> &v2,
	double &U1, double &U2, double &z, double &auc){

	size_t i, j, k;
	vector<AResultFloat> vv; //preparation
	int s1 = v1.size();
	int s2 = v2.size();
	vv.resize(v1.size() + v2.size());
	size_t kk = 0;
	for(i=0; i<v1.size(); i++){
		vv[kk].i = 0;
		vv[kk].f = v1[i];
		kk++;
	}
	for(i=0; i<v2.size(); i++){
		vv[kk].i = 1;
		vv[kk].f = v2[i];
		kk++;
	}

	sort(vv.begin(), vv.end());

	U1 = 0;
	U2 = 0;	
	for(j=0; j<vv.size(); j++){
		if(vv[j].i==0){
			U1+=j+1.0;
		}else{
			U2+=j+1.0;
		}
	}
	U1-=(s1*(s1+1.0)/2.0);
	U2-=(s2*(s2+1.0)/2.0);

	//normal approximation 		
	double mu = s1*s2 / 2.0;
	double sigma_u = sqrt((s1*s2)*(s1+s2+1.0)/12.0);
	z = (U1 - mu) / sigma_u;
	auc = U1/(s1*s2);
	return true;
}

vector<string> Do_Mann_Whitney_U_Test(
	gsl_rng *rnd,
	float **vall, 
	const vector<utype> &veciGenes, //process
	const vector<string> &vecstrProcess, //process

	const vector<utype> &veciallGenes, //background
	const vector<string> &vecstrGenes, //background
	CSeekIntIntMap *gmap
){

	unsigned int seed = gsl_rng_get(rnd);
	vector<string> outstr;

	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);

	size_t i, j, k;

	vector< vector<utype> > vi;
	vector< vector<string> > si;
	vi.resize(100);
	si.resize(100);
	int size = veciGenes.size();

	for(i=0; i<100; i++){
		vi[i] = vector<utype>();
		si[i] = vector<string>();
		GetRandomGenes(r, size, si[i], vi[i], vecstrGenes, veciallGenes, gmap);
	}

	vector<double> U1, U2;
	vector<float> z;
	vector<int> s1, s2;
	vector<float> auc;

	U1.resize(100); U2.resize(100); z.resize(100); 
	s1.resize(100); s2.resize(100); auc.resize(100);

	for(i=0; i<100; i++){
		vector<AResultFloat> vv; //preparation
		s1[i] = 0;
		s2[i] = 0;
		for(j=0; j<veciGenes.size(); j++){
			for(k=0; k<veciGenes.size(); k++){
				if(CMeta::IsNaN(vall[vi[i][j]][vi[i][k]])) continue;
				s1[i]++;
			}
		}
		for(j=0; j<veciGenes.size(); j++){
			for(k=0; k<veciGenes.size(); k++){
				if(CMeta::IsNaN(vall[veciGenes[j]][veciGenes[k]])) continue;
				s2[i]++;
			}
		}
		vv.resize(s1[i]+s2[i]);

		int kk = 0;
		for(j=0; j<veciGenes.size(); j++){
			for(k=0; k<veciGenes.size(); k++){
				if(CMeta::IsNaN(vall[vi[i][j]][vi[i][k]])) continue;
				vv[kk].f = vall[vi[i][j]][vi[i][k]];
				vv[kk].i = 0;
				kk++;
			}
		}
		for(j=0; j<veciGenes.size(); j++){
			for(k=0; k<veciGenes.size(); k++){
				if(CMeta::IsNaN(vall[veciGenes[j]][veciGenes[k]])) continue;
				vv[kk].f = vall[veciGenes[j]][veciGenes[k]];
				vv[kk].i = 1;
				kk++;
			}
		}
		sort(vv.begin(), vv.end());
	
		U1[i] = 0;
		U2[i] = 0;
		for(j=0; j<vv.size(); j++){
			if(vv[j].i==0){
				U1[i]+=j+1.0;
			}else{
				U2[i]+=j+1.0;
			}
		}
		U1[i]-=(s1[i]*(s1[i]+1.0)/2.0);
		U2[i]-=(s2[i]*(s2[i]+1.0)/2.0);
	
		/*double U = 0;
		if(U1<U2){
			U = U1;
		}else{
			U = U2;
		}*/

		//normal approximation 		
		double mu = s1[i]*s2[i] / 2.0;
		double sigma_u = sqrt((s1[i]*s2[i])*(s1[i]+s2[i]+1.0)/12.0);
		z[i] = (U1[i] - mu) / sigma_u;
		auc[i] = U1[i]/(s1[i]*s2[i]);
		//printf("%d %.2f %.2f %.2f %.2f\n", s1, (float) U1, (float) z, (float) U1/(s1*s2), (float) U2/(s1*s2));
	}
	
	for(i=0; i<100; i++){
		ostringstream ost;
		ost.setf(ios::fixed);
		ost << s1[i] << " " << U1[i] << " " << setprecision (2) << z[i] << " " << setprecision (2) << auc[i];
		outstr.push_back(ost.str());
	}

	return outstr;
}

vector<string> Do_Mann_Whitney_U_Test_By_Gene(gsl_rng *rnd, float **vall,
	const vector<string> &vecstrProcess,
	const vector<utype> &veciProcess,
	const vector<string> &vecstrGenes,
	const vector<utype> &veciallGenes,
	CSeekIntIntMap *gmap
){
	unsigned int seed = gsl_rng_get(rnd);

	vector<string> outstr;
	size_t i, j, jj, k;
	int size = veciProcess.size();
	int trials = 1;
	int perGeneTrials = 11;

	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);

	for(i=0; i<trials; i++){
		vector<utype> veciRandom;
		vector<string> vecstrRandom;
		GetRandomGenes(r, size, vecstrRandom, veciRandom, vecstrGenes, veciallGenes, gmap);

		outstr.push_back("Process");

		//Process
		for(j=0; j<size; j++){
			utype thisGene = veciProcess[j];
			vector<float> vProcess;
			for(k=0; k<size; k++){
				float x = vall[thisGene][veciProcess[k]];
				if(CMeta::IsNaN(x)) continue;
				vProcess.push_back(x);
			}

			ostringstream ost; 
			ost << "Gene " << vecstrProcess[j] << " size " << vProcess.size();
			outstr.push_back(ost.str());
			//fprintf(stdout, "Gene %s size %d\n", vecstrProcess[j].c_str(), vProcess.size());

			vector<AResultFloat> ar;
			vector<double> U1, U2, z, auc;

			ar.resize(perGeneTrials);
			U1.resize(perGeneTrials); 
			U2.resize(perGeneTrials); 
			z.resize(perGeneTrials); 
			auc.resize(perGeneTrials);

			for(jj=0; jj<perGeneTrials; jj++){
				vector<utype> veciRandomGene;
				vector<string> vecstrRandomGene;
				GetRandomGenes(r, size, vecstrRandomGene, veciRandomGene, vecstrGenes, veciallGenes, gmap);

				vector<float> vRandom;
				for(k=0; k<size; k++){
					float x = vall[thisGene][veciRandomGene[k]];
					if(CMeta::IsNaN(x)) continue;
					vRandom.push_back(x);
				}
				Get_Mann_Whitney_U_Statistics(vProcess, vRandom, U1[jj], U2[jj], z[jj], auc[jj]);

				ar[jj].i = jj;
				ar[jj].f = auc[jj];

				//fprintf(stdout, "%d %d %.2f %.2f\n", (int) U1, (int) U2, z, auc);
			}

			sort(ar.begin(), ar.end());
			int med = perGeneTrials/2;
			int ar_med = ar[med].i;
			//Show the median
			ostringstream ost2;
			ost2.setf(ios::fixed);
			ost2 << (int) U1[ar_med] << " " << (int) U2[ar_med] << " " << setprecision(2) << z[ar_med] << " " << setprecision(2) << auc[ar_med];
			outstr.push_back(ost2.str());
			//fprintf(stdout, "%d %d %.2f %.2f\n", (int) U1[ar_med], (int) U2[ar_med], z[ar_med], auc[ar_med]);

		}

		outstr.push_back("Random");
		//fprintf(stdout, "Random\n");
		//If process genes were just randomly selected geneset
		for(j=0; j<size; j++){
			utype thisGene = veciRandom[j];
			vector<float> vRan;
			for(k=0; k<size; k++){
				float x = vall[thisGene][veciRandom[k]];
				if(CMeta::IsNaN(x)) continue;
				vRan.push_back(x);
			}

			ostringstream ost; 
			ost << "Gene number " << j << " size " << vRan.size();
			outstr.push_back(ost.str());

			vector<AResultFloat> ar;
			vector<double> U1, U2, z, auc;

			ar.resize(perGeneTrials);
			U1.resize(perGeneTrials); 
			U2.resize(perGeneTrials); 
			z.resize(perGeneTrials); 
			auc.resize(perGeneTrials);

			for(jj=0; jj<perGeneTrials; jj++){
				vector<utype> veciRandomGene;
				vector<string> vecstrRandomGene;
				GetRandomGenes(r, size, vecstrRandomGene, veciRandomGene, vecstrGenes, veciallGenes, gmap);
				vector<float> vRandom;
				for(k=0; k<size; k++){
					float x = vall[thisGene][veciRandomGene[k]];
					if(CMeta::IsNaN(x)) continue;
					vRandom.push_back(x);
				}
				Get_Mann_Whitney_U_Statistics(vRan, vRandom, U1[jj], U2[jj], z[jj], auc[jj]);

				//fprintf(stdout, "%d %d %.2f %.2f\n", (int) U1, (int) U2, z, auc);
				ar[jj].i = jj;
				ar[jj].f = auc[jj];
			}

			sort(ar.begin(), ar.end());
			int med = perGeneTrials/2;
			int ar_med = ar[med].i;
			//Show the median
			//fprintf(stdout, "%d %d %.2f %.2f\n", (int) U1[ar_med], (int) U2[ar_med], z[ar_med], auc[ar_med]);
			ostringstream ost2;
			ost2.setf(ios::fixed);
			ost2 << (int) U1[ar_med] << " " << (int) U2[ar_med] << " " << setprecision (2) << z[ar_med] << " " << setprecision (2) << auc[ar_med];
			outstr.push_back(ost2.str());

		}


	}

	return outstr;
}

vector<string> Do_One(const char *file, gsl_rng *rnd, CSeekDataset *vcd, float **vall,
	map<string, utype> &mapstrintGene, vector<string> &vecstrGenes){

	size_t i, j, k;
	vector<string> ostr;

	vector<string> vecstrUserGenes;
	CSeekStrIntMap mapTmp;
	if(!CSeekTools::ReadListOneColumn(file, vecstrUserGenes, mapTmp))
		return ostr;

	vector<utype> userGenes;
	for(i=0; i<vecstrUserGenes.size(); i++){
		userGenes.push_back(mapstrintGene[vecstrUserGenes[i]]);
	}

	CSeekIntIntMap *gmap = vcd->GetGeneMap();

	vector<string> struserGenes;
	for(i=0; i<userGenes.size(); i++){
		if(CSeekTools::IsNaN(gmap->GetForward(userGenes[i]))){
			ostringstream ost;
			ost << "Error: gene " << vecstrUserGenes[i] << " is not found in this dataset";
			ostr.push_back(ost.str());
			continue;
		}
		struserGenes.push_back(vecstrGenes[userGenes[i]]);
	}

	if(struserGenes.size()==0 || struserGenes.size()<5){
		return ostr;
	}

	vector<utype> veciGenes;
	veciGenes.clear(); //only user gene set
	veciGenes.resize(struserGenes.size());
	for( i = 0; i < struserGenes.size( ); ++i ){
		//veciGenes[ i ] = Dat.GetGene( struserGenes[i] );
		veciGenes[i] = mapstrintGene[struserGenes[i]];
	}

	vector<utype> veciallGenes;
	veciallGenes.clear(); //all genes
	veciallGenes.resize(vecstrGenes.size());
	for( i = 0; i < vecstrGenes.size( ); ++i ){
		//veciallGenes[ i ] = Dat.GetGene( vecstrGenes[i] );
		veciallGenes[i] = mapstrintGene[vecstrGenes[i]];
	}

	//float** v2 = LoadGenes(struserGenes, veciGenes, vecstrGenes, 
	//	veciallGenes, gmap, vcd, vall);

	vector<string> outstr = 
	//	Do_Mann_Whitney_U_Test_By_Gene(rnd, vall, struserGenes, veciGenes, vecstrGenes, veciallGenes, gmap);
	Do_T_Test(rnd, vall, veciGenes, struserGenes, veciallGenes, vecstrGenes, gmap);
	//Do_Mann_Whitney_U_Test(rnd, vall, veciGenes, struserGenes, veciallGenes, vecstrGenes, gmap);


	for(i=0; i<outstr.size(); i++){
		ostr.push_back(outstr[i]);
	}
	
	//Do_T_Test(vall, veciGenes, struserGenes, veciallGenes, vecstrGenes, gmap);
	//ostr = Do_Mann_Whitney_U_Test(vall, veciGenes, struserGenes, veciallGenes, vecstrGenes, gmap);

	/*
	int *a1 = (int*)malloc(2*sizeof(int));
	int *a2 = (int*)malloc(size*sizeof(int));
	for(i=0; i<size; i++){
		a2[i] = i;
	}

	fprintf(stderr, "Random to random\n");
	for(i=0; i<100; i++){
		gsl_ran_choose(rnd, a1, 2, a2, size, sizeof(int));
		double pvalt, t;
		CalculateWelch(rand_mean[a1[0]], rand_stdev[a1[0]], veciGenes.size(), 
			rand_mean[a1[1]], rand_stdev[a1[1]], veciGenes.size(), pvalt, t);
		fprintf(stderr, "%.2f %.3E\n", t, pvalt);
	}

	free(a1);
	free(a2);
	*/
	//CSeekTools::Free2DArray(v2);

	//for(i=0; i<100; i++){
	//	CSeekTools::Free2DArray(randomGenes[i]);
	//}
	return ostr;
}


int main( int iArgs, char** aszArgs ) {
	static const size_t	c_iBuffer	= 1024;
	gengetopt_args_info	sArgs;
	ifstream			ifsm;
	istream*			pistm;
	vector<string>		vecstrLine, vecstrGenes, vecstrDatasets, vecstrQuery, vecstrUserDatasets;
	char				acBuffer[ c_iBuffer ];
	size_t				i, j, k;


	omp_set_num_threads(8);
	/*utype ii;
	omp_set_num_threads(8);
	#pragma omp parellel for \
	private(i) \
	schedule(static, 1)
	for(ii=0; ii<100; ii++){
		utype tid = omp_get_thread_num();
		fprintf(stderr, "Hello World %d\n", tid);
	}*/


	const gsl_rng_type *T;
	gsl_rng *rnd;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rnd = gsl_rng_alloc(T);

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }

	vector<string> vecstrGeneID;
	map<string, utype> mapstrintGene;
	if(!CSeekTools::ReadListTwoColumns(sArgs.input_arg, vecstrGeneID, vecstrGenes))
		return 1;
	for(i=0; i<vecstrGenes.size(); i++)
		mapstrintGene[vecstrGenes[i]] = i;

	if(sArgs.db_flag==1){
		vector<float> quant;
		CSeekTools::ReadQuantFile(sArgs.quant_arg, quant);

		vector<string> vecstrDataset, vDP;

		if(!CSeekTools::ReadListTwoColumns(sArgs.dataset_list_arg, 
		vecstrDataset, vDP))
			return false;

		CDatabase *DB = new CDatabase(false);
		DB->Open(sArgs.db_dir_arg, vecstrGenes, vecstrDataset.size(), sArgs.db_num_arg);

		string strPrepInputDirectory = sArgs.prep_arg;
		string strSinfoInputDirectory = sArgs.sinfo_arg;
		vector<CSeekDataset*> vc;
		vc.resize(vecstrDataset.size());
		size_t i, j, k;
		for(i=0; i<vecstrDataset.size(); i++){
			vc[i] = new CSeekDataset();
			string strFileStem = vecstrDataset[i];
			string strAvgPath = strPrepInputDirectory+"/"+
				strFileStem + ".gavg";
			string strPresencePath = strPrepInputDirectory+"/"+
				strFileStem + ".gpres";
			string strSinfoPath = strSinfoInputDirectory+"/"+
				strFileStem + ".sinfo";
			vc[i]->ReadGeneAverage(strAvgPath);
			vc[i]->ReadGenePresence(strPresencePath);
			vc[i]->ReadDatasetAverageStdev(strSinfoPath);
			vc[i]->InitializeGeneMap();
		}
	
		size_t iGenes = vecstrGenes.size();
		size_t iDatasets = vecstrDataset.size();
		
		vector<string> vecstrQuery;
		CSeekTools::ReadMultiGeneOneLine(sArgs.query_arg, vecstrQuery);
		//Need to load the query

		vector<char> cQuery;
		CSeekTools::InitVector(cQuery, iGenes, (char)0);
		for(i=0; i<vecstrQuery.size(); i++){
			if(mapstrintGene.find(vecstrQuery[i])==mapstrintGene.end()) continue;
			utype k = mapstrintGene.find(vecstrQuery[i])->second;
			cQuery[k] = 1;
		}
		vector<utype> allQ;
		for(i=0; i<cQuery.size(); i++)
			if(cQuery[i]==1)
				allQ.push_back(i);
		allQ.resize(allQ.size());

		for(i=0; i<iDatasets; i++)
			if(vc[i]->GetDBMap()!=NULL)
				vc[i]->DeleteQueryBlock();
	
		for(i=0; i<iDatasets; i++)
			vc[i]->InitializeQueryBlock(allQ);

		size_t m, d;
		for(i=0; i<allQ.size(); i++){
			m = allQ[i];
			vector<unsigned char> Qi;
			if(!DB->GetGene(m, Qi)){
				cerr << "Gene does not exist" << endl;
				continue;
			}
			utype db;
			CSeekIntIntMap *qu = NULL;
			unsigned char **r = NULL;
			for(j=0; j<vecstrDataset.size(); j++){
				if((qu=vc[j]->GetDBMap())==NULL)
					continue;
				if(CSeekTools::IsNaN(db=(qu->GetForward(m)))) 
					continue;
				for(r = vc[j]->GetMatrix(), k=0; k<iGenes; k++)
					r[db][k] = Qi[k*vecstrDataset.size()+j];
			}
			Qi.clear();
		}

		if(!!sArgs.count_pair_flag){
			map<float,vector<int> > countPairs;
			float point = 0.0;
			while(point<=5.0){
				countPairs[point] = vector<int>();
				CSeekTools::InitVector(countPairs[point], iDatasets, (int)0);
				point+=0.25;
			}
			for(k=0; k<iDatasets; k++){
				CSeekIntIntMap *mapQ = vc[k]->GetDBMap();
				CSeekIntIntMap *mapG = vc[k]->GetGeneMap();
				if(mapQ==NULL) continue;
				unsigned char **f = vc[k]->GetMatrix();
				size_t qi, qj;
				for(qi=0; qi<allQ.size(); qi++){
					utype gene_qi = allQ[qi];
					utype iQ = mapQ->GetForward(gene_qi);
					if(CSeekTools::IsNaN(iQ)) continue;
					for(qj=qi+1; qj<allQ.size(); qj++){
						utype gene_qj = allQ[qj];
						utype jQ = mapG->GetForward(gene_qj);
						if(CSeekTools::IsNaN(jQ)) continue;
						unsigned char uc = f[iQ][gene_qj];
						if(uc==255) continue;
						float vv = quant[uc];
						point = 0.0;
						while(point<=5.0){
							if(vv>point)
								countPairs[point][k]++;
							point+=0.25;
						}
					}
				}
			}
			for(i=0; i<iDatasets; i++)
				vc[i]->DeleteQueryBlock();
			point = 0.0;			
			while(point<=5.0){
				sort(countPairs[point].begin(), countPairs[point].end(), greater<int>());
				float tmp = 0;
				for(i=0; i<10; i++)
					tmp+=(float)countPairs[point][i];
				tmp/=10.0;
				fprintf(stderr, "%.2f\t%.1f pairs\n", point, tmp);
				point+=0.25;
			}
		}

		if(!!sArgs.histogram_flag){
			srand(unsigned(time(0)));
			vector<int> dID;
			for(k=0; k<iDatasets; k++)
				dID.push_back(k);
			random_shuffle(dID.begin(), dID.end());
			utype kk;
			for(kk=0; kk<100; kk++){
				k = dID[kk];
				CSeekIntIntMap *mapQ = vc[k]->GetDBMap();
				CSeekIntIntMap *mapG = vc[k]->GetGeneMap();
				if(mapQ==NULL) continue;
				unsigned char **f = vc[k]->GetMatrix();
				size_t qi, ii;
				for(qi=0; qi<allQ.size(); qi++){
					utype gene_qi = allQ[qi];
					utype iQ = mapQ->GetForward(gene_qi);
					if(CSeekTools::IsNaN(iQ)) continue;
					vector<float> z_score;
					const vector<utype> &allRGenes = mapG->GetAllReverse();
					for(ii=0; ii<mapG->GetNumSet(); ii++){
						size_t ij = allRGenes[ii];
						float a = vc[k]->GetGeneAverage(ij);
						unsigned char x = f[iQ][ij];
						if(x==255) continue;
						//float xnew = quant[x] - a;
						float xnew = quant[x];
						z_score.push_back(xnew);
					}
					sort(z_score.begin(), z_score.end());
					float pts[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
					fprintf(stderr, "Query %s\tDataset %d\t", vecstrGenes[gene_qi].c_str(), k);
					for(ii=0; ii<9; ii++)
						fprintf(stderr, "%.2f\t", z_score[(int)(pts[ii]*z_score.size())]);
					fprintf(stderr, "\n");
				}
			}
			return 0;
		}

	}
		
	CSeekStrIntMap mapTmp;
	vector<string> vecstrList;
	if(!CSeekTools::ReadListOneColumn(sArgs.gene_set_list_arg, vecstrList, mapTmp))
		return 1;

	if(sArgs.dab_flag==1){
		CSeekDataset *vcd = new CSeekDataset();

		string strAvg = sArgs.gavg_input_arg;
		string strPres = sArgs.gpres_input_arg;
		vcd->ReadGeneAverage(strAvg);
		vcd->ReadGenePresence(strPres);
		vcd->InitializeGeneMap();

		CDataPair Dat;
		if(!Dat.Open(sArgs.dabinput_arg, false, false)){
			cerr << "Error opening dab file" << endl;
			return 1;
		}
		//fprintf(stderr, "Read.\n");
		//getchar();
	
		CSeekIntIntMap *gmap = vcd->GetGeneMap();
		vector<utype> veciallGenes;
		veciallGenes.clear(); //all genes
		veciallGenes.resize(vecstrGenes.size());
		for( i = 0; i < vecstrGenes.size( ); ++i ){
			veciallGenes[i] = mapstrintGene[vecstrGenes[i]];
			//fprintf(stderr, "Gene %s is id %d %d\n", vecstrGenes[i].c_str(), i, Dat.GetGene(vecstrGenes[i]));
		}
		//getchar();

		float** vall = ReadAllGenes(vecstrGenes, veciallGenes, gmap, vcd, Dat);
		//fprintf(stderr, "Finished loading all\n");

		vector< vector<string> > vx;
		vx.resize(vecstrList.size());

		//getchar();	
		#pragma omp parallel for \
		private(i) \
		shared(vx, vcd, vall, mapstrintGene, vecstrGenes, vecstrList) \
		schedule(dynamic)
		for(i=0; i<vecstrList.size(); i++){
			fprintf(stderr, "Doing: %s\n", vecstrList[i].c_str());
			vx[i] = Do_One(vecstrList[i].c_str(), rnd, vcd, vall, mapstrintGene, vecstrGenes);
		}

		for(i=0; i<vecstrList.size(); i++){
			printf("File: %s\n", vecstrList[i].c_str());
			for(j=0; j<vx[i].size(); j++){
				cout << vx[i][j] << endl;
			}
		}
			

	

    }

		
//#ifdef WIN32
//	pthread_win32_process_detach_np( );
//#endif // WIN32
	return 0; 

}

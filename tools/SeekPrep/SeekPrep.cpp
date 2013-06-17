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

bool InitializeDataset(size_t &iDatasets, vector<string> &vecstrDatasets,
	string &strPrepInputDirectory, vector<CSeekDataset*> &vc){
	vc.clear();
	vc.resize(iDatasets);
	size_t i, j, k;

	for(i=0; i<iDatasets; i++){
		vc[i] = new CSeekDataset();
		string strFileStem = vecstrDatasets[i];
		string strAvgPath = strPrepInputDirectory + "/" + strFileStem
			+ ".gavg";
		string strPresencePath = strPrepInputDirectory + "/" + strFileStem
			+ ".gpres";
		vc[i]->ReadGeneAverage(strAvgPath);
		vc[i]->ReadGenePresence(strPresencePath);
		vc[i]->InitializeGeneMap();
	}
	return true;
}


bool InitializeDB(size_t &iDatasets, size_t &iGenes,
	vector<string> &vecstrGenes, vector<CSeekDataset*> &vc, CDatabaselet &DBL){

	utype i,j,k;
	vector<char> cQuery;
	CSeekTools::InitVector(cQuery, iGenes, (char) 0);
	vector<string> allQuery;
	map<string, size_t> mapstrGenes;
	for(i=0; i<vecstrGenes.size(); i++){
		mapstrGenes[vecstrGenes[i]] = i;
	}

	/* Databaselet mapping */
	map<string, size_t> dbmap;
	vector<utype> veciQuery;
	for(i=0; i<DBL.GetGenes(); i++){
		string strQuery = DBL.GetGene(i);
		dbmap[strQuery] = i;
		allQuery.push_back(strQuery);
		/* global mapping */
		cQuery[mapstrGenes[strQuery]] = 1;
		veciQuery.push_back((utype) mapstrGenes[strQuery]);
	}

	//fprintf(stderr, "Start initializing query map...\n");
	for(i=0; i<iDatasets; i++){
		vc[i]->InitializeQueryBlock(veciQuery);
	}

	//fprintf(stderr, "Finished initializing map\n");

	//fprintf(stderr, "Start making gene-centric...\n");
	for(i=0; i<DBL.GetGenes(); i++){
		vector<unsigned char> Q;
		/* expanded */
		DBL.Get(i, Q);
		utype m = mapstrGenes[DBL.GetGene(i)];
		CSeekIntIntMap *qu = NULL;
		utype db;
		for(j=0; j<iDatasets; j++){
			if((qu=vc[j]->GetDBMap())==NULL) continue;
			if(CSeekTools::IsNaN(db=qu->GetForward(m))) continue;
			unsigned char **r = vc[j]->GetMatrix();
			for(k=0; k<iGenes; k++){
				unsigned char c = Q[k*iDatasets + j];
				r[db][k] = c;
			}
		}
	}

	//fprintf(stderr, "Finished making gene-centric\n");

	return true;
}

bool OpenDBFiles(string &DBFile, vector<unsigned char *> &cc, bool &useNibble){
	CDatabaselet *CD;
	(CD = new CDatabaselet(useNibble))->Open(DBFile);
	cc.push_back(CD->GetCharImage());
	return true;
}


bool OpenDB(string &DBFile, bool &useNibble, size_t &iDatasets,
	size_t &m_iGenes, vector<string> &vecstrGenes,
	map<utype, utype> &mapiPlatform, const vector<float> &quant,
	vector<CSeekDataset*> &vc, CFullMatrix<float> &platform_avg,
	CFullMatrix<float> &platform_stdev, vector<string> &vecstrQuery,
	const bool &logit){

	size_t i, j, k;

	CDatabaselet CD(useNibble);
	CD.Open(DBFile);
	InitializeDB(iDatasets, m_iGenes, vecstrGenes, vc, CD);

	vector<string> presGenes;
	for(i=0; i<CD.GetGenes(); i++)
		presGenes.push_back(CD.GetGene(i));

	size_t numPlatforms = platform_avg.GetRows();
	map<string, size_t> mapstriGenes;
	for(i=0; i<vecstrGenes.size(); i++)
		mapstriGenes[vecstrGenes[i]] = i;

	//fprintf(stderr, "Start calculating platform average\n");
	for(i=0; i<CD.GetGenes(); i++){
		vector<float> sum, sq_sum, mean, stdev;
		vector<int> num;
		sum.resize(numPlatforms);
		sq_sum.resize(numPlatforms);
		mean.resize(numPlatforms);
		stdev.resize(numPlatforms);
		num.resize(numPlatforms);

		for(k=0; k<numPlatforms; k++){
			sum[k] = 0;
			sq_sum[k] = 0;
			mean[k] = 0;
			stdev[k] = 0;
			num[k] = 0;
		}

		string thisGene = CD.GetGene(i);
		size_t geneID = mapstriGenes[thisGene];
		vecstrQuery.push_back(thisGene);
		for(k=0; k<iDatasets; k++){
			CSeekIntIntMap *mapQ = vc[k]->GetDBMap();
			if(mapQ==NULL) continue;
			unsigned char **f = vc[k]->GetMatrix();
			utype iQ = mapQ->GetForward(geneID);
			if(CSeekTools::IsNaN(iQ)){
				continue;
			}
			utype platform_id = mapiPlatform[k];
			if(platform_id>=(utype) numPlatforms){
				printf("Error, platforms are equal %d %d",
					(int) platform_id, (int) numPlatforms); getchar();
			}
			for(j=0; j<m_iGenes; j++){
				unsigned char uc = f[iQ][j];
				float v = 0;
				if(uc==255) v = CMeta::GetNaN();
				else{
					float vv = -1;
					if(logit) vv = log(quant[uc]) - log((float)(1.0-quant[uc]));
					else vv = quant[uc];
					v = vv - vc[k]->GetGeneAverage(j);
					/*if(isnan(vv) || isinf(vv) || isnan(vc[k]->GetGeneAverage(j)) ||
						isinf(vc[k]->GetGeneAverage(j))){
						fprintf(stderr, "%d %.5f %.5f %.5f\n", (int) uc, quant[uc], vv, vc[k]->GetGeneAverage(j));
					}*/
					//v = quant[uc];
					sum[platform_id] += v;
					num[platform_id]++;
					sq_sum[platform_id] += v*v;
				}
			}
		}

		for(k=0; k<numPlatforms; k++){
			if(num[k]==0) continue;
			mean[k] = sum[k] / (float) num[k];
			stdev[k] = sq_sum[k] / (float) num[k] - mean[k] * mean[k];
			stdev[k] = sqrt(stdev[k]);
			fprintf(stderr, "%.5f %.5f\n", mean[k], stdev[k]);
			platform_avg.Set(k, geneID, mean[k]);
			platform_stdev.Set(k, geneID, stdev[k]);
		}
	}

	//fprintf(stderr, "Finished calculating platform average\n");
	//fprintf(stderr, "Start deleting\n");

	for(i=0; i<iDatasets; i++)
		vc[i]->DeleteQueryBlock();
	//fprintf(stderr, "Finished deleting\n");
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
	vector<string>		vecstrLine, vecstrGenes, vecstrDBs, vecstrQuery;
	char				acBuffer[ c_iBuffer ];
	utype				i, j, k, l;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }

	/* reading gene-mapping file */
	if( sArgs.input_arg ) {
		ifsm.open( sArgs.input_arg );
		pistm = &ifsm; }
	else
		pistm = &cin;

	map<string, size_t> mapstriGenes;
	while( !pistm->eof( ) ) {
		pistm->getline( acBuffer, c_iBuffer - 1 );
		acBuffer[ c_iBuffer - 1 ] = 0;
		vecstrLine.clear( );
		CMeta::Tokenize( acBuffer, vecstrLine );
		if( vecstrLine.size( ) < 2 ) {
			cerr << "Ignoring line: " << acBuffer << endl;
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

	if( sArgs.input_arg ) ifsm.close( );

	omp_set_num_threads(1);
	int numThreads = omp_get_max_threads();

	/* PCL mode */
	if(sArgs.pclbin_flag==1){

		if(sArgs.sinfo_flag==1){
			string fileName = CMeta::Basename(sArgs.pclinput_arg);
			string fileStem = CMeta::Deextension(fileName);
			char outFile[125];
			sprintf(outFile, "%s/%s.sinfo", sArgs.dir_out_arg, fileStem.c_str());

			string pclfile = sArgs.pclinput_arg;
			if(!CMeta::IsExtension(pclfile, ".bin")){
				fprintf(stderr, "Input file is not bin type!\n");
				return 1;
			}
			
			CPCL pcl;
			if(!pcl.Open(pclfile.c_str())){
				fprintf(stderr, "Error opening file\n");
			}
			CMeasurePearNorm pn;
			CDat Dat;
			int numG = pcl.GetGeneNames().size();

			vector<int> veciGenes;
			veciGenes.resize(numG);
			for (j = 0; j < veciGenes.size(); ++j)
				veciGenes[j] = j;
			int iOne, iTwo;
			Dat.Open(pcl.GetGeneNames());

			for (k = 0; k < Dat.GetGenes(); ++k)
				for (j = (k + 1); j < Dat.GetGenes(); ++j)
					Dat.Set(k, j, CMeta::GetNaN());
			for (k = 0; k < numG; ++k) {
				if ((iOne = veciGenes[k]) == -1)
					continue;
				float *adOne = &pcl.Get(iOne)[0];
				for (j = (k + 1); j < numG; ++j){
					if ((iTwo = veciGenes[j]) != -1){
						float x = (float) pn.Measure(adOne,
							pcl.GetExperiments(), &pcl.Get(iTwo)[0],
							pcl.GetExperiments());
						Dat.Set(k, j, x);
						//fprintf(stderr, "%d %d %.5f\n", k, j, x);
					}
				}
			}
			
			double gmean, gstdev;
			size_t iN;

			Dat.AveStd(gmean, gstdev, iN);
			fprintf(stderr, "%.5f %.5f\n", (float) gmean, (float) gstdev);
			vector<float> vv;
			vv.resize(2);
			vv[0] = (float) gmean;
			vv[1] = (float) gstdev;
			CSeekTools::WriteArray(outFile, vv);

		}

		//if calculating gene variance per dataset
		else if(sArgs.gexpvarmean_flag==1){
			string fileName = CMeta::Basename(sArgs.pclinput_arg);
			string fileStem = CMeta::Deextension(fileName);
			char outFile[125];

			string pclfile = sArgs.pclinput_arg;
			if(!CMeta::IsExtension(pclfile, ".bin")){
				fprintf(stderr, "Input file is not bin type!\n");
				return 1;
			}
			
			CPCL pcl;
			if(!pcl.Open(pclfile.c_str())){
				fprintf(stderr, "Error opening file\n");
				return 1;
			}

			vector<float> var;
			vector<float> avg;

			CSeekTools::InitVector(var, vecstrGenes.size(), (float) CMeta::GetNaN());
			CSeekTools::InitVector(avg, vecstrGenes.size(), (float) CMeta::GetNaN());

			int totNumExperiments = pcl.GetExperiments();
			if(totNumExperiments<=2){
				fprintf(stderr, "This dataset is skipped because it contains <=2 columns\n");
				fprintf(stderr, "An empty vector will be returned\n");
			}else{
				for(j=0; j<vecstrGenes.size(); j++){
					utype g = pcl.GetGene(vecstrGenes[j]);
					if(CSeekTools::IsNaN(g)) continue; //gene does not exist in the dataset
					float *val = pcl.Get(g);
					vector<float> rowVal;
					for(k=0; k<totNumExperiments; k++)
						rowVal.push_back(val[k]);
					float mean = 0;
					float variance = 0;
					for(k=0; k<rowVal.size(); k++)
						mean+=rowVal[k];
					mean/=rowVal.size();
					for(k=0; k<rowVal.size(); k++)
						variance += (rowVal[k] - mean) * (rowVal[k] - mean);
					variance /= rowVal.size();
					var[j] = variance;
					avg[j] = mean;
					//fprintf(stderr, "%.5f %.5f\n", mean, variance);
				}
				//fprintf(stderr, "done\n"); 
			}
			//fprintf(stderr, "G\n"); 
			sprintf(outFile, "%s/%s.gexpvar", sArgs.dir_out_arg, fileStem.c_str());
			CSeekTools::WriteArray(outFile, var);
			sprintf(outFile, "%s/%s.gexpmean", sArgs.dir_out_arg, fileStem.c_str());
			CSeekTools::WriteArray(outFile, avg);
		}

	}


	/* DB mode */
	if(sArgs.db_flag==1){

		if(!sArgs.quant_arg){
			fprintf(stderr, "Must give quant file\n");
			return -1;
		}

		vector<float> quant;
		string strQuantFile = sArgs.quant_arg;
		//fprintf(stderr, "%s\n", strQuantFile.c_str());
		CSeekTools::ReadQuantFile(strQuantFile, quant);

		//fprintf(stderr, "quant %d %.5f\n", 170, quant[170]);
		//getchar();

		if(sArgs.gplat_flag==1){
			bool logit = false;
			if(sArgs.logit_flag==1) logit = true;

			vector<CSeekPlatform> vp;
			map<string, utype> mapstriPlatform;
			map<utype, string> mapistrPlatform;
			map<utype, utype> mapiPlatform;
			vector<string> vecstrPlatforms;

			vector<string> vecstrDatasets;

			ifsm.open(sArgs.dset_arg);
			i = 0;
			while(!pistm->eof()){
				pistm->getline(acBuffer, c_iBuffer -1);
				if(acBuffer[0]==0) break;
				acBuffer[c_iBuffer-1] = 0;
				vecstrLine.clear();
				CMeta::Tokenize( acBuffer, vecstrLine );
				/* read dataset name */
				vecstrDatasets.push_back(vecstrLine[0]);
				/* just read the platform information */
				string pl = vecstrLine[1];
				map< string, utype >::const_iterator iter;
				iter = mapstriPlatform.find(pl);
				if(iter== mapstriPlatform.end()){
					utype s = mapstriPlatform.size();
					mapstriPlatform[pl] = s;
					mapistrPlatform[s] = pl;
					vecstrPlatforms.push_back(pl);
				}
				utype platform_id = mapstriPlatform[pl];
				mapiPlatform[i] = platform_id;
				i++;
			}
			ifsm.close();

			vector<string> dblist;
			ifsm.open(sArgs.dblist_arg);
			i = 0;
			while(!pistm->eof()){
				pistm->getline(acBuffer, c_iBuffer -1);
				if(acBuffer[0]==0) break;
				acBuffer[c_iBuffer-1] = 0;
				dblist.push_back(acBuffer);
			}
			dblist.resize(dblist.size());
			ifsm.close();

			bool useNibble = false;
			if(sArgs.useNibble_flag==1) useNibble = true;

			size_t numPlatforms = mapstriPlatform.size();
			size_t iDatasets = vecstrDatasets.size();
			size_t m_iGenes = vecstrGenes.size();
			CFullMatrix<float> platform_avg, platform_stdev;
			platform_avg.Initialize(numPlatforms, m_iGenes);
			platform_stdev.Initialize(numPlatforms, m_iGenes);

			//printf("Size: %d %d\n", numPlatforms, m_iGenes); getchar();

			for(i=0; i<numPlatforms; i++){
				for(j=0; j<m_iGenes; j++){
					platform_avg.Set(i, j, CMeta::GetNaN());
					platform_stdev.Set(i, j, CMeta::GetNaN());
				}
			}

			/*if(iDatasets<numThreads){
				numThreads = iDatasets;
				omp_set_num_threads(numThreads);
			}*/

			string strPrepInputDirectory = sArgs.dir_prep_in_arg;
			vector<CSeekDataset*> *vc = new vector<CSeekDataset*>[numThreads];
			CFullMatrix<float> *platform_avg_threads =
				new CFullMatrix<float>[numThreads];
			CFullMatrix<float> *platform_stdev_threads=
				new CFullMatrix<float>[numThreads];

			for(i=0; i<numThreads; i++){
				InitializeDataset(iDatasets, vecstrDatasets,
					strPrepInputDirectory, vc[i]);
				platform_avg_threads[i].Initialize(numPlatforms, m_iGenes);
				platform_stdev_threads[i].Initialize(numPlatforms, m_iGenes);
				for(j=0; j<numPlatforms; j++){
					for(k=0; k<m_iGenes; k++){
						platform_avg_threads[i].Set(j, k, CMeta::GetNaN());
						platform_stdev_threads[i].Set(j, k, CMeta::GetNaN());
					}
				}
			}

			//printf("Dataset initialized"); getchar();
			vector<string> vecstrQuery;

			//#pragma omp parallel for \
			shared(vc, dblist, iDatasets, m_iGenes, vecstrGenes, mapiPlatform, quant, \
			platform_avg_threads, platform_stdev_threads, vecstrQuery, logit) \
			private(i) firstprivate(useNibble) schedule(dynamic)
			for(i=0; i<dblist.size(); i++){
				int tid = omp_get_thread_num();
				string DBFile = dblist[i];
				fprintf(stderr, "opening db file %s\n", DBFile.c_str());
				OpenDB(DBFile, useNibble, iDatasets, m_iGenes,
				vecstrGenes, mapiPlatform, quant, vc[tid],
				platform_avg_threads[tid], platform_stdev_threads[tid],
				vecstrQuery, logit);
				fprintf(stderr, "finished opening db file %s\n",
					DBFile.c_str());
			}

			for(i=0; i<numThreads; i++){
				for(j=0; j<numPlatforms; j++){
					for(k=0; k<m_iGenes; k++){
						float ca = platform_avg_threads[i].Get(j, k);
						float cs = platform_stdev_threads[i].Get(j, k);

						if(ca==CMeta::GetNaN() || cs==CMeta::GetNaN()){
							continue;
						}
						platform_avg.Set(j, k, ca);
						platform_stdev.Set(j, k, cs);
					}
				}
			}

			//for(i=0; i<numPlatforms; i++){
				//printf("Platform %s\n", mapistrPlatform[i].c_str());
				/*for(j=0; j<vecstrQuery.size(); j++){
					size_t iGene = mapstriGenes[vecstrQuery[j]];
					printf("Gene %s %.5f %.5f\n", vecstrQuery[j].c_str(), platform_avg.Get(i, iGene),
						platform_stdev.Get(i,iGene));
				}*/
			//}

			char outFile[125];
			sprintf(outFile, "%s/all_platforms.gplatavg", sArgs.dir_out_arg);
			platform_avg.Save(outFile);
			sprintf(outFile, "%s/all_platforms.gplatstdev", sArgs.dir_out_arg);
			platform_stdev.Save(outFile);
			sprintf(outFile, "%s/all_platforms.gplatorder", sArgs.dir_out_arg);
			ofstream outfile;
			outfile.open(outFile);
			for(i=0; i<vecstrPlatforms.size(); i++){
				outfile << vecstrPlatforms[i] << "\n";
			}
			outfile.close();
		}

	} else if(sArgs.dab_flag==1){
		
		if(sArgs.norm_flag==1){
			CDataPair Dat;
			char outFile[1024];
			fprintf(stderr, "Opening file...\n");
			//if(!Dat.Open(sArgs.dabinput_arg, false, 2, false, false, true)){
			if(!Dat.Open(sArgs.dabinput_arg, false, false, 2, false, false)){
				cerr << "error opening file" << endl;
				return 1;
			}
			fprintf(stderr, "Finished opening file\n");
			string fileName = CMeta::Basename(sArgs.dabinput_arg);
			string fileStem = CMeta::Deextension(fileName);
			sprintf(outFile, "%s/%s.2.dab", sArgs.dir_out_arg,
				fileStem.c_str());
			int max_rank = 3000;
			float rbp_p = 0.99;
			//cutoff, expTransform, divideNorm, subtractNorm
			//CSeekWriter::NormalizeDAB(Dat, vecstrGenes, true, false, true, false);
			//CSeekWriter::RankNormalizeDAB(Dat, vecstrGenes, max_rank, rbp_p);
			//Dat.Save(outFile);
			vector<map<utype,unsigned short> > umat;
			CSeekWriter::GetSparseRankMatrix(Dat, umat, max_rank, vecstrGenes);
			CSeekWriter::WriteSparseMatrix(Dat, umat, max_rank, vecstrGenes, outFile);
			/*fprintf(stderr, "Begin\n");
			vector<unsigned short> l;
			vector<vector<float> > mat;
			CSeekWriter::ReadSparseMatrixAsArray(l, outFile);
			fprintf(stderr, "Begin 2\n");
			CSeekWriter::ReadSparseMatrix(l, mat, 0.99, vecstrGenes);*/
		}

		if(sArgs.gavg_flag==1){
			bool logit = false;
			if(sArgs.logit_flag==1) logit = true;

			CDataPair Dat;
			char outFile[125];
			if(!Dat.Open(sArgs.dabinput_arg, false, false)){
				cerr << "error opening file" << endl;
				return 1;
			}
			vector<float> vecGeneAvg;
			string fileName = CMeta::Basename(sArgs.dabinput_arg);
			string fileStem = CMeta::Deextension(fileName);
			sprintf(outFile, "%s/%s.gavg", sArgs.dir_out_arg,
				fileStem.c_str());
			CSeekWriter::GetGeneAverage(Dat, vecstrGenes, vecGeneAvg, logit, sArgs.top_avg_percent_arg);

			//DEBUGGING
			for(i=0; i<vecGeneAvg.size(); i++){
				fprintf(stderr, "%s\t%.3f\n", vecstrGenes[i].c_str(), vecGeneAvg[i]);
			}

			CSeekTools::WriteArray(outFile, vecGeneAvg);
		}

		else if(sArgs.gpres_flag==1){
			CDataPair Dat;
			char outFile[125];
			if(!Dat.Open(sArgs.dabinput_arg, false, false)){
				cerr << "error opening file" << endl;
				return 1;
			}
			vector<char> vecGenePresence;
			string fileName = CMeta::Basename(sArgs.dabinput_arg);
			string fileStem = CMeta::Deextension(fileName);
			sprintf(outFile, "%s/%s.gpres", sArgs.dir_out_arg,
				fileStem.c_str());
			CSeekWriter::GetGenePresence(Dat, vecstrGenes, vecGenePresence);

			//DEBUGGING
			for(i=0; i<vecGenePresence.size(); i++){
				fprintf(stderr, "%s\t%d\n", vecstrGenes[i].c_str(), vecGenePresence[i]);
			}

			CSeekTools::WriteArray(outFile, vecGenePresence);
		}

	}
#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; }

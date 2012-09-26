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
	size_t				i;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }

	vector<string> vecstrGeneID;
	map<string, ushort> mapstrintGene;
	if(!CSeekTools::ReadListTwoColumns(sArgs.input_arg, vecstrGeneID, vecstrGenes))
		return false;

	for(i=0; i<vecstrGenes.size(); i++)
		mapstrintGene[vecstrGenes[i]] = i;

	if(sArgs.dataset_flag==1){
		vector<string> vecstrDP, vecstrUserDP;
		//dataset-platform mapping (required)
		if(!CSeekTools::ReadListTwoColumns(sArgs.db_arg, vecstrDatasets, vecstrDP))
			return false;

		//dataset filter
		if(!CSeekTools::ReadListTwoColumns(sArgs.dset_list_arg, vecstrUserDatasets, vecstrUserDP))
			return false;

		map<string, ushort> mapstrintDataset;
		map<string, string> mapstrstrDatasetPlatform;
		ushort i;
		for(i=0; i<vecstrDatasets.size(); i++){
			mapstrstrDatasetPlatform[vecstrDatasets[i]] = vecstrDP[i];
			mapstrintDataset[vecstrDatasets[i]] = i;
		}

		CSeekIntIntMap *mapUserDatasets = new CSeekIntIntMap(vecstrDatasets.size());
		for(i=0; i<vecstrUserDatasets.size(); i++){
			mapUserDatasets->Add(mapstrintDataset[vecstrUserDatasets[i]]);
		}
		
		//fprintf(stderr, "Finished reading dataset\n");

		size_t iDatasets = vecstrDatasets.size();
		vector<string> vecstrPlatforms;
		vector<CSeekPlatform> vp;
		map<string, ushort> mapstriPlatform;
		CSeekTools::ReadPlatforms(sArgs.platform_dir_arg, vp, vecstrPlatforms,
			mapstriPlatform);
		//fprintf(stderr, "Finished reading platform\n");

		vector<CSeekDataset*> vc;
		vc.resize(iDatasets);
		string strPrepInputDirectory = sArgs.dir_prep_in_arg;
		string strGvarInputDirectory = sArgs.dir_gvar_in_arg;
		string strSinfoInputDirectory = sArgs.dir_sinfo_in_arg;

		for(i=0; i<iDatasets; i++){
			vc[i] = new CSeekDataset();
			string strFileStem = vecstrDatasets[i];
			string strAvgPath = strPrepInputDirectory + "/" +
				strFileStem + ".gavg";
			string strPresencePath = strPrepInputDirectory + "/" +
				strFileStem + ".gpres";

			if(strGvarInputDirectory!="NA"){
				string strGvarPath = strGvarInputDirectory + "/" + 
					strFileStem + ".gexpvar";
				vc[i]->ReadGeneVariance(strGvarPath);
			}
			if(strSinfoInputDirectory!="NA"){
				string strSinfoPath = strSinfoInputDirectory + "/" + 
					strFileStem + ".sinfo";
				vc[i]->ReadDatasetAverageStdev(strSinfoPath);
			}

			vc[i]->ReadGeneAverage(strAvgPath);
			vc[i]->ReadGenePresence(strPresencePath);
			string strPlatform =
				mapstrstrDatasetPlatform.find(strFileStem)->second;
			ushort platform_id = mapstriPlatform.find(strPlatform)->second;
			vc[i]->SetPlatform(vp[platform_id]);
		}

			
		//fprintf(stderr, "Finished reading prep\n");

		for(i=0; i<iDatasets; i++) vc[i]->InitializeGeneMap();

		vector<vector<string> > vecstrQueries;
		string multiQuery = sArgs.multi_query_arg;

		if(multiQuery=="NA"){
			vecstrQueries.resize(1);
			if(!CSeekTools::ReadMultiGeneOneLine(sArgs.single_query_arg, vecstrQueries[0]))
				return false;
		}else{	
			if(!CSeekTools::ReadMultipleQueries(multiQuery, vecstrQueries))
				return false;
		}

		bool toOutput = false;
		string output_file = sArgs.output_file_arg;
		FILE *out = NULL;

		if(output_file=="NA"){
			toOutput = false;
		}else{
			toOutput = true;
			out = fopen(output_file.c_str(), "w");
		}

		ushort ii=0;
		float cutoff = sArgs.gvar_cutoff_arg;
		bool toCutoff = false;
		if(cutoff<0){
			toCutoff = false;
		}else{
			toCutoff = true;
		}

		if(sArgs.order_stat_single_gene_query_flag==1){
			if(strGvarInputDirectory=="NA"){
				fprintf(stderr, "Order statistics mode, but need to provide gvar!\n");
				return false;
			}
			if(cutoff<0){
				fprintf(stderr, "Need to provide positive Gvar cutoff\n");
				return false;
			}
		}

		for(ii=0; ii<vecstrQueries.size(); ii++){

		vector<string> &vecstrQuery = vecstrQueries[ii];

		//fprintf(stderr, "Finished reading query\n");
		bool isFirst = true;
		vector<int> count;
		ushort j;
		CSeekTools::InitVector(count, vecstrQuery.size(), (int) 0);

		//only analyze and report on user datasets
		const vector<ushort> &allRDatasets = mapUserDatasets->GetAllReverse();
		ushort iUserDatasets = mapUserDatasets->GetNumSet();
		ushort dd;
		for(dd=0; dd<iUserDatasets; dd++){
			i = allRDatasets[dd];
			//fprintf(stdout, "%d %d %s\n", i, dd, vecstrDatasets[i].c_str());
			CSeekIntIntMap *si = vc[i]->GetGeneMap();

			if(toCutoff && sArgs.order_stat_single_gene_query_flag==1){
				if(mapstrintGene.find(vecstrQuery[0])==mapstrintGene.end())
					continue;	
				if(CSeekTools::IsNaN(si->GetForward(
					mapstrintGene[vecstrQuery[0]]))) continue;
				float gene_var = vc[i]->GetGeneVariance(mapstrintGene[vecstrQuery[0]]);
				if(gene_var < cutoff) continue;
				if(isFirst){
					isFirst = false;
					if(toOutput){
						fprintf(out, "%s", vecstrDatasets[i].c_str());
					}else{
						fprintf(stdout, "%s", vecstrDatasets[i].c_str());
					}
				}else{
					if(toOutput){
						fprintf(out, " %s", vecstrDatasets[i].c_str());
					}else{
						fprintf(stdout, " %s", vecstrDatasets[i].c_str());
					}
				}

			}else{
				ushort present = 0;
				for(j=0, present=0; j<vecstrQuery.size(); j++){
					if(mapstrintGene.find(vecstrQuery[j])==mapstrintGene.end())
						continue;
					if(CSeekTools::IsNaN(si->GetForward(
						mapstrintGene[vecstrQuery[j]]))) continue;
					count[j]++;
					present++;
				}
				if(present==vecstrQuery.size()){
					if(isFirst){
						isFirst = false;
						if(toOutput){
							fprintf(out, "%s", vecstrDatasets[i].c_str());
						}else{
							fprintf(stdout, "%s", vecstrDatasets[i].c_str());
						}
					}else{
						if(toOutput){
							fprintf(out, " %s", vecstrDatasets[i].c_str());
						}else{
							fprintf(stdout, " %s", vecstrDatasets[i].c_str());
						}
					}
					//fprintf(stderr, "%s\t%s\t%d\t%d\n", 
					//	vecstrDatasets[i].c_str(), vecstrDP[i].c_str(), 
					//	present, si->GetNumSet());
				}
			}
		}

		for(j=0; j<vecstrQuery.size(); j++)
			fprintf(stderr, "Gene %s: %d\n", vecstrQuery[j].c_str(), count[j]);

		if(toOutput){
			fprintf(out, "\n");
		}else{
			fprintf(stdout, "\n");
		}

		}

		if(toOutput){
			fclose(out);
		}		

	}
	else if(sArgs.databaselet_flag==1){
		bool useNibble = false;
		if(sArgs.is_nibble_flag==1) useNibble = true;
		CDatabase DB(useNibble);
		vector<string> vecstrDP;
		vector<string> vecstrQuery;
		if(!CSeekTools::ReadListTwoColumns(sArgs.db_arg, vecstrDatasets, vecstrDP))
			return false;
		if(!CSeekTools::ReadMultiGeneOneLine(sArgs.single_query_arg, vecstrQuery))
			return false;

		string strInputDirectory = sArgs.dir_in_arg;
		DB.Open(strInputDirectory);

		size_t iDatasets = DB.GetDatasets();
		size_t iGenes = DB.GetGenes();

		size_t j,k;
		vector<CSeekDataset*> vc;
		vc.clear();
		vc.resize(iDatasets);
		for(i=0; i<iDatasets; i++){
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
		vector<ushort> allQ;
		CSeekTools::InitVector(cQuery, iGenes, (char) 0);

		for(i=0; i<vecstrQuery.size(); i++){
			k = DB.GetGene(vecstrQuery[i]);
			if(k==-1) continue;
			cQuery[k] = 1;
		}

		for(k=0; k<iGenes; k++)
			if(cQuery[k]==1)
				allQ.push_back(k);

		for(i=0; i<iDatasets; i++){
			vc[i]->InitializeGeneMap();
			vc[i]->InitializeQueryBlock(allQ);
		}

		vector<unsigned char> *Q =
			new vector<unsigned char>[vecstrQuery.size()];

		for(i=0; i<vecstrQuery.size(); i++)
			if(!DB.GetGene(vecstrQuery[i], Q[i]))
				cerr << "Gene does not exist" << endl;

		//printf("Before"); getchar();
		for(i=0; i<vecstrQuery.size(); i++){
			if(DB.GetGene(vecstrQuery[i])==-1) continue;
			size_t m = DB.GetGene(vecstrQuery[i]);
			size_t l = 0;
			for(j=0; j<iDatasets; j++){
				CSeekIntIntMap *qu = vc[j]->GetDBMap();
				if(qu==NULL) continue;
				unsigned char **r = vc[j]->GetMatrix();
				ushort query = qu->GetForward(m);
				if(CSeekTools::IsNaN(query)) continue;
			    for(k=0; k<iGenes; k++){
			    	unsigned char c = Q[i][k*iDatasets + j];
			    	r[query][k] = c;
			    }
			}
		}

		for(i=0; i<iDatasets; i++){
			printf("Dataset %ld\n", i);
			unsigned char **r = vc[i]->GetMatrix();
			CSeekIntIntMap *mapG = vc[i]->GetGeneMap();
			CSeekIntIntMap *mapDB = vc[i]->GetDBMap();
			if(mapDB==NULL) continue;
			for(j=0; j<mapDB->GetNumSet(); j++){
				if(vecstrDatasets[i]=="GSE19470.GPL5175.pcl"){
				printf("Row %ld\n", j);
				for(k=0; k<mapG->GetNumSet(); k++){
					printf("%d ", r[j][mapG->GetReverse(k)]);
				}
				printf("\n");
				}
			}
		}


	}else{
		cerr << "Must give a db list." << endl;
		return 1;

	}

#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; }

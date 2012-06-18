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
	vector<string>		vecstrLine, vecstrGenes, vecstrDatasets, vecstrQuery;
	char				acBuffer[ c_iBuffer ];
	size_t				i;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }

	if( sArgs.input_arg ) {
		ifsm.open( sArgs.input_arg );
		pistm = &ifsm; }
	else
		pistm = &cin;
	while( !pistm->eof( ) ) {
		pistm->getline( acBuffer, c_iBuffer - 1 );
		acBuffer[ c_iBuffer - 1 ] = 0;
		vecstrLine.clear( );
		CMeta::Tokenize( acBuffer, vecstrLine );
		if( vecstrLine.size( ) < 2 ) {
			cerr << "Ignoring line: " << acBuffer << endl;
			continue; }
		if( !( i = atoi( vecstrLine[ 0 ].c_str( ) ) ) ) {
			cerr << "Illegal gene ID: " << vecstrLine[ 0 ] << " for " << vecstrLine[ 1 ] << endl;
			return 1; }
		i--;
		if( vecstrGenes.size( ) <= i )
			vecstrGenes.resize( i + 1 );
		vecstrGenes[ i ] = vecstrLine[ 1 ]; }
	if( sArgs.input_arg )
		ifsm.close( );

	bool useNibble = false;
	if(sArgs.is_nibble_flag==1){
		useNibble = true;
	}

	CDatabase DB(useNibble);

	if(sArgs.db_arg){
		ifsm.open(sArgs.db_arg);
		while(!pistm->eof()){
			pistm->getline(acBuffer, c_iBuffer -1);
			if(acBuffer[0]==0){
				break;
			}
			acBuffer[c_iBuffer-1] = 0;
			vecstrDatasets.push_back(acBuffer);
		}
		vecstrDatasets.resize(vecstrDatasets.size());
		ifsm.close();

		ifsm.open(sArgs.query_arg);
		while(!pistm->eof()){
			pistm->getline(acBuffer, c_iBuffer -1);
			if(acBuffer[0]==0){
				break;
			}
			acBuffer[c_iBuffer-1] = 0;
			vecstrQuery.push_back(acBuffer);
		}
		vecstrQuery.resize(vecstrQuery.size());
		ifsm.close();

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
			//string strFileStem = CMeta::Deextension(CMeta::Basename(vecstrDatasets[i].c_str()));
			string strAvgPath = strPrepInputDirectory + "/" + strFileStem + ".gavg";
			string strPresencePath = strPrepInputDirectory + "/" + strFileStem + ".gpres";
			vc[i]->ReadGeneAverage(strAvgPath);
			vc[i]->ReadGenePresence(strPresencePath);
		}

		vector<char> cQuery;
		CSeekTools::InitVector(cQuery, iGenes, (char) 0);

		for(i=0; i<vecstrQuery.size(); i++){
			k = DB.GetGene(vecstrQuery[i]);
			if(k==-1) continue;
			cQuery[k] = 1;
		}

		for(i=0; i<iDatasets; i++){
			vc[i]->InitializeQuery(cQuery);
		}

		vector<unsigned char> *Q =
			new vector<unsigned char>[vecstrQuery.size()];

		for(i=0; i<vecstrQuery.size(); i++){
			if(!DB.GetGene(vecstrQuery[i], Q[i])){
				cerr << "Gene does not exist" << endl;
			}
		}

		//printf("Before"); getchar();
		for(i=0; i<vecstrQuery.size(); i++){
			if(DB.GetGene(vecstrQuery[i])==-1){
				continue;
			}
			size_t m = DB.GetGene(vecstrQuery[i]);
			size_t l = 0;
			for(j=0; j<iDatasets; j++){
				CSeekIntIntMap *qu = vc[j]->GetQueryMap();
			    size_t query = qu->GetForward(m);
			    if(query==-1) continue;
			    for(k=0; k<iGenes; k++){
			    	unsigned char c = Q[i][k*iDatasets + j];
			    	vc[j]->SetQueryNoMapping(query, k, c);
			    }
			}
		}

		/*for(i=0; i<iDatasets; i++){
			printf("Dataset %ld\n", i);
			CSeekMatrix<unsigned char> *cm = vc[i]->GetMatrix();
			for(j=0; j<cm->GetNumRow(); j++){
				printf("Row %ld\n", j);
				for(k=0; k<1000; k++){
					printf("%d ", cm->Get(j, k));
				}
				printf("\n");
			}
		}*/
		/*size_t j;
		for(i=0; i<vecstrQuery.size(); i++){
			printf("Query: %s\n", vecstrQuery[i].c_str());
			for(j=0; j<Q[i].size(); j++){
				printf("%d ", (int) Q[i][j]);
			}
			printf("\n");
			getchar();
		}*/

		//printf("Done"); getchar();

	}else{
		cerr << "Must give a db list." << endl;
		return 1;

	}

#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; }

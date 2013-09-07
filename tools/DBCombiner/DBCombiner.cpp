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
	vector<string>		vecstrLine, vecstrGenes, vecstrDBs;
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
			cerr << "Illegal gene ID: " << vecstrLine[ 0 ] <<
				" for " << vecstrLine[ 1 ] << endl;
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

	if(sArgs.reorganize_flag==1){
		vector<string> vecstrDataset;
		ifstream ifsm2;
		ifsm2.open(sArgs.dataset_arg);
		while(!ifsm2.eof()){
			ifsm2.getline(acBuffer, c_iBuffer-1);
			if(acBuffer[0]==0) break;
			acBuffer[c_iBuffer-1] = 0;
			vector<string> vecstrLine;
			CMeta::Tokenize(acBuffer, vecstrLine);
			vecstrDataset.push_back(vecstrLine[0]);
		}
		ifsm2.close();

		CDatabase db(false);
		db.Open(sArgs.db_dir_arg, vecstrGenes, vecstrDataset.size(), 
			sArgs.src_db_num_arg);
		db.Reorganize(sArgs.dest_db_dir_arg, sArgs.dest_db_num_arg);
		return 0;
	}

	if(sArgs.combine_flag==1){
		CDatabase DB(useNibble);

		bool fSplit = false;
		if(sArgs.split_flag==1){
			fSplit = true;
		}

		if(sArgs.db_arg){
			ifsm.open(sArgs.db_arg);
			while(!pistm->eof()){
				pistm->getline(acBuffer, c_iBuffer -1);
				if(acBuffer[0]==0){
					break;
				}
				acBuffer[c_iBuffer-1] = 0;
				vecstrDBs.push_back(acBuffer);
			}
			vecstrDBs.resize(vecstrDBs.size());
			ifsm.close();

			//printf("Reading DBS"); getchar();
			vector<CDatabaselet*> DBS;
			DBS.resize(vecstrDBs.size());
			for(i=0; i<vecstrDBs.size(); i++){
				DBS[i] = new CDatabaselet(useNibble);
				DBS[i]->Open(vecstrDBs[i]);
			}
			//printf("Finished reading DBS"); getchar();

			CDatabaselet::Combine(DBS, sArgs.dir_out_arg, vecstrGenes, fSplit);
			for(i=0; i<vecstrDBs.size(); i++){
				free(DBS[i]);
			}

		}else{
			cerr << "Must give a db list." << endl;
			return 1;

		}
	}
#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; }

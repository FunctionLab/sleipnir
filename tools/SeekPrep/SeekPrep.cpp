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
	vector<string>		vecstrLine, vecstrGenes, vecstrDBs, vecstrQuery;
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

	if(sArgs.dab_arg){
		CDataPair Dat;
		char outFile[125];
		if(!Dat.Open(sArgs.dab_arg, false, false)){
			cerr << "error opening file" << endl;
			return 1;
		}

		if(sArgs.gavg_flag==1){
			vector<float> vecGeneAvg;
			string fileName = CMeta::Basename(sArgs.dab_arg);
			string fileStem = CMeta::Deextension(fileName);
			sprintf(outFile, "%s/%s.gavg", sArgs.dir_out_arg, fileStem.c_str());
			CSeekWriter::GetGeneAverage(Dat, vecstrGenes, vecGeneAvg);
			CSeekTools::WriteArray(outFile, vecGeneAvg);
		}

		if(sArgs.gpres_flag==1){
			vector<char> vecGenePresence;
			string fileName = CMeta::Basename(sArgs.dab_arg);
			string fileStem = CMeta::Deextension(fileName);
			sprintf(outFile, "%s/%s.gpres", sArgs.dir_out_arg, fileStem.c_str());
			CSeekWriter::GetGenePresence(Dat, vecstrGenes, vecGenePresence);
			CSeekTools::WriteArray(outFile, vecGenePresence);
		}

	}else{
		cerr << "Must give a dab." << endl;
		return 1;

	}

#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; }

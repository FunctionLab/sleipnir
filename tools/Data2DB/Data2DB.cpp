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
	vector<string>		vecstrLine, vecstrGenes;
	char				acBuffer[ c_iBuffer ];
	CDatabase			DB;
	CBayesNetSmile		BNSmile;
	size_t				i;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );
#ifdef SMILEXML_LIB
	EnableXdslFormat( );
#endif

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

	if( !BNSmile.Open( sArgs.network_arg ) ) {
		cerr << "Could not open: " << sArgs.network_arg << endl;
		return 1; }
	DB.SetMemmap( !!sArgs.memmap_flag );
	DB.SetBuffer( !!sArgs.buffer_flag );
	DB.SetBlockOut( sArgs.block_files_arg );
	DB.SetBlockIn( sArgs.block_datasets_arg );
	if( !DB.Open( vecstrGenes, sArgs.dir_in_arg, &BNSmile, sArgs.dir_out_arg, min((size_t)sArgs.files_arg,
		vecstrGenes.size( )) ) ) {
		cerr << "Could not open data" << endl;
		return 1; }

	CMeta::Shutdown( );
#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; }

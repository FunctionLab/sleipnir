#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CDatasetCompact		Data;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( sArgs.load_arg ) {
		CDatasetCompactMap	DataMap;

		if( !DataMap.Open( sArgs.load_arg ) ) {
			cerr << "Couldn't open: " << sArgs.load_arg << endl;
			return 1; }
#ifdef _MSC_VER
		Sleep( INFINITE );
#else // _MSC_VER
		sleep( -1 );
#endif // _MSC_VER
		return 0; }

	if( sArgs.input_arg ) {
		if( !Data.Open( ifstream( sArgs.input_arg, ios_base::binary ) ) ) {
			cerr << "Couldn't open: " << sArgs.input_arg << endl;
			return 1; } }
	else if( sArgs.network_arg ) {
		CDataPair		Answers;
		CBayesNetSmile	BNSmile;

		if( !sArgs.answers_arg ) {
			cmdline_parser_print_help( );
			return 1; }
		if( !BNSmile.Open( sArgs.network_arg ) ) {
			cerr << "Couldn't open: " << sArgs.network_arg << endl;
			return 1; }
		if( !Answers.Open( sArgs.answers_arg, false ) ) {
			cerr << "Couldn't open: " << sArgs.answers_arg << endl;
			return 1; }
		if( !Data.Open( Answers, sArgs.directory_arg, &BNSmile ) ) {
			cerr << "Couldn't open: " << sArgs.directory_arg << endl;
			return 1; } }
	else {
		vector<string>	vecstrFiles;

		vecstrFiles.resize( sArgs.inputs_num );
		copy( sArgs.inputs, sArgs.inputs + sArgs.inputs_num, vecstrFiles.begin( ) );
		if( !Data.Open( vecstrFiles ) ) {
			cerr << "Couldn't open inputs" << endl;
			return 1; } }

	if( sArgs.output_arg )
		Data.Save( ofstream( sArgs.output_arg, ios_base::binary ), true );
	else
		Data.Save( cout, false );

	CMeta::Shutdown( );
	return 0; }

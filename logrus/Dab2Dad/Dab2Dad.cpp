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
	}
	else if( sArgs.input_arg ) {
		ifstream	ifsm;

		ifsm.open( sArgs.input_arg, ios_base::binary );
		if( !Data.Open( ifsm ) ) {
			cerr << "Couldn't open: " << sArgs.input_arg << endl;
			return 1; }
		Data.Save( cout, false );
		cout.flush( ); }
	else {
		vector<string>	vecstrFiles;

		vecstrFiles.resize( sArgs.inputs_num );
		copy( sArgs.inputs, sArgs.inputs + sArgs.inputs_num, vecstrFiles.begin( ) );
		if( !Data.Open( vecstrFiles ) ) {
			cerr << "Couldn't open inputs" << endl;
			return 1; }
		if( sArgs.output_arg ) {
			ofstream	ofsm;

			ofsm.open( sArgs.output_arg, ios_base::binary );
			Data.Save( ofsm, true );
			ofsm.close( ); }
		else
			Data.Save( cout, false ); }

	CMeta::Shutdown( );
	return 0; }

#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	COrthology			Orthology;
	ifstream			ifsm;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	ifsm.open( sArgs.input_arg );
	if( !Orthology.Open( ifsm ) ) {
		cerr << "Could not open: " << sArgs.input_arg << endl;
		return 1; }
	ifsm.close( );

	CMeta::Shutdown( );
	return 0; }

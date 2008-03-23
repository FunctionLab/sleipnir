#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CBayesNetSmile		BNSmile;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );
#if !( defined(_MSC_VER) && defined(_DEBUG) )
#ifdef SMILEXML_LIB
	EnableXdslFormat( );
#endif
#endif // !( defined(_MSC_VER) && defined(_DEBUG) )

	if( !BNSmile.Open( sArgs.input_arg ) ) {
		cerr << "Couldn't open: " << sArgs.input_arg << endl;
		return 1; }
	BNSmile.Save( sArgs.output_arg );

	CMeta::Shutdown( );
	return 0; }

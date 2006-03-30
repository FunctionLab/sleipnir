#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	CDat				Dat;
	size_t				i, j;
	ofstream			ofsm;
	gengetopt_args_info	sArgs;
	float				d;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }

	CMeta::Startup( sArgs.verbosity_arg );
	if( !Dat.Open( sArgs.input_arg ) ) {
		cerr << "Could not open input: " << sArgs.input_arg << endl;
		return 1; }

	Dat.Normalize( !sArgs.zscore_flag );
	if( sArgs.flip_flag )
		for( i = 0; i < Dat.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
				if( !CMeta::IsNaN( d = Dat.Get( i, j ) ) )
					Dat.Set( i, j, 1 - d );

	ofsm.open( sArgs.output_arg ? sArgs.output_arg : sArgs.input_arg, ios_base::binary );
	Dat.Save( ofsm, true );
	ofsm.close( );

	CMeta::Shutdown( );
	return 0; }

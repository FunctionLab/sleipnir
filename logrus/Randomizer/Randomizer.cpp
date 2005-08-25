#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	CDat				Dat;
	size_t				i, j;
	ofstream			ofsm;
	gengetopt_args_info	sArgs;
	float				dFrac;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }

	CMeta::Startup( sArgs.verbosity_arg, sArgs.random_arg );
	if( !Dat.Open( sArgs.input_arg ) ) {
		cerr << "Could not open input: " << sArgs.input_arg << endl;
		return 1; }

	dFrac = 0;
	for( i = 0; i < Dat.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
			if( !CMeta::IsNaN( Dat.Get( i, j ) ) ) {
				dFrac++;
				Dat.Set( i, j, CMeta::GetNaN( ) ); }
	if( sArgs.count_arg )
		dFrac = (float)sArgs.count_arg;
	dFrac /= Dat.GetGenes( ) * ( Dat.GetGenes( ) - 1 ) / 2;
	for( i = 0; i < Dat.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
			if( ( (float)rand( ) / RAND_MAX ) < dFrac )
				Dat.Set( i, j, 1 );

	ofsm.open( sArgs.input_arg, ios_base::binary );
	Dat.Save( ofsm, true );
	ofsm.close( );

	CMeta::Shutdown( );
	return 0; }

#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	CDat				Dat;
	CPCL				PCL;
	size_t				i, j;
	ofstream			ofsm;
	gengetopt_args_info	sArgs;
	float				d;
	istream*			pistm;
	ifstream			ifsm;
	ostream*			postm;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( !strcmp( "pcl", sArgs.type_arg ) ) {
		if( sArgs.input_arg ) {
			ifsm.open( sArgs.input_arg );
			pistm = &ifsm; }
		else
			pistm = &cin;
		if( !PCL.Open( *pistm, sArgs.skip_arg ) ) {
			cerr << "Could not open input: " << ( sArgs.input_arg ? sArgs.input_arg : "standard input" ) << endl;
			return 1; }
		if( sArgs.input_arg )
			ifsm.close( );

		PCL.Normalize( );
		if( sArgs.output_arg ) {
			ofsm.open( sArgs.output_arg );
			postm = &ofsm; }
		else
			postm = &cout;
		PCL.Save( *postm );
		if( sArgs.output_arg )
			ofsm.close( ); }
	else if( !sArgs.input_arg ) {
		cmdline_parser_print_help( );
		return 1; }
	else {
		if( !Dat.Open( sArgs.input_arg ) ) {
			cerr << "Could not open input: " << sArgs.input_arg << endl;
			return 1; }

		Dat.Normalize( !sArgs.zscore_flag );
		if( sArgs.flip_flag )
			for( i = 0; i < Dat.GetGenes( ); ++i )
				for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
					if( !CMeta::IsNaN( d = Dat.Get( i, j ) ) )
						Dat.Set( i, j, 1 - d );

		Dat.Save( sArgs.output_arg ? sArgs.output_arg : sArgs.input_arg ); }

	return 0; }
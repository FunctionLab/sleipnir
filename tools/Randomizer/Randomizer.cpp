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
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, and Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	CDat				Dat;
	size_t				i, j;
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

	Dat.Save( sArgs.input_arg );

	return 0; }

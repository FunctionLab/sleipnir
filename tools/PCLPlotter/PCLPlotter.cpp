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
	gengetopt_args_info	sArgs;
	CPCL				PCL;
	vector<float>		vecdSumsIn, vecdSumSqsIn, vecdSumsOut, vecdSumSqsOut;
	vector<size_t>		veciCountsIn, veciCountsOut;
	size_t				i, j;
	float				d;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg );

	if( !PCL.Open( sArgs.input_arg, sArgs.skip_arg ) ) {
		cerr << "Could not open: " << ( sArgs.input_arg ? sArgs.input_arg : "stdin" ) << endl;
		return 1; }

	vecdSumsIn.resize( PCL.GetExperiments( ) );
	vecdSumsOut.resize( PCL.GetExperiments( ) );
	vecdSumSqsIn.resize( PCL.GetExperiments( ) );
	vecdSumSqsOut.resize( PCL.GetExperiments( ) );
	veciCountsIn.resize( PCL.GetExperiments( ) );
	veciCountsOut.resize( PCL.GetExperiments( ) );
	for( i = 0; i < PCL.GetGenes( ); ++i ) {
		bool			fIn			= PCL.GetFeatures( ) && ( PCL.GetFeature( i, 1 )[ 0 ] == '*' );
		vector<float>&	vecdSums	= fIn ? vecdSumsIn : vecdSumsOut;
		vector<float>&	vecdSumSqs	= fIn ? vecdSumSqsIn : vecdSumSqsOut;
		vector<size_t>&	veciCounts	= fIn ? veciCountsIn : veciCountsOut;

		for( j = 0; j < PCL.GetExperiments( ); ++j )
			if( !CMeta::IsNaN( d = PCL.Get( i, j ) ) ) {
				veciCounts[ j ]++;
				vecdSums[ j ] += d;
				vecdSumSqs[ j ] += d * d; } }
	for( i = 0; i < PCL.GetExperiments( ); ++i ) {
		if( j = veciCountsIn[ i ] ) {
			d = ( vecdSumsIn[ i ] /= j );
			vecdSumSqsIn[ i ] = sqrt( ( vecdSumSqsIn[ i ] / j ) - ( d * d ) ); }
		if( j = veciCountsOut[ i ] ) {
			d = ( vecdSumsOut[ i ] /= j );
			vecdSumSqsOut[ i ] = sqrt( ( vecdSumSqsOut[ i ] / j ) - ( d * d ) ); } }

	for( i = 0; i < PCL.GetExperiments( ); ++i )
		cout << ( i ? "\t" : "" ) << PCL.GetExperiment( i );
	cout << endl;
	for( i = 0; i < PCL.GetExperiments( ); ++i )
		cout << ( i ? "\t" : "" ) << vecdSumsIn[ i ];
	cout << endl;
	for( i = 0; i < PCL.GetExperiments( ); ++i )
		cout << ( i ? "\t" : "" ) << vecdSumSqsIn[ i ];
	cout << endl;
	for( i = 0; i < PCL.GetExperiments( ); ++i )
		cout << ( i ? "\t" : "" ) << vecdSumsOut[ i ];
	cout << endl;
	for( i = 0; i < PCL.GetExperiments( ); ++i )
		cout << ( i ? "\t" : "" ) << vecdSumSqsOut[ i ];
	cout << endl;
	

	return 0; }

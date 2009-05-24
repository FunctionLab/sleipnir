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
	CPCL				PCLIn, PCLOut;
	CDataMatrix			MatU, MatV;
	vector<float>		vecdS;
	size_t				i, iNonzero;
	float				dSum, dCur;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg );

	if( !PCLIn.Open( sArgs.input_arg, sArgs.skip_arg ) ) {
		cerr << "Could not open: " << ( sArgs.input_arg ? sArgs.input_arg : "stdin" ) << endl;
		return 1; }
	PCLIn.Get( ).SVD( MatU, MatV, vecdS );

	PCLOut.Open( PCLIn );
	if( sArgs.umatrix_arg ) {
		for( i = 0; i < PCLOut.GetGenes( ); ++i )
			PCLOut.Set( i, MatU.Get( i ) );
		PCLOut.Save( sArgs.umatrix_arg ); }

	for( dSum = 0,i = 0; i < vecdS.size( ); ++i )
		dSum += vecdS[ i ];
	for( dCur = 0,i = 0; i < vecdS.size( ); ++i )
		if( ( dCur += ( vecdS[ i ] / dSum ) ) > sArgs.reprojection_arg )
			break;
	i++;
	iNonzero = min( vecdS.size( ), i );
	for( ; i < vecdS.size( ); ++i ) {
		dSum -= vecdS[ i ];
		vecdS[ i ] = 0; }
	if( sArgs.signal_balance_flag ) {
		dSum /= iNonzero;
		fill( vecdS.begin( ), vecdS.end( ), dSum ); }
	else {
		for( dCur = 0,i = 0; i < iNonzero; ++i )
			dCur += ( vecdS[ i ] = pow( vecdS[ i ] / dSum, (float)sArgs.transform_arg ) );
		for( i = 0; i < iNonzero; ++i )
			vecdS[ i ] *= dSum / dCur; }

	MatU.Multiply( vecdS );
	MatU.Multiply( MatV, true );
	for( i = 0; i < PCLOut.GetGenes( ); ++i )
		PCLOut.Set( i, MatU.Get( i ) );
	PCLOut.Save( sArgs.output_arg );

	return 0; }

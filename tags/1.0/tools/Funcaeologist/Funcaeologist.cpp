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

static const char	c_szDab[]	= ".dab";

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CGenome				Genome;
	CGenes				GenesQuery( Genome );
	ifstream			ifsm;
	size_t				iFunction, i, j;
	vector<size_t>		veciQuery;
	float				d;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( sArgs.genes_arg )
		ifsm.open( sArgs.genes_arg );
	if( !GenesQuery.Open( sArgs.genes_arg ? ifsm : cin ) ) {
		cerr << "Could not open: " << ( sArgs.genes_arg ? sArgs.genes_arg : "genes" ) << endl;
		return 1; }
	if( sArgs.genes_arg )
		ifsm.close( );

	veciQuery.resize( GenesQuery.GetGenes( ) );
	cout << "Function	Score	QueryIn	QueryOut	FuncIn	FuncOut" << endl;
	for( iFunction = 0; iFunction < sArgs.inputs_num; ++iFunction ) {
		CGenes			GenesFunction( Genome );
		CDat			Dat;
		string			strDab;
		float			dFuncIn, dFuncOut, dQueryIn, dQueryOut;
		size_t			iFuncIn, iFuncOut, iQueryIn, iQueryOut, iOne, iTwo;
		vector<size_t>	veciFunction;

		if( !GenesFunction.Open( sArgs.inputs[ iFunction ] ) ) {
			cerr << "Could not open: " << sArgs.inputs[ iFunction ] << endl;
			return 1; }
		strDab = (string)sArgs.directory_arg + '/' + CMeta::Deextension( CMeta::Basename(
			sArgs.inputs[ iFunction ] ) ) + c_szDab;
		if( !Dat.Open( strDab.c_str( ), sArgs.memmap_flag && !sArgs.normalize_flag ) ) {
			cerr << "Could not open: " << strDab << endl;
			return 1; }
		if( sArgs.normalize_flag )
			Dat.Normalize( );

		for( i = 0; i < veciQuery.size( ); ++i )
			veciQuery[ i ] = Dat.GetGene( GenesQuery.GetGene( i ).GetName( ) );
		veciFunction.resize( GenesFunction.GetGenes( ) );
		for( i = 0; i < veciFunction.size( ); ++i )
			veciFunction[ i ] = Dat.GetGene( GenesFunction.GetGene( i ).GetName( ) );

		dFuncIn = dFuncOut = dQueryIn = dQueryOut = 0;
		for( iQueryIn = iFuncIn = iFuncOut = i = 0; i < veciFunction.size( ); ++i ) {
			if( ( iOne = veciFunction[ i ] ) == -1 )
				continue;
			for( j = ( i + 1 ); j < veciFunction.size( ); ++j )
				if( ( ( iTwo = veciFunction[ j ] ) != -1 ) &&
					!CMeta::IsNaN( d = Dat.Get( iOne, iTwo ) ) ) {
					iFuncIn++;
					dFuncIn += d; }
			for( j = 0; j < Dat.GetGenes( ); ++j )
				if( !CMeta::IsNaN( d = Dat.Get( iOne, j ) ) ) {
					iFuncOut++;
					dFuncOut += d; }
			for( j = 0; j < veciQuery.size( ); ++j )
				if( ( ( iTwo = veciQuery[ j ] ) != -1 ) &&
					!CMeta::IsNaN( d = Dat.Get( iOne, iTwo ) ) ) {
					iQueryIn++;
					dQueryIn += d; } }
		for( iQueryOut = i = 0; i < veciQuery.size( ); ++i ) {
			if( ( iOne = veciQuery[ i ] ) == -1 )
				continue;
			for( j = 0; j < Dat.GetGenes( ); ++j )
				if( !CMeta::IsNaN( d = Dat.Get( iOne, j ) ) ) {
					iQueryOut++;
					dQueryOut += d; } }
		dFuncIn /= iFuncIn;
		dFuncOut /= iFuncOut;
		dQueryIn /= iQueryIn;
		dQueryOut /= iQueryOut;

		cout << CMeta::Basename( sArgs.inputs[ iFunction ] ) << '\t' <<
			( dQueryIn / dQueryOut / dFuncIn * dFuncOut ) << '\t' << dQueryIn << '\t' << dQueryOut << '\t' <<
			dFuncIn << '\t' << dFuncOut << endl; }

	return 0; }

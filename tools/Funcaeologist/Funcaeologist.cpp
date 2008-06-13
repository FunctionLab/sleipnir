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

static int MainSet( const gengetopt_args_info& );
static int MainBackground( const gengetopt_args_info& );
static float In( const vector<size_t>&, const vector<size_t>&, const CDat& );

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta = CMeta( sArgs.verbosity_arg, sArgs.random_arg );

	return ( sArgs.sizes_arg ? MainBackground : MainSet )( sArgs ); }

int MainSet( const gengetopt_args_info& sArgs ) {
	CGenome			Genome;
	CGenes			GenesQuery( Genome );
	ifstream		ifsm;
	size_t			iFunction, i, j;
	vector<size_t>	veciQuery;
	float			d;

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

int MainBackground( const gengetopt_args_info& sArgs ) {
	CDat			Dat;
	CDataMatrix		MatAves, MatStds;
	vector<size_t>	veciSizes;
	size_t			i, j, iIndexOne, iIndexTwo, iCountOne, iCountTwo;
	float			d, dAve, dStd, dOutOne, dOutTwo;
	vector<float>	vecdOut;

	if( !Dat.Open( sArgs.input_arg, sArgs.memmap_flag && !sArgs.normalize_flag ) ) {
		cerr << "Could not open: " << ( sArgs.input_arg ? sArgs.input_arg : "stdin" ) << endl;
		return 1; }
	if( sArgs.normalize_flag )
		Dat.Normalize( );
	{
		CPCL	PCLSizes( false );

		if( !PCLSizes.Open( sArgs.sizes_arg, 0 ) ) {
			cerr << "Could not open: " << sArgs.sizes_arg << endl;
			return 1; }
		veciSizes.resize( PCLSizes.GetGenes( ) );
		for( i = 0; i < veciSizes.size( ); ++i )
			veciSizes[ i ] = atoi( PCLSizes.GetGene( i ).c_str( ) );
	}

	vecdOut.resize( Dat.GetGenes( ) );
	fill( vecdOut.begin( ), vecdOut.end( ), 0.0f );
	for( i = 0; i < Dat.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = Dat.Get( i, j ) ) ) {
				vecdOut[ i ] += d;
				vecdOut[ j ] += d; }
	for( i = 0; i < vecdOut.size( ); ++i )
		vecdOut[ i ] /= Dat.GetGenes( ) - 1;

	MatAves.Initialize( veciSizes.size( ), veciSizes.size( ) );
	MatAves.Clear( );
	MatStds.Initialize( veciSizes.size( ), veciSizes.size( ) );
	MatStds.Clear( );
	for( iIndexOne = 0; iIndexOne < veciSizes.size( ); ++iIndexOne ) {
		vector<size_t>	veciOne;

		veciOne.resize( veciSizes[ iIndexOne ] );
		for( iCountOne = 0; iCountOne < (size_t)sArgs.count_arg; ++iCountOne ) {
			for( dOutOne = 0,i = 0; i < veciOne.size( ); ++i ) {
				veciOne[ i ] = rand( ) % Dat.GetGenes( );
				dOutOne += vecdOut[ veciOne[ i ] ]; }
			dOutOne /= veciOne.size( );
			for( iIndexTwo = iIndexOne; iIndexTwo < veciSizes.size( ); ++iIndexTwo ) {
				vector<size_t>	veciTwo;

				if( !( iCountOne % 10 ) )
					cerr << veciSizes[ iIndexOne ] << ':' << veciSizes[ iIndexTwo ] << '\t' << iCountOne << '/'
						<< sArgs.count_arg << endl;
				veciTwo.resize( veciSizes[ iIndexTwo ] );
				for( iCountTwo = 0; iCountTwo < (size_t)sArgs.count_arg; ++iCountTwo ) {
					for( dOutTwo = 0,i = 0; i < veciTwo.size( ); ++i ) {
						veciTwo[ i ] = rand( ) % Dat.GetGenes( );
						dOutTwo += vecdOut[ veciTwo[ i ] ]; }
					dOutTwo /= veciTwo.size( );

					d = ( ( veciOne.size( ) * dOutOne ) + ( veciTwo.size( ) * dOutTwo ) ) /
						( veciOne.size( ) + veciTwo.size( ) );
					d = In( veciOne, veciTwo, Dat ) / d;
					MatAves.Get( iIndexOne, iIndexTwo ) += d;
					MatStds.Get( iIndexOne, iIndexTwo ) += d * d; } } } }

	for( i = 0; i < veciSizes.size( ); ++i )
		cout << '\t' << veciSizes[ i ];
	cout << endl;
	for( i = 0; i < MatAves.GetRows( ); ++i ) {
		cout << veciSizes[ i ];
		for( j = 0; j < i; ++j )
			cout << '\t';
		for( ; j < MatAves.GetColumns( ); ++j ) {
			iCountOne = sArgs.count_arg * sArgs.count_arg;
			dAve = ( MatAves.Get( i, j ) /= iCountOne );
			MatStds.Set( i, j, dStd = sqrt( ( MatStds.Get( i, j ) / iCountOne ) - ( dAve * dAve ) ) );
			cout << '\t' << dStd; }
		cout << endl; }

	return 0; }

float In( const vector<size_t>& veciOne, const vector<size_t>& veciTwo, const CDat& Dat ) {
	size_t	i, j, iIn;
	float	d, dIn;

	for( dIn = 0,iIn = i = 0; i < veciOne.size( ); ++i )
		for( j = 0; j < veciTwo.size( ); ++j )
			if( !CMeta::IsNaN( d = Dat.Get( veciOne[ i ], veciTwo[ j ] ) ) ) {
				iIn++;
				dIn += d; }

	return ( iIn ? ( dIn / iIn ) : 0 ); }

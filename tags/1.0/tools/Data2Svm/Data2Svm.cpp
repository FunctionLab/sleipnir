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

static const char	c_szRBF[]			= "rbf";
static const char	c_szPolynomial[]	= "poly";

int main( int iArgs, char** aszArgs ) {
	CPCL				Data;
	CSVM				SVM;
	ifstream			ifsm;
	ofstream			ofsm;
	CGenome				Genome;
	CGenes				Genes( Genome ), GenesEx( Genome );
	gengetopt_args_info	sArgs;
	vector<float>		vecdResults;
	size_t				i, j;
	float				dAve, dStd;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg, sArgs.random_arg );

	ifsm.open( sArgs.input_arg );
	if( !Data.Open( ifsm, sArgs.skip_arg ) ) {
		cerr << "Could not open: " << sArgs.input_arg << endl;
		return 1; }
	ifsm.close( );
	if( sArgs.normalize_flag )
		Data.Normalize( );
	if( sArgs.random_features_flag )
		Data.Randomize( );

	if( sArgs.genes_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.genes_arg ); }
	if( !Genes.Open( sArgs.genes_arg ? ifsm : cin ) ) {
		cerr << "Could not open: " << ( sArgs.genes_arg ? sArgs.genes_arg : "gene input" ) << endl;
		return 1; }
	if( sArgs.genes_arg )
		ifsm.close( );

	if( sArgs.genex_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.genex_arg );
		if( !GenesEx.Open( ifsm ) ) {
			cerr << "Could not open: " << sArgs.genex_arg << endl;
			return 1; }
		ifsm.close( );

		for( i = 0; i < GenesEx.GetGenes( ); ++i )
			if( ( j = Data.GetGene( GenesEx.GetGene( i ).GetName( ) ) ) != -1 )
				Data.MaskGene( j ); }

	if( sArgs.alphas_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.alphas_arg );
		if( !SVM.OpenAlphas( ifsm ) ) {
			cerr << "Could not open: " << sArgs.alphas_arg << endl;
			return 1; }
		ifsm.close( ); }

	if( !strcmp( sArgs.kernel_arg, c_szRBF ) )
		SVM.SetKernel( CSVM::EKernelRBF );
	else if( !strcmp( sArgs.kernel_arg, c_szPolynomial ) )
		SVM.SetKernel( CSVM::EKernelPolynomial );
	else
		SVM.SetKernel( CSVM::EKernelLinear );

	SVM.SetCache( sArgs.cache_arg );
	SVM.SetIterations( sArgs.iterations_arg );
	SVM.SetGamma( sArgs.gamma_arg );
	SVM.SetDegree( sArgs.degree_arg );
	if( sArgs.tradeoff_given )
		SVM.SetTradeoff( sArgs.tradeoff_arg );
	SVM.SetVerbosity( 0 );

	SVM.Learn( Data, Genes );
	if( sArgs.model_arg ) {
		ofsm.open( sArgs.model_arg );
		SVM.Save( sArgs.model_arg ? (ostream&)ofsm : cout );
		ofsm.close( ); }

	if( sArgs.heldout_flag )
		for( i = 0; i < Data.GetGenes( ); ++i )
			Data.MaskGene( i, !Data.IsMasked( i ) );
	SVM.Evaluate( Data, vecdResults );
	if( sArgs.random_output_flag )
		random_shuffle( vecdResults.begin( ), vecdResults.end( ) );

	dAve = (float)CStatistics::Average( vecdResults );
	dStd = (float)sqrt( CStatistics::Variance( vecdResults, dAve ) );
	for( i = 0; i < vecdResults.size( ); ++i )
		vecdResults[ i ] = ( vecdResults[ i ] - dAve ) / dStd;
	for( i = j = 0; i < Data.GetGenes( ); ++i )
		if( !Data.IsMasked( i ) )
			cout << Data.GetGene( i ) << '\t' << vecdResults[ j++ ] << endl;

	return 0; }

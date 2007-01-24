#include "stdafx.h"
#include "cmdline.h"

static const char	c_szRBF[]			= "rbf";
static const char	c_szPolynomial[]	= "poly";

int main( int iArgs, char** aszArgs ) {
	CPCLSet				PCLs;
	CDataset			Data;
	CSVM				SVM;
	CDataPair			Answers;
	CDat				Dat;
	vector<string>		vecstrInputs;
	size_t				i;
	ifstream			ifsm;
	ofstream			ofsm;
	CGenome				Genome;
	CGenes				GenesIn( Genome ), GenesEx( Genome );
	gengetopt_args_info	sArgs;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg, sArgs.random_arg );

	vecstrInputs.resize( sArgs.inputs_num );
	for( i = 0; i < vecstrInputs.size( ); ++i )
		vecstrInputs[ i ] = sArgs.inputs[ i ];
	if( sArgs.pcl_flag ) {
		if( !PCLs.Open( vecstrInputs, sArgs.skip_arg ) ) {
			cerr << "Could not open PCLs" << endl;
			return 1; } }
	else {
		if( !Data.Open( vecstrInputs ) ) {
			cerr << "Could not open DATs" << endl;
			return 1; } }

	if( sArgs.alphas_arg ) {
		ifsm.open( sArgs.alphas_arg );
		if( !SVM.OpenAlphas( ifsm ) ) {
			cerr << "Could not open: " << sArgs.alphas_arg << endl;
			return 1; }
		ifsm.close( ); }

	if( sArgs.genes_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.genes_arg );
		if( !GenesIn.Open( ifsm ) ) {
			cerr << "Could not open: " << sArgs.genes_arg << endl;
			return 1; }
		ifsm.close( ); }
	if( sArgs.genex_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.genex_arg );
		if( !GenesEx.Open( ifsm ) ) {
			cerr << "Could not open: " << sArgs.genex_arg << endl;
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
	if( sArgs.binary_arg && !sArgs.output_arg ) {
		SVM.Learn( sArgs.binary_arg );
		if( sArgs.model_arg )
			ofsm.open( sArgs.model_arg );
		SVM.Save( sArgs.model_arg ? (ostream&)ofsm : cout );
		if( sArgs.model_arg )
			ofsm.close( );
		else
			cout.flush( ); }
	else if( sArgs.input_arg ) {
		if( !Answers.Open( sArgs.input_arg, false ) ) {
			cerr << "Could not open: " << sArgs.input_arg << endl;
			return 1; }
		if( GenesIn.GetGenes( ) )
			Answers.FilterGenes( GenesIn, CDat::EFilterInclude );
		if( GenesEx.GetGenes( ) )
			Answers.FilterGenes( GenesEx, CDat::EFilterExclude );
		if( sArgs.pcl_flag )
			SVM.Learn( PCLs, Answers );
		else
			SVM.Learn( &Data, Answers );
		if( sArgs.model_arg )
			ofsm.open( sArgs.model_arg );
		SVM.Save( sArgs.model_arg ? (ostream&)ofsm : cout );
		if( sArgs.model_arg )
			ofsm.close( );
		else
			cout.flush( ); }
	else if( sArgs.model_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.model_arg );
		if( !SVM.Open( ifsm ) ) {
			cerr << "Could not open: " << sArgs.model_arg << endl;
			return 1; }

		if( sArgs.binary_arg )
			SVM.Evaluate( sArgs.binary_arg, Dat );
		else {
			const vector<string>&	vecstrGenes	= sArgs.pcl_flag ? PCLs.GetGeneNames( ) :
				Data.GetGeneNames( );

			Dat.Open( vecstrGenes );
			if( GenesIn.GetGenes( ) ) {
				if( sArgs.pcl_flag )
					SVM.Evaluate( PCLs, GenesIn, Dat );
				else
					SVM.Evaluate( &Data, GenesIn, Dat ); }
			else {
				if( sArgs.pcl_flag )
					SVM.Evaluate( PCLs, Dat );
				else
					SVM.Evaluate( &Data, Dat ); } }
		Dat.Normalize( );
		if( sArgs.output_arg )
			Dat.Save( sArgs.output_arg ); }
	else
		cmdline_parser_print_help( );

	CMeta::Shutdown( );
	return 0; }

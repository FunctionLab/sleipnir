#include "stdafx.h"
#include "cmdline.h"

static const char	c_szRBF[]			= "rbf";
static const char	c_szPolynomial[]	= "poly";

int main( int iArgs, char** aszArgs ) {
	CPCLSet				PCLs;
	CSVM				SVM;
	CDataPair			Answers;
	CDat				Dat;
	vector<string>		vecstrPCLs;
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

	vecstrPCLs.resize( sArgs.inputs_num );
	for( i = 0; i < vecstrPCLs.size( ); ++i )
		vecstrPCLs[ i ] = sArgs.inputs[ i ];
	if( !PCLs.Open( vecstrPCLs, sArgs.skip_arg ) ) {
		cerr << "Could not open PCLs" << endl;
		return 1; }

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
	if( sArgs.input_arg ) {
		if( !Answers.Open( sArgs.input_arg, false ) ) {
			cerr << "Could not open: " << sArgs.input_arg << endl;
			return 1; }
		if( GenesIn.GetGenes( ) )
			Answers.FilterGenes( GenesIn, CDat::EFilterInclude );
		if( GenesEx.GetGenes( ) )
			Answers.FilterGenes( GenesEx, CDat::EFilterExclude );
		SVM.Learn( PCLs, Answers );
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
		Dat.Open( PCLs.GetGeneNames( ) );
		if( GenesIn.GetGenes( ) )
			SVM.Evaluate( PCLs, GenesIn, Dat );
		else
			SVM.Evaluate( PCLs, Dat );
		Dat.Normalize( );
		if( sArgs.output_arg ) {
			ofsm.open( sArgs.output_arg, ios_base::binary );
			Dat.Save( ofsm, true );
			ofsm.close( ); } }
	else
		cmdline_parser_print_help( );

	CMeta::Shutdown( );
	return 0; }

#include "stdafx.h"
#include "cmdline.h"

static const char	c_szRBF[]			= "rbf";
static const char	c_szPolynomial[]	= "poly";

int init_svm( const gengetopt_args_info&, CSVM& );
int main_one( const gengetopt_args_info&, const CPCLSet&, const CDataset&, const CGenes&,
	const CGenes& );
int main_many( const gengetopt_args_info&, const CPCLSet&, const CGenes&, const CGenes& );

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CPCLSet				PCLs;
	CDataset			Data;
	vector<string>		vecstrInputs;
	ifstream			ifsm;
	CGenome				Genome;
	CGenes				GenesIn( Genome ), GenesEx( Genome );
	int					iRet;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg, sArgs.random_arg );

	vecstrInputs.resize( sArgs.inputs_num );
	copy( sArgs.inputs, sArgs.inputs + sArgs.inputs_num, vecstrInputs.begin( ) );
	if( sArgs.pcl_flag ) {
		if( !PCLs.Open( vecstrInputs, sArgs.skip_arg ) ) {
			cerr << "Could not open PCLs" << endl;
			return 1; } }
	else {
		if( !Data.Open( vecstrInputs ) ) {
			cerr << "Could not open DATs" << endl;
			return 1; } }

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

	iRet = sArgs.genewise_flag ? main_many( sArgs, PCLs, GenesIn, GenesEx ) :
		main_one( sArgs, PCLs, Data, GenesIn, GenesEx );

	CMeta::Shutdown( );
	return iRet; }

int init_svm( const gengetopt_args_info& sArgs, CSVM& SVM ) {
	ifstream	ifsm;

	if( sArgs.alphas_arg ) {
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

	return 0; }

int main_one( const gengetopt_args_info& sArgs, const CPCLSet& PCLs, const CDataset& Data,
	const CGenes& GenesIn, const CGenes& GenesEx ) {
	CSVM				SVM;
	CDataPair			Answers;
	CDat				Dat;
	vector<string>		vecstrInputs;
	ofstream			ofsm;
	ifstream			ifsm;
	int					iRet;

	if( iRet = init_svm( sArgs, SVM ) )
		return iRet;
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
	else {
		cmdline_parser_print_help( );
		return 1; }

	return 0; }

int main_many( const gengetopt_args_info& sArgs, const CPCLSet& PCLs, const CGenes& GenesIn,
	const CGenes& GenesEx ) {
	vector<CSVM>	vecSVMs;
	size_t			i, j, iGene;
	int				iRet;
	CDataPair		Answers;
	ofstream		ofsm;
	CPCL			PCL;

	if( PCLs.GetPCLs( ) == 0 )
		PCL.Open( cin, sArgs.skip_arg );
	else {
		vector<string>	vecstrExperiments, vecstrFeatures;
		size_t			iPCL, iExp;

		for( iPCL = 0; iPCL < PCLs.GetPCLs( ); ++iPCL )
			for( iExp = 0; iExp < PCLs.Get( iPCL ).GetExperiments( ); ++iExp )
				vecstrExperiments.push_back( PCLs.Get( iPCL ).GetExperiment( iExp ) );
		vecstrFeatures.resize( PCLs.Get( 0 ).GetFeatures( ) - 1 );
		for( i = 0; i < vecstrFeatures.size( ); ++i )
			vecstrFeatures[ i ] = PCLs.Get( 0 ).GetFeature( i + 1 );
		PCL.Open( PCLs.GetGeneNames( ), vecstrExperiments, vecstrFeatures );
		for( iGene = 0; iGene < PCL.GetGenes( ); ++iGene )
			for( i = iPCL = 0; iPCL < PCLs.GetPCLs( ); ++iPCL )
				for( iExp = 0; iExp < PCLs.Get( iPCL ).GetExperiments( ); ++iExp )
					PCL.Set( iGene, i++, PCLs.Get( iPCL, iGene, iExp ) ); }

	vecSVMs.resize( PCL.GetGenes( ) );
	for( i = 0; i < vecSVMs.size( ); ++i )
		if( iRet = init_svm( sArgs, vecSVMs[ i ] ) )
			return iRet;

	if( sArgs.input_arg ) {
		vector<size_t>	veciGenes;
		size_t			iTwo;

		if( !Answers.Open( sArgs.input_arg, false ) ) {
			cerr << "Could not open: " << sArgs.input_arg << endl;
			return 1; }
		if( GenesIn.GetGenes( ) )
			Answers.FilterGenes( GenesIn, CDat::EFilterInclude );
		if( GenesEx.GetGenes( ) )
			Answers.FilterGenes( GenesEx, CDat::EFilterExclude );
		veciGenes.resize( PCL.GetGenes( ) );
		for( i = 0; i < veciGenes.size( ); ++i )
			veciGenes[ i ] = Answers.GetGene( PCL.GetGene( i ) );
		for( i = 0; i < vecSVMs.size( ); ++i ) {
			CGenome			Genome;
			CGenes			GenesPos( Genome ), GenesNeg( Genome );
			vector<string>	vecstrPos, vecstrNeg;
			float			d;

			if( ( iGene = veciGenes[ i ] ) == -1 )
				continue;
			for( j = 0; j < PCL.GetGenes( ); ++j )
				if( ( i != j ) && ( ( iTwo = veciGenes[ j ] ) != -1 ) &&
					!CMeta::IsNaN( d = Answers.Get( iGene, iTwo ) ) ) {
					if( d > 0 )
						vecstrPos.push_back( PCL.GetGene( j ) );
					else
						vecstrNeg.push_back( PCL.GetGene( j ) ); }
			GenesPos.Open( vecstrPos );
			GenesNeg.Open( vecstrNeg );

			vecSVMs[ i ].Learn( PCL, GenesPos, GenesNeg );
			if( sArgs.model_arg ) {
				ofsm.clear( );
				ofsm.open( ( (string)sArgs.model_arg + "/" + CMeta::Filename( PCL.GetGene( i ) ) +
					".svm" ).c_str( ) ); }
			vecSVMs[ i ].Save( sArgs.model_arg ? (ostream&)ofsm : cout );
			if( sArgs.model_arg )
				ofsm.close( );
			else
				cout.flush( ); } }
	else {
		cmdline_parser_print_help( );
		return 1; }

	return 0; }

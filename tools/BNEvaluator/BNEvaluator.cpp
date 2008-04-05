#include "stdafx.h"
#include "cmdline.h"

static int Genes( const char*, CGenes& );

int main( int iArgs, char** aszArgs ) {
	IBayesNet*			pBN;
	gengetopt_args_info	sArgs;
	CPCL				PCLOut;
	CPCLPair			PCLIn;
	CGenome				Genome;
	CGenes				GenesIn( Genome ), GenesEx( Genome );
	size_t				i, j;
	int					iRet;
	vector<string>		vecstrGenes, vecstrExperiments, vecstrNodes, vecstrFeatures;
	vector<size_t>		veciGenes;
	ofstream			ofsm;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );
#ifdef SMILEXML_LIB
	EnableXdslFormat( );
#endif

	CBayesNetSmile	BNSmile( !!sArgs.group_flag );
#ifdef PNL_ENABLED
	CBayesNetPNL	BNPNL( !!sArgs.group_flag );
#endif // PNL_ENABLED
	CBayesNetFN		BNFN;

	if( sArgs.function_flag ) {
		if( !BNFN.Open( sArgs.input_arg ) ) {
			cerr << "Couldn't open: " << sArgs.input_arg << endl;
			return 1; }
		pBN = &BNFN; }
	else {
		if( !BNSmile.Open( sArgs.input_arg ) ) {
			cerr << "Couldn't open: " << sArgs.input_arg << endl;
			return 1; }
#ifdef PNL_ENABLED
		if( sArgs.pnl_flag ) {
			BNSmile.Convert( BNPNL );
			pBN = &BNPNL; }
		else
#endif // PNL_ENABLED
			pBN = &BNSmile; }

	if( ( iRet = Genes( sArgs.genes_arg, GenesIn ) ) ||
		( iRet = Genes( sArgs.genex_arg, GenesEx ) ) )
		return iRet;
	if( !PCLIn.Open( sArgs.data_arg, sArgs.skip_arg ) ) {
		cerr << "Couldn't open: " << ( sArgs.data_arg ? sArgs.data_arg : "PCL" ) << endl;
		return 1; }

	for( i = 0; i < PCLIn.GetGenes( ); ++i ) {
		if( GenesIn.GetGenes( ) && !GenesIn.IsGene( PCLIn.GetGene( i ) ) )
			continue;
		if( GenesEx.GetGenes( ) && GenesEx.IsGene( PCLIn.GetGene( i ) ) )
			continue;
		veciGenes.push_back( i );
		vecstrGenes.push_back( PCLIn.GetGene( i ) ); }

	pBN->GetNodes( vecstrNodes );
	for( i = 0; i < vecstrNodes.size( ); ++i ) {
		for( j = 0; j < PCLIn.GetExperiments( ); ++j )
			if( vecstrNodes[ i ] == PCLIn.GetExperiment( j ) )
				break;
		if( j < PCLIn.GetExperiments( ) )
			continue;
		for( j = 0; j < pBN->GetValues( i ); ++j ) {
			char	acTmp[ 1024 ];

#pragma warning( disable : 4996 )
			sprintf( acTmp, "%s:%d", vecstrNodes[ i ].c_str( ), j );
#pragma warning( default : 4996 )
			vecstrExperiments.push_back( acTmp ); } }

	PCLOut.Open( vecstrGenes, vecstrExperiments, vecstrFeatures );
	pBN->Evaluate( PCLIn, PCLOut, !!sArgs.zero_flag, sArgs.algorithm_arg );
	if( sArgs.output_arg )
		ofsm.open( sArgs.output_arg );
	PCLOut.Save( sArgs.output_arg ? (ostream&)ofsm : cout );
	if( sArgs.output_arg )
		ofsm.close( );

	return 0; }

static int Genes( const char* szGenes, CGenes& Genes ) {
	ifstream	ifsm;

	if( !szGenes )
		return 0;

	ifsm.open( szGenes );
	if( !Genes.Open( ifsm ) ) {
		cerr << "Couldn't open: " << szGenes << endl;
		return 1; }
	ifsm.close( );
	return 0; }

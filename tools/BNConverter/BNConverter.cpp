#include "stdafx.h"
#include "cmdline.h"

static void Evaluate( const IDataset*, const IBayesNet*, bool, ostream& );
static void DebugDataset( const IDataset* );

int main( int iArgs, char** aszArgs ) {
	CDatasetCompact		Data;
	CDatasetCompactMap	DataMap;
	CDataset			DataFull;
	CDataPair			Answers;
	CDataMask			Train, Test;
	IDataset*			pData;
	gengetopt_args_info	sArgs;
	ofstream			ofsm;
	IBayesNet*			pNet;
	size_t				i;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg, sArgs.random_arg );
#if SMILEXML_LIB
	EnableXdslFormat( );
#endif

	CBayesNetSmile	BNSmile( !!sArgs.group_flag ), BNDefault( !!sArgs.group_flag );
#ifdef PNL_ENABLED
	CBayesNetPNL	BNPNL( !!sArgs.group_flag );
#endif // PNL_ENABLED
	CBayesNetFN		BNFN;

	if( sArgs.function_flag ) {
		if( !BNFN.Open( sArgs.input_arg ) ) {
			cerr << "Couldn't open: " << sArgs.input_arg << endl;
			return 1; }
		pNet = &BNFN; }
	else {
		if( !BNSmile.Open( sArgs.input_arg ) ) {
			cerr << "Couldn't open: " << sArgs.input_arg << endl;
			return 1; }
		if( sArgs.default_arg ) {
			if( !BNDefault.Open( sArgs.default_arg ) ) {
				cerr << "Couldn't open: " << sArgs.default_arg << endl;
				return 1; }
			BNSmile.SetDefault( BNDefault ); }
#ifdef PNL_ENABLED
		if( sArgs.pnl_flag ) {
			BNSmile.Convert( BNPNL );
			pNet = &BNPNL; }
		else
#endif // PNL_ENABLED
			pNet = &BNSmile; }
	if( sArgs.randomize_flag )
		pNet->Randomize( );

	if( sArgs.dataset_arg ) {
		if( !DataMap.Open( sArgs.dataset_arg ) ) {
			cerr << "Couldn't open: " << sArgs.dataset_arg << endl;
			return 1; }

		if( sArgs.genes_arg && !DataMap.FilterGenes( sArgs.genes_arg, CDat::EFilterInclude ) ) {
			cerr << "Couldn't open: " << sArgs.genes_arg << endl;
			return 1; }
		if( sArgs.genet_arg && !DataMap.FilterGenes( sArgs.genet_arg, CDat::EFilterTerm ) ) {
			cerr << "Couldn't open: " << sArgs.genet_arg << endl;
			return 1; }
		if( sArgs.genex_arg && !DataMap.FilterGenes( sArgs.genex_arg, CDat::EFilterExclude ) ) {
			cerr << "Couldn't open: " << sArgs.genex_arg << endl;
			return 1; }
		DataMap.FilterAnswers( );
		pData = &DataMap; }
	else {
		if( !Answers.Open( sArgs.answers_arg, pNet->IsContinuous( 0 ) ) ) {
			cerr << "Couldn't open: " << sArgs.answers_arg << endl;
			return 1; }
		if( sArgs.genes_arg && !Answers.FilterGenes( sArgs.genes_arg, CDat::EFilterInclude ) ) {
			cerr << "Couldn't open: " << sArgs.genes_arg << endl;
			return 1; }
		if( sArgs.genet_arg && !Answers.FilterGenes( sArgs.genet_arg, CDat::EFilterTerm ) ) {
			cerr << "Couldn't open: " << sArgs.genet_arg << endl;
			return 1; }
		if( sArgs.genex_arg && !Answers.FilterGenes( sArgs.genex_arg, CDat::EFilterExclude ) ) {
			cerr << "Couldn't open: " << sArgs.genex_arg << endl;
			return 1; }

		if( pNet->IsContinuous( ) ) {
			if( !DataFull.Open( Answers, sArgs.datadir_arg, pNet ) ) {
				cerr << "Couldn't open: " << sArgs.datadir_arg << endl;
				return 1; }
			pData = &DataFull; }
		else {
			if( !Data.Open( Answers, sArgs.datadir_arg, pNet ) ) {
				cerr << "Couldn't open: " << sArgs.datadir_arg << endl;
				return 1; }
			Data.FilterAnswers( );
			pData = &Data; } }

	if( sArgs.test_arg ) {
		Train.AttachRandom( pData, (float)( 1 - sArgs.test_arg ) );
		Test.AttachComplement( Train );
		pData = &Train; }
	if( sArgs.test_arg < 1 ) {
		if( sArgs.checkpoint_flag )
			for( i = 0; i < (size_t)sArgs.iterations_arg; ++i ) {
				pNet->Learn( pData, 1, !!sArgs.zero_flag, !!sArgs.elr_flag );
				cerr << "Iteration " << i << '/' << sArgs.iterations_arg << " complete" << endl;
				pNet->Save( sArgs.output_arg ); }
		else
			pNet->Learn( pData, sArgs.iterations_arg, !!sArgs.zero_flag, !!sArgs.elr_flag ); }

	if( sArgs.murder_given )
		pNet->Randomize( sArgs.murder_arg );
	pNet->Save( sArgs.output_arg );

	if( sArgs.eval_train_arg && ( sArgs.test_arg < 1 ) ) {
		ofsm.open( sArgs.eval_train_arg, ios_base::binary );
		Evaluate( pData, pNet, !!sArgs.zero_flag, ofsm );
		ofsm.close( ); }
	if( sArgs.eval_test_arg && sArgs.test_arg ) {
		ofsm.clear( );
		ofsm.open( sArgs.eval_test_arg, ios_base::binary );
		Evaluate( &Test, pNet, !!sArgs.zero_flag, ofsm );
		ofsm.close( ); }

	CMeta::Shutdown( );
	return 0; }

static void Evaluate( const IDataset* pData, const IBayesNet* pNet, bool fZero, ostream& ostm ) {
	size_t	i, j;
	float	d;
	CDat	Dat;

	Dat.Open( pData->GetGeneNames( ) );
	pNet->Evaluate( pData, Dat, fZero );
	if( !pNet->IsContinuous( ) )
		for( i = 0; i < Dat.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
				if( !CMeta::IsNaN( d = Dat.Get( i, j ) ) )
					Dat.Set( i, j, 1 - d );
	Dat.Save( ostm ); }

static void DebugDataset( const IDataset* pData ) {
	size_t	i, j;

	for( i = 0; i < pData->GetGenes( ); ++i ) {
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j )
			cerr << ( pData->IsExample( i, j ) ? 1 : 0 );
		if( ( i + 1 ) < pData->GetGenes( ) )
			cerr << endl; } }

#include "stdafx.h"
#include "cmdline.h"

void evaluate( const IDataset*, const IBayesNet*, bool, ostream& );
void debug_dataset( const IDataset* );

int main( int iArgs, char** aszArgs ) {
	CDatasetCompact			Data;
	CDataPair				Answers;
	CDataMask				Train, Test;
	IDataset*				pData;
	gengetopt_args_info		sArgs;
	ifstream				ifsm;
	ofstream				ofsm;
	IBayesNet*				pNet;
	CGenome					Genome;
	CGenes					GenesIn( Genome ), GenesEx( Genome );
	size_t					i;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg, sArgs.random_arg );

	CBayesNetSmile	BNSmile( !!sArgs.group_flag );
	CBayesNetPNL	BNPNL( !!sArgs.group_flag );

	if( !BNSmile.Open( sArgs.input_arg ) ) {
		cerr << "Couldn't open: " << sArgs.input_arg << endl;
		return 1; }
	if( sArgs.randomize_flag )
		BNSmile.Randomize( );

	if( !Answers.Open( sArgs.answers_arg, BNSmile.IsContinuous( ) ) ) {
		cerr << "Couldn't open: " << sArgs.answers_arg << endl;
		return 1; }
	if( sArgs.genes_arg ) {
		ifsm.open( sArgs.genes_arg );
		if( !GenesIn.Open( ifsm ) ) {
			cerr << "Couldn't open: " << sArgs.genes_arg << endl;
			return 1; }
		Answers.FilterGenes( GenesIn, CDat::EFilterInclude );
		ifsm.close( ); }
	if( sArgs.genex_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.genex_arg );
		if( !GenesEx.Open( ifsm ) ) {
			cerr << "Couldn't open: " << sArgs.genex_arg << endl;
			return 1; }
		Answers.FilterGenes( GenesEx, CDat::EFilterExclude );
		ifsm.close( ); }
	if( !Data.Open( Answers, sArgs.datadir_arg, &BNSmile ) ) {
		cerr << "Couldn't open: " << sArgs.datadir_arg << endl;
		return 1; }

	pData = &Data;
	if( sArgs.test_arg ) {
		Train.AttachRandom( pData, (float)( 1 - sArgs.test_arg ) );
		Test.AttachComplement( Train );
		pData = &Train; }

	if( sArgs.pnl_flag ) {
		BNSmile.Convert( BNPNL );
		pNet = &BNPNL; }
	else
		pNet = &BNSmile;
	if( sArgs.test_arg < 1 ) {
		if( sArgs.checkpoint_flag )
			for( i = 0; i < sArgs.iterations_arg; ++i ) {
				pNet->Learn( pData, 1, !!sArgs.zero_flag );
				cerr << "Iteration " << (unsigned int)i << '/' << sArgs.iterations_arg <<
					" complete" << endl;
				pNet->Save( sArgs.output_arg ); }
		else
			pNet->Learn( pData, sArgs.iterations_arg, !!sArgs.zero_flag ); }

	if( sArgs.murder_given ) {
		pNet->Reverse( sArgs.murder_arg );
		pNet->Save( sArgs.output_arg ); }

	if( sArgs.eval_train_arg && ( sArgs.test_arg < 1 ) ) {
		ofsm.open( sArgs.eval_train_arg, ios_base::binary );
		evaluate( pData, pNet, !!sArgs.zero_flag, ofsm );
		ofsm.close( ); }
	if( sArgs.eval_test_arg && sArgs.test_arg ) {
		ofsm.clear( );
		ofsm.open( sArgs.eval_test_arg, ios_base::binary );
		evaluate( &Test, pNet, !!sArgs.zero_flag, ofsm );
		ofsm.close( ); }
	if( !sArgs.checkpoint_flag )
		pNet->Save( sArgs.output_arg );

	CMeta::Shutdown( );
	return 0; }

void evaluate( const IDataset* pData, const IBayesNet* pNet, bool fZero, ostream& ostm ) {
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
	Dat.Save( ostm, true ); }

void debug_dataset( const IDataset* pData ) {
	size_t	i, j;

	for( i = 0; i < pData->GetGenes( ); ++i ) {
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j )
			cerr << ( pData->IsExample( i, j ) ? 1 : 0 );
		if( ( i + 1 ) < pData->GetGenes( ) )
			cerr << endl; } }

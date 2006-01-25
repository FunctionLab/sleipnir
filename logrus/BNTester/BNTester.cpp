#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	IBayesNet*				pBN;
	CDatasetCompact			Data;
	CDat					Dat;
	CGenome					Genome;
	CGenes					GenesIn( Genome ), GenesEx( Genome );
	ifstream				ifsm;
	ofstream				ofsm;
	vector<vector<float> >	vecvecdResults;
	float					d;
	size_t					i, j, k;
	gengetopt_args_info		sArgs;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );
#if defined(_MSC_VER) && !defined(_DEBUG)
	EnableXdslFormat( );
#endif // defined(_MSC_VER) && !defined(_DEBUG)

	CBayesNetSmile	BNSmile( !!sArgs.group_flag );
	CBayesNetPNL	BNPNL( !!sArgs.group_flag );

	if( !BNSmile.Open( sArgs.input_arg ) ) {
		cerr << "Couldn't open: " << sArgs.input_arg << endl;
		return 1; }
	if( sArgs.genes_arg ) {
		ifsm.open( sArgs.genes_arg );
		if( !GenesIn.Open( ifsm ) ) {
			cerr << "Couldn't open: " << sArgs.genes_arg << endl;
			return 1; }
		ifsm.close( ); }
	if( sArgs.genex_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.genex_arg );
		if( !GenesEx.Open( ifsm ) ) {
			cerr << "Couldn't open: " << sArgs.genex_arg << endl;
			return 1; }
		ifsm.close( ); }
	if( !Data.Open( sArgs.datadir_arg, &BNSmile, GenesIn, GenesEx ) ) {
		cerr << "Couldn't open: " << sArgs.datadir_arg << endl;
		return 1; }

	if( sArgs.output_arg )
		Dat.Open( Data.GetGeneNames( ) );
	vecvecdResults.clear( );
	cerr << "Evaluating..." << endl;
	if( sArgs.pnl_flag ) {
		BNSmile.Convert( BNPNL );
		pBN = &BNPNL; }
	else
		pBN = &BNSmile;
	if( sArgs.output_arg )
		pBN->Evaluate( &Data, Dat, !!sArgs.zero_flag );
	else
		pBN->Evaluate( &Data, vecvecdResults, !!sArgs.zero_flag );

	if( sArgs.output_arg ) {
		cerr << "Saving..." << endl;
		for( i = 0; i < Dat.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
				if( !CMeta::IsNaN( d = Dat.Get( i, j ) ) )
					Dat.Set( i, j, 1 - d );
		ofsm.open( sArgs.output_arg, ios_base::binary );
		Dat.Save( ofsm, true );
		ofsm.close( ); }
	else {
		cerr << "Storing..." << endl;
		for( k = i = 0; i < Data.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < Data.GetGenes( ); ++j ) {
				if( !Data.IsExample( i, j ) )
					continue;
				d = vecvecdResults[ k++ ][ 0 ];
				if( !pBN->IsContinuous( ) )
					d = 1 - d;
				cout << Data.GetGene( i ) << '\t' << Data.GetGene( j ) << '\t' << d <<
					endl; }
		cout.flush( ); }

	CMeta::Shutdown( );
	return 0; }

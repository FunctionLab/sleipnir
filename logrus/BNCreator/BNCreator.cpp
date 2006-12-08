#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info		sArgs;
	CDataPair				Answers;
	size_t					iArg, i;
	vector<CBayesNetSmile*>	vecpBNs;
	CBayesNetSmile			BNIn, BNOut, BNDefault;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );
	EnableXdslFormat( );

	if( sArgs.default_arg && !BNDefault.Open( sArgs.default_arg ) ) {
		cerr << "Couldn't open: " << sArgs.default_arg << endl;
		return 1; }

	if( !Answers.Open( sArgs.answers_arg, false ) ) {
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

	{
		CDatasetCompact	Data;
		vector<string>	vecstrDummy;

		if( !Data.Open( Answers, vecstrDummy ) ) {
			cerr << "Couldn't open answer set" << endl;
			return 1; }
		vecstrDummy.push_back( "FR" );
		if( !BNIn.Open( &Data, vecstrDummy ) ) {
			cerr << "Couldn't create base network" << endl;
			return 1; }
		if( sArgs.default_arg )
			BNIn.SetDefault( BNDefault );
		if( !BNIn.Learn( &Data, 1 ) ) {
			cerr << "Couldn't learn base network" << endl;
			return 1; }
	}

	for( iArg = 0; iArg < sArgs.inputs_num; ++iArg ) {
		CDatasetCompact	Data;
		vector<string>	vecstrNames;
		CBayesNetSmile*	pBN;

		vecstrNames.push_back( sArgs.inputs[ iArg ] );
		if( !Data.Open( Answers, vecstrNames ) ) {
			cerr << "Couldn't open: " << sArgs.inputs[ iArg ] << endl;
			return 1; }
		vecstrNames.insert( vecstrNames.begin( ), sArgs.answers_arg );
		for( i = 0; i < vecstrNames.size( ); ++i )
			vecstrNames[ i ] = CMeta::Filename( CMeta::Deextension( CMeta::Basename( vecstrNames[ i ].c_str( ) ) ) );
		vecpBNs.push_back( pBN = new CBayesNetSmile( ) );
		if( !pBN->Open( &Data, vecstrNames ) ) {
			cerr << "Couldn't create network for: " << sArgs.inputs[ iArg ] << endl;
			return 1; }

		if( sArgs.default_arg )
			pBN->SetDefault( BNDefault );
		if( !pBN->Learn( &Data, 1, !!sArgs.zero_flag ) ) {
			cerr << "Couldn't learn network for: " << sArgs.inputs[ iArg ] << endl;
			return 1; } }

	if( !BNOut.Open( BNIn, vecpBNs ) ) {
		cerr << "Couldn't merge networks" << endl;
		return 1; }
	BNOut.Save( sArgs.output_arg );

	for( i = 0; i < vecpBNs.size( ); ++i )
		delete vecpBNs[ i ];

	CMeta::Shutdown( );
	return 0; }

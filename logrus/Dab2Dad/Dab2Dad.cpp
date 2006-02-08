#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CDatasetCompact		Data;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );
#if !( defined(_MSC_VER) && defined(_DEBUG) )
    EnableXdslFormat( );
#endif // !( defined(_MSC_VER) && defined(_DEBUG) )

	if( sArgs.load_arg ) {
		CDatasetCompactMap	DataMap;

		if( !DataMap.Open( sArgs.load_arg ) ) {
			cerr << "Couldn't open: " << sArgs.load_arg << endl;
			return 1; }
#ifdef _MSC_VER
		Sleep( INFINITE );
#else // _MSC_VER
		sleep( -1 );
#endif // _MSC_VER
		return 0; }

	if( sArgs.input_arg ) {
		ifstream	ifsm;

		ifsm.open( sArgs.input_arg, ios_base::binary );
		if( !Data.Open( ifsm ) ) {
			cerr << "Couldn't open: " << sArgs.input_arg << endl;
			return 1; }
		ifsm.close( ); }
	else if( sArgs.network_arg ) {
		CDataPair		Answers;
		CBayesNetSmile	BNSmile;

		if( !sArgs.answers_arg ) {
			cmdline_parser_print_help( );
			return 1; }
		if( !BNSmile.Open( sArgs.network_arg ) ) {
			cerr << "Couldn't open: " << sArgs.network_arg << endl;
			return 1; }
		if( !Answers.Open( sArgs.answers_arg, false ) ) {
			cerr << "Couldn't open: " << sArgs.answers_arg << endl;
			return 1; }
		if( !Data.Open( Answers, sArgs.directory_arg, &BNSmile ) ) {
			cerr << "Couldn't open: " << sArgs.directory_arg << endl;
			return 1; } }
	else {
		vector<string>	vecstrFiles;

		vecstrFiles.resize( sArgs.inputs_num );
		copy( sArgs.inputs, sArgs.inputs + sArgs.inputs_num, vecstrFiles.begin( ) );
		if( !Data.Open( vecstrFiles ) ) {
			cerr << "Couldn't open inputs" << endl;
			return 1; } }

	if( sArgs.mask_arg ) {
		CDat		Mask;
		vector<int>	veciGenes;
		size_t		i, j;
		int			iOne, iTwo;
		float		d;

		if( !Mask.Open( sArgs.mask_arg ) ) {
			cerr << "Couldn't open: " << sArgs.mask_arg << endl;
			return 1; }
		veciGenes.resize( Data.GetGenes( ) );
		for( i = 0; i < Data.GetGenes( ); ++i )
			veciGenes[ i ] = (int)Mask.GetGene( Data.GetGene( i ) );
		for( i = 0; i < Data.GetGenes( ); ++i ) {
			if( ( iOne = veciGenes[ i ] ) == -1 ) {
				for( j = ( i + 1 ); j < Data.GetGenes( ); ++j )
					Data.Remove( i, j );
				continue; }
			for( j = ( i + 1 ); j < Data.GetGenes( ); ++j )
				if( ( ( iTwo = veciGenes[ j ] ) == -1 ) ||
					CMeta::IsNaN( d = Mask.Get( iOne, iTwo ) ) || !d )
					Data.Remove( i, j ); } }

	if( sArgs.output_arg ) {
		ofstream	ofsm;

		ofsm.open( sArgs.output_arg, ios_base::binary );
		Data.Save( ofsm, true );
		ofsm.close( ); }
	else
		Data.Save( cout, false );

	CMeta::Shutdown( );
	return 0; }

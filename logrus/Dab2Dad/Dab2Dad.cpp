#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CDatasetCompact		Data;
	CGenome				Genome;
	CGenes				GenesIn( Genome ), GenesEx( Genome );
	ifstream			ifsm;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );
    EnableXdslFormat( );

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

	if( sArgs.input_arg ) {
		ifsm.clear( );
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
		if( !Data.Open( Answers, sArgs.directory_arg, &BNSmile, GenesIn, GenesEx, !!sArgs.everything_flag ) ) {
			cerr << "Couldn't open: " << sArgs.directory_arg << endl;
			return 1; } }
	else if( sArgs.lookup1_arg && sArgs.lookup2_arg ) {
		size_t	i;

		for( i = 0; i < sArgs.inputs_num; ++i ) {
			CDataPair	Dat;
			size_t		iOne, iTwo;
			float		d;

			if( !Dat.Open( sArgs.inputs[ i ], false, !!sArgs.memmap_flag ) ) {
				cerr << "Couldn't open: " << sArgs.inputs[ i ] << endl;
				return 1; }
			cout << sArgs.inputs[ i ];
			if( ( ( iOne = Dat.GetGene( sArgs.lookup1_arg ) ) != -1 ) &&
				( ( iTwo = Dat.GetGene( sArgs.lookup2_arg ) ) != -1 ) &&
				!CMeta::IsNaN( d = Dat.Get( iOne, iTwo ) ) )
				cout << '\t' << d;
			cout << endl; }
		return 0; }
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

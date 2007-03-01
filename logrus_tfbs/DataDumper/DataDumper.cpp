#include "stdafx.h"
#include "cmdline.h"

static const char	c_acDab[]	= ".dab";

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	size_t						i, j, k, iPairs, iPair, iArg;
	map<string,size_t>			mapZeros;
	vector<string>				vecstrNames;
	CDataPair					Answers;
	CFullMatrix<unsigned char>	MatData;
	vector<size_t>				veciGenes;
	float						d;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( sArgs.zeros_arg ) {
		ifstream		ifsm;
		vector<string>	vecstrZeros;
		char			acLine[ 1024 ];

		ifsm.open( sArgs.zeros_arg );
		if( !ifsm.is_open( ) ) {
			cerr << "Couldn't open: " << sArgs.zeros_arg << endl;
			return 1; }
		while( !ifsm.eof( ) ) {
			ifsm.getline( acLine, ARRAYSIZE(acLine) - 1 );
			acLine[ ARRAYSIZE(acLine) - 1 ] = 0;
			vecstrZeros.clear( );
			CMeta::Tokenize( acLine, vecstrZeros );
			if( vecstrZeros.empty( ) )
				continue;
			mapZeros[ vecstrZeros[ 0 ] ] = atoi( vecstrZeros[ 1 ].c_str( ) ); } }

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

	for( iPairs = i = 0; i < Answers.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < Answers.GetGenes( ); ++j )
			if( !CMeta::IsNaN( Answers.Get( i, j ) ) )
				iPairs++;
	MatData.Initialize( iPairs, sArgs.inputs_num );
	MatData.Clear( );

	veciGenes.resize( Answers.GetGenes( ) );
	for( iArg = 0; iArg < sArgs.inputs_num; ++iArg ) {
		CDatasetCompact						Data;
		size_t								iOne, iTwo, iZero, iVal;
		map<string,size_t>::const_iterator	iterZero;

		vecstrNames.clear( );
		vecstrNames.push_back( sArgs.inputs[ iArg ] );
		if( !Data.Open( Answers, vecstrNames, true ) ) {
			cerr << "Couldn't open: " << sArgs.inputs[ iArg ] << endl;
			return 1; }
		vecstrNames[ 0 ] = CMeta::Filename( CMeta::Deextension( CMeta::Basename( vecstrNames[ 0 ].c_str( ) ) ) );
		iZero = ( ( iterZero = mapZeros.find( vecstrNames[ 0 ] ) ) == mapZeros.end( ) ) ? -1 :
			iterZero->second;
		for( i = 0; i < veciGenes.size( ); ++i )
			veciGenes[ i ] = Data.GetGene( Answers.GetGene( i ) );
		for( iPair = i = 0; i < veciGenes.size( ); ++i ) {
			iOne = veciGenes[ i ];
			for( j = ( i + 1 ); j < veciGenes.size( ); ++j ) {
				if( CMeta::IsNaN( Answers.Get( i, j ) ) )
					continue;
				if( ( iOne != -1 ) && ( ( iTwo = veciGenes[ j ] ) != -1 ) ) {
					iVal = Data.GetDiscrete( iOne, iTwo, 1 );
					if( ( iVal != -1 ) || ( ( iVal = iZero ) != -1 ) || ( sArgs.zero_flag && !( iVal = 0 ) ) )
						MatData.Set( iPair, iArg, (unsigned char)( iVal + 1 ) ); }
				iPair++; } } }

	cout << "Gene 1	Gene 2	Answer";
	for( i = 0; i < sArgs.inputs_num; ++i )
		cout << '\t' << sArgs.inputs[ i ];
	cout << endl;
	for( iPair = i = 0; i < Answers.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < Answers.GetGenes( ); ++j ) {
			if( CMeta::IsNaN( d = Answers.Get( i, j ) ) )
				continue;
			cout << Answers.GetGene( i ) << '\t' << Answers.GetGene( j ) << '\t' << d;
			for( k = 0; k < MatData.GetColumns( ); ++k )
				cout << '\t' << ( (int)MatData.Get( iPair, k ) - 1 );
			cout << endl;
			iPair++; }

	CMeta::Shutdown( );
	return 0; }

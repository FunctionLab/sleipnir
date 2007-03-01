#include "stdafx.h"
#include "cmdline.h"

static const char	c_acDab[]	= ".dab";

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	size_t				i;
	map<string,size_t>	mapZeros;
	CBayesNetSmile		BNIn;
	vector<string>		vecstrNames;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );
	EnableXdslFormat( );

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

	if( sArgs.input_arg ) {
		vector<string>	vecstrFiles, vecstrGenes;
		CDat			DatYes, DatNo;
		CDatasetCompact	Data;
		vector<size_t>	veciGenes;
		size_t			j, k, iOne, iTwo, iBin, iZero, iIndex;
		float			d;
		CDataMatrix		MatCPT;
		ofstream		ofsm;
		char			acTemp[ L_tmpnam + 1 ];
		const char*		szTemp;
		float*			adYes;
		float*			adNo;
		CGenome			Genome;
		CGenes			GenesIn( Genome ), GenesEx( Genome );

		if( !BNIn.Open( sArgs.input_arg ) ) {
			cerr << "Couldn't open: " << sArgs.input_arg << endl;
			return 1; }

		BNIn.GetNodes( vecstrFiles );
		vecstrFiles.erase( vecstrFiles.begin( ) );
		for( i = 0; i < vecstrFiles.size( ); ++i )
			vecstrFiles[ i ] = (string)sArgs.directory_arg + '/' + vecstrFiles[ i ] + c_acDab;
		if( !Data.OpenGenes( vecstrFiles ) ) {
			cerr << "Couldn't open: " << sArgs.directory_arg << endl;
			return 1; }
		if( sArgs.genes_arg ) {
			ifstream	ifsm;

			ifsm.open( sArgs.genes_arg );
			if( !GenesIn.Open( ifsm ) ) {
				cerr << "Couldn't open: " << sArgs.genes_arg << endl;
				return 1; }
			ifsm.close( );
			for( i = 0; i < GenesIn.GetGenes( ); ++i )
				vecstrGenes.push_back( GenesIn.GetGene( i ).GetName( ) ); }

#pragma warning( disable : 4996 )
		if( !( szTemp = tmpnam( acTemp ) ) ) {
			cerr << "Couldn't create temp file: " << acTemp << endl;
			return 1; }
#pragma warning( default : 4996 )
		DatYes.Open( vecstrGenes.size( ) ? vecstrGenes : Data.GetGeneNames( ), false, sArgs.output_arg );
		DatNo.Open( DatYes.GetGeneNames( ), false ); // , szTemp );
		adYes = new float[ DatYes.GetGenes( ) ];
		adNo = new float[ DatNo.GetGenes( ) ];
		BNIn.GetCPT( 0, MatCPT );
		d = log( MatCPT.Get( 0, 0 ) );
		for( i = 0; i < DatNo.GetGenes( ); ++i )
			adNo[ i ] = d;
		for( i = 0; i < DatNo.GetGenes( ); ++i )
			DatNo.Set( i, adNo );
		d = log( MatCPT.Get( 1, 0 ) );
		for( i = 0; i < DatYes.GetGenes( ); ++i )
			adYes[ i ] = d;
		for( i = 0; i < DatYes.GetGenes( ); ++i )
			DatYes.Set( i, adYes );

		BNIn.GetNodes( vecstrNames );
		veciGenes.resize( DatYes.GetGenes( ) );
		for( i = 0; i < vecstrFiles.size( ); ++i ) {
			vector<string>						vecstrDatum;
			map<string,size_t>::const_iterator	iterZero;

			vecstrDatum.push_back( vecstrFiles[ i ] );
			if( !Data.Open( vecstrDatum ) ) {
				cerr << "Couldn't open: " << vecstrFiles[ i ] << endl;
				return 1; }
			iZero = ( ( iterZero = mapZeros.find( vecstrNames[ i + 1 ] ) ) == mapZeros.end( ) ) ? -1 :
					iterZero->second;
			BNIn.GetCPT( i + 1, MatCPT );
			for( j = 0; j < veciGenes.size( ); ++j )
				veciGenes[ j ] = Data.GetGene( DatYes.GetGene( j ) );
			for( j = 0; j < DatYes.GetGenes( ); ++j ) {
				if( ( iOne = veciGenes[ j ] ) == -1 )
					continue;
				memcpy( adYes, DatYes.Get( j ), ( DatYes.GetGenes( ) - j - 1 ) * sizeof(*adYes) );
				memcpy( adNo, DatNo.Get( j ), ( DatNo.GetGenes( ) - j - 1 ) * sizeof(*adNo) );
				for( k = ( j + 1 ); k < DatYes.GetGenes( ); ++k ) {
					if( ( ( iTwo = veciGenes[ k ] ) == -1 ) || !Data.IsExample( iOne, iTwo ) )
						continue;
					if( ( iBin = Data.GetDiscrete( iOne, iTwo, 0 ) ) == -1 )
						iBin = iZero;
					if( iBin == -1 )
						continue;
					adNo[ iIndex = ( k - j - 1 ) ] += log( MatCPT.Get( iBin, 0 ) );
					adYes[ iIndex ] += log( MatCPT.Get( iBin, 1 ) ); }
				DatYes.Set( j, adYes );
				DatNo.Set( j, adNo ); } }
		for( i = 0; i < DatYes.GetGenes( ); ++i ) {
			memcpy( adYes, DatYes.Get( i ), ( DatYes.GetGenes( ) - i - 1 ) * sizeof(*adYes) );
			memcpy( adNo, DatNo.Get( i ), ( DatNo.GetGenes( ) - i - 1 ) * sizeof(*adNo) );
			for( j = 0; j < ( DatYes.GetGenes( ) - i - 1 ); ++j ) {
				d = exp( adYes[ j ] );
				adYes[ j ] = d / ( d + exp( adNo[ j ] ) ); }
			DatYes.Set( i, adYes ); }
		_unlink( szTemp ); }
	else {
		CDataPair				Answers;
		size_t					iArg;
		vector<CBayesNetSmile*>	vecpBNs;
		CBayesNetSmile			BNDefault, BNOut;
		vector<size_t>			veciZeros;

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

			if( !Data.Open( Answers, vecstrDummy, true ) ) {
				cerr << "Couldn't open answer set" << endl;
				return 1; }
			vecstrDummy.push_back( "FR" );
			if( !BNIn.Open( &Data, vecstrDummy, veciZeros ) ) {
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
			CBayesNetSmile*	pBN;

			vecstrNames.clear( );
			vecstrNames.push_back( sArgs.inputs[ iArg ] );
			if( !Data.Open( Answers, vecstrNames, true ) ) {
				cerr << "Couldn't open: " << sArgs.inputs[ iArg ] << endl;
				return 1; }
			vecstrNames.insert( vecstrNames.begin( ), sArgs.answers_arg );
			for( i = 0; i < vecstrNames.size( ); ++i )
				vecstrNames[ i ] = CMeta::Filename( CMeta::Deextension( CMeta::Basename( vecstrNames[ i ].c_str( ) ) ) );
			vecpBNs.push_back( pBN = new CBayesNetSmile( ) );
			veciZeros.resize( vecstrNames.size( ) );
			for( i = 0; i < veciZeros.size( ); ++i ) {
				map<string,size_t>::const_iterator	iterZero;

				veciZeros[ i ] = ( ( iterZero = mapZeros.find( vecstrNames[ i ] ) ) == mapZeros.end( ) ) ? -1 :
					iterZero->second; }
			if( !pBN->Open( &Data, vecstrNames, veciZeros ) ) {
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
			delete vecpBNs[ i ]; }

	CMeta::Shutdown( );
	return 0; }

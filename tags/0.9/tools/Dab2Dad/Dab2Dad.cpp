/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CDatasetCompact		Data;
	CGenome				Genome;
	CGenes				GenesIn( Genome ), GenesEx( Genome );
	ifstream			ifsm;
	size_t				i, j, k, iOne, iTwo;
	float				d;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );
#ifdef SMILEXML_LIB
    EnableXdslFormat( );
#endif

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
	else if( sArgs.lookupp_arg ) {
		CDat			DatLookup;
		vector<size_t>	veciGenes;
		CPCL			PCLLookup;
		vector<string>	vecstrNames, vecstrExperiments, vecstrDummy;

		ifsm.clear( );
		ifsm.open( sArgs.lookupp_arg );
		if( !DatLookup.Open( ifsm, CDat::EFormatText, 1 ) ) {
			cerr << "Couldn't open: " << sArgs.lookupp_arg << endl;
			return 1; }
		ifsm.close( );

		for( i = 0; i < DatLookup.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < DatLookup.GetGenes( ); ++j )
				if( DatLookup.Get( i, j ) == 1 )
					vecstrNames.push_back( DatLookup.GetGene( i ) + " - " + DatLookup.GetGene( j ) );
		vecstrExperiments.resize( sArgs.inputs_num );
		copy( sArgs.inputs, sArgs.inputs + sArgs.inputs_num, vecstrExperiments.begin( ) );
		PCLLookup.Open( vecstrNames, vecstrExperiments, vecstrDummy );
		
		veciGenes.resize( DatLookup.GetGenes( ) );
		for( i = 0; i < sArgs.inputs_num; ++i ) {
			CDataPair	Dat;
			size_t		iGene;

			if( !Dat.Open( sArgs.inputs[ i ], false, !!sArgs.memmap_flag ) ) {
				cerr << "Couldn't open: " << sArgs.inputs[ i ] << endl;
				return 1; }
			for( j = 0; j < veciGenes.size( ); ++j )
				veciGenes[ j ] = Dat.GetGene( DatLookup.GetGene( j ) );
			for( iGene = j = 0; j < DatLookup.GetGenes( ); ++j ) {
				iOne = veciGenes[ j ];
				for( k = ( j + 1 ); k < DatLookup.GetGenes( ); ++k ) {
					if( DatLookup.Get( j, k ) != 1 )
						continue;
					iTwo = veciGenes[ k ];
					if( ( iOne != -1 ) && ( iTwo != -1 ) )
						PCLLookup.Set( iGene, i, Dat.Get( iOne, iTwo ) );
					iGene++; } } }
		PCLLookup.Save( cout );
		return 0; }
	else if( sArgs.lookups_arg ) {
		CGenes			GenesLk( Genome );
		vector<size_t>	veciGenes;

		ifsm.clear( );
		ifsm.open( sArgs.lookups_arg );
		if( !GenesLk.Open( ifsm ) ) {
			cerr << "Couldn't open: " << sArgs.lookups_arg << endl;
			return 1; }
		ifsm.close( );
		veciGenes.resize( GenesLk.GetGenes( ) );
		for( i = 0; i < sArgs.inputs_num; ++i ) {
			CDataPair	Dat;

			if( !Dat.Open( sArgs.inputs[ i ], false, !!sArgs.memmap_flag ) ) {
				cerr << "Couldn't open: " << sArgs.inputs[ i ] << endl;
				return 1; }
			for( j = 0; j < GenesLk.GetGenes( ); ++j )
				veciGenes[ j ] = Dat.GetGene( GenesLk.GetGene( j ).GetName( ) );
			cout << sArgs.inputs[ i ];
			if( sArgs.lookup1_arg ) {
				if( ( iOne = Dat.GetGene( sArgs.lookup1_arg ) ) != -1 )
					for( j = 0; j < veciGenes.size( ); ++j )
						if( ( ( iTwo = veciGenes[ j ] ) != -1 ) &&
							!CMeta::IsNaN( d = Dat.Get( iOne, iTwo ) ) )
							cout << '\t' << d; }
			else
				for( j = 0; j < veciGenes.size( ); ++j ) {
					if( ( iOne = veciGenes[ j ] ) == -1 )
						continue;
					for( k = ( j + 1 ); k < veciGenes.size( ); ++k )
						if( ( ( iTwo = veciGenes[ k ] ) != -1 ) &&
							!CMeta::IsNaN( d = Dat.Get( iOne, iTwo ) ) )
							cout << '\t' << Dat.Get( iOne, iTwo ); }
			cout << endl; }
		return 0; }
	else if( sArgs.lookup1_arg && sArgs.lookup2_arg ) {
		for( i = 0; i < sArgs.inputs_num; ++i ) {
			CDataPair	Dat;

			if( !Dat.Open( sArgs.inputs[ i ], false, !!sArgs.memmap_flag ) ) {
				cerr << "Couldn't open: " << sArgs.inputs[ i ] << endl;
				return 1; }
			cout << sArgs.inputs[ i ];
			if( ( ( iOne = Dat.GetGene( sArgs.lookup1_arg ) ) != -1 ) &&
				( ( iTwo = Dat.GetGene( sArgs.lookup2_arg ) ) != -1 ) &&
				!CMeta::IsNaN( d = Dat.Get( iOne, iTwo ) ) )
				cout << '\t' << ( sArgs.quantize_flag ? Dat.Quantize( d ) : d );
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
					CMeta::IsNaN( d = Mask.Get( iOne, iTwo ) ) || ( d <= 0 ) )
					Data.Remove( i, j ); } }

	if( sArgs.output_arg ) {
		ofstream	ofsm;

		ofsm.open( sArgs.output_arg, ios_base::binary );
		Data.Save( ofsm, true );
		ofsm.close( ); }
	else
		Data.Save( cout, false );

	return 0; }

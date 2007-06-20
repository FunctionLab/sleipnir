#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info			sArgs;
	CDataPair					DatOne, DatTwo;
	size_t						i, j, iDatOne, iDatTwo, iGeneOne, iGeneTwo;
	vector<size_t>				veciGenes, veciOne, veciTwo;
	CMeasureMutualInformation	MeasureMI;
	float*						adOne;
	float*						adTwo;
	double						d;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( sArgs.table_flag ) {
		for( i = 0; i < sArgs.inputs_num; ++i )
			cout << '\t' << sArgs.inputs[ i ];
		cout << endl; }
	for( iDatOne = 0; iDatOne < sArgs.inputs_num; ++iDatOne ) {
		if( sArgs.table_flag ) {
			cout << sArgs.inputs[ iDatOne ];
			for( i = 0; i <= iDatOne; ++i )
				cout << '\t'; }
		if( ( iDatOne + 1 ) == sArgs.inputs_num )
			break;
		DatOne.Open( sArgs.inputs[ iDatOne ], false, !!sArgs.memmap_flag );
		veciGenes.resize( DatOne.GetGenes( ) );
		for( iDatTwo = ( iDatOne + 1 ); iDatTwo < sArgs.inputs_num; ++iDatTwo ) {
			DatTwo.Open( sArgs.inputs[ iDatTwo ], false, !!sArgs.memmap_flag );
			for( i = 0; i < veciGenes.size( ); ++i )
				veciGenes[ i ] = DatTwo.GetGene( DatOne.GetGene( i ) );
			veciOne.clear( );
			veciTwo.clear( );
			for( i = 0; i < veciGenes.size( ); ++i ) {
				if( ( iGeneOne = veciGenes[ i ] ) == -1 )
					continue;
				for( j = ( i + 1 ); j < veciGenes.size( ); ++j ) {
					if( ( iGeneTwo = veciGenes[ j ] ) == -1 )
						continue;
					veciOne.push_back( DatOne.Quantize( DatOne.Get( i, j ) ) );
					veciTwo.push_back( DatTwo.Quantize( DatTwo.Get( iGeneOne, iGeneTwo ) ) ); } }
			adOne = new float[ veciOne.size( ) ];
			adTwo = new float[ veciTwo.size( ) ];
			for( i = 0; i < veciOne.size( ); ++i ) {
				adOne[ i ] = (float)veciOne[ i ];
				adTwo[ i ] = (float)veciTwo[ i ]; }
			d = MeasureMI.Measure( adOne, veciOne.size( ), adTwo, veciTwo.size( ), IMeasure::EMapNone,
				NULL, NULL );
			if( sArgs.table_flag )
				cout << '\t' << d;
			else
				cout << sArgs.inputs[ iDatOne ] << '\t' << sArgs.inputs[ iDatTwo ] << '\t' << d <<
					endl; }
		if( sArgs.table_flag )
			cout << endl; }

	CMeta::Shutdown( );
	return 0; }

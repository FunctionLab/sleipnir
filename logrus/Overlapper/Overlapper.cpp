#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CDataPair			DatOne, DatTwo;
	size_t				i, j, iOne, iTwo;
	float				dOne, dTwo;
	vector<size_t>		veciGenes;
	CFullMatrix<size_t>	MatConfusion;
	size_t*				ai;
	const char*			szOne;
	const char*			szTwo;
	bool				fOneSmall;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( !DatOne.Open( sArgs.first_arg, false, !!sArgs.memmap_flag ) ) {
		cerr << "Couldn't open: " << sArgs.first_arg << endl;
		return 1; }
	if( !DatTwo.Open( sArgs.second_arg, false, !!sArgs.memmap_flag ) ) {
		cerr << "Couldn't open: " << sArgs.second_arg << endl;
		return 1; }

	cout << sArgs.first_arg << ": " << DatOne.GetGenes( ) << " genes" << endl;
	cout << sArgs.second_arg << ": " << DatTwo.GetGenes( ) << " genes" << endl << endl;

	fOneSmall = DatOne.GetGenes( ) < DatTwo.GetGenes( );
	szOne = fOneSmall ? sArgs.first_arg : sArgs.second_arg;
	szTwo = fOneSmall ? sArgs.second_arg : sArgs.first_arg;
	{
		const CDataPair&	DatSmall	= fOneSmall ? DatOne : DatTwo;
		CDataPair&			DatBig		= fOneSmall ? DatTwo : DatOne;
		vector<bool>		vecfShared;

		MatConfusion.Initialize( DatSmall.GetValues( ) + 1, DatBig.GetValues( ) + 1 );
		MatConfusion.Clear( );
		veciGenes.resize( DatSmall.GetGenes( ) );
		for( i = 0; i < veciGenes.size( ); ++i )
			veciGenes[ i ] = DatBig.GetGene( DatSmall.GetGene( i ) );
		for( i = 0; i < DatSmall.GetGenes( ); ++i ) {
			if( ( iOne = veciGenes[ i ] ) == -1 ) {
				for( j = ( i + 1 ); j < DatSmall.GetGenes( ); ++j )
					if( !CMeta::IsNaN( dOne = DatSmall.Get( i, j ) ) )
						MatConfusion.Get( DatSmall.Quantize( dOne ), DatBig.GetValues( ) )++;
				continue; }
			for( j = ( i + 1 ); j < DatSmall.GetGenes( ); ++j )
				if( !CMeta::IsNaN( dOne = DatSmall.Get( i, j ) ) )
					MatConfusion.Get( DatSmall.Quantize( dOne ), ( ( iTwo = veciGenes[ j ] ) == -1 ) || CMeta::IsNaN( dTwo = DatBig.Get( iOne, iTwo ) ) ?
						DatBig.GetValues( ) : DatBig.Quantize( dTwo ) )++; }
		vecfShared.resize( DatBig.GetGenes( ) );
		for( i = 0; i < vecfShared.size( ); ++i )
			vecfShared[ i ] = ( DatSmall.GetGene( DatBig.GetGene( i ) ) != -1 );
		for( i = 0; i < DatBig.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < DatBig.GetGenes( ); ++j )
				if( !( vecfShared[ i ] && vecfShared[ j ] ) && !CMeta::IsNaN( dTwo = DatBig.Get( i, j ) ) )
					MatConfusion.Get( DatSmall.GetValues( ), DatBig.Quantize( dTwo ) )++;
	}

	cout << '\t' << szTwo << endl << szOne;
	for( i = 0; i < MatConfusion.GetRows( ); ++i ) {
		for( j = 0; j < MatConfusion.GetColumns( ); ++j )
			cout << '\t' << MatConfusion.Get( i, j );
		cout << endl; }

	CMeta::Shutdown( );
	return 0; }

#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CDat				DatOne, DatTwo;
	size_t				i, j, iOne, iTwo;
	float				dOne, dTwo;
	vector<size_t>		veciGenes;
	size_t				aiOne[ 2 ], aiTwo[ 2 ], aiAgree[ 2 ], aiDisagree[ 2 ];
	size_t*				ai;
	const char*			szOne;
	const char*			szTwo;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( !DatOne.Open( sArgs.first_arg ) ) {
		cerr << "Couldn't open: " << sArgs.first_arg << endl;
		return 1; }
	if( !DatTwo.Open( sArgs.second_arg ) ) {
		cerr << "Couldn't open: " << sArgs.second_arg << endl;
		return 1; }

	memset( aiOne, 0, sizeof(aiOne) );
	memset( aiTwo, 0, sizeof(aiTwo) );
	memset( aiAgree, 0, sizeof(aiAgree) );
	memset( aiDisagree, 0, sizeof(aiDisagree) );
	szOne = ( DatOne.GetGenes( ) < DatTwo.GetGenes( ) ) ? sArgs.first_arg : sArgs.second_arg;
	szTwo = ( DatOne.GetGenes( ) < DatTwo.GetGenes( ) ) ? sArgs.second_arg : sArgs.first_arg;
	{
		const CDat&	DatSmall	= ( DatOne.GetGenes( ) < DatTwo.GetGenes( ) ) ? DatOne : DatTwo;
		CDat&	DatBig			= ( DatOne.GetGenes( ) < DatTwo.GetGenes( ) ) ? DatTwo : DatOne;

		veciGenes.resize( DatSmall.GetGenes( ) );
		for( i = 0; i < veciGenes.size( ); ++i )
			veciGenes[ i ] = DatBig.GetGene( DatSmall.GetGene( i ) );
		for( i = 0; i < DatSmall.GetGenes( ); ++i ) {
			if( ( iOne = veciGenes[ i ] ) == -1 ) {
				for( j = ( i + 1 ); j < DatSmall.GetGenes( ); ++j )
					if( !CMeta::IsNaN( dOne = DatSmall.Get( i, j ) ) )
						aiOne[ ( dOne > 0 ) ? 1 : 0 ]++;
				continue; }
			for( j = ( i + 1 ); j < DatSmall.GetGenes( ); ++j )
				if( !CMeta::IsNaN( dOne = DatSmall.Get( i, j ) ) ) {
					if( ( ( iTwo = veciGenes[ j ] ) == -1 ) || CMeta::IsNaN( dTwo = DatBig.Get( iOne, iTwo ) ) )
						ai = aiOne;
					else {
						DatBig.Set( iOne, iTwo, CMeta::GetNaN( ) );
						ai = ( ( dOne > 0 ) == ( dTwo > 0 ) ) ? aiAgree : aiDisagree; }
					ai[ ( dOne > 0 ) ? 1 : 0 ]++; } }
		for( i = 0; i < DatBig.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < DatBig.GetGenes( ); ++j )
				if( !CMeta::IsNaN( dTwo = DatBig.Get( i, j ) ) )
					aiTwo[ ( dTwo > 0 ) ? 1 : 0 ]++;
	}

	cout << szOne << " +" << endl;
	cout << aiAgree[ 1 ] << '\t' << aiOne[ 1 ] << '\t' << aiDisagree[ 1 ] << endl;
	cout << aiTwo[ 1 ] << "\t\t" << aiTwo[ 0 ] << '\t' << szTwo << " -" << endl;
	cout << aiDisagree[ 0 ] << '\t' << aiOne[ 0 ] << '\t' << aiAgree[ 0 ] << endl;
	cout << szOne << " -" << endl;

	cout << endl;
/*
	dOne = (float)aiAgree[ 1 ] / ( aiAgree[ 1 ] + aiOne[ 1 ] + aiDisagree[ 1 ] );
	dTwo = (float)aiAgree[ 1 ] / ( aiAgree[ 1 ] + aiTwo[ 1 ] + aiDisagree[ 0 ] );
	cout << ( ( dOne + dTwo ) / 2 ) << '\t';
	dOne = (float)aiAgree[ 0 ] / ( aiAgree[ 0 ] + aiOne[ 0 ] + aiDisagree[ 0 ] );
	dTwo = (float)aiAgree[ 0 ] / ( aiAgree[ 0 ] + aiTwo[ 0 ] + aiDisagree[ 1 ] );
	cout << ( ( dOne + dTwo ) / 2 ) << endl;
*/
	cout << ( (float)aiAgree[ 1 ] / ( aiAgree[ 1 ] + aiDisagree[ 0 ] + aiDisagree[ 1 ] ) ) << '\t';
	cout << ( (float)aiAgree[ 0 ] / ( aiAgree[ 0 ] + aiDisagree[ 0 ] + aiDisagree[ 1 ] ) ) << endl;

	CMeta::Shutdown( );
	return 0; }

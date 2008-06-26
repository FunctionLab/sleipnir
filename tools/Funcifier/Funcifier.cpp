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

enum EShared {
	ESharedIgnore	= 0,
	ESharedDiscard	= ESharedIgnore + 1,
	ESharedOneOnly	= ESharedDiscard + 1
};

static const char*	c_aszShared[]	= {"ignore", "discard", "oneonly"};

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info		sArgs;
	CDat					DatIn;
	size_t					iF1, iF2, i, j, iOne, iTwo, iCountIn, iCountOut, iSharedOne, iSharedTwo;
	CGenome					Genome;
	vector<CGenes*>			vecpGenes;
	vector<string>			vecstrNames;
	vector<vector<size_t> >	vecveciGenes;
	float					d, dAveIn, dAveOut;
	ofstream				ofsm;
	EShared					eShared;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta = CMeta( sArgs.verbosity_arg );

	for( i = 0; i < ARRAYSIZE(c_aszShared); ++i )
		if( !strcmp( sArgs.shared_arg, c_aszShared[ i ] ) )
			break;
	if( i >= ARRAYSIZE(c_aszShared) ) {
		cmdline_parser_print_help( );
		return 1; }
	eShared = (EShared)i;

	if( !DatIn.Open( sArgs.input_arg, sArgs.memmap_flag && !sArgs.normalize_flag ) ) {
		cerr << "Could not open: " << sArgs.input_arg << endl;
		return 1; }
	if( sArgs.normalize_flag )
		DatIn.Normalize( CDat::ENormalizeSigmoid );

	vecpGenes.resize( sArgs.inputs_num );
	vecstrNames.resize( vecpGenes.size( ) );
	vecveciGenes.resize( vecpGenes.size( ) );
	for( i = 0; i < vecpGenes.size( ); ++i ) {
		vecstrNames[ i ] = CMeta::Basename( sArgs.inputs[ i ] );
		while( ( j = vecstrNames[ i ].rfind( '.' ) ) != string::npos )
			vecstrNames[ i ] = vecstrNames[ i ].substr( 0, j );
		vecpGenes[ i ] = new CGenes( Genome );
		if( !vecpGenes[ i ]->Open( sArgs.inputs[ i ] ) ) {
			cerr << "Could not open: " << sArgs.inputs[ i ] << endl;
			return 1; }

		vecveciGenes[ i ].resize( vecpGenes[ i ]->GetGenes( ) );
		for( j = 0; j < vecpGenes[ i ]->GetGenes( ); ++j )
			vecveciGenes[ i ][ j ] = DatIn.GetGene( vecpGenes[ i ]->GetGene( j ).GetName( ) ); }

	{
		CDat	DatOut;

		DatOut.Open( vecstrNames );
		for( iF1 = 0; iF1 < vecveciGenes.size( ); ++iF1 )
			for( iF2 = ( iF1 + 1 ); iF2 < vecveciGenes.size( ); ++iF2 ) {
				map<const CGene*, size_t>			mappiGenes;
				map<const CGene*, size_t>::iterator	iterGene;

				for( i = 0; i < vecpGenes[ iF1 ]->GetGenes( ); ++i )
					mappiGenes[ &vecpGenes[ iF1 ]->GetGene( i ) ] = 1;
				for( i = 0; i < vecpGenes[ iF2 ]->GetGenes( ); ++i ) {
					const CGene*	pGene	= &vecpGenes[ iF2 ]->GetGene( i );

					if( ( iterGene = mappiGenes.find( pGene ) ) == mappiGenes.end( ) )
						mappiGenes[ pGene ] = 1;
					else
						iterGene->second++; }

				dAveIn = 0;
				for( iCountIn = i = 0; i < vecveciGenes[ iF1 ].size( ); ++i ) {
					if( ( iOne = vecveciGenes[ iF1 ][ i ] ) == -1 )
						continue;
					iSharedOne = mappiGenes[ &vecpGenes[ iF1 ]->GetGene( i ) ];
					if( ( eShared == ESharedDiscard ) && ( iSharedOne > 1 ) )
						continue;
					for( j = 0; j < vecveciGenes[ iF2 ].size( ); ++j ) {
						iSharedTwo = mappiGenes[ &vecpGenes[ iF2 ]->GetGene( j ) ];
						switch( eShared ) {
							case ESharedDiscard:
								if( iSharedTwo > 1 )
									continue;
								break;

							case ESharedOneOnly:
								if( ( iSharedOne > 1 ) && ( iSharedTwo > 1 ) )
									continue;
								break; }

						if( ( ( iTwo = vecveciGenes[ iF2 ][ j ] ) != -1 ) &&
							!CMeta::IsNaN( d = DatIn.Get( iOne, iTwo ) ) ) {
							iCountIn++;
							dAveIn += d; } } }
				DatOut.Set( iF1, iF2, dAveIn / iCountIn ); }
		DatOut.Save( sArgs.output_arg );
	}

	if( sArgs.colors_arg ) {
		ofsm.open( sArgs.colors_arg );
		if( !ofsm.is_open( ) ) {
			cerr << "Could not open: " << sArgs.colors_arg << endl;
			return 1; }
		for( i = 0; i < vecpGenes.size( ); ++i ) {
			dAveIn = dAveOut = 0;
			for( iCountIn = iCountOut = iF1 = 0; iF1 < vecpGenes[ i ]->GetGenes( ); ++iF1 ) {
				if( ( iOne = vecveciGenes[ i ][ iF1 ] ) == -1 )
					continue;
				for( iF2 = ( iF1 + 1 ); iF2 < vecpGenes[ i ]->GetGenes( ); ++iF2 )
					if( ( ( iTwo = vecveciGenes[ i ][ iF2 ] ) != -1 ) &&
						!CMeta::IsNaN( d = DatIn.Get( iOne, iTwo ) ) ) {
//cerr << DatIn.GetGene( iOne ) << '\t' << DatIn.GetGene( iTwo ) << '\t' << d << endl;
						iCountIn++;
						dAveIn += d; }
				for( j = 0; j < DatIn.GetGenes( ); ++j )
					if( !vecpGenes[ i ]->IsGene( DatIn.GetGene( j ) ) &&
						!CMeta::IsNaN( d = DatIn.Get( iOne, j ) ) ) {
						iCountOut++;
						dAveOut += d; } }
			dAveIn /= iCountIn;
			dAveOut /= iCountOut;
			ofsm << ( dAveIn - dAveOut ) << '\t' << dAveIn << '\t' << dAveOut << '\t' <<
				vecstrNames[ i ] << endl; }
		ofsm.close( ); }

	for( i = 0; i < vecpGenes.size( ); ++i )
		delete vecpGenes[ i ];

	return 0; }

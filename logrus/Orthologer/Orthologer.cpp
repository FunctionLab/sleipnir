#include "stdafx.h"
#include "cmdline.h"

static const size_t	c_iConfirmation	= 512;

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info		sArgs;
	COrthology				Orthology;
	ifstream				ifsm;
	size_t					i, j, iOrganism, iCluster, iOne, iTwo;
	vector<size_t>			veciHits, veciMap;
	vector<CDatasetCompact>	vecData;
	vector<string>			vecstrData;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( sArgs.input_arg )
		ifsm.open( sArgs.input_arg );
	if( !Orthology.Open( sArgs.input_arg ? ifsm : cin ) ) {
		cerr << "Could not open: " << ( sArgs.input_arg || "standard input" ) << endl;
		return 1; }
	if( sArgs.input_arg )
		ifsm.close( );

	veciHits.resize( Orthology.GetGenomes( ) );
	veciMap.resize( sArgs.inputs_num );
	for( iOrganism = 0; iOrganism < sArgs.inputs_num; ++iOrganism ) {
		CDat	Answers;

		ifsm.clear( );
		ifsm.open( sArgs.inputs[ iOrganism ] );
		if( !Answers.OpenGenes( ifsm, true ) ) {
			cerr << "Could not open: " << sArgs.input_arg << endl;
			return 1; }
		ifsm.close( );

		for( i = 0; i < Orthology.GetGenomes( ); ++i ) {
			const CGenome*	pGenome	= Orthology.GetGenome( i );

			veciHits[ i ] = 0;
			for( j = 0; j < min( Answers.GetGenes( ), c_iConfirmation ); ++j )
				if( pGenome->FindGene( Answers.GetGene( j ) ) != -1 )
					veciHits[ i ]++; }

		veciMap[ iOrganism ] = 0;
		for( i = 1; i < veciHits.size( ); ++i )
			if( veciHits[ i ] > veciHits[ veciMap[ iOrganism ] ] )
				veciMap[ iOrganism ] = i; }

//	for( i = 0; i < sArgs.inputs_num; ++i )
//		cerr << sArgs.inputs[ i ] << ": " << Orthology.GetOrganism( veciMap[ i ] ) << endl;

	vecData.resize( sArgs.inputs_num );
	vecstrData.resize( 1 );
	for( i = 0; i < vecData.size( ); ++i ) {
		vecstrData[ 0 ] = sArgs.inputs[ i ];
		if( !vecData[ i ].Open( vecstrData ) ) {
			cerr << "Could not open inputs" << endl;
			return 1; } }

	for( iOrganism = 0; iOrganism < sArgs.inputs_num; ++iOrganism ) {
		CDat						Dat;
		vector<string>				vecstrGenes;
		map<const CGene*,size_t>	mapiGenes;
		size_t						iA, iGenome;

		for( i = 0; i < vecData[ iOrganism ].GetGenes( ); ++i )
			vecstrGenes.push_back( vecData[ iOrganism ].GetGene( i ) );
		for( iCluster = 0; iCluster < Orthology.GetClusters( ); ++iCluster )
			for( iOne = 0; iOne < Orthology.GetGenes( iCluster ); ++iOne ) {
				const CGene&	Gene	= Orthology.GetGene( iCluster, iOne );

				if( Orthology.GetOrganism( Gene ) == veciMap[ iOrganism ] ) {
					if( ( i = vecData[ iOrganism ].GetGene( Gene.GetName( ) ) ) == -1 ) {
						i = vecstrGenes.size( );
						vecstrGenes.push_back( Gene.GetName( ) ); }
					mapiGenes[ &Gene ] = i; } }

		Dat.Open( vecstrGenes, false );
		for( iOne = 0; iOne < vecData[ iOrganism ].GetGenes( ); ++iOne )
			for( iTwo = ( iOne + 1 ); iTwo < vecData[ iOrganism ].GetGenes( ); ++iTwo )
				Dat.Set( iOne, iTwo, ( ( ( i = vecData[ iOrganism ].GetDiscrete( iOne, iTwo, 0 ) ) == -1 ) ||
					!i ) ? CMeta::GetNaN( ) : i );

		for( iCluster = 0; iCluster < Orthology.GetClusters( ); ++iCluster )
			for( iOne = 0; iOne < Orthology.GetGenes( iCluster ); ++iOne ) {
				const CGene&	GeneOne	= Orthology.GetGene( iCluster, iOne );

				if( Orthology.GetOrganism( GeneOne ) != veciMap[ iOrganism ] )
					continue;
				iA = mapiGenes[ &GeneOne ];
				for( iTwo = ( iOne + 1 ); iTwo < Orthology.GetGenes( iCluster ); ++iTwo ) {
					const CGene&	GeneTwo	= Orthology.GetGene( iCluster, iTwo );

					if( ( Orthology.GetOrganism( GeneTwo ) ) == veciMap[ iOrganism ] )
						Dat.Set( iA, mapiGenes[ &GeneTwo ], 1 ); } }

		for( iGenome = 0; iGenome < Orthology.GetGenomes( ); ++iGenome ) {
			if( veciMap[ iGenome ] == iOrganism )
				continue;
			{
				const CDatasetCompact&	Datum	= vecData[ veciMap[ iGenome ] ];

				for( iCluster = 0; 

/*
	for( iOrganism = 0; iOrganism < sArgs.inputs_num; ++iOrganism ) {
		CDat						DatIn, DatOut;
		vector<string>				vecstrGenes;
		map<const CGene*,size_t>	mapiGenes;
		size_t						iOrgOne, iOrgTwo;

		if( !DatIn.Open( sArgs.inputs[ iOrganism ], !!sArgs.memmap_flag ) ) {
			cerr << "Could not open: " << sArgs.inputs[ iOrganism ] << endl;
			return 1; }
		for( i = 0; i < DatIn.GetGenes( ); ++i )
			vecstrGenes.push_back( DatIn.GetGene( i ) );
		for( iCluster = 0; iCluster < Orthology.GetClusters( ); ++iCluster )
			for( iOne = 0; iOne < Orthology.GetGenes( iCluster ); ++iOne ) {
				const CGene&	Gene	= Orthology.GetGene( iCluster, iOne );

				if( Orthology.GetOrganism( Gene ) == veciMap[ iOrganism ] ) {
					if( ( i = DatIn.GetGene( Gene.GetName( ) ) ) == -1 ) {
						i = vecstrGenes.size( );
						vecstrGenes.push_back( Gene.GetName( ) ); }
					mapiGenes[ &Gene ] = i; } }

		DatOut.Open( vecstrGenes, false );
		for( iOne = 0; iOne < DatIn.GetGenes( ); ++iOne )
			for( iTwo = ( iOne + 1 ); iTwo < DatIn.GetGenes( ); ++iTwo )
				DatOut.Set( iOne, iTwo, ( CMeta::IsNaN( d = DatIn.Get( iOne, iTwo ) ) || ( d <= 0 ) ) ?
					CMeta::GetNaN( ) : d );

		for( iCluster = 0; iCluster < Orthology.GetClusters( ); ++iCluster )
			for( iOne = 0; iOne < Orthology.GetGenes( iCluster ); ++iOne ) {
				const CGene&	GeneOne	= Orthology.GetGene( iCluster, iOne );

				iOrgOne = Orthology.GetOrganism( GeneOne );
				iA = mapiGenes[ &GeneOne ];
				for( iTwo = ( iOne + 1 ); iTwo < Orthology.GetGenes( iCluster ); ++iTwo ) {
					const CGene&	GeneTwo	= Orthology.GetGene( iCluster, iTwo );

					if( ( iOrgTwo = Orthology.GetOrganism( GeneTwo ) ) != iOrgOne )
						continue;
					if( iOrgOne == veciMap[ iOrganism ] )
						DatOut.Set( iA, mapiGenes[ &GeneTwo ], 1 ); } } }
*/

	CMeta::Shutdown( );
	return 0; }

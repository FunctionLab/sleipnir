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

static const size_t	c_iConfirmation		= 512;
static const char	c_szOrthologized[]	= "_ortho.dab";

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info		sArgs;
	COrthology				Orthology;
	ifstream				ifsm;
	size_t					i, j, iOrganism, iCluster, iOne, iTwo, iOrgOne, iIndexOne, iIndexTwo;
	size_t					iGeneOne, iGeneTwo;
	vector<size_t>			veciHits, veciOrgsIn2Orth, veciOrgsOrth2In;
	vector<CDatasetCompact>	vecData;
	vector<string>			vecstrData;
	CBinaryMatrix			MatRelated;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta = CMeta( sArgs.verbosity_arg );

	if( sArgs.input_arg )
		ifsm.open( sArgs.input_arg );
	if( !Orthology.Open( sArgs.input_arg ? ifsm : cin ) ) {
		cerr << "Could not open: " << ( sArgs.input_arg || "standard input" ) << endl;
		return 1; }
	if( sArgs.input_arg )
		ifsm.close( );

	veciHits.resize( Orthology.GetGenomes( ) );
	veciOrgsIn2Orth.resize( sArgs.inputs_num );
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

		veciOrgsIn2Orth[ iOrganism ] = 0;
		for( i = 1; i < veciHits.size( ); ++i )
			if( veciHits[ i ] > veciHits[ veciOrgsIn2Orth[ iOrganism ] ] )
				veciOrgsIn2Orth[ iOrganism ] = i; }
	veciOrgsOrth2In.resize( Orthology.GetGenomes( ) );
	for( i = 0; i < veciOrgsIn2Orth.size( ); ++i )
		veciOrgsOrth2In[ veciOrgsIn2Orth[ i ] ] = i;

	for( i = 0; i < sArgs.inputs_num; ++i )
		cerr << sArgs.inputs[ i ] << ": " << Orthology.GetOrganism( veciOrgsIn2Orth[ i ] ) << endl;

	vecData.resize( sArgs.inputs_num );
	vecstrData.resize( 1 );
	for( i = 0; i < vecData.size( ); ++i ) {
		vecstrData[ 0 ] = sArgs.inputs[ i ];
		if( !vecData[ i ].Open( vecstrData ) ) {
			cerr << "Could not open inputs" << endl;
			return 1; } }

	vector<vector<pair<size_t, size_t> > >	vecvecpriiOrthology;
	vecvecpriiOrthology.resize( Orthology.GetClusters( ) );
	for( i = 0; i < vecvecpriiOrthology.size( ); ++i ) {
		vecvecpriiOrthology[ i ].resize( Orthology.GetGenes( i ) );
		for( j = 0; j < vecvecpriiOrthology[ i ].size( ); ++j ) {
			const CGene&			GeneOne	= Orthology.GetGene( i, j );
			const CDatasetCompact&	Dataset	= vecData[ iOrgOne = veciOrgsOrth2In[
				Orthology.GetOrganism( GeneOne ) ] ];

			vecvecpriiOrthology[ i ][ j ].first = Dataset.GetGene( GeneOne.GetName( ) );
			vecvecpriiOrthology[ i ][ j ].second = iOrgOne; } }

	MatRelated.Initialize( vecvecpriiOrthology.size( ) );
	for( iOne = 0; iOne < MatRelated.GetSize( ); ++iOne ) {
		if( !( iOne % 1000 ) )
			cerr << "Linking group " << iOne << '/' << MatRelated.GetSize( ) << endl;
		for( iTwo = ( iOne + 1 ); iTwo < MatRelated.GetSize( ); ++iTwo )
			for( iGeneOne = 0; iGeneOne < vecvecpriiOrthology[ iOne ].size( ); ++iGeneOne ) {
				const pair<size_t, size_t>&	priiOne	= vecvecpriiOrthology[ iOne ][ iGeneOne ];
				const CDatasetCompact&		Dataset	= vecData[ priiOne.second ];

				if( priiOne.first == -1 )
					continue;
				for( iGeneTwo = 0; iGeneTwo < Orthology.GetGenes( iTwo ); ++iGeneTwo ) {
					const pair<size_t, size_t>&	priiTwo	= vecvecpriiOrthology[ iTwo ][ iGeneTwo ];

					if( ( priiOne.second == priiTwo.second ) && ( priiTwo.first != -1 ) &&
						( ( i = Dataset.GetDiscrete( priiOne.first, priiTwo.first, 0 ) ) != -1 ) && i ) {
						if( sArgs.verbosity_arg > 6 )
							cerr << "Linked " << iOne << ':' << iTwo << " by " <<
								sArgs.inputs[ priiOne.second ] << ' ' << Dataset.GetGene(
								priiOne.first ) << ':' << Dataset.GetGene( priiTwo.first ) << endl;
						MatRelated.Set( iOne, iTwo, true );
						break; } }
				if( iGeneTwo < Orthology.GetGenes( iTwo ) )
					break; } }

	for( iOrganism = 0; iOrganism < sArgs.inputs_num; ++iOrganism ) {
		CDat						Dat;
		vector<string>				vecstrGenes;
		map<const CGene*,size_t>	mapGenes;

		cerr << "Processing " << sArgs.inputs[ iOrganism ] << endl;
		vecstrGenes.resize( vecData[ iOrganism ].GetGenes( ) );
		for( i = 0; i < vecstrGenes.size( ); ++i )
			vecstrGenes[ i ] = vecData[ iOrganism ].GetGene( i );
		for( iCluster = 0; iCluster < Orthology.GetClusters( ); ++iCluster )
			for( iOne = 0; iOne < Orthology.GetGenes( iCluster ); ++iOne ) {
				const CGene&	Gene	= Orthology.GetGene( iCluster, iOne );

				if( ( Orthology.GetOrganism( Gene ) == veciOrgsIn2Orth[ iOrganism ] ) &&
					( vecData[ iOrganism ].GetGene( Gene.GetName( ) ) == -1 ) ) {
					mapGenes[ &Gene ] = vecstrGenes.size( );
					vecstrGenes.push_back( Gene.GetName( ) ); } }

		Dat.Open( vecstrGenes, false );
		cerr << "Integrating orthology clusters" << endl;
		for( iOne = 0; iOne < vecData[ iOrganism ].GetGenes( ); ++iOne )
			for( iTwo = ( iOne + 1 ); iTwo < vecData[ iOrganism ].GetGenes( ); ++iTwo )
				Dat.Set( iOne, iTwo, ( ( ( i = vecData[ iOrganism ].GetDiscrete( iOne, iTwo, 0 ) ) ==
					-1 ) || !i ) ? CMeta::GetNaN( ) : i );
		for( iOne = 0; iOne < vecvecpriiOrthology.size( ); ++iOne )
			for( iGeneOne = 0; iGeneOne < vecvecpriiOrthology[ iOne ].size( ); ++iGeneOne ) {
				const pair<size_t, size_t>&	priiOne	= vecvecpriiOrthology[ iOne ][ iGeneOne ];

				if( priiOne.second != iOrganism )
					continue;
				iIndexOne = ( priiOne.first == -1 ) ? mapGenes[ &Orthology.GetGene( iOne,
					iGeneOne ) ] : priiOne.first;
				for( iGeneTwo = ( iGeneOne + 1 ); iGeneTwo < vecvecpriiOrthology[ iOne ].size( );
					++iGeneTwo ) {
					const pair<size_t, size_t>&	priiTwo	= vecvecpriiOrthology[ iOne ][ iGeneTwo ];

					if( priiOne.second == priiTwo.second ) {
						iIndexTwo = ( priiTwo.first == -1 ) ? mapGenes[ &Orthology.GetGene( iOne,
							iGeneTwo ) ] : priiTwo.first;
						if( CMeta::IsNaN( Dat.Get( iIndexOne, iIndexTwo ) ) )
							Dat.Set( iIndexOne, iIndexTwo, (float)sArgs.weight1_arg ); } } }

		cerr << "Integrating orthology links" << endl;
		for( iOne = 0; iOne < MatRelated.GetSize( ); ++iOne ) {
			if( !( iOne % 1000 ) )
				cerr << "Linking group " << iOne << '/' << MatRelated.GetSize( ) << endl;
			for( iTwo = ( iOne + 1 ); iTwo < MatRelated.GetSize( ); ++iTwo )
				if( MatRelated.Get( iOne, iTwo ) )
					for( iGeneOne = 0; iGeneOne < vecvecpriiOrthology[ iOne ].size( ); ++iGeneOne ) {
						const pair<size_t, size_t>&	priiOne	= vecvecpriiOrthology[ iOne ][ iGeneOne ];

						if( priiOne.second != iOrganism )
							continue;
						iIndexOne = ( priiOne.first == -1 ) ? mapGenes[ &Orthology.GetGene( iOne,
							iGeneOne ) ] : priiOne.first;
						for( iGeneTwo = 0; iGeneTwo < vecvecpriiOrthology[ iTwo ].size( ); ++iGeneTwo ) {
							const pair<size_t, size_t>&	priiTwo	= vecvecpriiOrthology[ iTwo ][ iGeneTwo ];

							if( priiOne.second == priiTwo.second ) {
								iIndexTwo = ( priiTwo.first == -1 ) ? mapGenes[ &Orthology.GetGene(
									iTwo, iGeneTwo ) ] : priiTwo.first;
								if( CMeta::IsNaN( Dat.Get( iIndexOne, iIndexTwo ) ) )
									Dat.Set( iIndexOne, iIndexTwo, (float)sArgs.weight2_arg ); } } } }

		Dat.Save( ( CMeta::Deextension( sArgs.inputs[ iOrganism ] ) + c_szOrthologized ).c_str( ) ); }

	return 0; }

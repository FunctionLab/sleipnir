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

static const char	c_acQTC[]	= "qtc";

int main( int iArgs, char** aszArgs ) {
	size_t						i, j;
	ifstream					ifsm;
	CDat						DatIn, DatOut;
	vector<uint16_t>			vecsClusters;
	const CPCL*					pPCL;
	uint16_t					sClusters;
	CPCL						PCL, Ranks, Weights;
	gengetopt_args_info			sArgs;
	IMeasure*					pMeasure;
	CMeasurePearson				Pearson;
	CMeasureEuclidean			Euclidean;
	CMeasureKendallsTau			KendallsTau;
	CMeasureKolmogorovSmirnov	KolmSmir;
	CMeasureSpearman			Spearman( true );
	CMeasureQuickPearson		QuickPear;
	CMeasureNegate				EuclideanNeg( &Euclidean, false );
	IMeasure*					apMeasures[]	= { &Pearson, &EuclideanNeg, &KendallsTau,
		&KolmSmir, &Spearman, &QuickPear, NULL };

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta = CMeta( sArgs.verbosity_arg, sArgs.random_arg );

	pMeasure = NULL;
	for( i = 0; apMeasures[ i ]; ++i )
		if( !strcmp( apMeasures[ i ]->GetName( ), sArgs.distance_arg ) ) {
			if( ( pMeasure = apMeasures[ i ] ) == &EuclideanNeg )
				sArgs.normalize_flag = true;
			break; }
	if( !pMeasure ) {
		cmdline_parser_print_help( );
		return 1; }

	CMeasureAutocorrelate	Autocorrelate( pMeasure, false );
	if( sArgs.autocorrelate_flag )
		pMeasure = &Autocorrelate;

	if( sArgs.input_arg ) {
		ifsm.open( sArgs.input_arg );
		if( PCL.Open( ifsm, sArgs.skip_arg ) ) {
			ifsm.close( );
			cerr << "Opened PCL: " << sArgs.input_arg << endl; }
		else {
			ifsm.close( );
			if( !DatIn.Open( sArgs.input_arg ) ) {
				cerr << "Could not open input: " << sArgs.input_arg << endl;
				return 1; }
			if( !PCL.Open( cin, sArgs.skip_arg ) ) {
				cerr << "Could not open PCL" << endl;
				return 1; } } }
	else if( PCL.Open( cin, sArgs.skip_arg ) )
		cerr << "Opened PCL" << endl;
	else {
		cerr << "Could not open PCL" << endl;
		return 1; }

	if( sArgs.weights_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.weights_arg );
		if( !Weights.Open( ifsm, sArgs.skip_arg ) ) {
			cerr << "Could not open: " << sArgs.weights_arg << endl;
			return 1; }
		ifsm.close( );

		if( ( Weights.GetExperiments( ) != PCL.GetExperiments( ) ) ||
			( Weights.GetGenes( ) != PCL.GetGenes( ) ) ) {
			cerr << "Illegal data sizes: " << PCL.GetExperiments( ) << 'x' << PCL.GetGenes( ) << ", " <<
				Weights.GetExperiments( ) << 'x' << Weights.GetGenes( ) << endl;
			return 1; } }

	if( pMeasure->IsRank( ) ) {
		Ranks.Open( PCL );
		Ranks.RankTransform( );
		pPCL = &Ranks; }
	else
		pPCL = &PCL;

	if( sArgs.delta_arg ) {
		DatOut.Open( pPCL->GetGeneNames( ) );
		if( DatIn.GetGenes( ) )
			CClustQTC::Cluster( DatIn.Get( ), (float)sArgs.diamineter_arg, (float)sArgs.diameter_arg,
				(float)sArgs.delta_arg, sArgs.size_arg, DatOut.Get( ) );
		else
			CClustQTC::Cluster( pPCL->Get( ), pMeasure, (float)sArgs.diamineter_arg,
				(float)sArgs.diameter_arg, (float)sArgs.delta_arg, sArgs.size_arg, DatOut.Get( ),
				sArgs.weights_arg ? &Weights.Get( ) : NULL );
		DatOut.Save( sArgs.output_arg ); }
	else {
		if( !strcmp( sArgs.algorithm_arg, c_acQTC ) )
			sClusters = DatIn.GetGenes( ) ?
				CClustQTC::Cluster( DatIn.Get( ), (float)sArgs.diameter_arg, sArgs.size_arg, vecsClusters ) :
				CClustQTC::Cluster( pPCL->Get( ), pMeasure, (float)sArgs.diameter_arg, sArgs.size_arg,
				vecsClusters, sArgs.weights_arg ? &Weights.Get( ) : NULL );
		else
			sClusters = CClustKMeans::Cluster( pPCL->Get( ), pMeasure, sArgs.size_arg, vecsClusters,
				sArgs.weights_arg ? &Weights.Get( ) : NULL ) ? sArgs.size_arg : 0;

		if( !sClusters )
			cerr << "No clusters found" <<endl;
		else if( sArgs.output_arg ) {
			DatOut.Open( pPCL->GetGeneNames( ) );
			for( i = 0; i < vecsClusters.size( ); ++i ) {
				if( ( vecsClusters[ i ] + 1 ) == sClusters )
					continue;
				for( j = ( i + 1 ); j < vecsClusters.size( ); ++j ) {
					if( ( vecsClusters[ j ] + 1 ) == sClusters )
						continue;
					DatOut.Set( i, j, ( vecsClusters[ i ] == vecsClusters[ j ] ) ? 1.0f : 0.0f ); } }
			DatOut.Save( sArgs.output_arg ); }
		else
			for( i = 0; i < vecsClusters.size( ); ++i )
				cout << pPCL->GetGene( i ) << '\t' << vecsClusters[ i ] << endl; }

	return 0; }

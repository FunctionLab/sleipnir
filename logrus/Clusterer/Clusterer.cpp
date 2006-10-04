#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	size_t						i, j;
	ofstream					ofsm;
	ifstream					ifsm;
	istream*					pistm;
	CDat						Dat;
	vector<uint16_t>			vecsClusters;
	const CPCL*					pPCL;
	uint16_t					sClusters;
	CPCL						PCL, Ranks;
	gengetopt_args_info			sArgs;
	IMeasure*					pMeasure;
	CMeasurePearson				Pearson;
	CMeasureEuclidean			Euclidean;
	CMeasureKendallsTau			KendallsTau;
	CMeasureKolmogorovSmirnov	KolmSmir;
	CMeasureSpearman			Spearman;
	CMeasureQuickPearson		QuickPear;
	CMeasureNegate				EuclideanNeg( &Euclidean );
	IMeasure*					apMeasures[]	= { &Pearson, &EuclideanNeg, &KendallsTau,
		&KolmSmir, &Spearman, &QuickPear, NULL };

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg, sArgs.random_arg );

	pMeasure = NULL;
	for( i = 0; apMeasures[ i ]; ++i )
		if( !strcmp( apMeasures[ i ]->GetName( ), sArgs.distance_arg ) ) {
			if( ( pMeasure = apMeasures[ i ] ) == &EuclideanNeg )
				sArgs.normalize_flag = true;
			break; }
	if( !pMeasure ) {
		cmdline_parser_print_help( );
		return 1; }

	if( sArgs.input_arg ) {
		ifsm.open( sArgs.input_arg );
		pistm = &ifsm; }
	else
		pistm = &cin;
	if( !PCL.Open( *pistm, sArgs.skip_arg ) ) {
		cerr << "Could not open PCL" << endl;
		return 1; }
	if( sArgs.input_arg )
		ifsm.close( );

	if( pMeasure->IsRank( ) ) {
		Ranks.Open( PCL );
		Ranks.RankTransform( );
		pPCL = &Ranks; }
	else
		pPCL = &PCL;

	if( sArgs.delta_arg ) {
		Dat.Open( pPCL->GetGeneNames( ) );
		CClustQTC::Cluster( pPCL->Get( ), pMeasure, (float)sArgs.diamineter_arg, (float)sArgs.diameter_arg,
			(float)sArgs.delta_arg, sArgs.size_arg, !!sArgs.autocorrelate_flag, Dat.Get( ) );
		if( sArgs.output_arg ) {
			ofsm.open( sArgs.output_arg, ios_base::binary );
			Dat.Save( ofsm, true );
			ofsm.close( ); }
		else {
			Dat.Save( cout, false );
			cout.flush( ); } }
	else {
		sClusters = CClustQTC::Cluster( pPCL->Get( ), pMeasure, (float)sArgs.diameter_arg, sArgs.size_arg,
			!!sArgs.autocorrelate_flag, vecsClusters );

		if( sArgs.output_arg ) {
			Dat.Open( pPCL->GetGeneNames( ) );
			for( i = 0; i < vecsClusters.size( ); ++i ) {
				if( ( vecsClusters[ i ] + 1 ) == sClusters )
					continue;
				for( j = ( i + 1 ); j < vecsClusters.size( ); ++j ) {
					if( ( vecsClusters[ j ] + 1 ) == sClusters )
						continue;
					Dat.Set( i, j, ( vecsClusters[ i ] == vecsClusters[ j ] ) ? 1.0f : 0.0f ); } }
			ofsm.open( sArgs.output_arg, ios_base::binary );
			Dat.Save( ofsm, true );
			ofsm.close( ); }
		else
			for( i = 0; i < vecsClusters.size( ); ++i )
				cout << pPCL->GetGene( i ) << '\t' << vecsClusters[ i ] << endl; }

	CMeta::Shutdown( );
	return 0; }

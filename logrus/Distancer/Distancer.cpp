#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	CPCL						PCL;
	CDat						Dat;
	CDistanceMatrix				Dist;
	size_t						i, j;
	ofstream					ofsm;
	ifstream					ifsm;
	vector<string>				vecstrGenes;
	gengetopt_args_info			sArgs;
	IMeasure*					pMeasure;
	CMeasurePearson				Pearson;
	CMeasureEuclidean			Euclidean;
	CMeasureKendallsTau			KendallsTau;
	CMeasureKolmogorovSmirnov	KolmSmir;
	CMeasureSpearman			Spearman( true );
	CMeasureNegate				EuclideanNeg( &Euclidean );
	IMeasure*					apMeasures[]	= { &Pearson, &EuclideanNeg, &KendallsTau,
		&KolmSmir, &Spearman, NULL };

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }

	CMeta::Startup( sArgs.verbosity_arg );

	pMeasure = NULL;
	for( i = 0; apMeasures[ i ]; ++i )
		if( !strcmp( apMeasures[ i ]->GetName( ), sArgs.distance_arg ) ) {
			if( ( pMeasure = apMeasures[ i ] ) == &EuclideanNeg )
				sArgs.normalize_flag = true;
			break; }
	if( !pMeasure ) {
		cmdline_parser_print_help( );
		return 1; }

	CMeasureAutocorrelate		Autocorrelate( pMeasure );
	if( sArgs.autocorrelate_flag )
		pMeasure = &Autocorrelate;

	if( sArgs.input_arg ) {
		ifsm.open( sArgs.input_arg );
		if( !PCL.Open( ifsm, sArgs.skip_arg ) ) {
			cerr << "Could not open input: " << sArgs.input_arg << endl;
			return 1; }
		ifsm.close( );
		cerr << "Opened PCL: " << sArgs.input_arg << endl; }
	else if( !PCL.Open( cin, sArgs.skip_arg ) ) {
		cerr << "Could not open PCL" << endl;
		return 1; }

	if( ( pMeasure == &Spearman ) || ( pMeasure == &KendallsTau ) )
		PCL.RankTransform( );
	Dist.Initialize( PCL.GetGenes( ) );
	for( i = 0; i < PCL.GetGenes( ); ++i ) {
		if( !( i % 100 ) )
			cerr << "Processing gene " << (unsigned int)i << '/' <<
				(unsigned int)PCL.GetGenes( ) << endl;
		for( j = ( i + 1 ); j < PCL.GetGenes( ); ++j )
			Dist.Set( i, j, (float)pMeasure->Measure( PCL.Get( i ), PCL.GetExperiments( ),
				PCL.Get( j ), PCL.GetExperiments( ) ) ); }

/*
	if( sArgs.zscore_flag ) {
		size_t	iN;
		float	d, dAve, dDev;

		dAve = dDev = 0;
		for( iN = i = 0; i < Dist.GetSize( ); ++i )
			for( j = ( i + 1 ); j < Dist.GetSize( ); ++j )
				if( !CMeta::IsNaN( d = Dist.Get( i, j ) ) ) {
					iN++;
					dAve += d;
					dDev += d * d; }
		dAve /= iN;
		dDev = sqrt( ( dDev / iN ) - ( dAve * dAve ) );
		for( i = 0; i < Dist.GetSize( ); ++i )
			for( j = ( i + 1 ); j < Dist.GetSize( ); ++j )
				if( !CMeta::IsNaN( d = Dist.Get( i, j ) ) )
					Dist.Set( i, j, ( d - dAve ) / dDev ); }
*/

	Dat.Open( PCL.GetGeneNames( ), Dist );
	if( sArgs.normalize_flag )
		Dat.Normalize( );
	ofsm.open( sArgs.output_arg, ios_base::binary );
	Dat.Save( ofsm, true );
	ofsm.close( );

	CMeta::Shutdown( );
	return 0; }

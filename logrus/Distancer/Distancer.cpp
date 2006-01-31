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
	CMeasurePearNorm			PearNorm;
	IMeasure*					apMeasures[]	= { &Pearson, &EuclideanNeg, &KendallsTau,
		&KolmSmir, &Spearman, &PearNorm, NULL };

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

	if( pMeasure->IsRank( ) )
		PCL.RankTransform( );
	Dist.Initialize( PCL.GetGenes( ) );
	for( i = 0; i < PCL.GetGenes( ); ++i ) {
		if( !( i % 100 ) )
			cerr << "Processing gene " << (unsigned int)i << '/' <<
				(unsigned int)PCL.GetGenes( ) << endl;
		for( j = ( i + 1 ); j < PCL.GetGenes( ); ++j )
			Dist.Set( i, j, (float)pMeasure->Measure( PCL.Get( i ), PCL.GetExperiments( ),
				PCL.Get( j ), PCL.GetExperiments( ) ) ); }

	Dat.Open( PCL.GetGeneNames( ), Dist );
	if( sArgs.normalize_flag || sArgs.zscore_flag )
		Dat.Normalize( !!sArgs.normalize_flag );
	ofsm.open( sArgs.output_arg, ios_base::binary );
	Dat.Save( ofsm, true );
	ofsm.close( );

	CMeta::Shutdown( );
	return 0; }

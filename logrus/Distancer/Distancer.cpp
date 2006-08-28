#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	CPCL						PCL;
	CDat						Dat;
	CDistanceMatrix				Dist;
	size_t						i, j, iOne, iTwo;
	float						d;
	ofstream					ofsm;
	ifstream					ifsm;
	vector<string>				vecstrGenes;
	CGenome						Genome;
	CGenes						GenesIn( Genome );
	vector<size_t>				veciGenes;
	const float*				adOne;
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

	if( sArgs.genes_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.genes_arg );
		if( !GenesIn.Open( ifsm ) ) {
			cerr << "Couldn't open: " << sArgs.genes_arg << endl;
			return 1; }
		ifsm.close( ); }
	else
		GenesIn.Open( PCL.GetGeneNames( ) );
	veciGenes.resize( GenesIn.GetGenes( ) );
	for( i = 0; i < veciGenes.size( ); ++i )
		veciGenes[ i ] = PCL.GetGene( GenesIn.GetGene( i ).GetName( ) );

	if( pMeasure->IsRank( ) )
		PCL.RankTransform( );
	Dist.Initialize( GenesIn.GetGenes( ) );
	for( i = 0; i < Dist.GetSize( ); ++i )
		for( j = ( i + 1 ); j < Dist.GetSize( ); ++j )
			Dist.Set( i, j, CMeta::GetNaN( ) );
	for( i = 0; i < GenesIn.GetGenes( ); ++i ) {
		if( !( i % 100 ) )
			cerr << "Processing gene " << i << '/' << GenesIn.GetGenes( ) << endl;
		if( ( iOne = veciGenes[ i ] ) == -1 )
			continue;
		adOne = PCL.Get( iOne );
		for( j = ( i + 1 ); j < GenesIn.GetGenes( ); ++j )
			if( ( iTwo = veciGenes[ j ] ) != -1 )
				Dist.Set( i, j, (float)pMeasure->Measure(
					adOne, PCL.GetExperiments( ), PCL.Get( iTwo ), PCL.GetExperiments( ) ) ); }

	Dat.Open( GenesIn.GetGeneNames( ), Dist );
	if( sArgs.normalize_flag || sArgs.zscore_flag )
		Dat.Normalize( !!sArgs.normalize_flag );
	if( sArgs.cutoff_given )
		for( i = 0; i < Dat.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
				if( !CMeta::IsNaN( d = Dat.Get( i, j ) ) && ( d < sArgs.cutoff_arg ) )
					Dat.Set( i, j, CMeta::GetNaN( ) );
	ofsm.open( sArgs.output_arg, ios_base::binary );
	Dat.Save( ofsm, true );
	ofsm.close( );

	CMeta::Shutdown( );
	return 0; }

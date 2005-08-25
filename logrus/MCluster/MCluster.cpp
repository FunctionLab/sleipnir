#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	CHierarchy*					pHier;
	CPCL						PCL, Weights;
	const CPCL*					pPCL;
	CDat						Dat;
	size_t						i, j, iGene;
	ofstream					ofsm;
	ifstream					ifsm;
	vector<size_t>				veciGenes, veciPCL;
	vector<float>				vecdPCL;
	vector<string>				vecstrGenes;
	gengetopt_args_info			sArgs;
	IMeasure*					pMeasure;
	CMeasurePearson				Pearson;
	CMeasureEuclidean			Euclidean;
	CMeasureKendallsTau			KendallsTau;
	CMeasureKolmogorovSmirnov	KolmSmir;
	CMeasureSpearman			Spearman;
	CMeasureNegate				EuclideanNeg( &Euclidean );
	IMeasure*					apMeasures[]	= { &Pearson, &EuclideanNeg, &KendallsTau,
		&KolmSmir, &Spearman, NULL };

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }

	CMeta::Startup( sArgs.verbosity_arg );
	if( sArgs.input_arg ) {
		if( Dat.Open( sArgs.input_arg ) ) {
			cerr << "Opened dat: " << sArgs.input_arg << endl;
			if( !PCL.Open( cin, sArgs.skip_arg ) ) {
				cerr << "Could not open PCL" << endl;
				return 1; } }
		else {
			ifsm.open( sArgs.input_arg );
			if( !PCL.Open( ifsm, sArgs.skip_arg ) ) {
				cerr << "Could not open input: " << sArgs.input_arg << endl;
				return 1; }
			ifsm.close( );
			cerr << "Opened PCL: " << sArgs.input_arg << endl; } }
	else if( !PCL.Open( cin, sArgs.skip_arg ) ) {
		cerr << "Could not open PCL" << endl;
		return 1; }

	if( !Dat.GetGenes( ) ) {
		pMeasure = NULL;
		for( i = 0; apMeasures[ i ]; ++i )
			if( !strcmp( apMeasures[ i ]->GetName( ), sArgs.distance_arg ) ) {
				if( ( pMeasure = apMeasures[ i ] ) == &EuclideanNeg )
					sArgs.normalize_flag = true;
				break; }
		if( !pMeasure ) {
			cmdline_parser_print_help( );
			return 1; }

		if( sArgs.weights_arg ) {
			ifsm.clear( );
			ifsm.open( sArgs.weights_arg );
			if( !Weights.Open( ifsm, sArgs.skip_arg ) ) {
				cerr << "Couldn't open: " << sArgs.weights_arg << endl;
				return 1; }
			ifsm.close( );

			if( ( Weights.GetExperiments( ) != PCL.GetExperiments( ) ) ||
				( Weights.GetGenes( ) != PCL.GetGenes( ) ) ) {
				cerr << "Illegal data sizes: " << (unsigned int)PCL.GetExperiments( ) <<
					'x' << (unsigned int)PCL.GetGenes( ) << ", " <<
					(unsigned int)Weights.GetExperiments( ) << 'x' <<
					(unsigned int)Weights.GetGenes( ) << endl;
				return 1; } }

		{
			CPCL	Ranks;

			if( ( pMeasure == &Spearman ) || ( pMeasure == &KendallsTau ) ) {
				Ranks.Open( PCL );
				Ranks.RankTransform( );
				pPCL = &Ranks; }
			else
				pPCL = &PCL;
			Dat.Open( PCL.GetGeneNames( ) );
			for( i = 0; i < PCL.GetGenes( ); ++i )
				for( j = ( i + 1 ); j < PCL.GetGenes( ); ++j )
					Dat.Set( i, j, (float)pMeasure->Measure( pPCL->Get( i ),
						pPCL->GetExperiments( ), pPCL->Get( j ), pPCL->GetExperiments( ),
						IMeasure::EMapCenter, sArgs.weights_arg ? Weights.Get( i ) : NULL,
						sArgs.weights_arg ? Weights.Get( j ) : NULL ) );
		} }

	for( i = 0; i < Dat.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
			if( CMeta::IsNaN( Dat.Get( i, j ) ) )
				Dat.Set( i, j, 0 );
	if( sArgs.normalize_flag )
		Dat.Normalize( );
	pHier = CClustHierarchical::Cluster( Dat.Get( ) );

	vecdPCL.resize( Dat.GetGenes( ) );
	for( i = 0; i < vecdPCL.size( ); ++i ) {
		if( ( iGene = PCL.GetGene( Dat.GetGene( i ) ) ) == -1 ) {
			vecdPCL[ i ] = CMeta::GetNaN( );
			continue; }
		vecdPCL[ i ] = 0;
		for( j = 0; j < PCL.GetExperiments( ); ++j )
			vecdPCL[ i ] += PCL.Get( iGene, j );
		vecdPCL[ i ] /= PCL.GetExperiments( ); }
	pHier->SortChildren( vecdPCL );

	pHier->GetGenes( veciGenes );
	veciPCL.resize( PCL.GetGenes( ) );
	for( i = 0; i < veciGenes.size( ); ++i )
		veciPCL[ i ] = PCL.GetGene( Dat.GetGene( veciGenes[ i ] ) );
	for( j = 0; j < PCL.GetGenes( ); ++j )
		if( Dat.GetGene( PCL.GetGene( j ) ) == -1 )
			veciPCL[ i++ ] = j;
	PCL.SortGenes( veciPCL );
	for( i = 0; i < veciGenes.size( ); ++i )
		veciPCL[ i ] = veciGenes[ i ];
	for( ; i < veciPCL.size( ); ++i )
		veciPCL[ i ] = i;

	ofsm.open( sArgs.output_arg );
	pHier->Save( ofsm, Dat.GetGenes( ) );
	ofsm.close( );
	pHier->Destroy( );
	PCL.Save( cout, &veciPCL );
	cout.flush( );

	CMeta::Shutdown( );
	return 0; }

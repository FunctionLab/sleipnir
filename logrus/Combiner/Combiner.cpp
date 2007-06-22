#include "stdafx.h"

static int MainDABs( const gengetopt_args_info& );
static int MainDATs( const gengetopt_args_info& );
static int MainPCLs( const gengetopt_args_info& );

static const TPFnCombiner	c_apfnCombiners[]	= { MainPCLs, MainDATs, MainDABs, NULL };
static const char*			c_aszCombiners[]	= { "pcl", "dat", "dab", NULL };
static const char			c_szMean[]			= "mean";
static const char			c_szGMean[]			= "gmean";
static const char			c_szHMean[]			= "hmean";
static const char			c_szMax[]			= "max";
static const char			c_szMin[]			= "min";
static const char			c_szVote[]			= "vote";

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	int					iRet;
	size_t				i;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	for( i = 0; c_aszCombiners[ i ]; ++i )
		if( !strcmp( c_aszCombiners[ i ], sArgs.type_arg ) ) {
			iRet = c_apfnCombiners[ i ]( sArgs );
			break; }

	CMeta::Shutdown( );
	return iRet; }

static int MainPCLs( const gengetopt_args_info& sArgs ) {
	CPCL						PCL, PCLNew;
	size_t						i, j, iArg, iExp, iGene;
	vector<string>				vecstrGenes, vecstrExps, vecstrFeatures;
	set<string>					setstrGenes;
	set<string>::const_iterator	iterGenes;
	ifstream					ifsm;
	ofstream					ofsm;

	for( iArg = 0; iArg < sArgs.inputs_num; ++iArg ) {
		ifsm.clear( );
		ifsm.open( sArgs.inputs[ iArg ] );
		PCL.Open( ifsm, sArgs.skip_arg );
		if( !iArg )
			for( i = 0; i < PCL.GetFeatures( ); ++i )
				vecstrFeatures.push_back( PCL.GetFeature( i ) );
		for( i = 0; i < PCL.GetExperiments( ); ++i )
			vecstrExps.push_back( PCL.GetExperiment( i ) );
		for( i = 0; i < PCL.GetGenes( ); ++i )
			setstrGenes.insert( PCL.GetGene( i ) );
		ifsm.close( ); }
	vecstrGenes.resize( setstrGenes.size( ) );
	copy( setstrGenes.begin( ), setstrGenes.end( ), vecstrGenes.begin( ) );

	PCLNew.Open( vecstrGenes, vecstrExps, vecstrFeatures );
	iExp = 0;
	for( iArg = 0; iArg < sArgs.inputs_num; ++iArg ) {
		cerr << "Processing " << sArgs.inputs[ iArg ] << "..." << endl;
		ifsm.clear( );
		ifsm.open( sArgs.inputs[ iArg ] );
		PCL.Open( ifsm, sArgs.skip_arg );
		for( i = 0; i < PCLNew.GetGenes( ); ++i )
			if( ( iGene = PCL.GetGene( vecstrGenes[ i ] ) ) != -1 ) {
				if( !iArg )
					for( j = 1; j < PCLNew.GetFeatures( ); ++j )
						PCLNew.SetFeature( i, j, PCL.GetFeature( iGene, j ) );
				for( j = 0; j < PCL.GetExperiments( ); ++j )
					PCLNew.Set( i, iExp + j, PCL.Get( iGene, j ) ); }
		iExp += PCL.GetExperiments( );
		ifsm.close( ); }

	if( sArgs.output_arg ) {
		ofsm.open( sArgs.output_arg );
		PCLNew.Save( ofsm );
		ofsm.close( ); }
	else {
		PCLNew.Save( cout );
		cout.flush( ); }

	return 0; }

int MainDATs( const gengetopt_args_info& sArgs ) {
	CDataset					Dataset;
	CDat						DatOut, DatCur;
	CHalfMatrix<unsigned char>	HMatCounts;
	size_t						i, j, k, iOne, iTwo;
	vector<size_t>				veciGenes;
	float						d;
	vector<string>				vecstrFiles;

	if( !sArgs.inputs_num )
		return 1;

	vecstrFiles.resize( sArgs.inputs_num );
	copy( sArgs.inputs, sArgs.inputs + sArgs.inputs_num, vecstrFiles.begin( ) );
	if( !Dataset.OpenGenes( vecstrFiles ) ) {
		cerr << "Couldn't open: " << vecstrFiles[ 0 ];
		for( i = 1; i < vecstrFiles.size( ); ++i )
			cerr << ", " << vecstrFiles[ i ];
		cerr << endl;
		return 1; }

	DatOut.Open( Dataset.GetGeneNames( ), false, sArgs.memmap_flag ? sArgs.output_arg : NULL );
	if( !strcmp( c_szMax, sArgs.method_arg ) )
		d = -FLT_MAX;
	else if( !strcmp( c_szMin, sArgs.method_arg ) )
		d = FLT_MAX;
	else
		d = 0;
	for( i = 0; i < DatOut.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < DatOut.GetGenes( ); ++j )
			DatOut.Set( i, j, d );
	if( !d ) {
		HMatCounts.Initialize( DatOut.GetGenes( ) );
		HMatCounts.Clear( ); }
	veciGenes.resize( DatOut.GetGenes( ) );
	for( i = 0; i < sArgs.inputs_num; ++i ) {
		if( !DatCur.Open( sArgs.inputs[ i ], !!sArgs.memmap_flag && !sArgs.normalize_flag ) ) {
			cerr << "Couldn't open: " << sArgs.inputs[ i ] << endl;
			return 1; }
		cerr << "Opened: " << sArgs.inputs[ i ] << endl;
		if( sArgs.normalize_flag )
			DatCur.Normalize( false );
		for( j = 0; j < veciGenes.size( ); ++j )
			veciGenes[ j ] = DatCur.GetGene( DatOut.GetGene( j ) );
		for( j = 0; j < veciGenes.size( ); ++j ) {
			if( ( iOne = veciGenes[ j ] ) == -1 )
				continue;
			for( k = ( j + 1 ); k < veciGenes.size( ); ++k ) {
				if( ( ( iTwo = veciGenes[ k ] ) == -1 ) || CMeta::IsNaN( d = DatCur.Get( iOne, iTwo ) ) )
					continue;
				if( !strcmp( c_szMean, sArgs.method_arg ) ) {
					DatOut.Get( j, k ) += d;
					HMatCounts.Get( j, k )++; }
				else if( !strcmp( c_szGMean, sArgs.method_arg ) ) {
					DatOut.Get( j, k ) *= d;
					HMatCounts.Get( j, k )++; }
				else if( !strcmp( c_szHMean, sArgs.method_arg ) ) {
					DatOut.Get( j, k ) += 1 / d;
					HMatCounts.Get( j, k )++; }
				else if( !strcmp( c_szMax, sArgs.method_arg ) ) {
					if( d > DatOut.Get( j, k ) )
						DatOut.Set( j, k, d ); }
				else if( !strcmp( c_szMin, sArgs.method_arg ) ) {
					if( d < DatOut.Get( j, k ) )
						DatOut.Set( j, k, d ); } } } }
	for( i = 0; i < DatOut.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < DatOut.GetGenes( ); ++j )
			if( !strcmp( c_szMean, sArgs.method_arg ) )
				DatOut.Set( i, j, ( k = HMatCounts.Get( i, j ) ) ? ( DatOut.Get( i, j ) / k ) :
					CMeta::GetNaN( ) );
			else if( !strcmp( c_szGMean, sArgs.method_arg ) )
				DatOut.Set( i, j, ( k = HMatCounts.Get( i, j ) ) ?
					(float)pow( (double)DatOut.Get( i, j ), 1.0 / k ) : CMeta::GetNaN( ) );
			else if( !strcmp( c_szHMean, sArgs.method_arg ) )
				DatOut.Set( i, j, ( k = HMatCounts.Get( i, j ) ) ? ( k / DatOut.Get( i, j ) ) :
					CMeta::GetNaN( ) );
			else if( !strcmp( c_szMax, sArgs.method_arg ) ) {
				if( DatOut.Get( i, j ) == -FLT_MAX )
					DatOut.Set( i, j, CMeta::GetNaN( ) ); }
			else if( !strcmp( c_szMin, sArgs.method_arg ) ) {
				if( DatOut.Get( i, j ) == FLT_MAX )
					DatOut.Set( i, j, CMeta::GetNaN( ) ); }

	if( !sArgs.memmap_flag )
		DatOut.Save( sArgs.output_arg );

	return 0; }

static int MainDABs( const gengetopt_args_info& sArgs ) {
	CDatasetCompact	Dataset;
	size_t			i;
	vector<string>	vecstrFiles;
	ofstream		ofsm;

	if( !sArgs.inputs_num )
		return 1;

	vecstrFiles.resize( sArgs.inputs_num );
	copy( sArgs.inputs, sArgs.inputs + sArgs.inputs_num, vecstrFiles.begin( ) );
	if( !Dataset.Open( vecstrFiles ) ) {
		cerr << "Couldn't open: " << vecstrFiles[ 0 ];
		for( i = 1; i < vecstrFiles.size( ); ++i )
			cerr << ", " << vecstrFiles[ i ];
		cerr << endl;
		return 1; }

	if( sArgs.output_arg ) {
		ofsm.open( sArgs.output_arg );
		Dataset.Save( ofsm, true );
		ofsm.close( ); }
	else {
		Dataset.Save( cout, false );
		cout.flush( ); }

	return 0; }

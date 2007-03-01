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

static int MainDATs( const gengetopt_args_info& sArgs ) {
	CDataset		Dataset;
	CDataSubset		DataSubset;
	IDataset*		pData;
	CDat			Dat;
	size_t			i, j, k, iN, iSubset, iSubsets;
	float			d, dMax;
	vector<string>	vecstrFiles;

	if( !sArgs.inputs_num )
		return 1;

	vecstrFiles.resize( sArgs.inputs_num );
	copy( sArgs.inputs, sArgs.inputs + sArgs.inputs_num, vecstrFiles.begin( ) );
	if( sArgs.subset_arg ) {
		DataSubset.Initialize( vecstrFiles, sArgs.subset_arg );
		iSubsets = ( DataSubset.GetGenes( ) / sArgs.subset_arg ) +
			( ( DataSubset.GetGenes( ) % sArgs.subset_arg ) ? 1 : 0 );
		pData = &DataSubset; }
	else {
		if( !Dataset.Open( vecstrFiles ) ) {
			cerr << "Couldn't open: " << vecstrFiles[ 0 ];
			for( i = 1; i < vecstrFiles.size( ); ++i )
				cerr << ", " << vecstrFiles[ i ];
			cerr << endl;
			return 1; }
		iSubsets = 1;
		pData = &Dataset; }

	Dat.Open( pData->GetGeneNames( ) );
	for( iSubset = 0; iSubset < iSubsets; ++iSubset ) {
		if( sArgs.subset_arg ) {
			cerr << "Processing subset " << iSubset << '/' << iSubsets << endl;
			DataSubset.Open( iSubset * sArgs.subset_arg ); }
		for( i = 0; i < Dat.GetGenes( ); ++i ) {
			for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j ) {
				if( !pData->IsExample( i, j ) )
					continue;
				dMax = CMeta::GetNaN( );
				if( !strcmp( c_szMean, sArgs.method_arg ) ) {
					for( iN = k = 0; k < pData->GetExperiments( ); ++k )
						if( !CMeta::IsNaN( d = pData->GetContinuous( i, j, k ) ) ) {
							dMax = ( CMeta::IsNaN( dMax ) ? 0 : dMax ) + d;
							iN++; }
					dMax /= iN; }
				else if( !strcmp( c_szGMean, sArgs.method_arg ) ) {
					for( iN = k = 0; k < pData->GetExperiments( ); ++k )
						if( !CMeta::IsNaN( d = pData->GetContinuous( i, j, k ) ) ) {
							dMax = ( CMeta::IsNaN( dMax ) ? 1 : dMax ) * d;
							iN++; }
					dMax = pow( dMax, 1.0f / iN ); }
				else if( !strcmp( c_szHMean, sArgs.method_arg ) ) {
					for( iN = k = 0; k < pData->GetExperiments( ); ++k )
						if( !CMeta::IsNaN( d = pData->GetContinuous( i, j, k ) ) ) {
							dMax = ( CMeta::IsNaN( dMax ) ? 0 : dMax ) + ( 1 / d );
							iN++; }
					dMax = iN / dMax; }
				else if( !strcmp( c_szMax, sArgs.method_arg ) ) {
					for( k = 0; k < pData->GetExperiments( ); ++k )
						if( !CMeta::IsNaN( d = pData->GetContinuous( i, j, k ) ) &&
							( CMeta::IsNaN( dMax ) || ( d > dMax ) ) )
							dMax = d; }
				else if( !strcmp( c_szMin, sArgs.method_arg ) ) {
					for( k = 0; k < pData->GetExperiments( ); ++k )
						if( !CMeta::IsNaN( d = pData->GetContinuous( i, j, k ) ) &&
							( CMeta::IsNaN( dMax ) || ( d < dMax ) ) )
							dMax = d; }
				else if( !strcmp( c_szVote, sArgs.method_arg ) ) {
					size_t	iPos, iNeg;

					for( iPos = iNeg = k = 0; k < pData->GetExperiments( ); ++k )
						if( !CMeta::IsNaN( d = pData->GetContinuous( i, j, k ) ) ) {
							if( d > 0 )
								iPos++;
							else
								iNeg++; }
					dMax = ( iPos == iNeg ) ? CMeta::GetNaN( ) : ( ( iPos > iNeg ) ? 1 : 0 ); }
				Dat.Set( i, j, dMax ); } } }

	Dat.Save( sArgs.output_arg );

	return 0; }

int MainDATs2( const gengetopt_args_info& sArgs ) {
	CDatasetCompact	Dataset;
	CDat			Dat;
	size_t			i, j, k;
	vector<string>	vecstrFiles;

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

	Dat.Open( Dataset.GetGeneNames( ) );
	for( i = 0; i < Dat.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
			for( k = 0; k < Dataset.GetExperiments( ); ++k )
				if( Dataset.GetDiscrete( i, j, k ) != -1 ) {
					Dat.Set( i, j, 1 );
					break; }

	Dat.Save( sArgs.output_arg );

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

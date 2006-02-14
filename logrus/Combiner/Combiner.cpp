#include "stdafx.h"

static int MainDABs( const gengetopt_args_info& );
static int MainDATs( const gengetopt_args_info& );
static int MainPCLs( const gengetopt_args_info& );

static const TPFnCombiner	c_apfnCombiners[]	= { MainPCLs, MainDATs, MainDABs, NULL };
static const char*			c_aszCombiners[]	= { "pcl", "dat", "dab", NULL };

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
	vector<string>				vecstrGenes, vecstrExps;
	set<string>					setstrGenes;
	set<string>::const_iterator	iterGenes;
	ifstream					ifsm;
	ofstream					ofsm;

	for( iArg = 0; iArg < sArgs.inputs_num; ++iArg ) {
		ifsm.clear( );
		ifsm.open( sArgs.inputs[ iArg ] );
		PCL.Open( ifsm, 2 );
		for( i = 0; i < PCL.GetExperiments( ); ++i )
			vecstrExps.push_back( PCL.GetExperiment( i ) );
		for( i = 0; i < PCL.GetGenes( ); ++i )
			setstrGenes.insert( PCL.GetGene( i ) );
		ifsm.close( ); }
	vecstrGenes.resize( setstrGenes.size( ) );
	i = 0;
	for( iterGenes = setstrGenes.begin( ); iterGenes != setstrGenes.end( ); ++iterGenes )
		vecstrGenes[ i++ ] = *iterGenes;

	PCLNew.Open( vecstrGenes, vecstrExps );
	iExp = 0;
	for( iArg = 0; iArg < sArgs.inputs_num; ++iArg ) {
		ifsm.clear( );
		ifsm.open( sArgs.inputs[ iArg ] );
		PCL.Open( ifsm, 2 );
		for( i = 0; i < PCL.GetGenes( ); ++i )
			if( ( iGene = PCL.GetGene( vecstrGenes[ i ] ) ) != -1 )
				for( j = 0; j < PCL.GetExperiments( ); ++j )
					PCLNew.Set( i, iExp + j, PCL.Get( iGene, j ) );
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
	CDataSubset		Dataset;
	CDat			Dat;
	size_t			i, j, k, iSubset, iSubsets, iN;
	float			d, dMax;
	vector<string>	vecstrFiles;
	ofstream		ofsm;

	if( !sArgs.inputs_num )
		return 1;

	vecstrFiles.resize( sArgs.inputs_num );
	copy( sArgs.inputs, sArgs.inputs + sArgs.inputs_num, vecstrFiles.begin( ) );
	if( !Dataset.Initialize( vecstrFiles, sArgs.subset_arg ) ) {
		cerr << "Couldn't open: " << vecstrFiles[ 0 ];
		for( i = 1; i < vecstrFiles.size( ); ++i )
			cerr << ", " << vecstrFiles[ i ];
		cerr << endl;
		return 1; }

	Dat.Open( Dataset.GetGeneNames( ) );
	iSubsets = ( Dataset.GetGenes( ) / sArgs.subset_arg ) + ( ( Dataset.GetGenes( ) %
		sArgs.subset_arg ) ? 1 : 0 );
	for( iSubset = 0; iSubset < iSubsets; ++iSubset ) {
		cerr << "Opening subset " << (unsigned int)iSubset << "/" <<
			(unsigned int)iSubsets << endl;
		if( !Dataset.Open( iSubset * sArgs.subset_arg ) ) {
			cerr << "Couldn't open subset " << (unsigned int)iSubset << endl;
			return 1; }
		for( i = 0; i < Dat.GetGenes( ); ++i ) {
			for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j ) {
				if( !Dataset.IsExample( i, j ) )
					continue;
				dMax = CMeta::GetNaN( );
				for( iN = k = 0; k < Dataset.GetExperiments( ); ++k )
					if( !CMeta::IsNaN( d = Dataset.GetContinuous( i, j, k ) ) ) {
						dMax += d;
						iN++; }
				dMax /= iN;
/* Max
					if( !CMeta::IsNaN( d = Dataset.GetContinuous( i, j, k ) ) &&
						( CMeta::IsNaN( dMax ) || ( d > dMax ) ) )
						dMax = d;
//*/
				Dat.Set( i, j, dMax ); } } }

	if( sArgs.output_arg ) {
		ofsm.open( sArgs.output_arg, ios_base::binary );
		Dat.Save( ofsm, true );
		ofsm.close( ); }
	else {
		Dat.Save( cout, false );
		cout.flush( ); }

	return 0; }

int MainDATs2( const gengetopt_args_info& sArgs ) {
	CDatasetCompact	Dataset;
	CDat			Dat;
	size_t			i, j, k;
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

	Dat.Open( Dataset.GetGeneNames( ) );
	for( i = 0; i < Dat.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
			for( k = 0; k < Dataset.GetExperiments( ); ++k )
				if( Dataset.GetDiscrete( i, j, k ) != -1 ) {
					Dat.Set( i, j, 1 );
					break; }

	if( sArgs.output_arg ) {
		ofsm.open( sArgs.output_arg );
		Dat.Save( ofsm, true );
		ofsm.close( ); }
	else {
		Dat.Save( cout, false );
		cout.flush( ); }

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

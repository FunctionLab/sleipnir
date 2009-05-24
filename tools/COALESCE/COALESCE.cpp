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

enum EFile {
	EFilePCL,
	EFileWIG,
	EFileError
};

const char*					c_aszMatchTypes[]	= {"pvalue", "rmse", "js", NULL};
const SMotifMatch::EType	c_aeMatchTypes[]	=
	{SMotifMatch::ETypePValue, SMotifMatch::ETypeRMSE, SMotifMatch::ETypeJensenShannon, SMotifMatch::ETypePValue};

int main_postprocess( const gengetopt_args_info&, CCoalesceMotifLibrary& );
bool recluster( const gengetopt_args_info&, size_t, CCoalesceMotifLibrary&, const CHierarchy&,
	const vector<CCoalesceCluster>&, const vector<string>&, const CPCL&, size_t& );
bool trim( const gengetopt_args_info&, const CPCL&, vector<CCoalesceCluster>& );

EFile open_pclwig( const char* szFile, size_t iSkip, CFASTA& FASTA, CPCL& PCL ) {

	if( FASTA.Open( szFile ) )
		return EFilePCL;
	if( PCL.Open( szFile, iSkip ) )
		return EFileWIG;

	cerr << "Could not open: " << szFile << endl;
	return EFileError; }

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info			sArgs;
	CFASTA						FASTA;
	CPCL						PCL;
	CCoalesce					Coalesce;
	vector<CCoalesceCluster>	vecClusters;
	size_t						i, j;
	set<string>					setstrTypes;
	vector<CFASTA*>				vecpFASTAs;

#ifdef WIN32
	pthread_win32_process_attach_np( );
#endif // WIN32

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg, sArgs.random_arg );

	CCoalesceMotifLibrary		Motifs( sArgs.k_arg );
	Motifs.SetPenaltyGap( (float)sArgs.penalty_gap_arg );
	Motifs.SetPenaltyMismatch( (float)sArgs.penalty_mismatch_arg );

	if( sArgs.postprocess_arg )
		return main_postprocess( sArgs, Motifs );

	if( sArgs.sequences_arg ) {
		vector<string>	vecstrTypes;

		CMeta::Tokenize( sArgs.sequences_arg, vecstrTypes, "," );
		for( i = 0; i < vecstrTypes.size( ); ++i )
			setstrTypes.insert( vecstrTypes[ i ] ); }

	if( !PCL.Open( sArgs.input_arg, sArgs.skip_arg ) ) {
		cerr << "Could not open: " << ( sArgs.input_arg ? sArgs.input_arg : "stdin" ) << endl;
		return 1; }
	if( sArgs.fasta_arg && !FASTA.Open( sArgs.fasta_arg, setstrTypes ) ) {
		cerr << "Could not open: " << sArgs.fasta_arg << endl;
		return 1; }

	if( sArgs.datasets_arg ) {
		static const size_t	c_iBuffer	= 131072;
		ifstream			ifsm;
		char				acBuffer[ c_iBuffer ];

		ifsm.open( sArgs.datasets_arg );
		if( !ifsm.is_open( ) ) {
			cerr << "Could not open: " << sArgs.datasets_arg << endl;
			return 1; }
		while( !ifsm.eof( ) ) {
			vector<string>	vecstrLine;
			set<size_t>		setiConditions;

			ifsm.getline( acBuffer, c_iBuffer - 1 );
			acBuffer[ c_iBuffer - 1 ] = 0;
			CMeta::Tokenize( acBuffer, vecstrLine );
			for( i = 0; i < vecstrLine.size( ); ++i )
				for( j = 0; j < PCL.GetExperiments( ); ++j )
					if( vecstrLine[ i ] == PCL.GetExperiment( j ) ) {
						setiConditions.insert( j );
						break; }
			Coalesce.AddDataset( setiConditions ); } }

	Coalesce.SetMotifs( Motifs );
	Coalesce.SetProbabilityGene( (float)sArgs.prob_gene_arg );
	Coalesce.SetPValueCondition( (float)sArgs.pvalue_cond_arg );
	Coalesce.SetZScoreCondition( (float)sArgs.zscore_cond_arg );
	Coalesce.SetPValueMotif( (float)sArgs.pvalue_motif_arg );
	Coalesce.SetZScoreMotif( (float)sArgs.zscore_motif_arg );
	Coalesce.SetPValueCorrelation( (float)sArgs.pvalue_correl_arg );
	Coalesce.SetNumberCorrelation( sArgs.number_correl_arg );
	Coalesce.SetPValueMerge( (float)sArgs.pvalue_merge_arg );
	Coalesce.SetCutoffMerge( (float)sArgs.cutoff_merge_arg );
	Coalesce.SetBasesPerMatch( sArgs.bases_arg );
	Coalesce.SetSizeMinimum( sArgs.size_minimum_arg );
	Coalesce.SetSizeMerge( sArgs.size_merge_arg );
	Coalesce.SetSizeMaximum( sArgs.size_maximum_arg );
	Coalesce.SetThreads( sArgs.threads_arg );
	if( sArgs.output_arg )
		Coalesce.SetDirectoryIntermediate( sArgs.output_arg );
	if( sArgs.cache_arg )
		Coalesce.SetSequenceCache( sArgs.cache_arg );

	vecpFASTAs.resize( sArgs.inputs_num );
	for( i = 0; i < vecpFASTAs.size( ); ++i ) {
		vecpFASTAs[ i ] = new CFASTA( );
		if( !vecpFASTAs[ i ]->Open( sArgs.inputs[ i ] ) ) {
			cerr << "Could not open: " << sArgs.inputs[ i ] << endl;
			return 1; }
		Coalesce.AddWiggle( *vecpFASTAs[ i ] ); }

	if( sArgs.progressive_flag )
		Coalesce.AddOutputIntermediate( cout );
	if( !Coalesce.Cluster( PCL, FASTA, vecClusters ) ) {
		cerr << "Clustering failed" << endl;
		return 1; }

	for( i = 0; i < vecpFASTAs.size( ); ++i )
		delete vecpFASTAs[ i ];
	if( !sArgs.progressive_flag )
		for( i = 0; i < vecClusters.size( ); ++i )
			vecClusters[ i ].Save( cout, i, PCL, &Motifs );

#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; }

int main_postprocess( const gengetopt_args_info& sArgs, CCoalesceMotifLibrary& Motifs ) {
	string						strDir, strFile, strBase;
	vector<CCoalesceCluster>	vecClustersFrom, vecClustersTo;
	bool						fFailed;
	CDistanceMatrix				MatSim;
	size_t						i, j, iDatasets;
	CHierarchy*					pHier;
	vector<string>				vecstrClusters;
	CPCL						PCL;
	ifstream					ifsm;

	if( sArgs.input_arg ) {
		if( !PCL.Open( sArgs.input_arg, sArgs.skip_arg ) ) {
			cerr << "Could not open: " << sArgs.input_arg << endl;
			return 1; } }
	else if( !PCL.Open( cin, sArgs.skip_arg ) ) {
		cerr << "Could not open: stdin" << endl;
		return 1; }
	if( sArgs.known_motifs_arg ) {
		ifsm.open( sArgs.known_motifs_arg );
		if( !Motifs.OpenKnown( ifsm ) ) {
			cerr << "Could not open: " << sArgs.known_motifs_arg << endl;
			return 1; }
		ifsm.close( ); }
	fFailed = false;
	strDir = sArgs.postprocess_arg;
	FOR_EACH_DIRECTORY_FILE(strDir, strFile)
		if( !CMeta::IsExtension( strFile, CPCL::GetExtension( ) ) )
			continue;
		strBase = CMeta::Deextension( strFile );
		strFile = strDir + '/' + strFile;

		if( !fFailed )
			vecClustersFrom.push_back( CCoalesceCluster( ) );
		if( fFailed = ( ( iDatasets = vecClustersFrom.back( ).Open( strFile, sArgs.skip_arg, PCL,
			&Motifs ) ) == -1 ) ) {
			cerr << "Could not open: " << strFile << endl;
			continue; }
		if( !( vecstrClusters.size( ) % 50 ) )
			cerr << "Opened cluster " << vecstrClusters.size( ) << "..." << endl;
		vecstrClusters.push_back( strBase ); }
	if( fFailed )
		vecClustersFrom.pop_back( );
	if( vecClustersFrom.empty( ) ) {
		char	szBuffer[ 16 ];

		fFailed = false;
		ifsm.clear( );
		ifsm.open( sArgs.postprocess_arg );
		while( !ifsm.eof( ) ) {
			if( !fFailed )
				vecClustersFrom.push_back( CCoalesceCluster( ) );
			if( fFailed = ( ( iDatasets = vecClustersFrom.back( ).Open( ifsm, PCL, &Motifs ) ) == -1 ) )
				continue;
			if( !( vecstrClusters.size( ) % 50 ) )
				cerr << "Opened cluster " << vecstrClusters.size( ) << "..." << endl;
			sprintf_s( szBuffer, "c%06d", vecstrClusters.size( ) );
			vecstrClusters.push_back( szBuffer ); }
		if( fFailed )
			vecClustersFrom.pop_back( ); }

	cerr << "Trimming clusters..." << endl;
	if( !trim( sArgs, PCL, vecClustersFrom ) )
		return 1;

	cerr << "Calculating cluster similarities..." << endl;
	MatSim.Initialize( vecClustersFrom.size( ) );
	for( i = 0; i < MatSim.GetSize( ); ++i )
		for( j = ( i + 1 ); j < MatSim.GetSize( ); ++j ) {
			if( sArgs.verbosity_arg >= 7 )
				cerr << "Comparing clusters:	" << i << '\t' << j << endl;
			MatSim.Set( i, j, vecClustersFrom[ i ].GetSimilarity( vecClustersFrom[ j ], PCL.GetGenes( ),
				iDatasets ) ); }

	i = 0;
	return ( ( ( pHier = CClustHierarchical::Cluster( MatSim ) ) &&
		recluster( sArgs, MatSim.GetSize( ) * ( MatSim.GetSize( ) - 1 ) / 2, Motifs, *pHier, vecClustersFrom,
		vecstrClusters, PCL, i ) ) ? 0 : 1 ); }

bool recluster( const gengetopt_args_info& sArgs, size_t iPairs, CCoalesceMotifLibrary& Motifs,
	const CHierarchy& Hier, const vector<CCoalesceCluster>& vecClustersFrom,
	const vector<string>& vecstrClustersFrom, const CPCL& PCL, size_t& iID ) {
	SMotifMatch::EType	eMatchType;
	size_t				i;

	for( i = 0; c_aszMatchTypes[ i ]; ++i )
		if( !strcmp( c_aszMatchTypes[ i ], sArgs.known_type_arg ) )
			break;
	eMatchType = c_aeMatchTypes[ i ];

	if( Hier.IsGene( ) || ( Hier.GetSimilarity( ) >= sArgs.cutoff_postprocess_arg ) ) {
		CCoalesceCluster	Cluster;

		cerr << "Creating output cluster " << iID << endl;
		if( !Cluster.Open( Hier, vecClustersFrom, vecstrClustersFrom,
			(float)sArgs.fraction_postprocess_arg, (float)sArgs.cutoff_merge_arg, sArgs.max_motifs_arg,
			&Motifs ) )
			return false;
		if( Cluster.GetGenes( ).size( ) < (size_t)sArgs.size_minimum_arg ) {
			cerr << "Cluster too small: " << Cluster.GetGenes( ).size( ) << endl;
			return true; }

		if( !Cluster.LabelMotifs( Motifs, eMatchType, (float)sArgs.penalty_gap_arg,
			(float)sArgs.penalty_mismatch_arg, (float)sArgs.known_cutoff_arg ) )
			return false;
		if( sArgs.output_arg )
			Cluster.Save( sArgs.output_arg, iID, PCL, &Motifs );
		Cluster.Save( cout, iID, PCL, &Motifs, (float)sArgs.min_info_arg, (float)sArgs.penalty_gap_arg,
			(float)sArgs.penalty_mismatch_arg, !!sArgs.remove_rcs_flag );
		iID++;
		return true; }

	return ( recluster( sArgs, iPairs, Motifs, Hier.Get( false ), vecClustersFrom, vecstrClustersFrom, PCL,
		iID ) && recluster( sArgs, iPairs, Motifs, Hier.Get( true ), vecClustersFrom, vecstrClustersFrom,
		PCL, iID ) ); }

bool trim( const gengetopt_args_info& sArgs, const CPCL& PCL, vector<CCoalesceCluster>& vecClusters ) {
	vector<size_t>				veciCounts, veciGenes, veciRemove;
	CHalfMatrix<size_t>			MatCounts;
	CDistanceMatrix				MatScores;
	size_t						i, j, iOne, iCluster;
	float						dAve, dStd, d, dCutoff;
	vector<float>				vecdScores;

	veciCounts.resize( PCL.GetGenes( ) );
	MatCounts.Initialize( PCL.GetGenes( ) );
	MatCounts.Clear( );
	for( iCluster = 0; iCluster < vecClusters.size( ); ++iCluster ) {
		const CCoalesceCluster&	Cluster	= vecClusters[ iCluster ];

		veciGenes.resize( Cluster.GetGenes( ).size( ) );
		copy( Cluster.GetGenes( ).begin( ), Cluster.GetGenes( ).end( ), veciGenes.begin( ) );
		for( i = 0; i < veciGenes.size( ); ++i ) {
			veciCounts[ iOne = veciGenes[ i ] ]++;
			for( j = ( i + 1 ); j < veciGenes.size( ); ++j )
				MatCounts.Get( iOne, veciGenes[ j ] )++; } }
	MatScores.Initialize( MatCounts.GetSize( ) );
	for( i = 0; i < MatScores.GetSize( ); ++i )
		for( j = ( i + 1 ); j < MatScores.GetSize( ); ++j ) {
			iOne = MatCounts.Get( i, j );
			iCluster = veciCounts[ i ] + veciCounts[ j ] - iOne;
			MatScores.Set( i, j, iCluster ? ( (float)iOne / iCluster ) : CMeta::GetNaN( ) ); }

	dAve = dStd = 0;
	for( iOne = i = 0; i < MatScores.GetSize( ); ++i )
		for( j = ( i + 1 ); j < MatScores.GetSize( ); ++j )
			if( !CMeta::IsNaN( d = MatScores.Get( i, j ) ) ) {
				iOne++;
				dAve += d;
				dStd += d * d; }
	dAve /= iOne;
	dStd = sqrt( ( dStd / iOne ) - ( dAve * dAve ) );
	dCutoff = dAve + ( (float)sArgs.cutoff_trim_arg * dStd );
	if( sArgs.verbosity_arg >= 6 )
		cerr << "Cutoff " << dCutoff << " (" << dAve << ',' << dStd << ')' << endl;

	for( iCluster = 0; iCluster < vecClusters.size( ); ++iCluster ) {
		CCoalesceCluster&	Cluster	= vecClusters[ iCluster ];

		if( sArgs.verbosity_arg >= 6 )
			cerr << "Trimming cluster " << iCluster << endl;
		veciGenes.resize( Cluster.GetGenes( ).size( ) );
		copy( Cluster.GetGenes( ).begin( ), Cluster.GetGenes( ).end( ), veciGenes.begin( ) );
		vecdScores.resize( veciGenes.size( ) );
		for( i = 0; i < vecdScores.size( ); ++i ) {
			iOne = veciGenes[ i ];
			for( j = ( i + 1 ); j < vecdScores.size( ); ++j ) {
				d = MatScores.Get( iOne, veciGenes[ j ] );
				vecdScores[ i ] += d;
				vecdScores[ j ] += d; } }
		for( i = 0; i < vecdScores.size( ); ++i )
			vecdScores[ i ] /= veciGenes.size( ) - 1;
		if( sArgs.verbosity_arg >= 7 )
			for( i = 0; i < vecdScores.size( ); ++i )
				cerr << "Gene " << i << " (" << PCL.GetGene( veciGenes[ i ] ) << ") " << vecdScores[ i ] <<
					endl;

		veciRemove.clear( );
		for( i = 0; i < vecdScores.size( ); ++i )
			if( vecdScores[ i ] < dCutoff ) {
				if( sArgs.verbosity_arg >= 6 )
					cerr << "Removing " << PCL.GetGene( veciGenes[ i ] ) << " at " << vecdScores[ i ] << endl;
				veciRemove.push_back( veciGenes[ i ] ); }
		Cluster.RemoveGenes( veciRemove ); }

	return true; }

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
#include "pst.h"

enum EFile {
	EFilePCL,
	EFileWIG,
	EFileError
};

int main_postprocess( const gengetopt_args_info&, CCoalesceMotifLibrary& );
bool recluster( const gengetopt_args_info&, size_t, CCoalesceMotifLibrary&, const CHierarchy&,
	const vector<CCoalesceCluster>&, const vector<string>&, vector<CCoalesceCluster>& );

EFile open_pclwig( const char* szFile, size_t iSkip, CFASTA& FASTA, CPCL& PCL ) {

	if( FASTA.Open( szFile ) )
		return EFilePCL;
	if( PCL.Open( szFile, iSkip ) )
		return EFileWIG;

	cerr << "Could not open: " << szFile << endl;
	return EFileError; }

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info			sArgs;
	CFASTA						FASTA, FASTANucleosomes, FASTAConservation;
	CPCL						PCL;
	CCoalesce					Coalesce;
	vector<CCoalesceCluster>	vecClusters;
	size_t						i, j;
	set<string>					setstrTypes;

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
	Coalesce.SetPValueMotif( (float)sArgs.pvalue_motif_arg );
	Coalesce.SetPValueCorrelation( (float)sArgs.pvalue_correl_arg );
	Coalesce.SetNumberCorrelation( sArgs.number_correl_arg );
	Coalesce.SetPValueMerge( (float)sArgs.pvalue_merge_arg );
	Coalesce.SetCutoffMerge( (float)sArgs.cutoff_merge_arg );
	Coalesce.SetBasesPerMatch( sArgs.bases_arg );
	Coalesce.SetSizeMinimum( sArgs.size_minimum_arg );
	Coalesce.SetSizeMaximum( sArgs.size_maximum_arg );
	Coalesce.SetThreads( sArgs.threads_arg );
	if( sArgs.output_arg )
		Coalesce.SetDirectoryIntermediate( sArgs.output_arg );
	if( sArgs.cache_arg )
		Coalesce.SetSequenceCache( sArgs.cache_arg );
	if( sArgs.nucleosomes_arg ) {
		if( !FASTANucleosomes.Open( sArgs.nucleosomes_arg ) ) {
			cerr << "Could not open: " << sArgs.nucleosomes_arg << endl;
			return 1; }
		Coalesce.AddWiggle( FASTANucleosomes ); }
	if( sArgs.conservation_arg ) {
		if( !FASTAConservation.Open( sArgs.conservation_arg ) ) {
			cerr << "Could not open: " << sArgs.conservation_arg << endl;
			return 1; }
		Coalesce.AddWiggle( FASTAConservation ); }
	if( !Coalesce.Cluster( PCL, FASTA, vecClusters ) ) {
		cerr << "Clustering failed" << endl;
		return 1; }

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

	if( sArgs.input_arg ) {
		if( !PCL.Open( sArgs.input_arg, sArgs.skip_arg ) ) {
			cerr << "Could not open: " << sArgs.input_arg << endl;
			return 1; } }
	else if( !PCL.Open( cin, sArgs.skip_arg ) ) {
		cerr << "Could not open: stdin" << endl;
		return 1; }
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
		vecstrClusters.push_back( strBase ); }
	if( fFailed )
		vecClustersFrom.pop_back( );

	MatSim.Initialize( vecClustersFrom.size( ) );
	for( i = 0; i < MatSim.GetSize( ); ++i )
		for( j = ( i + 1 ); j < MatSim.GetSize( ); ++j ) {
			if( sArgs.verbosity_arg >= 7 )
				cerr << "Comparing clusters:	" << i << '\t' << j << endl;
			MatSim.Set( i, j, vecClustersFrom[ i ].GetSimilarity( vecClustersFrom[ j ], PCL.GetGenes( ),
				iDatasets ) ); }
	if( !( ( pHier = CClustHierarchical::Cluster( MatSim ) ) &&
		recluster( sArgs, MatSim.GetSize( ) * ( MatSim.GetSize( ) - 1 ) / 2, Motifs, *pHier, vecClustersFrom,
		vecstrClusters, vecClustersTo ) ) )
		return 1;

	for( i = 0; i < vecClustersTo.size( ); ++i ) {
		vecClustersTo[ i ].RemoveMotifs( (float)sArgs.min_zscore_arg );
		if( sArgs.output_arg )
			vecClustersTo[ i ].Save( sArgs.output_arg, i, PCL, &Motifs );
		vecClustersTo[ i ].Save( cout, i, PCL, &Motifs, (float)sArgs.min_info_arg,
			!!sArgs.remove_rcs_flag ); }

	return 0; }

bool recluster( const gengetopt_args_info& sArgs, size_t iPairs, CCoalesceMotifLibrary& Motifs,
	const CHierarchy& Hier, const vector<CCoalesceCluster>& vecClustersFrom,
	const vector<string>& vecstrClustersFrom, vector<CCoalesceCluster>& vecClustersTo ) {
	bool	fRet;

	if( Hier.IsGene( ) || ( Hier.GetSimilarity( ) >= sArgs.cutoff_postprocess_arg ) ) {
		cerr << "Creating output cluster " << vecClustersTo.size( ) << endl;
		vecClustersTo.push_back( CCoalesceCluster( ) );
		fRet = vecClustersTo.back( ).Open( Hier, vecClustersFrom, vecstrClustersFrom,
			(float)sArgs.fraction_postprocess_arg, (float)sArgs.cutoff_merge_arg, &Motifs );
		if( fRet && ( vecClustersTo.back( ).GetGenes( ).size( ) < (size_t)sArgs.size_minimum_arg ) ) {
			cerr << "Cluster too small: " << vecClustersTo.back( ).GetGenes( ).size( ) << endl;
			vecClustersTo.pop_back( ); }
		return fRet; }

	return ( recluster( sArgs, iPairs, Motifs, Hier.Get( false ), vecClustersFrom, vecstrClustersFrom,
		vecClustersTo ) && recluster( sArgs, iPairs, Motifs, Hier.Get( true ), vecClustersFrom,
		vecstrClustersFrom, vecClustersTo ) ); }

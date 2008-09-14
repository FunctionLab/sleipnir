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

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info			sArgs;
	CFASTA						FASTA;
	CPCL						PCL;
	CCoalesce					Coalesce;
	vector<CCoalesceCluster>	vecClusters;
	size_t						i, j;
	set<string>					setstrTypes;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg, sArgs.random_arg );

	CCoalesceMotifLibrary		Motifs( sArgs.k_arg );

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
	Coalesce.SetFractionCorrelation( (float)sArgs.frac_correl_arg );
	Coalesce.SetPValueMerge( (float)sArgs.pvalue_merge_arg );
	Coalesce.SetCutoffMerge( (float)sArgs.cutoff_merge_arg );
	Coalesce.SetPenaltyGap( (float)sArgs.penalty_gap_arg );
	Coalesce.SetPenaltyMismatch( (float)sArgs.penalty_mismatch_arg );
	Coalesce.SetBasesPerMatch( sArgs.bases_arg );
	Coalesce.SetSizeMinimum( sArgs.size_minimum_arg );
	Coalesce.SetSizeMaximum( sArgs.size_maximum_arg );
	if( sArgs.intermediate_arg )
		Coalesce.SetDirectoryIntermediate( sArgs.intermediate_arg );
	if( sArgs.cache_arg )
		Coalesce.SetSequenceCache( sArgs.cache_arg );
	if( !Coalesce.Cluster( PCL, FASTA, vecClusters ) ) {
		cerr << "Clustering failed" << endl;
		return 1; }

	for( i = 0; i < vecClusters.size( ); ++i )
		vecClusters[ i ].Save( cout, i, PCL, &Motifs );

	return 0; }

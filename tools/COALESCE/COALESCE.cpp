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

#include "measure.h"
#include "statistics.h"

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info			sArgs;
	CFASTA						FASTA;
	CPCL						PCL;
	CCoalesce					Coalesce;
	vector<CCoalesceCluster>	vecClusters;
	size_t						i;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg, sArgs.random_arg );

	if( !PCL.Open( sArgs.input_arg, sArgs.skip_arg ) ) {
		cerr << "Could not open: " << ( sArgs.input_arg ? sArgs.input_arg : "stdin" ) << endl;
		return 1; }
	if( sArgs.fasta_arg && !FASTA.Open( sArgs.fasta_arg ) ) {
		cerr << "Could not open: " << sArgs.fasta_arg << endl;
		return 1; }

	Coalesce.SetProbabilityGene( (float)sArgs.prob_gene_arg );
	Coalesce.SetPValueCondition( (float)sArgs.pvalue_cond_arg );
	Coalesce.SetPValueMotif( (float)sArgs.pvalue_motif_arg );
	Coalesce.SetPValueCorrelation( (float)sArgs.pvalue_correl_arg );
	Coalesce.SetOutputIntermediate( !!sArgs.intermediate_flag );
	if( !Coalesce.Cluster( PCL, FASTA, vecClusters ) ) {
		cerr << "Clustering failed" << endl;
		return 1; }

	for( i = 0; i < vecClusters.size( ); ++i ) {
		const CCoalesceCluster&		Cluster	= vecClusters[ i ];
		set<size_t>::const_iterator	iter;

		cout << "Cluster\t" << i << endl;
		cout << "Genes";
		for( iter = Cluster.GetGenes( ).begin( ); iter != Cluster.GetGenes( ).end( ); ++iter )
			cout << '\t' << PCL.GetGene( *iter );
		cout << endl << "Conditions";
		for( iter = Cluster.GetConditions( ).begin( ); iter != Cluster.GetConditions( ).end( ); ++iter )
			cout << '\t' << *iter;
		cout << endl;
		if( sArgs.fasta_arg ) {
			cout << "Motifs";
			for( iter = Cluster.GetMotifs( ).begin( ); iter != Cluster.GetMotifs( ).end( ); ++iter )
				cout << '\t' << *iter;
			cout << endl; } }
	
	return 0; }

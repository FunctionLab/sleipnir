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

struct SSorter {
	const vector<float>&	m_vecdScores;

	SSorter( const vector<float>& vecdScores ) : m_vecdScores(vecdScores) { }

	bool operator()( size_t iOne, size_t iTwo ) {

		return ( m_vecdScores[iOne] > m_vecdScores[iTwo] ); }
};

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CPCL				PCL;
	CDat				Dat;
	size_t				i, j, k, iGene;
	vector<bool>		vecfClinical;
	vector<size_t>		veciGenes2PCL, veciPCL2Genes, veciIndices, veciScores;
	vector<float>		vecdScores;
	CGenome				Genome;
	float				d;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg );

	if( !PCL.Open( sArgs.input_arg, sArgs.skip_arg ) ) {
		cerr << "Could not open: " << ( sArgs.input_arg ? sArgs.input_arg : "stdin" ) << endl;
		return 1; }
	if( PCL.GetFeatures( ) < 2 ) {
		cerr << "PCL requires at least one clinical variable feature" << endl;
		return 1; }
	if( sArgs.spearman_flag )
		PCL.RankTransform( );
	if( sArgs.global_arg && !Dat.Open( sArgs.global_arg, !!sArgs.memmap_flag ) ) {
		cerr << "Could not open: " << sArgs.global_arg << endl;
		return 1; }

	vecfClinical.resize( PCL.GetGenes( ) );
	veciPCL2Genes.resize( PCL.GetGenes( ) );
	fill( veciPCL2Genes.begin( ), veciPCL2Genes.end( ), -1 );
	for( iGene = i = 0; i < vecfClinical.size( ); ++i )
		if( !( vecfClinical[i] = !!atoi( PCL.GetFeature( i, 1 ).c_str( ) ) ) )
			veciPCL2Genes[i] = iGene++;
	veciGenes2PCL.resize( iGene );
	vecdScores.resize( veciGenes2PCL.size( ) );
	veciScores.resize( vecdScores.size( ) );
	veciIndices.resize( iGene );
	for( iGene = i = 0; i < PCL.GetGenes( ); ++i )
		if( !vecfClinical[i] )
			veciGenes2PCL[iGene++] = i;
	for( i = 0; i < PCL.GetGenes( ); ++i ) {
		CGenes			GenesFinal( Genome );
		vector<size_t>	veciFinal;

		if( !vecfClinical[i] )
			continue;
		cerr << "Processing: " << PCL.GetGene( i ) << endl;
		for( iGene = j = 0; j < PCL.GetGenes( ); ++j )
			if( !vecfClinical[j] ) {
				vecdScores[iGene] = (float)CMeasurePearson::Pearson( PCL.Get( i ), PCL.GetExperiments( ), PCL.Get( j ), PCL.GetExperiments( ),
					IMeasure::EMapNone, NULL, NULL, &veciScores[iGene] );
				iGene++; }
		for( j = 0; j < veciIndices.size( ); ++j )
			veciIndices[j] = j;
		sort( veciIndices.begin( ), veciIndices.end( ), SSorter( vecdScores ) );

		if( sArgs.global_arg && sArgs.initial_arg ) {
			CDat	DatCopy;
			CGenome	Genome;
			CGenes	GenesInitial( Genome );

			for( j = 0; ( j < (size_t)sArgs.initial_arg ) && ( j < veciIndices.size( ) ); ++j )
				GenesInitial.AddGene( PCL.GetGene( veciGenes2PCL[veciIndices[j]] ) );
			DatCopy.Open( Dat );
			cerr << "Filtering...";
			DatCopy.FilterGenes( GenesInitial, sArgs.hefalmp_flag ? CDat::EFilterHefalmp : CDat::EFilterPixie, sArgs.final_arg, CMeta::GetNaN( ) );
			cerr << "  done!" << endl;
			for( j = 0; j < DatCopy.GetGenes( ); ++j )
				for( k = ( j + 1 ); k < DatCopy.GetGenes( ); ++k )
					if( !CMeta::IsNaN( DatCopy.Get( j, k ) ) ) {
						GenesFinal.AddGene( Dat.GetGene( j ) );
						GenesFinal.AddGene( Dat.GetGene( k ) ); }
			veciFinal.reserve( GenesFinal.GetGenes( ) );
			for( j = 0; j < GenesFinal.GetGenes( ); ++j )
				if( ( k = PCL.GetGene( GenesFinal.GetGene( j ).GetName( ) ) ) != -1 )
					veciFinal.push_back( k );
			iGene = veciFinal.size( ); }
		else {
			veciFinal.reserve( sArgs.final_arg );
			for( j = 0; ( j < (size_t)sArgs.final_arg ) && ( j < veciIndices.size( ) ); ++j )
				veciFinal.push_back( veciGenes2PCL[veciIndices[j]] );
			iGene = veciGenes2PCL.size( ); }

		for( j = 0; j < veciFinal.size( ); ++j ) {
			k = veciPCL2Genes[veciFinal[j]];
			d = (float)( sArgs.spearman_flag ? CStatistics::PValueSpearman : CStatistics::PValuePearson )( vecdScores[k], veciScores[k] );
			cout << PCL.GetGene( i ) << '\t' << PCL.GetGene( veciFinal[j] ) << '\t' << vecdScores[k] << '\t' << veciScores[k] << '\t' <<
				( d * iGene ) << endl; } }

	return 0; }

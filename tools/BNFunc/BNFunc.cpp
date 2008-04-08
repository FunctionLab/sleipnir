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
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, and Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "stdafx.h"
#include "cmdline.h"

const char	c_szBP[]	= "bp";
const char	c_szCC[]	= "cc";

double overlap( const CGenome&, const CSlim&, size_t, size_t );

int main( int iArgs, char** aszArgs ) {
	CGenome								Genome;
	gengetopt_args_info					sArgs;
	ifstream							ifsmAnno, ifsmOnto, ifsmTerm;
	COntologyGO							GO;
	COntologyGO::ENamespace				eName;
	COntologyKEGG						KEGG;
	COntologyMIPS						MIPS;
	const IOntology*					pOnto;
	CDat								Dat;
	CSlim								Slim, SlimNeg;
	ofstream							ofsm;
	size_t								i, j, k;
	map<const CGene*,bool>				mapGenes;
	map<const CGene*,bool>::iterator	iterGene;
	int									iRet;

	iRet = cmdline_parser2( iArgs, aszArgs, &sArgs, 0, 1, 0 );
	if( sArgs.config_arg )
		iRet = cmdline_parser_configfile( sArgs.config_arg, &sArgs, 0, 0, 1 ) && iRet;
	if( iRet ) {
		cmdline_parser_print_help( );
		return 1; }

	CMeta::Startup( sArgs.verbosity_arg, sArgs.random_arg );
	if( sArgs.go_onto_arg ) {
		ifsmOnto.open( sArgs.go_onto_arg );
		if( sArgs.go_anno_arg )
			ifsmAnno.open( sArgs.go_anno_arg );
		if( !strcmp( sArgs.go_name_arg, c_szBP ) )
			eName = COntologyGO::ENamespaceBP;
		else if( !strcmp( sArgs.go_name_arg, c_szCC ) )
			eName = COntologyGO::ENamespaceCC;
		else
			eName = COntologyGO::ENamespaceMF;
		if( !GO.Open( ifsmOnto, ifsmAnno, Genome, eName, !!sArgs.dbids_flag ) ) {
			cerr << "Couldn't open: ";
			if( sArgs.go_anno_arg )
				cerr << sArgs.go_anno_arg << ", ";
			cerr << sArgs.go_onto_arg << endl;
			return 1; }
		ifsmOnto.close( );
		if( sArgs.go_anno_arg )
			ifsmAnno.close( );
		pOnto = &GO; }
	else if( sArgs.mips_onto_arg ) {
		ifsmOnto.open( sArgs.mips_onto_arg );
		if( sArgs.mips_anno_arg )
			ifsmAnno.open( sArgs.mips_anno_arg );
		if( !MIPS.Open( ifsmOnto, ifsmAnno, Genome ) ) {
			cerr << "Couldn't open: " << sArgs.mips_anno_arg << ", " <<
				sArgs.mips_onto_arg << endl;
			return 1; }
		ifsmOnto.close( );
		if( sArgs.mips_anno_arg )
			ifsmAnno.close( );
		pOnto = &MIPS; }
	else if( sArgs.kegg_arg ) {
		ifsmOnto.open( sArgs.kegg_arg );
		if( !KEGG.Open( ifsmOnto, Genome, sArgs.kegg_org_arg ) ) {
			cerr << "Couldn't open: " << sArgs.kegg_arg << endl;
			return 1; }
		ifsmOnto.close( );
		pOnto = &KEGG; }
	else {
		cerr << "No ontology found." << endl;
		return 1; }

	ifsmOnto.clear( );
	ifsmOnto.open( sArgs.input_arg );
	if( !Slim.Open( ifsmOnto, pOnto ) ) {
		cerr << "Couldn't open: " << sArgs.input_arg << endl;
		return 1; }
	ifsmOnto.close( );

	if( sArgs.sql_arg ) {
		ofsm.open( ( (string)sArgs.sql_arg + "/functions.txt" ).c_str( ) );
		for( i = 0; i < Slim.GetSlims( ); ++i )
			ofsm << ( i + 1 ) << '\t' << Slim.GetSlim( i ) << endl;
		ofsm.close( );

		ofsm.open( ( (string)sArgs.sql_arg + "/genes.txt" ).c_str( ) );
		for( i = 0; i < Genome.GetGenes( ); ++i )
			ofsm << ( i + 1 ) << '\t' << Genome.GetGene( i ).GetName( ) << endl;
		ofsm.close( );

		ofsm.open( ( (string)sArgs.sql_arg + "/gene_names.txt" ).c_str( ) );
		for( k = i = 0; i < Genome.GetGenes( ); ++i ) {
			const CGene&	Gene	= Genome.GetGene( i );

			for( j = 0; j < Gene.GetSynonyms( ); ++j )
				ofsm << ( k++ + 1 ) << '\t' << ( i + 1 ) << '\t' << Gene.GetSynonym( j ) << endl; }
		ofsm.close( );

		ofsm.open( ( (string)sArgs.sql_arg + "/functions_genes.txt" ).c_str( ) );
		for( i = 0; i < Slim.GetSlims( ); ++i )
			for( j = 0; j < Slim.GetGenes( i ); ++j )
				ofsm << ( i + 1 ) << '\t' << ( Genome.GetGene( Slim.GetGene( i, j ).GetName( ) ) + 1 ) << endl;
		ofsm.close( );

		return 0; }

	if( sArgs.nsets_flag ) {
		for( i = 0; i < Slim.GetSlims( ); ++i ) {
			map<const CGene*,bool>				mapGenes;
			map<const CGene*,bool>::iterator	iterGene;

			for( j = 0; j < Slim.GetSlims( ); ++j ) {
				if( i == j )
					continue;
				if( overlap( Genome, Slim, i, j ) < sArgs.nsetlap_arg )
					for( k = 0; k < Slim.GetGenes( j ); ++k )
						mapGenes[ &Slim.GetGene( j, k ) ] = false;
				else
					for( k = 0; k < Slim.GetGenes( j ); ++k ) {
						const CGene*	pGene	= &Slim.GetGene( j, k );

						if( mapGenes.find( pGene ) == mapGenes.end( ) )
							mapGenes[ pGene ] = true; } }

			ofsm.open( ( (string)sArgs.directory_arg + '/' + CMeta::Filename( Slim.GetSlim( i ) ) ).c_str( ) );
			for( iterGene = mapGenes.begin( ); iterGene != mapGenes.end( ); ++iterGene )
				if( iterGene->second )
					ofsm << iterGene->first->GetName( ) << endl;
			ofsm.close( ); }
		return 0; }

	if( sArgs.output_arg ) {
		if( sArgs.negatives_arg ) {
			ifsmOnto.clear( );
			ifsmOnto.open( sArgs.negatives_arg );
			if( !( SlimNeg.Open( ifsmOnto, pOnto ) && Dat.Open( Slim, SlimNeg ) ) ) {
				cerr << "Couldn't open: " << sArgs.negatives_arg << endl;
				return 1; }
			ifsmOnto.close( ); }
		else if( !Dat.Open( Slim ) ) {
			cerr << "Couldn't open: " << sArgs.input_arg << endl;
			return 1; }

		Dat.Save( sArgs.output_arg ); }

	if( sArgs.test_arg ) {
		for( i = 0; i < Slim.GetSlims( ); ++i )
			for( j = 0; j < Slim.GetGenes( i ); ++j )
				mapGenes[ &Slim.GetGene( i, j ) ] = false;
		for( i = 0; i < Slim.GetSlims( ); ++i )
			mapGenes[ &Slim.GetGene( i, rand( ) % Slim.GetGenes( i ) ) ] = true;
		for( iterGene = mapGenes.begin( ); iterGene != mapGenes.end( ); ++iterGene )
			if( ( (float)rand( ) / RAND_MAX ) < sArgs.test_arg )
				iterGene->second = true; }

	if( sArgs.directory_arg )
		for( i = 0; i < Slim.GetSlims( ); ++i ) {
			ofsm.clear( );
			ofsm.open( ( (string)sArgs.directory_arg + '\\' +
				CMeta::Filename( Slim.GetSlim( i ) ) ).c_str( ) );
			for( j = 0; j < Slim.GetGenes( i ); ++j ) {
				const CGene&	Gene	= Slim.GetGene( i, j );
				const string&	strName	= ( sArgs.synonyms_flag && Gene.GetSynonyms( ) ) ?
					Gene.GetSynonym( 0 ) : Gene.GetName( );
			
				ofsm << strName;
				if( sArgs.allids_flag )
					for( k = 0; k < Gene.GetSynonyms( ); ++k )
						ofsm << '\t' << ( ( !k && sArgs.synonyms_flag ) ? Gene.GetName( ) :
							Gene.GetSynonym( k ) );
				ofsm << endl; }
			ofsm.close( ); }
	if( sArgs.test_arg )
		for( iterGene = mapGenes.begin( ); iterGene != mapGenes.end( ); ++iterGene )
			if( iterGene->second ) {
				const string&	strName	= ( sArgs.synonyms_flag && iterGene->first->GetSynonyms( ) ) ?
					iterGene->first->GetSynonym( 0 ) : iterGene->first->GetName( );

				cout << strName << endl; }

	return 0; }

double overlap( const CGenome& Genome, const CSlim& Slim, size_t iOne, size_t iTwo ) {
	size_t				i, iOverlap;
	set<const CGene*>	setOne;

	for( i = 0; i < Slim.GetGenes( iOne ); ++i )
		setOne.insert( &Slim.GetGene( iOne, i ) );
	for( iOverlap = i = 0; i < Slim.GetGenes( iTwo ); ++i )
		if( setOne.find( &Slim.GetGene( iTwo, i ) ) != setOne.end( ) )
			iOverlap++;

	return CStatistics::HypergeometricCDF( iOverlap, Slim.GetGenes( iOne ), Slim.GetGenes( iTwo ),
		Genome.GetGenes( ) ); }

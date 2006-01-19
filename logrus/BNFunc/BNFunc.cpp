#include "stdafx.h"
#include "cmdline.h"

const char	c_szBP[]	= "bp";
const char	c_szCC[]	= "cc";

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
	size_t								i, j;
	map<const CGene*,bool>				mapGenes;
	map<const CGene*,bool>::iterator	iterGene;

	if( cmdline_parser2( iArgs, aszArgs, &sArgs, 0, 1, 0 ) || ( sArgs.config_arg &&
		cmdline_parser_configfile( sArgs.config_arg, &sArgs, 0, 0, 1 ) ) ) {
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
		if( !GO.Open( ifsmOnto, ifsmAnno, Genome, eName ) ) {
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
		if( !KEGG.Open( ifsmOnto, Genome ) ) {
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

	ofsm.open( sArgs.output_arg, ios_base::binary );
	Dat.Save( ofsm, true );
	ofsm.close( );

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
			for( j = 0; j < Slim.GetGenes( i ); ++j )
				ofsm << Slim.GetGene( i, j ).GetName( ) << endl;
			ofsm.close( ); }
	if( sArgs.test_arg )
		for( iterGene = mapGenes.begin( ); iterGene != mapGenes.end( ); ++iterGene )
			if( iterGene->second )
				cout << iterGene->first->GetName( ) << endl;

	CMeta::Shutdown( );
	return 0; }

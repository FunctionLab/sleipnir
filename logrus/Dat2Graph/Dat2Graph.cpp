#include "stdafx.h"
#include "cmdline.h"

int open_genes( const char* szFile, CGenes& Genes ) {
	ifstream	ifsm;

	if( szFile ) {
		ifsm.open( szFile );
		if( !Genes.Open( ifsm ) ) {
			cerr << "Could not open: " << szFile << endl;
			return 1; } }

	return 0; }

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	ifstream			ifsm;
	CDat				Dat;
	float				d, dCutoff;
	CGenome				Genome;
	CGenes				GenesIn( Genome ), GenesQr( Genome );
	int					iRet;
	size_t				i, j;

	if( cmdline_parser2( iArgs, aszArgs, &sArgs, 0, 1, 0 ) && ( sArgs.config_arg &&
		cmdline_parser_configfile( sArgs.config_arg, &sArgs, 0, 0, 1 ) ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( sArgs.features_arg ) {
		ifsm.open( sArgs.features_arg );
		if( !Genome.Open( ifsm ) ) {
			cerr << "Could not open: " << sArgs.features_arg << endl;
			return 1; }
		ifsm.close( ); }

	if( iRet = open_genes( sArgs.genes_arg, GenesIn ) )
		return iRet;
	if( iRet = open_genes( sArgs.geneq_arg, GenesQr ) )
		return iRet;

	if( sArgs.input_arg ) {
		if( !Dat.Open( sArgs.input_arg, sArgs.memmap_flag && !( sArgs.normalize_flag ||
			sArgs.genes_arg || sArgs.geneq_arg ) ) ) {
			cerr << "Couldn't open: " << sArgs.input_arg << endl;
			return 1; } }
	else if( !Dat.Open( cin, CDat::EFormatText ) ) {
		cerr << "Couldn't open input" << endl;
		return 1; }
	if( GenesIn.GetGenes( ) )
		Dat.FilterGenes( GenesIn, CDat::EFilterInclude );
	if( GenesQr.GetGenes( ) ) {
		if( sArgs.cutoff_given )
			for( i = 0; i < Dat.GetGenes( ); ++i )
				for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
					if( !CMeta::IsNaN( d = Dat.Get( i, j ) ) && ( d < sArgs.cutoff_arg ) )
						Dat.Set( i, j, CMeta::GetNaN( ) );
		Dat.FilterGenes( GenesQr, CDat::EFilterPixie, sArgs.neighbors_arg ); }
	if( sArgs.normalize_flag )
		Dat.Normalize( );

	dCutoff = (float)( sArgs.cutoff_given ? sArgs.cutoff_arg : HUGE_VAL );
	if( !strcmp( sArgs.format_arg, "dot" ) )
		Dat.SaveDOT( cout, dCutoff, &Genome, true );
	else if( !strcmp( sArgs.format_arg, "gdf" ) )
		Dat.SaveGDF( cout, dCutoff );
	else if( !strcmp( sArgs.format_arg, "net" ) )
		Dat.SaveNET( cout, dCutoff );
	else if( !strcmp( sArgs.format_arg, "matisse" ) )
		Dat.SaveMATISSE( cout, dCutoff, &Genome );

	CMeta::Shutdown( );
	return 0; }

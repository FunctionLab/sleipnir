#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	ifstream			ifsm;
	CDat				Dat;
	float				dCutoff;
	CGenome				Genome;

	if( cmdline_parser2( iArgs, aszArgs, &sArgs, 0, 1, 0 ) || ( sArgs.config_arg &&
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

	if( sArgs.input_arg ) {
		if( !Dat.Open( sArgs.input_arg ) ) {
			cerr << "Couldn't open: " << sArgs.input_arg << endl;
			return 1; } }
	else if( !Dat.Open( cin, CDat::EFormatText ) ) {
		cerr << "Couldn't open input" << endl;
		return 1; }

	dCutoff = (float)( sArgs.cutoff_given ? sArgs.cutoff_arg : HUGE_VAL );
	if( !strcmp( sArgs.format_arg, "dot" ) )
		Dat.SaveDOT( cout, dCutoff, &Genome, true );
	else if( !strcmp( sArgs.format_arg, "gdf" ) )
		Dat.SaveGDF( cout, dCutoff );
	else if( !strcmp( sArgs.format_arg, "net" ) )
		Dat.SaveNET( cout, dCutoff );

	CMeta::Shutdown( );
	return 0; }

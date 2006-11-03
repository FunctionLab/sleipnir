#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CPCL				PCL;
	CDat				Dat;
	ofstream			ofsm;
	int					iRet;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( iRet = CPCL::Distance( sArgs.input_arg, sArgs.skip_arg, sArgs.distance_arg, !!sArgs.normalize_flag,
		!!sArgs.zscore_flag, !!sArgs.autocorrelate_flag, sArgs.genes_arg, sArgs.cutoff_given ? (float)sArgs.cutoff_arg :
		CMeta::GetNaN( ), sArgs.limit_arg, PCL, Dat ) ) {
		cmdline_parser_print_help( );
		return iRet; }

	if( sArgs.output_arg ) {
		ofsm.open( sArgs.output_arg, ios_base::binary );
		Dat.Save( ofsm, true );
		ofsm.close( ); }
	else {
		Dat.Save( cout, false );
		cout.flush( ); }

	CMeta::Shutdown( );
	return 0; }

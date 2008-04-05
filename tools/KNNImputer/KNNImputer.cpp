#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CPCL				PCL;
	CDat				Dat;
	int					iRet;
	ofstream			ofsm;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( iRet = CPCL::Distance( sArgs.input_arg, sArgs.skip_arg, sArgs.distance_arg, false, false,
		!!sArgs.autocorrelate_flag, sArgs.genes_arg, CMeta::GetNaN( ), sArgs.limit_arg, PCL, Dat ) ) {
		cmdline_parser_print_help( );
		return iRet; }

	PCL.Impute( sArgs.neighbors_arg, (float)sArgs.missing_arg, Dat );

	if( sArgs.output_arg ) {
		ofsm.open( sArgs.output_arg );
		PCL.Save( ofsm );
		ofsm.close( ); }
	else {
		PCL.Save( cout );
		cout.flush( ); }

	return 0; }

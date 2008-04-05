#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CPCL				PCL;
	CDat				Dat;
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
	if( sArgs.flip_flag )
		Dat.Invert( );

	Dat.Save( sArgs.output_arg );

	return 0; }

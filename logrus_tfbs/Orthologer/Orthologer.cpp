#include "stdafx.h"
#include "cmdline.h"

static const size_t	c_iConfirmation	= 1024;

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	COrthology			Orthology;
	ifstream			ifsm;
	size_t				i, j, k;
	vector<size_t>		veciHits, veciMap;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( sArgs.input_arg )
		ifsm.open( sArgs.input_arg );
	if( !Orthology.Open( sArgs.input_arg ? ifsm : cin ) ) {
		cerr << "Could not open: " << ( sArgs.input_arg || " standard input" ) << endl;
		return 1; }
	if( sArgs.input_arg )
		ifsm.close( );

	veciHits.resize( Orthology.GetGenomes( ) );
	veciMap.resize( sArgs.inputs_num );
	for( i = 0; i < sArgs.inputs_num; ++i ) {
		CDat	Answers;

		ifsm.clear( );
		ifsm.open( sArgs.inputs[ i ] );
		if( !Answers.OpenGenes( ifsm, true ) ) {
			cerr << "Could not open: " << sArgs.input_arg << endl;
			return 1; }
		ifsm.close( );

		for( j = 0; j < Orthology.GetGenomes( ); ++j ) {
			const CGenome*	pGenome	= Orthology.GetGenome( j );

			veciHits[ j ] = 0;
			for( k = 0; k < min( Answers.GetGenes( ), c_iConfirmation ); ++k )
				if( pGenome->FindGene( Answers.GetGene( k ) ) != -1 )
					veciHits[ j ]++; }

		veciMap[ i ] = 0;
		for( j = 1; j < veciHits.size( ); ++j )
			if( veciHits[ j ] > veciHits[ veciMap[ i ] ] )
				veciMap[ i ] = j; }

	for( i = 0; i < sArgs.inputs_num; ++i )
		cerr << sArgs.inputs[ i ] << ": " << Orthology.GetOrganism( veciMap[ i ] ) << endl;

	CMeta::Shutdown( );
	return 0; }

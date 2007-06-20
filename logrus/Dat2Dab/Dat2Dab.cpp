#include "stdafx.h"
#include "cmdline.h"
#include "mathb.h"

int Dat2Dab( const CGenes&, const gengetopt_args_info& );

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CGenome				Genome;
	CGenes				Genes( Genome );
	ifstream			ifsm;
	CDat				Dat;
	size_t				i, j;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( sArgs.genes_arg ) {
		ifsm.open( sArgs.genes_arg );
		if( !Genes.Open( ifsm ) ) {
			cerr << "Could not open: " << sArgs.genes_arg << endl;
			return 1; }
		ifsm.close( ); }

	if( sArgs.input_arg ) {
		if( !Dat.Open( sArgs.input_arg, !!sArgs.memmap_flag ) ) {
			cerr << "Could not open: " << sArgs.input_arg << endl;
			return 1; } }
	else if( !Dat.Open( cin, CDat::EFormatText, (float)HUGE_VAL, !!sArgs.duplicates_flag ) ) {
		cerr << "Could not open input" << endl;
		return 1; }

/*
float d;
for( i = 0; i < Dat.GetGenes( ); ++i )
for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
if( !CMeta::IsNaN( d = Dat.Get( i, j ) ) )
Dat.Set( i, j, (float)CMath::Sigmoid(
//1,4.548,0.4906,0 // pixie
//1,16.77,0.1213,0 // gerstein
//1,0.8945,0.5995,0 // lee
//1,3.406,0.5154,0 // albert
//1,9.033,0.2307,0 // mefit
//1,4,1.5,0 // tan
//1,0,1,0
,d ) );
Dat.Save( sArgs.input_arg, true );
return 0;
//*/

	if( sArgs.rank_flag )
		Dat.Rank( );
	if( sArgs.normalize_flag || sArgs.zscore_flag )
		Dat.Normalize( !!sArgs.normalize_flag );
	if( sArgs.zero_flag )
		for( i = 0; i < Dat.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
				if( CMeta::IsNaN( Dat.Get( i, j ) ) )
					Dat.Set( i, j, 0 );
	if( sArgs.flip_flag )
		Dat.Invert( );
	if( Genes.GetGenes( ) )
		Dat.FilterGenes( Genes, CDat::EFilterInclude );

	if( sArgs.cutoff_given )
		for( i = 0; i < Dat.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
				if( Dat.Get( i, j ) < sArgs.cutoff_arg )
					Dat.Set( i, j, CMeta::GetNaN( ) );

	if( sArgs.output_arg )
		Dat.Save( sArgs.output_arg );
	else {
		Dat.Save( cout, CDat::EFormatText );
		cout.flush( ); }

	CMeta::Shutdown( );
	return 0; }

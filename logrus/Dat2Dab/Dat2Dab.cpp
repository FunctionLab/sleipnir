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

	if( sArgs.remap_arg ) {
		static const size_t	c_iBuffer	= 1024;
		char								acBuffer[ c_iBuffer ];
		vector<string>						vecstrTokens;
		map<string,string>					mapNames;
		map<string,string>::const_iterator	iterName;

		ifsm.clear( );
		ifsm.open( sArgs.remap_arg );
		if( !ifsm.is_open( ) ) {
			cerr << "Could not open: " << sArgs.remap_arg << endl;
			return 1; }
		while( !ifsm.eof( ) ) {
			ifsm.getline( acBuffer, c_iBuffer - 1 );
			acBuffer[ c_iBuffer - 1 ] = 0;
			vecstrTokens.clear( );
			CMeta::Tokenize( acBuffer, vecstrTokens );
			if( vecstrTokens.empty( ) )
				continue;
			if( vecstrTokens.size( ) != 2 ) {
				cerr << "Illegal remap line (" << vecstrTokens.size( ) << "): " << acBuffer << endl;
				return 1; }
			if( vecstrTokens[ 0 ] == vecstrTokens[ 1 ] )
				continue;
			if( ( iterName = mapNames.find( vecstrTokens[ 0 ] ) ) == mapNames.end( ) )
				mapNames[ vecstrTokens[ 0 ] ] = vecstrTokens[ 1 ];
			else if( iterName->second != vecstrTokens[ 1 ] ) {
				cerr << "Ambiguous mapping: " << vecstrTokens[ 0 ] << " to " << iterName->second <<
					", " << vecstrTokens[ 1 ] << endl;
				return 1; } }

		for( i = 0; i < Dat.GetGenes( ); ++i )
			if( ( iterName = mapNames.find( Dat.GetGene( i ) ) ) != mapNames.end( ) ) {
				if( Dat.GetGene( iterName->second ) ) {
					cerr << "Duplicate mapping: " << Dat.GetGene( i ) << " to " <<
						iterName->second << endl;
					return 1; }
				Dat.SetGene( i, iterName->second ); } }

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

	if( sArgs.lookup1_given && sArgs.lookup2_given ) {
		if( ( i = Dat.GetGene( sArgs.lookup1_arg ) ) == -1 ) {
			cerr << "Unknown gene: " << sArgs.lookup1_arg << endl;
			return 1; }
		if( ( j = Dat.GetGene( sArgs.lookup2_arg ) ) == -1 ) {
			cerr << "Unknown gene: " << sArgs.lookup2_arg << endl;
			return 1; }
		cout << Dat.GetGene( i ) << '\t' << Dat.GetGene( j ) << '\t' << Dat.Get( i, j ) << endl; }
	else if( sArgs.output_arg )
		Dat.Save( sArgs.output_arg );
	else {
		Dat.Save( cout, CDat::EFormatText );
		cout.flush( ); }

	CMeta::Shutdown( );
	return 0; }

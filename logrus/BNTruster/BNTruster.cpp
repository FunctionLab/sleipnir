#include "stdafx.h"
#include "cmdline.h"

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info		sArgs;
	size_t					i, iDSL;
	vector<unsigned char>	vecbDatum;
	vector<float>			vecdOut;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );
#if !( defined(_MSC_VER) && defined(_DEBUG) )
	EnableXdslFormat( );
#endif // !( defined(_MSC_VER) && defined(_DEBUG) )

	for( iDSL = 0; iDSL < sArgs.inputs_num; ++iDSL ) {
		CBayesNetSmile	BNSmile;
		size_t			iNode;
		unsigned char	bValue;

		if( !BNSmile.Open( sArgs.inputs[ iDSL ] ) ) {
			cerr << "Couldn't open: " << sArgs.inputs[ iDSL ] << endl;
			return 1; }
		vecbDatum.resize( BNSmile.GetNodes( ).size( ) );

		if( !iDSL ) {
			for( iNode = 1; iNode < BNSmile.GetNodes( ).size( ); ++iNode )
				cout << '\t' << BNSmile.GetNodes( )[ iNode ];
			cout << endl; }

		cout << sArgs.inputs[ iDSL ];
		for( iNode = 1; iNode < BNSmile.GetNodes( ).size( ); ++iNode ) {
			float	d;

			for( i = 0; i < vecbDatum.size( ); ++i )
				vecbDatum[ i ] = 0;
			vecdOut.clear( );
			for( bValue = 0; bValue < BNSmile.GetValues( iNode ); ++bValue ) {
				vecbDatum[ iNode ] = bValue + 1;
				BNSmile.Evaluate( vecbDatum, vecdOut, false ); }

			d = 0;
			for( i = 0; i < vecdOut.size( ); ++i )
				d += fabs( 0.5f - vecdOut[ i ] );
			cout << '\t' << ( d / BNSmile.GetValues( iNode ) ); }
		cout << endl; }

	CMeta::Shutdown( );
	return 0; }

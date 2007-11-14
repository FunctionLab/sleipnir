#include "stdafx.h"

static int MainPosteriors( const gengetopt_args_info& );
static int MainSums( const gengetopt_args_info& );
static int MainRatios( const gengetopt_args_info& );

static const TPFnTruster	c_apfnTrusters[]	= { MainPosteriors, MainSums, MainRatios, NULL };
static const char*			c_aszTrusters[]		= { "posteriors", "sums", "ratios", NULL };

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	size_t				i;
	int					iRet;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );
	EnableXdslFormat( );

	for( iRet = 1,i = 0; c_aszTrusters[ i ]; ++i )
		if( !strcmp( c_aszTrusters[ i ], sArgs.type_arg ) ) {
			iRet = c_apfnTrusters[ i ]( sArgs );
			break; }

	CMeta::Shutdown( );
	return iRet; }

int MainPosteriors( const gengetopt_args_info& sArgs ) {
	size_t					i, iDSL;
	vector<unsigned char>	vecbDatum;
	vector<float>			vecdOut, vecdCur;
	vector<string>			vecstrNodes;
	map<string,size_t>		mapstriNodes;

	for( iDSL = 0; iDSL < sArgs.inputs_num; ++iDSL ) {
		CBayesNetSmile	BNSmile;
		vector<string>	vecstrCur;

		if( !BNSmile.Open( sArgs.inputs[ iDSL ] ) ) {
			cerr << "Couldn't open: " << sArgs.inputs[ iDSL ] << endl;
			return 1; }
		BNSmile.GetNodes( vecstrCur );
		for( i = 1; i < vecstrCur.size( ); ++i )
			if( mapstriNodes.find( vecstrCur[ i ] ) == mapstriNodes.end( ) ) {
				mapstriNodes[ vecstrCur[ i ] ] = vecstrNodes.size( );
				vecstrNodes.push_back( vecstrCur[ i ] ); } }
	if( !sArgs.bins_flag ) {
		for( i = 0; i < vecstrNodes.size( ); ++i )
			cout << '\t' << vecstrNodes[ i ];
		cout << endl; }

	vecdCur.resize( vecstrNodes.size( ) );
	for( iDSL = 0; iDSL < sArgs.inputs_num; ++iDSL ) {
		CBayesNetSmile	BNSmile;
		size_t			iNode;
		unsigned char	bValue;
		float			dPrior;
		vector<string>	vecstrCur;
		CDataMatrix		MatFR;

		if( !BNSmile.Open( sArgs.inputs[ iDSL ] ) ) {
			cerr << "Couldn't open: " << sArgs.inputs[ iDSL ] << endl;
			return 1; }
		cerr << "Opened: " << sArgs.inputs[ iDSL ] << endl;
		for( i = 0; i < vecdCur.size( ); ++i )
			vecdCur[ i ] = CMeta::GetNaN( );
		BNSmile.GetNodes( vecstrCur );
		vecbDatum.resize( vecstrCur.size( ) );
		for( i = 0; i < vecbDatum.size( ); ++i )
			vecbDatum[ i ] = 0;

		vecdOut.clear( );
		BNSmile.Evaluate( vecbDatum, vecdOut, false );
		dPrior = 1 - vecdOut[ 0 ];
		for( iNode = 1; iNode < vecstrCur.size( ); ++iNode ) {
			float			d, dSum;
			CDataMatrix		MatCPT;
			vector<float>	vecdProbs;

			vecbDatum[ iNode - 1 ] = 0;
			vecdProbs.clear( );
			BNSmile.Evaluate( vecbDatum, vecdProbs, false, iNode );
			for( dSum = 0,i = 0; i < vecdProbs.size( ); ++i )
				dSum += vecdProbs[ i ];
			vecdProbs.push_back( 1 - dSum );
			vecdOut.clear( );
			for( bValue = 0; bValue < BNSmile.GetValues( iNode ); ++bValue ) {
				vecbDatum[ iNode ] = bValue + 1;
				BNSmile.Evaluate( vecbDatum, vecdOut, false ); }

			for( dSum = 0,i = 0; i < vecdOut.size( ); ++i ) {
				vecdOut[ i ] = 1 - vecdOut[ i ];
				d = fabs( dPrior - vecdOut[ i ] );
				if( sArgs.bins_flag )
					cout << iDSL << '\t' << ( iNode - 1 ) << '\t' << i << '\t' << ( ( vecdOut[ i ] > dPrior ) ?
						( d / ( 1 - dPrior ) ) : -( d / dPrior ) ) << endl;
				dSum += d * vecdProbs[ i ]; }
			vecdCur[ mapstriNodes[ vecstrCur[ iNode ] ] ] = dSum / BNSmile.GetValues( iNode ); }

		if( !sArgs.bins_flag ) {
			cout << sArgs.inputs[ iDSL ];
			for( i = 0; i < vecdCur.size( ); ++i ) {
				cout << '\t';
				if( !CMeta::IsNaN( vecdCur[ i ] ) )
					cout << vecdCur[ i ]; }
			cout << endl; } }

	return 0; }

int MainSums( const gengetopt_args_info& sArgs ) {
	size_t	iDSL;

	for( iDSL = 0; iDSL < sArgs.inputs_num; ++iDSL ) {
		CBayesNetSmile	BNSmile;
		size_t			iNode;
		unsigned char	bValue;
		float			dSum;
		vector<string>	vecstrNodes;

		if( !BNSmile.Open( sArgs.inputs[ iDSL ] ) ) {
			cerr << "Couldn't open: " << sArgs.inputs[ iDSL ] << endl;
			return 1; }
		BNSmile.GetNodes( vecstrNodes );
		if( !iDSL ) {
			for( iNode = 1; iNode < vecstrNodes.size( ); ++iNode )
				cout << '\t' << vecstrNodes[ iNode ];
			cout << endl; }

		cout << sArgs.inputs[ iDSL ];
		for( iNode = 1; iNode < vecstrNodes.size( ); ++iNode ) {
			CDataMatrix	MatCPT;

			dSum = 0;
			BNSmile.GetCPT( iNode, MatCPT );
			for( bValue = 0; bValue < MatCPT.GetRows( ); ++bValue )
				dSum += fabs( MatCPT.Get( bValue, 0 ) - MatCPT.Get( bValue, 1 ) );
			cout << '\t' << dSum; }
		cout << endl; }

	return 0; }

int MainRatios( const gengetopt_args_info& sArgs ) {
	size_t	iDSL;

	for( iDSL = 0; iDSL < sArgs.inputs_num; ++iDSL ) {
		CBayesNetSmile	BNSmile;
		size_t			iNode;
		unsigned char	bValue;
		float			dProd;
		vector<string>	vecstrNodes;

		if( !BNSmile.Open( sArgs.inputs[ iDSL ] ) ) {
			cerr << "Couldn't open: " << sArgs.inputs[ iDSL ] << endl;
			return 1; }
		BNSmile.GetNodes( vecstrNodes );
		if( !iDSL ) {
			for( iNode = 1; iNode < vecstrNodes.size( ); ++iNode )
				cout << '\t' << vecstrNodes[ iNode ];
			cout << endl; }

		cout << sArgs.inputs[ iDSL ];
		for( iNode = 1; iNode < vecstrNodes.size( ); ++iNode ) {
			CDataMatrix	MatCPT;
			float		dMin, dMax;

			dProd = 1;
			BNSmile.GetCPT( iNode, MatCPT );
			for( bValue = 0; bValue < MatCPT.GetRows( ); ++bValue ) {
				dMin = MatCPT.Get( bValue, 0 );
				dMax = MatCPT.Get( bValue, 1 );
				if( dMin > dMax )
					swap( dMin, dMax );
				dProd *= dMax / dMin; }
			cout << '\t' << log( dProd ); }
		cout << endl; }

	return 0; }

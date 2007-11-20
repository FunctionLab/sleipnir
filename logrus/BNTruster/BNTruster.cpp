#include "stdafx.h"

static int MainPosteriors( const gengetopt_args_info& );
static int MainSums( const gengetopt_args_info& );
static int MainRatios( const gengetopt_args_info& );

static const TPFnTruster	c_apfnTrusters[]	= { MainPosteriors, MainSums, MainRatios, NULL };
static const char*			c_aszTrusters[]		= { "posteriors", "sums", "ratios", NULL };

int main( int iArgs, char** aszArgs ) {
#ifdef WIN32
	pthread_win32_process_attach_np( );
#endif // WIN32
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
	pthread_exit( NULL );
#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return iRet; }

struct SData {
	const char*					m_szDSL;
	const map<string, size_t>*	m_pmapstriNodes;
	bool						m_fBins;
	vector<float>				m_vecdResults;
	vector<size_t>				m_veciBins;
};

void* posteriorsOneDSL( void* pData ) {
	SData*					psData;
	CBayesNetSmile			BNSmile;
	size_t					iNode;
	unsigned char			bValue;
	float					dPrior;
	vector<string>			vecstrNodes;
	CDataMatrix				MatFR;
	vector<unsigned char>	vecbDatum;
	vector<float>			vecdOut;
	size_t					i, j, iBins, iResult;

	psData = (SData*)pData;
	if( !BNSmile.Open( psData->m_szDSL ) ) {
		cerr << "Couldn't open: " << psData->m_szDSL << endl;
		return (void*)1; }
	cerr << "Opened: " << psData->m_szDSL << endl;
	BNSmile.GetNodes( vecstrNodes );
	vecbDatum.resize( vecstrNodes.size( ) );
	fill( vecbDatum.begin( ), vecbDatum.end( ), 0 );
	if( psData->m_fBins ) {
		psData->m_veciBins.resize( vecstrNodes.size( ) - 1 );
		for( iBins = 0,i = 1; i < vecstrNodes.size( ); ++i ) {
			psData->m_veciBins[ i - 1 ] = j = BNSmile.GetValues( i );
			iBins += j; }
		psData->m_vecdResults.resize( iBins ); }
	else {
		psData->m_vecdResults.resize( vecstrNodes.size( ) );
		fill( psData->m_vecdResults.begin( ), psData->m_vecdResults.end( ), CMeta::GetNaN( ) ); }

	vecdOut.clear( );
	BNSmile.Evaluate( vecbDatum, vecdOut, false );
	dPrior = 1 - vecdOut[ 0 ];
	for( iResult = 0,iNode = 1; iNode < vecstrNodes.size( ); ++iNode ) {
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
			if( psData->m_fBins )
				psData->m_vecdResults[ iResult++ ] = ( vecdOut[ i ] > dPrior ) ? ( d / ( 1 - dPrior ) ) :
					-( d / dPrior );
			else
				dSum += d * vecdProbs[ i ]; }
		if( !psData->m_fBins )
			psData->m_vecdResults[ psData->m_pmapstriNodes->find( vecstrNodes[ iNode ] )->second ] =
				dSum / BNSmile.GetValues( iNode ); }

	return NULL; }

int MainPosteriors( const gengetopt_args_info& sArgs ) {
	size_t					i, j, k, iDSLBase, iDSLOffset;
	vector<string>			vecstrNodes;
	map<string,size_t>		mapstriNodes;
	vector<pthread_t>		vecpthdThreads;
	vector<SData>			vecsData;

	for( iDSLOffset = 0; iDSLOffset < sArgs.inputs_num; ++iDSLOffset ) {
		CBayesNetSmile	BNSmile;
		vector<string>	vecstrCur;

		if( !BNSmile.Open( sArgs.inputs[ iDSLOffset ] ) ) {
			cerr << "Couldn't open: " << sArgs.inputs[ iDSLOffset ] << endl;
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

	vecpthdThreads.resize( sArgs.threads_arg );
	vecsData.resize( vecpthdThreads.size( ) );
	for( iDSLBase = 0; iDSLBase < sArgs.inputs_num; iDSLBase += sArgs.threads_arg ) {
		for( iDSLOffset = 0; ( iDSLOffset < (size_t)sArgs.threads_arg ) && ( ( iDSLBase + iDSLOffset ) <
			sArgs.inputs_num ); ++iDSLOffset ) {
			vecsData[ iDSLOffset ].m_fBins = !!sArgs.bins_flag;
			vecsData[ iDSLOffset ].m_pmapstriNodes = &mapstriNodes;
			vecsData[ iDSLOffset ].m_szDSL = sArgs.inputs[ iDSLBase + iDSLOffset ];
			if( pthread_create( &vecpthdThreads[ iDSLOffset ], NULL, posteriorsOneDSL,
				&vecsData[ iDSLOffset ] ) ) {
				cerr << "Couldn't create thread: " << sArgs.inputs[ iDSLBase + iDSLOffset ] << endl;
				return 1; } }
		for( i = 0; i < iDSLOffset; ++i ) {
			void*			pValue;
			float			d;
			size_t			iResult;
			const SData&	sDatum	= vecsData[ i ];

			pthread_join( vecpthdThreads[ i ], &pValue );
			if( pValue )
				return 1;
			if( sArgs.bins_flag )
				for( iResult = j = 0; j < sDatum.m_veciBins.size( ); ++j )
					for( k = 0; k < sDatum.m_veciBins[ j ]; ++k )
						cout << ( iDSLBase + i ) << '\t' << j << '\t' << k << '\t' <<
							sDatum.m_vecdResults[ iResult++ ] << endl;
			else {
				cout << sArgs.inputs[ iDSLBase + i ];
				for( j = 0; j < sDatum.m_vecdResults.size( ); ++j ) {
					cout << '\t';
					if( !CMeta::IsNaN( d = sDatum.m_vecdResults[ j ] ) )
						cout << d; }
				cout << endl; } } }

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

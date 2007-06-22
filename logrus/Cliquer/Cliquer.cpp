#include "stdafx.h"
#include "cmdline.h"

int cliques( const gengetopt_args_info&, const CDat&, const CDat&, const vector<size_t>& );
int heavy( const gengetopt_args_info&, CDat&, const CDat&, const vector<size_t>& );
bool connectivity( size_t, const vector<size_t>&, const vector<float>&, const vector<size_t>&,
	float, size_t, float, size_t, const CDat&, float&, size_t&, float&, size_t& );
void max_connectivity( const vector<bool>&, const vector<size_t>&, const vector<float>&,
	const vector<size_t>&, float, size_t, float, size_t, size_t, const CDat&, float&, size_t&, float&,
	size_t&, float&, size_t& );

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CDat				Dat, DatKnowns;
	vector<size_t>		veciGenes, veciKnowns;
	size_t				i;
	int					iRet;
	
	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( sArgs.input_arg ) {
		if( !Dat.Open( sArgs.input_arg, sArgs.memmap_flag && !sArgs.normalize_flag ) ) {
			cerr << "Could not open: " << sArgs.input_arg << endl;
			return 1; } }
	else if( !Dat.Open( cin, CDat::EFormatText ) ) {
		cerr << "Could not open input" << endl;
		return 1; }
	if( sArgs.normalize_flag )
		Dat.Normalize( );
	if( sArgs.knowns_arg ) {
		if( !DatKnowns.Open( sArgs.knowns_arg, !!sArgs.memmap_flag ) ) {
			cerr << "Could not open: " << sArgs.knowns_arg << endl;
			return 1; }
		veciKnowns.resize( Dat.GetGenes( ) );
		for( i = 0; i < veciKnowns.size( ); ++i )
			veciKnowns[ i ] = DatKnowns.GetGene( Dat.GetGene( i ) ); }

	iRet = sArgs.heavy_arg ? heavy( sArgs, Dat, DatKnowns, veciKnowns ) :
		cliques( sArgs, Dat, DatKnowns, veciKnowns );

	CMeta::Shutdown( );
	return iRet; }

int cliques( const gengetopt_args_info& sArgs, const CDat& Dat, const CDat& DatKnowns,
	const vector<size_t>& veciKnowns ) {
	vector<pair<vector<size_t>,float> >	vecprCliques;
	vector<size_t>						veciGenes;
	float								d, dCur;
	size_t								i, j, iOne, iTwo, iMin;
	int									iCol;

	vecprCliques.resize( sArgs.subgraphs_arg );
	for( iMin = i = 0; i < vecprCliques.size( ); ++i ) {
		vecprCliques[ i ].first.resize( sArgs.size_arg );
		vecprCliques[ i ].second = -FLT_MAX; }
	veciGenes.resize( sArgs.size_arg );
	for( i = 0; i < veciGenes.size( ); ++i )
		veciGenes[ i ] = i;
	while( veciGenes[ 0 ] <= ( Dat.GetGenes( ) - veciGenes.size( ) ) ) {
		if( !rand( ) && ( ( (float)rand( ) / RAND_MAX ) < 1e-3 ) ) {
			for( i = 0; i < veciGenes.size( ); ++i )
				cerr << veciGenes[ i ] << ',';
			cerr << endl; }
		for( dCur = 0,i = 0; i < veciGenes.size( ); ++i ) {
			iOne = sArgs.knowns_arg ? veciKnowns[ veciGenes[ i ] ] : -1;
			for( j = ( i + 1 ); j < veciGenes.size( ); ++j ) {
				if( ( iOne != -1 ) && ( ( iTwo = veciKnowns[ veciGenes[ j ] ] ) != -1 ) &&
					!CMeta::IsNaN( d = DatKnowns.Get( iOne, iTwo ) ) && ( d > 0 ) )
					continue;
				if( !CMeta::IsNaN( d = Dat.Get( veciGenes[ i ], veciGenes[ j ] ) ) )
					dCur += d; } }
		if( dCur > vecprCliques[ iMin ].second ) {

			copy( veciGenes.begin( ), veciGenes.end( ), vecprCliques[ iMin ].first.begin( ) );
			vecprCliques[ iMin ].second = dCur;
			for( iMin = i = 0; i < vecprCliques.size( ); ++i )
				if( vecprCliques[ i ].second < vecprCliques[ iMin ].second )
					iMin = i; }
		for( i = 0; i < veciGenes.size( ); ++i )
			if( ++veciGenes[ veciGenes.size( ) - 1 - i ] < ( Dat.GetGenes( ) - i ) ) {
				for( iCol = ( (int)veciGenes.size( ) - (int)i ); ( iCol > 0 ) &&
					( (size_t)iCol < veciGenes.size( ) ); ++iCol )
					veciGenes[ iCol ] = veciGenes[ iCol - 1 ] + 1;
				break; } }

	for( i = 0; i < vecprCliques.size( ); ++i ) {
		cout << vecprCliques[ i ].second;
		for( j = 0; j < vecprCliques[ i ].first.size( ); ++j )
			cout << '\t' << Dat.GetGene( vecprCliques[ i ].first[ j ] );
		cout << endl; }

	return 0; }

struct SSeed {
	size_t	m_iOne;
	size_t	m_iTwo;
	float	m_dIn;
	size_t	m_iIn;
	float	m_dOut;
	size_t	m_iOut;
	float	m_dRatio;

	SSeed( size_t iOne, size_t iTwo, float dIn, size_t iIn, float dOut, size_t iOut, float dRatio ) :
		m_iOne(iOne), m_iTwo(iTwo), m_dIn(dIn), m_iIn(iIn), m_dOut(dOut), m_iOut(iOut),
		m_dRatio(dRatio) { }

	bool operator<( const SSeed& sSeed ) const {

		return ( m_dRatio < sSeed.m_dRatio ); }
};

int heavy( const gengetopt_args_info& sArgs, CDat& Dat, const CDat& DatKnowns,
	const vector<size_t>& veciKnowns ) {
	vector<bool>	vecfCluster;
	vector<float>	vecdConnectivity;
	vector<size_t>	veciConnectivity, veciCluster;
	size_t			i, j, iIn, iOut, iMax, iMaxIn, iMaxOut;
	size_t			iClusters;
	float			d, dMaxIn, dMaxOut, dIn, dOut, dRatio;

	veciConnectivity.resize( Dat.GetGenes( ) );
	vecdConnectivity.resize( Dat.GetGenes( ) );
	for( i = 0; i < veciConnectivity.size( ); ++i )
		for( j = ( i + 1 ); j < veciConnectivity.size( ); ++j )
			if( !CMeta::IsNaN( d = Dat.Get( i, j ) ) ) {
				veciConnectivity[ i ]++;
				veciConnectivity[ j ]++;
				vecdConnectivity[ i ] += d;
				vecdConnectivity[ j ] += d; }

	iClusters = 0;
	vecfCluster.resize( Dat.GetGenes( ) );
	vecdConnectivity.resize( Dat.GetGenes( ) );
	while( ( sArgs.subgraphs_arg == -1 ) || ( iClusters < sArgs.subgraphs_arg ) ) {
		priority_queue<SSeed>	pqueSeeds;

		d = 4;
		veciCluster.resize( 1 );
		for( i = 0; i < Dat.GetGenes( ); ++i ) {
			veciCluster[ 0 ] = i;
			for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
				if( connectivity( j, veciCluster, vecdConnectivity, veciConnectivity, 0, 0,
					vecdConnectivity[ i ], veciConnectivity[ i ], Dat, dIn, iIn, dOut, iOut ) &&
					( ( dRatio = ( ( dIn / iIn ) / ( dOut / iOut ) ) ) > ( d / 2 ) ) ) {
					d = max( dRatio, 4.0f );
					pqueSeeds.push( SSeed( i, j, dIn, iIn, dOut, iOut, dRatio ) ); } }
		if( pqueSeeds.empty( ) )
			break;
		cerr << "Seeds remaining: " << pqueSeeds.size( ) << endl;

		do {
			const SSeed&	sSeed	= pqueSeeds.top( );

			for( i = 0; i < veciCluster.size( ); ++i )
				vecfCluster[ veciCluster[ i ] ] = false;
			vecfCluster[ sSeed.m_iOne ] = vecfCluster[ sSeed.m_iTwo ] = true;
			veciCluster.resize( 2 );
			veciCluster[ 0 ] = sSeed.m_iOne;
			veciCluster[ 1 ] = sSeed.m_iTwo;

			cerr << "Cluster " << iClusters << " seed: " << Dat.GetGene( sSeed.m_iOne ) << ", " <<
				Dat.GetGene( sSeed.m_iTwo ) << ", " << sSeed.m_dRatio << endl;
			dIn = sSeed.m_dIn;
			iIn = sSeed.m_iIn;
			dOut = sSeed.m_dOut;
			iOut = sSeed.m_iOut;
			while( true ) {
				cerr << "Cluster " << iClusters << ", " << veciCluster.size( ) << " genes" << endl;
				max_connectivity( vecfCluster, veciCluster, vecdConnectivity, veciConnectivity, dIn, iIn,
					dOut, iOut, -1, Dat, dRatio, iMax, dMaxIn, iMaxIn, dMaxOut, iMaxOut );
				if( dRatio >= ( sArgs.heavy_arg * sSeed.m_dRatio ) ) {
					vecfCluster[ iMax ] = true;
					veciCluster.push_back( iMax );
					dIn = dMaxIn;
					iIn = iMaxIn;
					dOut = dMaxOut;
					iOut = iMaxOut; }
				else
					break; }
			pqueSeeds.pop( ); }
		while( !pqueSeeds.empty( ) && ( veciCluster.size( ) < 3 ) );
		if( veciCluster.size( ) < 3 )
			break;

		cerr << "Found cluster: " << ( ( dIn / iIn ) / ( dOut / iOut ) ) << " (" << dIn << '/' << iIn <<
			", " << dOut << '/' << iOut << ')' << endl;
		iClusters++;
		cout << ( ( dIn / iIn ) / ( dOut / iOut ) );
		for( i = 0; i < veciCluster.size( ); ++i )
			cout << '\t' << Dat.GetGene( veciCluster[ i ] );
		cout << endl;
		cout.flush( );

		dRatio = dIn / iIn;
		for( i = 0; i < veciCluster.size( ); ++i )
			for( j = ( i + 1 ); j < veciCluster.size( ); ++j )
				if( !CMeta::IsNaN( d = Dat.Get( veciCluster[ i ], veciCluster[ j ] ) ) ) {
					dIn = min( dRatio, d );
					Dat.Set( veciCluster[ i ], veciCluster[ j ], d - dIn );
					vecdConnectivity[ veciCluster[ i ] ] -= dIn;
					vecdConnectivity[ veciCluster[ j ] ] -= dIn; } }

	return 0; }

bool connectivity( size_t iNew, const vector<size_t>& veciCluster,
	const vector<float>& vecdConnectivity, const vector<size_t>& veciConnectivity, float dIn,
	size_t iIn, float dOut, size_t iOut, const CDat& Dat, float& dSumIn, size_t& iEdgesIn,
	float& dSumOut, size_t& iEdgesOut ) {
	size_t	i;
	float	d;
	bool	fRet;

	iEdgesOut = iOut + veciConnectivity[ iNew ];
	dSumOut = dOut + vecdConnectivity[ iNew ];
	iEdgesIn = iIn;
	dSumIn = dIn;
	fRet = false;
	for( i = 0; i < veciCluster.size( ); ++i )
		if( !CMeta::IsNaN( d = Dat.Get( iNew, veciCluster[ i ] ) ) ) {
			fRet = true;
			iEdgesOut--;
			dSumOut -= d;
			iEdgesIn++;
			dSumIn += d; }

	return fRet; }

void max_connectivity( const vector<bool>& vecfCluster, const vector<size_t>& veciCluster,
	const vector<float>& vecdConnectivity, const vector<size_t>& veciConnectivity, float dIn,
	size_t iIn, float dOut, size_t iOut, size_t iStart, const CDat& Dat, float& dMaxRatio, size_t& iMax,
	float& dMaxIn, size_t& iMaxIn, float& dMaxOut, size_t& iMaxOut ) {
	size_t	i, iEdgesIn, iEdgesOut;
	float	dRatio, dSumIn, dSumOut;

	dMaxRatio = -FLT_MAX;
	for( iMax = 0,i = ( iStart == -1 ) ? 0 : ( iStart + 1 ); i < vecfCluster.size( ); ++i ) {
		if( vecfCluster[ i ] || !connectivity( i, veciCluster, vecdConnectivity, veciConnectivity, dIn,
			iIn, dOut, iOut, Dat, dSumIn, iEdgesIn, dSumOut, iEdgesOut ) )
			continue;
		dRatio = ( ( dSumIn / iEdgesIn ) / ( dSumOut / iEdgesOut ) );
		if( dRatio > dMaxRatio ) {
			dMaxRatio = dRatio;
			iMax = i;
			dMaxIn = dSumIn;
			iMaxIn = iEdgesIn;
			dMaxOut = dSumOut;
			iMaxOut = iEdgesOut; } } }

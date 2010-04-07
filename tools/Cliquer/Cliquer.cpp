/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "stdafx.h"
#include "cmdline.h"

int cliques( const gengetopt_args_info&, const CDat&, const CDat&, const vector<size_t>& );
int heavy( const gengetopt_args_info&, CDat&, const CDat&, const vector<size_t>& );
int heavy2( const gengetopt_args_info&, CDat&, const CDat&, const vector<size_t>& );
int motifs( const gengetopt_args_info&, CDat& );
bool connectivity( size_t, const vector<size_t>&, const vector<float>&, const vector<size_t>&,
	float, size_t, float, size_t, const CDat&, float&, size_t&, float&, size_t& );
void max_connectivity( const vector<bool>&, const vector<size_t>&, const vector<float>&,
	const vector<size_t>&, float, size_t, float, size_t, size_t, const CDat&, float&, size_t&, float&,
	size_t&, float&, size_t& );

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CDat				Dat, DatKnowns;
	vector<size_t>		veciGenes, veciKnowns;
	size_t				i, j;
	float				d;
	int					iRet;
	
	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg );

	if( sArgs.input_arg ) {
		if( !Dat.Open( sArgs.input_arg, sArgs.memmap_flag && !sArgs.normalize_flag &&
			!sArgs.heavy_arg && !sArgs.cutoff_arg ) ) {
			cerr << "Could not open: " << sArgs.input_arg << endl;
			return 1; } }
	else if( !Dat.Open( cin, CDat::EFormatText ) ) {
		cerr << "Could not open input" << endl;
		return 1; }
	if( sArgs.normalize_flag )
		Dat.Normalize( CDat::ENormalizeSigmoid );
	if( sArgs.cutoff_arg )
		for( i = 0; i < Dat.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
				if( !CMeta::IsNaN( d = Dat.Get( i, j ) ) && ( d < sArgs.cutoff_arg ) )
					Dat.Set( i, j, CMeta::GetNaN( ) );
	if( sArgs.knowns_arg ) {
		if( !DatKnowns.Open( sArgs.knowns_arg, !!sArgs.memmap_flag ) ) {
			cerr << "Could not open: " << sArgs.knowns_arg << endl;
			return 1; }
		veciKnowns.resize( Dat.GetGenes( ) );
		for( i = 0; i < veciKnowns.size( ); ++i )
			veciKnowns[ i ] = DatKnowns.GetGene( Dat.GetGene( i ) ); }

	iRet = sArgs.motifs_arg ? motifs( sArgs, Dat ) : ( sArgs.heavy_arg ?
		heavy2( sArgs, Dat, DatKnowns, veciKnowns ) :
		cliques( sArgs, Dat, DatKnowns, veciKnowns ) );

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
					!CMeta::IsNaN( d = DatKnowns.Get( iOne, iTwo ) ) )
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
	float	m_dTotal;
	size_t	m_iTotal;
	float	m_dRatio;

	SSeed( size_t iOne, size_t iTwo, float dIn, size_t iIn, float dTotal, size_t iTotal,
		float dRatio ) : m_iOne(iOne), m_iTwo(iTwo), m_dIn(dIn), m_iIn(iIn), m_dTotal(dTotal),
		m_iTotal(iTotal), m_dRatio(dRatio) { }

	bool operator<( const SSeed& sSeed ) const {

		return ( m_dRatio < sSeed.m_dRatio ); }
};

int heavy( const gengetopt_args_info& sArgs, CDat& Dat, const CDat& DatKnowns,
	const vector<size_t>& veciKnowns ) {
	vector<bool>	vecfCluster;
	vector<float>	vecdConnectivity;
	vector<size_t>	veciConnectivity, veciCluster;
	size_t			i, j, iIn, iTotal, iMax, iMaxIn, iMaxTotal;
	size_t			iClusters;
	float			d, dMaxIn, dMaxTotal, dIn, dTotal, dRatio;

	// veciConnectivity[ i ] contains the number of edges out of gene i
	// vecdConnectivity[ i ] contains the sum of edge weights out of gene i
	veciConnectivity.resize( Dat.GetGenes( ) );
	vecdConnectivity.resize( Dat.GetGenes( ) );
	for( i = 0; i < veciConnectivity.size( ); ++i )
		for( j = ( i + 1 ); j < veciConnectivity.size( ); ++j )
			if( !CMeta::IsNaN( d = Dat.Get( i, j ) ) ) {
				veciConnectivity[ i ]++;
				veciConnectivity[ j ]++;
				vecdConnectivity[ i ] += d;
				vecdConnectivity[ j ] += d; }

	cout << sArgs.heavy_arg << endl;
	iClusters = 0;
	vecfCluster.resize( Dat.GetGenes( ) );
	while( ( sArgs.subgraphs_arg == -1 ) || ( iClusters < (size_t)sArgs.subgraphs_arg ) ) {
		priority_queue<SSeed>	pqueSeeds;

		veciCluster.resize( 1 );
		for( i = 0; i < Dat.GetGenes( ); ++i ) {
			veciCluster[ 0 ] = i;
			for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
				if( connectivity( j, veciCluster, vecdConnectivity, veciConnectivity, 0, 0,
					vecdConnectivity[ i ], veciConnectivity[ i ], Dat, dIn, iIn, dTotal, iTotal ) &&
					( ( dRatio = ( ( dIn / iIn ) / ( dTotal / iTotal ) ) ) >= sArgs.specificity_arg ) )
					pqueSeeds.push( SSeed( i, j, dIn, iIn, dTotal, iTotal, dRatio ) ); }
		if( pqueSeeds.empty( ) )
			break;
		cerr << "Seeds remaining: " << pqueSeeds.size( ) << endl;

		do {
			const SSeed&	sSeed	= pqueSeeds.top( );

			fill( vecfCluster.begin( ), vecfCluster.end( ), false );
			vecfCluster[ sSeed.m_iOne ] = vecfCluster[ sSeed.m_iTwo ] = true;
			veciCluster.resize( 2 );
			veciCluster[ 0 ] = sSeed.m_iOne;
			veciCluster[ 1 ] = sSeed.m_iTwo;

			cerr << "Cluster " << iClusters << " seed: " << Dat.GetGene( sSeed.m_iOne ) << ", " <<
				Dat.GetGene( sSeed.m_iTwo ) << ", " << sSeed.m_dRatio << endl;
			dIn = sSeed.m_dIn;
			iIn = sSeed.m_iIn;
			dTotal = sSeed.m_dTotal;
			iTotal = sSeed.m_iTotal;
			while( true ) {
				cerr << "Cluster " << iClusters << ", " << veciCluster.size( ) << " genes" << endl;
				max_connectivity( vecfCluster, veciCluster, vecdConnectivity, veciConnectivity, dIn,
					iIn, dTotal, iTotal, -1, Dat, dRatio, iMax, dMaxIn, iMaxIn, dMaxTotal, iMaxTotal );
				if( ( dRatio > 1 ) && ( dRatio >= ( sArgs.heavy_arg * sSeed.m_dRatio ) ) ) {
					vecfCluster[ iMax ] = true;
					veciCluster.push_back( iMax );
					dIn = dMaxIn;
					iIn = iMaxIn;
					dTotal = dMaxTotal;
					iTotal = iMaxTotal; }
				else
					break; }
			pqueSeeds.pop( ); }
		while( !pqueSeeds.empty( ) && ( veciCluster.size( ) < 3 ) );
		if( veciCluster.size( ) < 3 )
			break;

		cerr << "Found cluster: " << ( ( dIn / iIn ) / ( dTotal / iTotal ) ) << " (" << dIn << '/' <<
			iIn << ", " << dTotal << '/' << iTotal << ')' << endl;
		iClusters++;
		cout << ( ( dIn / iIn ) / ( dTotal / iTotal ) );
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
	size_t iIn, float dTotal, size_t iTotal, const CDat& Dat, float& dSumIn, size_t& iEdgesIn,
	float& dSumTotal, size_t& iEdgesTotal ) {
	size_t	i;
	float	d;
	bool	fRet;

	iEdgesTotal = iTotal + veciConnectivity[ iNew ];
	dSumTotal = dTotal + vecdConnectivity[ iNew ];
	iEdgesIn = iIn;
	dSumIn = dIn;
	fRet = false;
	for( i = 0; i < veciCluster.size( ); ++i )
		if( !CMeta::IsNaN( d = Dat.Get( iNew, veciCluster[ i ] ) ) ) {
			fRet = true;
			iEdgesTotal--;
			dSumTotal -= d;
			iEdgesIn++;
			dSumIn += d; }

	return fRet; }

void max_connectivity( const vector<bool>& vecfCluster, const vector<size_t>& veciCluster,
	const vector<float>& vecdConnectivity, const vector<size_t>& veciConnectivity, float dIn,
	size_t iIn, float dTotal, size_t iTotal, size_t iStart, const CDat& Dat, float& dMaxRatio,
	size_t& iMax, float& dMaxIn, size_t& iMaxIn, float& dMaxTotal, size_t& iMaxTotal ) {
	size_t	i, iEdgesIn, iEdgesTotal;
	float	dRatio, dSumIn, dSumTotal;

	dMaxRatio = -FLT_MAX;
	for( iMax = 0,i = ( iStart == -1 ) ? 0 : ( iStart + 1 ); i < vecfCluster.size( ); ++i ) {
		if( vecfCluster[ i ] || !connectivity( i, veciCluster, vecdConnectivity, veciConnectivity, dIn,
			iIn, dTotal, iTotal, Dat, dSumIn, iEdgesIn, dSumTotal, iEdgesTotal ) )
			continue;
		dRatio = ( ( dSumIn / iEdgesIn ) / ( dSumTotal / iEdgesTotal ) );
		if( dRatio > dMaxRatio ) {
			dMaxRatio = dRatio;
			iMax = i;
			dMaxIn = dSumIn;
			iMaxIn = iEdgesIn;
			dMaxTotal = dSumTotal;
			iMaxTotal = iEdgesTotal; } } }

int heavy2( const gengetopt_args_info& sArgs, CDat& Dat, const CDat& DatKnowns,
	const vector<size_t>& veciKnowns ) {
	size_t							i, j, iClusters, iMax;
	float							d, dCur, dMax;
	vector<size_t>					veciCluster, veciScores;
	vector<pair<size_t, size_t> >	vecpriiSeeds;
	vector<float>					vecdScores;
	vector<bool>					vecfCluster;

	for( i = 0; i < Dat.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = Dat.Get( i, j ) ) && ( d >= sArgs.specificity_arg ) )
				vecpriiSeeds.push_back( pair<size_t, size_t>( i, j ) );

	cout << sArgs.heavy_arg << endl;
	iClusters = 0;
	veciScores.resize( Dat.GetGenes( ) );
	vecdScores.resize( Dat.GetGenes( ) );
	vecfCluster.resize( Dat.GetGenes( ) );
	while( ( sArgs.subgraphs_arg == -1 ) || ( iClusters < (size_t)sArgs.subgraphs_arg ) ) {
		priority_queue<SSeed>	pqueSeeds;
		bool					fHit;

		for( i = 0; i < vecpriiSeeds.size( ); ++i )
			if( ( d = Dat.Get( vecpriiSeeds[ i ].first, vecpriiSeeds[ i ].second ) ) >= sArgs.specificity_arg )
					pqueSeeds.push( SSeed( vecpriiSeeds[ i ].first, vecpriiSeeds[ i ].second, d, 1, 0, 0, d ) );
		cerr << "Seeds remaining: " << pqueSeeds.size( ) << endl;

		for( fHit = false; !pqueSeeds.empty( ); ) {
			const SSeed&	sSeed	= pqueSeeds.top( );

			fill( vecfCluster.begin( ), vecfCluster.end( ), false );
			vecfCluster[ sSeed.m_iOne ] = vecfCluster [ sSeed.m_iTwo ] = true;
			veciCluster.resize( 2 );
			veciCluster[ 0 ] = sSeed.m_iOne;
			veciCluster[ 1 ] = sSeed.m_iTwo;

			cerr << "Cluster " << iClusters << " seed: " << Dat.GetGene( sSeed.m_iOne ) << ", " <<
				Dat.GetGene( sSeed.m_iTwo ) << ", " << sSeed.m_dRatio << endl;
			for( i = 0; i < Dat.GetGenes( ); ++i ) {
				vecdScores[ i ] = 0;
				for( veciScores[ i ] = j = 0; j < veciCluster.size( ); ++j )
					if( !CMeta::IsNaN( d = Dat.Get( i, veciCluster[ j ] ) ) ) {
						vecdScores[ i ] += d;
						veciScores[ i ]++; } }
			while( true ) {
				cerr << "Cluster " << iClusters << ", " << veciCluster.size( ) << " genes" << endl;
				for( dMax = 0,iMax = i = 0; i < Dat.GetGenes( ); ++i ) {
					if( vecfCluster[ i ] || !veciScores[ i ] )
						continue;
					if( ( dCur = ( vecdScores[ i ] / veciScores[ i ] ) ) > dMax ) {
						dMax = dCur;
						iMax = i; } }
				if( dMax < ( sArgs.heavy_arg * sSeed.m_dRatio ) )
					break;
				for( i = 0; i < Dat.GetGenes( ); ++i )
					if( !CMeta::IsNaN( d = Dat.Get( i, iMax ) ) ) {
						vecdScores[ i ] += d;
						veciScores[ i ]++; }
				vecfCluster[ iMax ] = true;
				veciCluster.push_back( iMax ); }
			if( veciCluster.size( ) < 3 ) {
				pqueSeeds.pop( );
				continue; }

			fHit = true;
			for( dMax = 0, iMax = i = 0; i < veciCluster.size( ); ++i )
				for( j = ( i + 1 ); j < veciCluster.size( ); ++j )
					if( !CMeta::IsNaN( d = Dat.Get( veciCluster[ i ], veciCluster[ j ] ) ) ) {
						iMax++;
						dMax += d; }

			cerr << "Found cluster: " << ( dMax /= iMax ) << endl;
			iClusters++;
			cout << dMax;
			for( i = 0; i < veciCluster.size( ); ++i )
				cout << '\t' << Dat.GetGene( veciCluster[ i ] );
			cout << endl;
			cout.flush( );

			for( i = 0; i < veciCluster.size( ); ++i )
				for( j = ( i + 1 ); j < veciCluster.size( ); ++j )
					if( !CMeta::IsNaN( d = Dat.Get( veciCluster[ i ], veciCluster[ j ] ) ) ) {
						dCur = min( dMax, d );
						Dat.Set( veciCluster[ i ], veciCluster[ j ], d - dCur ); }
			break; }
		if( !fHit )
			break; }

	return 0; }

int motifs( const gengetopt_args_info& sArgs, CDat& Dat ) {
	size_t			i, j, k, iOne, iTwo, iThree;
	vector<size_t>	veciOne, veciTwo;
	vector<bool>	vecfSign;
	float			dOne, dTwo, dThree;

	for( i = 0; i < Dat.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
			if( !CMeta::IsNaN( dOne = Dat.Get( i, j ) ) && ( fabs( dOne ) >= sArgs.motifs_arg ) ) {
				iOne = ( dOne >= 0 ) ? 1 : 0;
				for( k = ( j + 1 ); k < Dat.GetGenes( ); ++k )
					if( !CMeta::IsNaN( dTwo = Dat.Get( j, k ) ) && ( fabs( dTwo ) >= sArgs.motifs_arg ) &&
						!CMeta::IsNaN( dThree = Dat.Get( i, k ) ) && ( fabs( dThree ) >= sArgs.motifs_arg ) ) {
						iTwo = ( dTwo >= 0 ) ? 1 : 0;
						iThree = ( dThree >= 0 ) ? 1 : 0;
						if( ( iOne + iTwo + iThree ) == 1 ) {
							cout << Dat.GetGene( i ) << '\t' << dOne << '\t' <<
								Dat.GetGene( j ) << '\t' << dTwo << '\t' <<
								Dat.GetGene( k ) << '\t' << dThree << endl; } } }

	return 0; }

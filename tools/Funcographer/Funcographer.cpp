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

int output_cograph( const char*, const CDat&, const CDat&, const CPCL&, float );
int explore_graphs( const gengetopt_args_info&, const CDat&, const CDat&, const CPCL& );
bool connectivity( size_t, const vector<size_t>&, const vector<float>&, const vector<size_t>&,
	float, size_t, float, size_t, const CDat&, float&, size_t&, float&, size_t& );
void max_connectivity( const vector<bool>&, const vector<size_t>&, const vector<float>&,
	const vector<size_t>&, float, size_t, float, size_t, size_t, const CDat&, float&, size_t&, float&,
	size_t&, float&, size_t& );
void output_datasets( const gengetopt_args_info&, const CDat&, const CDat&, const CPCL&,
	const vector<size_t>&, const vector<size_t>& );

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CDat				DatFunctions, DatData;
	CPCL				PCLTrusts;
	int					iRet;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta = CMeta( sArgs.verbosity_arg );

	if( !DatFunctions.Open( sArgs.functions_arg ) ) {
		cerr << "Could not open: " << sArgs.functions_arg << endl;
		return 1; }
	if( !DatData.Open( sArgs.datasets_arg ) ) {
		cerr << "Could not open: " << sArgs.datasets_arg << endl;
		return 1; }
	if( !PCLTrusts.Open( sArgs.trusts_arg, sArgs.skip_arg ) ) {
		cerr << "Could not open: " << sArgs.trusts_arg << endl;
		return 1; }
	DatFunctions.Normalize( !sArgs.output_arg );
	DatData.Normalize( !sArgs.output_arg );
	PCLTrusts.Normalize( sArgs.output_arg ? CPCL::ENormalizeZScore : CPCL::ENormalizeMinMax );

	iRet = sArgs.output_arg ? output_cograph( sArgs.output_arg, DatData, DatFunctions, PCLTrusts,
			(float)sArgs.adjust_data_arg ) :
			explore_graphs( sArgs, DatData, DatFunctions, PCLTrusts );

	return iRet; }

int output_cograph( const char* szOutput, const CDat& DatData, const CDat& DatFunctions,
	const CPCL& PCLTrusts, float dAdjustData ) {
	vector<string>	vecstrNames;
	size_t			i, j;
	vector<size_t>	veciData, veciFunctions;
	CDat			Dat;

	vecstrNames.resize( DatFunctions.GetGenes( ) + DatData.GetGenes( ) );
	for( i = 0; i < DatFunctions.GetGenes( ); ++i )
		vecstrNames[ i ] = DatFunctions.GetGene( i );
	for( i = 0; i < DatData.GetGenes( ); ++i )
		vecstrNames[ DatFunctions.GetGenes( ) + i ] = DatData.GetGene( i );
	Dat.Open( vecstrNames );
	for( i = 0; i < DatFunctions.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < DatFunctions.GetGenes( ); ++j )
			Dat.Set( i, j, DatFunctions.Get( i, j ) );
	for( i = 0; i < DatData.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < DatData.GetGenes( ); ++j )
			Dat.Set( DatFunctions.GetGenes( ) + i, DatFunctions.GetGenes( ) + j, DatData.Get( i, j ) +
				dAdjustData );

	veciFunctions.resize( PCLTrusts.GetGenes( ) );
	veciData.resize( PCLTrusts.GetExperiments( ) );
	for( i = 0; i < veciFunctions.size( ); ++i ) {
		veciFunctions[ i ] = Dat.GetGene( PCLTrusts.GetGene( i ) );
		if( veciFunctions[ i ] == -1 ) {
			cerr << "Unknown function: " << PCLTrusts.GetGene( i ) << endl;
			return 1; } }
	for( i = 0; i < veciData.size( ); ++i ) {
		veciData[ i ] = Dat.GetGene( PCLTrusts.GetExperiment( i ) );
		if( veciData[ i ] == -1 )
			cerr << "Unknown dataset: " << PCLTrusts.GetExperiment( i ) << endl; }
	for( i = 0; i < veciFunctions.size( ); ++i ) {
		if( veciFunctions[ i ] == -1 )
			continue;
		for( j = 0; j < veciData.size( ); ++j )
			if( veciData[ j ] != -1 )
				Dat.Set( veciFunctions[ i ], veciData[ j ], PCLTrusts.Get( i, j ) ); }
	Dat.Save( szOutput );

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

int explore_graphs( const gengetopt_args_info& sArgs, const CDat& DatData, const CDat& DatFunctions,
	const CPCL& PCLTrusts ) {
	vector<bool>	vecfCluster;
	vector<float>	vecdConnectivity;
	vector<size_t>	veciConnectivity, veciCluster, veciFunctions2Trusts;
	size_t			i, j, iIn, iTotal, iMax, iMaxIn, iMaxTotal;
	size_t			iClusters;
	float			d, dMaxIn, dMaxTotal, dIn, dTotal, dRatio, dHeavy;
	CDat			Dat;

	Dat.Open( DatFunctions );
	dHeavy = (float)( sArgs.heavy_arg ? sArgs.heavy_arg : ( 1 / sArgs.specificity_arg ) );
	veciFunctions2Trusts.resize( DatFunctions.GetGenes( ) );
	for( i = 0; i < veciFunctions2Trusts.size( ); ++i )
		veciFunctions2Trusts[ i ] = PCLTrusts.GetGene( DatFunctions.GetGene( i ) );

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
			dTotal = sSeed.m_dTotal;
			iTotal = sSeed.m_iTotal;
			while( true ) {
				cerr << "Cluster " << iClusters << ", " << veciCluster.size( ) << " genes" << endl;
				max_connectivity( vecfCluster, veciCluster, vecdConnectivity, veciConnectivity, dIn,
					iIn, dTotal, iTotal, -1, Dat, dRatio, iMax, dMaxIn, iMaxIn, dMaxTotal, iMaxTotal );
				if( dRatio >= ( dHeavy * sSeed.m_dRatio ) ) {
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
		if( veciCluster.size( ) >= (size_t)sArgs.size_functions_arg ) {
			iClusters++;
			cout << ( ( dIn / iIn ) / ( dTotal / iTotal ) );
			for( i = 0; i < veciCluster.size( ); ++i )
				cout << '\t' << Dat.GetGene( veciCluster[ i ] );
			cout << endl;
			output_datasets( sArgs, DatData, DatFunctions, PCLTrusts, veciCluster, veciFunctions2Trusts );
			cout.flush( ); }

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

struct SColink {
	size_t	m_iDataset;
	float	m_dWeight;

	SColink( size_t iDataset, float dWeight ) : m_iDataset(iDataset), m_dWeight(dWeight) { }

	bool operator<( const SColink& sColink ) const {

		return ( m_dWeight < sColink.m_dWeight ); }
};

void output_datasets( const gengetopt_args_info& sArgs, const CDat& DatData, const CDat& DatFunctions,
	const CPCL& PCLTrusts, const vector<size_t>& veciCluster, const vector<size_t>& veciFunctions2Trusts ) {
	size_t			i, j, iFunction;
	float			d;
	vector<size_t>	veciData, veciTrusts2Data;

	{
		priority_queue<SColink>	pqueLinks;
		vector<float>			vecdData;

		vecdData.resize( PCLTrusts.GetExperiments( ) );
		for( i = 0; i < veciCluster.size( ); ++i ) {
			if( ( iFunction = veciFunctions2Trusts[ veciCluster[ i ] ] ) == -1 )
				continue;
			for( j = 0; j < PCLTrusts.GetExperiments( ); ++j )
				if( !CMeta::IsNaN( d = PCLTrusts.Get( iFunction, j ) ) )
					vecdData[ j ] += d; }

		for( i = 0; i < vecdData.size( ); ++i )
			pqueLinks.push( SColink( i, vecdData[ i ] ) );

		for( i = 0; ( i < (size_t)sArgs.size_datasets_arg ) && !pqueLinks.empty( ); ++i ) {
			cout << '\t' << pqueLinks.top( ).m_dWeight << '\t' <<
				PCLTrusts.GetExperiment( pqueLinks.top( ).m_iDataset ) << endl;
			veciData.push_back( pqueLinks.top( ).m_iDataset );
			pqueLinks.pop( ); }
	}

	for( i = 0; i < veciCluster.size( ); ++i ) {
		if( veciFunctions2Trusts[ veciCluster[ i ] ] == -1 )
			continue;
		for( j = ( i + 1 ); j < veciCluster.size( ); ++j )
			if( veciFunctions2Trusts[ veciCluster[ j ] ] != -1 )
				cout << "\t\t" << DatFunctions.GetGene( veciCluster[ i ] ) << '\t' <<
					DatFunctions.GetGene( veciCluster[ j ] ) << '\t' <<
					DatFunctions.Get( veciCluster[ i ], veciCluster[ j ] ) << endl; }
	for( i = 0; i < veciCluster.size( ); ++i ) {
		if( veciFunctions2Trusts[ veciCluster[ i ] ] == -1 )
			continue;
		for( j = 0; j < veciData.size( ); ++j )
			cout << "\t\t" << DatFunctions.GetGene( veciCluster[ i ] ) << '\t' <<
				PCLTrusts.GetExperiment( veciData[ j ] ) << '\t' <<
				PCLTrusts.Get( veciFunctions2Trusts[ veciCluster[ i ] ], veciData[ j ] ) << endl; }
	veciTrusts2Data.resize( PCLTrusts.GetExperiments( ) );
	for( i = 0; i < veciTrusts2Data.size( ); ++i )
		veciTrusts2Data[ i ] = DatData.GetGene( PCLTrusts.GetExperiment( i ) );
	for( i = 0; i < veciData.size( ); ++i )
		for( j = ( i + 1 ); j < veciData.size( ); ++j )
			cout << "\t\t" << PCLTrusts.GetExperiment( veciData[ i ] ) << '\t' <<
				PCLTrusts.GetExperiment( veciData[ j ] ) << '\t' <<
				DatData.Get( veciTrusts2Data[ veciData[ i ] ], veciTrusts2Data[ veciData[ j ] ] ) << endl; }

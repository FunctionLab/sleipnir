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

static const char	c_szDab[]	= ".dab";

struct SBackground {
	size_t					m_iCountOne;
	size_t					m_iCountTwo;
	size_t					m_iSizeOne;
	size_t					m_iSizeTwo;
	const vector<float>*	m_pvecdOut;
	float					m_dPrior;
	const CDat*				m_pDat;
	const CDat*				m_pDatWithin;
	vector<float>			m_vecdValues;
	const vector<size_t>*	m_pveciGenes;
};

static int MainSet( const gengetopt_args_info& );
static int MainBackground( const gengetopt_args_info& );
static void* BackgroundThread( void* );
static void* BackgroundThreadSingle( void* );

int main( int iArgs, char** aszArgs ) {
#ifdef WIN32
	pthread_win32_process_attach_np( );
#endif // WIN32
	gengetopt_args_info	sArgs;
	int					iRet;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta = CMeta( sArgs.verbosity_arg, sArgs.random_arg );

	iRet = ( sArgs.sizes_arg ? MainBackground : MainSet )( sArgs );
	pthread_exit( NULL );
#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return iRet; }

int MainSet( const gengetopt_args_info& sArgs ) {
	CGenome			Genome;
	CGenes			GenesQuery( Genome );
	ifstream		ifsm;
	size_t			iFunction, i, j;
	vector<size_t>	veciQuery;
	float			d;

	if( sArgs.genes_arg )
		ifsm.open( sArgs.genes_arg );
	if( !GenesQuery.Open( sArgs.genes_arg ? ifsm : cin ) ) {
		cerr << "Could not open: " << ( sArgs.genes_arg ? sArgs.genes_arg : "genes" ) << endl;
		return 1; }
	if( sArgs.genes_arg )
		ifsm.close( );

	veciQuery.resize( GenesQuery.GetGenes( ) );
	cout << "Function	Score	QueryIn	QueryOut	FuncIn	FuncOut" << endl;
	for( iFunction = 0; iFunction < sArgs.inputs_num; ++iFunction ) {
		CGenes			GenesFunction( Genome );
		CDat			Dat;
		string			strDab;
		float			dFuncIn, dFuncOut, dQueryIn, dQueryOut;
		size_t			iFuncIn, iFuncOut, iQueryIn, iQueryOut, iOne, iTwo;
		vector<size_t>	veciFunction;

		if( !GenesFunction.Open( sArgs.inputs[ iFunction ] ) ) {
			cerr << "Could not open: " << sArgs.inputs[ iFunction ] << endl;
			return 1; }
		strDab = (string)sArgs.directory_arg + '/' + CMeta::Deextension( CMeta::Basename(
			sArgs.inputs[ iFunction ] ) ) + c_szDab;
		if( !Dat.Open( strDab.c_str( ), sArgs.memmap_flag && !sArgs.normalize_flag ) ) {
			cerr << "Could not open: " << strDab << endl;
			return 1; }
		if( sArgs.normalize_flag )
			Dat.Normalize( CDat::ENormalizeSigmoid );

		for( i = 0; i < veciQuery.size( ); ++i )
			veciQuery[ i ] = Dat.GetGene( GenesQuery.GetGene( i ).GetName( ) );
		veciFunction.resize( GenesFunction.GetGenes( ) );
		for( i = 0; i < veciFunction.size( ); ++i )
			veciFunction[ i ] = Dat.GetGene( GenesFunction.GetGene( i ).GetName( ) );

		dFuncIn = dFuncOut = dQueryIn = dQueryOut = 0;
		for( iQueryIn = iFuncIn = iFuncOut = i = 0; i < veciFunction.size( ); ++i ) {
			if( ( iOne = veciFunction[ i ] ) == -1 )
				continue;
			for( j = ( i + 1 ); j < veciFunction.size( ); ++j )
				if( ( ( iTwo = veciFunction[ j ] ) != -1 ) &&
					!CMeta::IsNaN( d = Dat.Get( iOne, iTwo ) ) ) {
					iFuncIn++;
					dFuncIn += d; }
			for( j = 0; j < Dat.GetGenes( ); ++j )
				if( !CMeta::IsNaN( d = Dat.Get( iOne, j ) ) ) {
					iFuncOut++;
					dFuncOut += d; }
			for( j = 0; j < veciQuery.size( ); ++j )
				if( ( ( iTwo = veciQuery[ j ] ) != -1 ) &&
					!CMeta::IsNaN( d = Dat.Get( iOne, iTwo ) ) ) {
					iQueryIn++;
					dQueryIn += d; } }
		for( iQueryOut = i = 0; i < veciQuery.size( ); ++i ) {
			if( ( iOne = veciQuery[ i ] ) == -1 )
				continue;
			for( j = 0; j < Dat.GetGenes( ); ++j )
				if( !CMeta::IsNaN( d = Dat.Get( iOne, j ) ) ) {
					iQueryOut++;
					dQueryOut += d; } }
		dFuncIn /= iFuncIn;
		dFuncOut /= iFuncOut;
		dQueryIn /= iQueryIn;
		dQueryOut /= iQueryOut;

		cout << CMeta::Basename( sArgs.inputs[ iFunction ] ) << '\t' <<
			( dQueryIn / dQueryOut / dFuncIn * dFuncOut ) << '\t' << dQueryIn << '\t' << dQueryOut << '\t' <<
			dFuncIn << '\t' << dFuncOut << endl; }

	return 0; }

int MainBackground( const gengetopt_args_info& sArgs ) {
	CDat				Dat, DatWithin;
	const CDat*			pDatWithin;
	CDataMatrix			MatAves, MatStds, MatPValues;
	vector<size_t>		veciSizes, veciGenes;
	size_t				i, j, iIndexOne, iIndexTwo, iValue, iChunk;
	float				d, dPrior;
	vector<float>		vecdOut, vecdValues;
	vector<pthread_t>	vecpthdThreads;
	vector<SBackground>	vecsBackground;
	long double			dTmp;

	if( !Dat.Open( sArgs.input_arg, sArgs.memmap_flag && !sArgs.normalize_flag ) ) {
		cerr << "Could not open: " << ( sArgs.input_arg ? sArgs.input_arg : "stdin" ) << endl;
		return 1; }
	if( sArgs.normalize_flag )
		Dat.Normalize( CDat::ENormalizeSigmoid );

	{
		CPCL	PCLSizes( false );

		if( !PCLSizes.Open( sArgs.sizes_arg, 0 ) ) {
			cerr << "Could not open: " << sArgs.sizes_arg << endl;
			return 1; }
		veciSizes.resize( PCLSizes.GetGenes( ) );
		for( i = 0; i < veciSizes.size( ); ++i )
			veciSizes[ i ] = atoi( PCLSizes.GetGene( i ).c_str( ) );
	}

	if( sArgs.genes_arg ) {
		CGenome	Genome;
		CGenes	Genes( Genome );

		if( !Genes.Open( sArgs.genes_arg ) ) {
			cerr << "Could not open: " << sArgs.genes_arg << endl;
			return 1; }
		for( i = 0; i < Genes.GetGenes( ); ++i )
			if( ( j = Dat.GetGene( Genes.GetGene( i ).GetName( ) ) ) != -1 )
				veciGenes.push_back( j ); }
	else {
		veciGenes.resize( Dat.GetGenes( ) );
		for( i = 0; i < veciGenes.size( ); ++i )
			veciGenes[ i ] = i; }

	{
		vector<long double>	vecdTmp;
		vector<size_t>		veciTmp;

		vecdTmp.resize( Dat.GetGenes( ) );
		fill( vecdTmp.begin( ), vecdTmp.end( ), 0 );
		veciTmp.resize( Dat.GetGenes( ) );
		fill( veciTmp.begin( ), veciTmp.end( ), 0 );
		for( i = 0; i < Dat.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
				if( !CMeta::IsNaN( d = Dat.Get( i, j ) ) ) {
					veciTmp[ i ]++;
					veciTmp[ j ]++;
					vecdTmp[ i ] += d;
					vecdTmp[ j ] += d; }
		vecdOut.resize( vecdTmp.size( ) );
		for( i = 0; i < vecdOut.size( ); ++i )
			vecdOut[ i ] = (float)( vecdTmp[ i ] / ( veciTmp[ i ] ? veciTmp[ i ] : 1 ) );
	}

	if( sArgs.input_within_arg ) {
		if( !DatWithin.Open( sArgs.input_within_arg ) ) {
			cerr << "Could not open: " << sArgs.input_within_arg << endl;
			return 1; }
		pDatWithin = &DatWithin; }
	else
		pDatWithin = &Dat;

	for( dTmp = 0,iChunk = i = 0; i < pDatWithin->GetGenes( ); ++i )
		for( j = ( i + 1 ); j < pDatWithin->GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = pDatWithin->Get( i, j ) ) ) {
				iChunk++;
				dTmp += d; }
	dPrior = (float)( dTmp / iChunk );

	vecpthdThreads.resize( min( sArgs.threads_arg, sArgs.count_arg ) );
	vecsBackground.resize( vecpthdThreads.size( ) );
	iChunk = 1 + ( ( sArgs.count_arg - 1 ) / vecpthdThreads.size( ) );
	vecdValues.resize( ( sArgs.singles_flag ? 1 : sArgs.count_arg ) * iChunk * vecpthdThreads.size( ) );
	MatAves.Initialize( veciSizes.size( ), veciSizes.size( ) );
	MatAves.Clear( );
	MatStds.Initialize( veciSizes.size( ), veciSizes.size( ) );
	MatStds.Clear( );
	MatPValues.Initialize( veciSizes.size( ), veciSizes.size( ) );
	for( iIndexOne = 0; iIndexOne < veciSizes.size( ); ++iIndexOne )
		for( iIndexTwo = iIndexOne; iIndexTwo < ( sArgs.singles_flag ? ( iIndexOne + 1 ) : veciSizes.size( ) );
			++iIndexTwo ) {
			cerr << veciSizes[ iIndexOne ] << ':' << veciSizes[ iIndexTwo ] << endl;
			for( i = 0; i < vecpthdThreads.size( ); ++i ) {
				vecsBackground[ i ].m_dPrior = dPrior;
				vecsBackground[ i ].m_iCountOne = iChunk;
				vecsBackground[ i ].m_iCountTwo = sArgs.count_arg;
				vecsBackground[ i ].m_iSizeOne = veciSizes[ iIndexOne ];
				vecsBackground[ i ].m_iSizeTwo = veciSizes[ iIndexTwo ];
				vecsBackground[ i ].m_pDat = &Dat;
				vecsBackground[ i ].m_pDatWithin = pDatWithin;
				vecsBackground[ i ].m_pvecdOut = &vecdOut;
				vecsBackground[ i ].m_pveciGenes = &veciGenes;
				if( pthread_create( &vecpthdThreads[ i ], NULL, sArgs.singles_flag ? BackgroundThreadSingle:
					BackgroundThread, &vecsBackground[ i ] ) ) {
					cerr << "Couldn't create thread for group: " << i << endl;
					return 1; } }
			for( iValue = i = 0; i < vecpthdThreads.size( ); ++i ) {
				pthread_join( vecpthdThreads[ i ], NULL );
				for( j = 0; j < vecsBackground[ i ].m_vecdValues.size( ); ++j )
					vecdValues[ iValue++ ] = vecsBackground[ i ].m_vecdValues[ j ]; }
			MatAves.Set( iIndexOne, iIndexTwo, (float)CStatistics::Average( vecdValues ) );
			if( sArgs.invgauss_flag ) {
				for( i = 0; i < vecdValues.size( ); ++i )
					MatStds.Get( iIndexOne, iIndexTwo ) += 1 / vecdValues[ i ];
				MatStds.Set( iIndexOne, iIndexTwo, vecdValues.size( ) / ( MatStds.Get( iIndexOne, iIndexTwo ) -
					( vecdValues.size( ) / MatAves.Get( iIndexOne, iIndexTwo ) ) ) ); }
			else
				MatStds.Set( iIndexOne, iIndexTwo, (float)sqrt( CStatistics::Variance( vecdValues,
					MatAves.Get( iIndexOne, iIndexTwo ) ) ) );
			MatPValues.Set( iIndexOne, iIndexTwo, (float)CStatistics::Percentile( vecdValues.begin( ),
				vecdValues.end( ), sArgs.percentile_arg ) ); }

	if( sArgs.singles_flag )
		for( i = 0; i < veciSizes.size( ); ++i )
			cout << veciSizes[ i ] << '\t' << MatAves.Get( i, i ) << '|' << MatStds.Get( i, i ) << '|' <<
				MatPValues.Get( i, i ) << endl;
	else {
		for( i = 0; i < veciSizes.size( ); ++i )
			cout << '\t' << veciSizes[ i ];
		cout << endl;
		for( i = 0; i < MatAves.GetRows( ); ++i ) {
			cout << veciSizes[ i ];
			for( j = 0; j < i; ++j )
				cout << '\t';
			for( ; j < MatAves.GetColumns( ); ++j ) {
				cout << '\t' << MatAves.Get( i, j ) << '|' << MatStds.Get( i, j ) << '|' <<
					MatPValues.Get( i, j ); }
			cout << endl; } }

	return 0; }

void* BackgroundThread( void* pData ) {
	size_t				i, j, iCountOne, iCountTwo, iEdgesOne, iEdgesTwo, iValue;
	SBackground*		psData;
	vector<size_t>		veciOne, veciTwo;
	long double			dBackgroundOne, dBackgroundTwo, dWithinOne, dWithinTwo, dWithin, dBackground, dBetween;
	vector<long double>	vecdBetweenOne, vecdBetweenTwo;
	float				d;

	psData = (SBackground*)pData;
	psData->m_vecdValues.resize( psData->m_iCountOne * psData->m_iCountTwo );
	veciOne.resize( psData->m_iSizeOne );
	veciTwo.resize( psData->m_iSizeTwo );
	iEdgesOne = veciOne.size( ) * ( veciOne.size( ) + 1 ) / 2;
	iEdgesTwo = veciTwo.size( ) * ( veciTwo.size( ) + 1 ) / 2;
	vecdBetweenOne.resize( veciOne.size( ) );
	vecdBetweenTwo.resize( veciTwo.size( ) );
	for( iValue = iCountOne = 0; iCountOne < psData->m_iCountOne; ++iCountOne ) {
		const vector<size_t>&	veciGenes	= *psData->m_pveciGenes;
		const vector<float>&	vecdOut		= *psData->m_pvecdOut;
		set<size_t>				setiOne;

		while( setiOne.size( ) < veciOne.size( ) )
			setiOne.insert( veciGenes[ rand( ) % veciGenes.size( ) ] );
		copy( setiOne.begin( ), setiOne.end( ), veciOne.begin( ) );

		dBackgroundOne = dWithinOne = 0;
		for( i = 0; i < veciOne.size( ); ++i ) {
			dBackgroundOne += vecdOut[ veciOne[ i ] ];
			dWithinOne += psData->m_dPrior;
			for( j = ( i + 1 ); j < veciOne.size( ); ++j )
				dWithinOne += psData->m_pDatWithin->Get( veciOne[ i ], veciOne[ j ] ); }

		for( iCountTwo = 0; iCountTwo < psData->m_iCountTwo; ++iCountTwo ) {
			set<size_t>				setiTwo;

			while( setiTwo.size( ) < veciTwo.size( ) )
				setiTwo.insert( veciGenes[ rand( ) % veciGenes.size( ) ] );
			copy( setiTwo.begin( ), setiTwo.end( ), veciTwo.begin( ) );

			dBackgroundTwo = dWithinTwo = 0;
			for( i = 0; i < veciTwo.size( ); ++i ) {
				dBackgroundTwo += (*psData->m_pvecdOut)[ veciTwo[ i ] ];
				dWithinTwo += psData->m_dPrior;
				for( j = ( i + 1 ); j < veciTwo.size( ); ++j )
					dWithinTwo += psData->m_pDatWithin->Get( veciTwo[ i ], veciTwo[ j ] ); }

// Alternative: calculate an average over all edges (weighted) rather than between sets (unweighted)
//			dWithin = ( dWithinOne + dWithinTwo ) / ( iEdgesOne + iEdgesTwo );
			dWithin = ( ( dWithinOne / iEdgesOne ) + ( dWithinTwo / iEdgesTwo ) ) / 2;
			dBackground = ( dBackgroundOne + dBackgroundTwo ) / ( veciOne.size( ) + veciTwo.size( ) );
// Alternative: calculate between scores by edge rather than node
/*
			dBetween = 0;
			for( i = 0; i < veciOne.size( ); ++i )
				for( j = 0; j < veciTwo.size( ); ++j )
					dBetween += ( veciOne[ i ] == veciTwo[ j ] ) ? 1 :
						psData->m_pDat->Get( veciOne[ i ], veciTwo[ j ] );
			dBetween /= veciOne.size( ) * veciTwo.size( );
*/
			fill( vecdBetweenOne.begin( ), vecdBetweenOne.end( ), 0 );
			fill( vecdBetweenTwo.begin( ), vecdBetweenTwo.end( ), 0 );
			for( i = 0; i < veciOne.size( ); ++i )
				for( j = 0; j < veciTwo.size( ); ++j ) {
					d = ( veciOne[ i ] == veciTwo[ j ] ) ? 1 :
						psData->m_pDat->Get( veciOne[ i ], veciTwo[ j ] );
					vecdBetweenOne[ i ] += d;
					vecdBetweenTwo[ j ] += d; }
			for( dBetween = 0,i = 0; i < vecdBetweenOne.size( ); ++i )
				dBetween += vecdBetweenOne[ i ] / veciTwo.size( );
			for( i = 0; i < vecdBetweenTwo.size( ); ++i )
				dBetween += vecdBetweenTwo[ i ] / veciOne.size( );
			dBetween /= veciOne.size( ) + veciTwo.size( );

			psData->m_vecdValues[ iValue++ ] = (float)( psData->m_dPrior * dBetween / dWithin /
				dBackground ); } }

	return NULL; }

void* BackgroundThreadSingle( void* pData ) {
	size_t			i, j, iCount, iEdges;
	SBackground*	psData;
	vector<size_t>	veciOne;
	long double		dWithin, dBackground;
	float			d;

	psData = (SBackground*)pData;
	psData->m_vecdValues.resize( psData->m_iCountOne );
	veciOne.resize( psData->m_iSizeOne );
	iEdges = veciOne.size( ) * ( veciOne.size( ) + 1 ) / 2;
	for( iCount = 0; iCount < psData->m_iCountOne; ++iCount ) {
		const vector<size_t>&	veciGenes	= *psData->m_pveciGenes;
		const vector<float>&	vecdOut		= *psData->m_pvecdOut;
		set<size_t>				setiOne;

		while( setiOne.size( ) < veciOne.size( ) )
			setiOne.insert( veciGenes[ rand( ) % veciGenes.size( ) ] );
		copy( setiOne.begin( ), setiOne.end( ), veciOne.begin( ) );

		dBackground = dWithin = 0;
		for( i = 0; i < veciOne.size( ); ++i ) {
			dBackground += vecdOut[ veciOne[ i ] ];
			dWithin += psData->m_dPrior;
			for( j = ( i + 1 ); j < veciOne.size( ); ++j )
				if( !CMeta::IsNaN( d = psData->m_pDat->Get( veciOne[ i ], veciOne[ j ] ) ) )
					dWithin += d; }
		dBackground /= veciOne.size( );
		dWithin /= iEdges;

		psData->m_vecdValues[ iCount ] = (float)( dWithin / dBackground ); }

	return NULL; }

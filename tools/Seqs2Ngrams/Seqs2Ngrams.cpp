#include "stdafx.h"
#include "cmdline.h"

static const size_t	c_iBuffer	= 1048576;
static const char	c_acBases[]	= {'T', 'A', 'C', 'G', 'N'};

static size_t ngram( const string&, size_t, size_t );

static inline size_t ngram( char c ) {
	size_t	i;

	for( i = 0; i < ARRAYSIZE(c_acBases); ++i )
		if( c == c_acBases[ i ] )
			return i;

	cerr << c << endl;
	exit( 1 ); }

static inline size_t ipow( size_t iBase, size_t iExp ) {
	size_t	i, iRet;

	for( iRet = 1,i = 0; i < iExp; ++i )
		iRet *= iBase;

	return iRet; }

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info		sArgs;
	char*					acLine;
	ifstream				ifsm;
	CDat					Dat;
	size_t					i, j, k, m, iGene, iNgrams, iNgram, iPlaces;
	vector<vector<string> >	vecvecstrSeqs;
	vector<size_t>			veciCounts;
	char*					szSeq;
	float*					adOne;
	float*					adTwo;
	CGenome					Genome;
	CGenes					GenesIn( Genome );
	CMeasureEuclidean		Euclidean;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	iPlaces = ipow( ARRAYSIZE(c_acBases), sArgs.n_arg - 1 );
	iNgrams = iPlaces * ARRAYSIZE(c_acBases);

	CMeasureSigmoid			EuclideanSig( &Euclidean, false, 1.0f / iNgrams );

	if( sArgs.genes_arg ) {
		ifsm.open( sArgs.genes_arg );
		if( !GenesIn.Open( ifsm ) ) {
			cerr << "Could not open: " << sArgs.genes_arg << endl;
			return 1; }
		ifsm.close( ); }

	if( sArgs.input_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.input_arg ); }
	{
		istream&							istm	= sArgs.input_arg ? (istream&)ifsm : cin;
		map<string,size_t>					mapGenes;
		map<string,size_t>::const_iterator	iterGene;
		vector<string>						vecstrGenes;

		acLine = new char[ c_iBuffer ];
		while( !istm.eof( ) ) {
			istm.getline( acLine, c_iBuffer - 1 );
			if( !( szSeq = strchr( acLine, '\t' ) ) )
				continue;
			*szSeq = 0;

			if( GenesIn.GetGenes( ) && !GenesIn.IsGene( acLine ) )
				continue;

			if( ( iterGene = mapGenes.find( acLine ) ) == mapGenes.end( ) ) {
				iGene = mapGenes.size( );
				mapGenes[ acLine ] = iGene;
				vecstrGenes.push_back( acLine );
				veciCounts.push_back( 1 );
				vecvecstrSeqs.push_back( vector<string>( ) ); }
			else
				veciCounts[ iGene = iterGene->second ]++;
			vecvecstrSeqs[ iGene ].push_back( ++szSeq ); }
		delete[] acLine;
		Dat.Open( vecstrGenes );
	}
	if( sArgs.input_arg )
		ifsm.close( );

	adOne = new float[ iNgrams ];
	adTwo = new float[ iNgrams ];
	for( i = 0; i < Dat.GetGenes( ); ++i ) {
		cerr << "Gene " << i << '/' << Dat.GetGenes( ) << endl;

		memset( adOne, 0, iNgrams * sizeof(*adOne) );
		for( j = 0; j < vecvecstrSeqs[ i ].size( ); ++j ) {
			const string&	strSeq	= vecvecstrSeqs[ i ][ j ];

			adOne[ iNgram = ngram( strSeq, 0, sArgs.n_arg ) ]++;
			for( k = 1; k < ( strSeq.length( ) - sArgs.n_arg ); ++k ) {
				iNgram = ( ( iNgram - ( ngram( strSeq[ k - 1 ] ) * iPlaces ) ) * ARRAYSIZE(c_acBases) ) +
					ngram( strSeq[ k + sArgs.n_arg - 1 ] );
				adOne[ iNgram ]++; } }
		for( j = 0; j < iNgrams; ++j )
			if( adOne[ j ] )
				adOne[ j ] /= veciCounts[ i ];

		for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j ) {
			memset( adTwo, 0, iNgrams * sizeof(*adTwo) );
			for( k = 0; k < vecvecstrSeqs[ j ].size( ); ++k ) {
				const string&	strSeq	= vecvecstrSeqs[ j ][ k ];

				adTwo[ iNgram = ngram( strSeq, 0, sArgs.n_arg ) ]++;
				for( m = 1; m < ( strSeq.length( ) - sArgs.n_arg ); ++m ) {
					iNgram = ( ( iNgram - ( ngram( strSeq[ m - 1 ] ) * iPlaces ) ) * ARRAYSIZE(c_acBases) ) +
						ngram( strSeq[ m + sArgs.n_arg - 1 ] );
					adTwo[ iNgram ]++; } }
			for( k = 0; k < iNgrams; ++k )
				if( adTwo[ k ] )
					adTwo[ k ] /= veciCounts[ j ];

			Dat.Set( i, j, (float)EuclideanSig.Measure( adOne, iNgrams, adTwo, iNgrams, IMeasure::EMapNone, NULL, NULL ) );
		} }
	delete[] adOne;
	delete[] adTwo;

	if( sArgs.normalize_flag || sArgs.zscore_flag )
		Dat.Normalize( !!sArgs.normalize_flag );

	Dat.Save( sArgs.output_arg );

	return 0; }

size_t ngram( const string& strSeq, size_t iOff, size_t iLen ) {
	size_t	i, iRet;

	for( iRet = i = 0; i < iLen; ++i )
		iRet = ( iRet * ARRAYSIZE(c_acBases) ) + ngram( strSeq[ iOff + i ] );

	return iRet; }

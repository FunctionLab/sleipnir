#include "stdafx.h"

enum ETFPN {
	ETFPN_TP	= 0,
	ETFPN_FP	= ETFPN_TP + 1,
	ETFPN_TN	= ETFPN_FP + 1,
	ETFPN_FN	= ETFPN_TN + 1
};

void Balance( CDat&, const CDat&, const vector<size_t>& );
int MainPR( const CDat&, const CDat&, bool, double, double, double );
int MainLLS( const CDat&, const CDat&, bool, double, double, double );

int main( int iArgs, char** aszArgs ) {
	const char*	szAnswers;
	const char*	szData;
	CDat		DatAnswers, DatData;
	double		dMin, dDelta, dMax;
	bool		fInvert;

	srand( 0 );
	if( iArgs < 3 ) {
		cerr << "Usage: " << aszArgs[ 0 ] << " <answers.txt> <data.dat> [invert] [min] [delta] [max]" <<
			endl;
		return 1; }
	szAnswers = aszArgs[ 1 ];
	szData = aszArgs[ 2 ];
	fInvert = ( iArgs > 3 ) ? !!atoi( aszArgs[ 3 ] ) : false;
	dMin = ( iArgs > 4 ) ? atof( aszArgs[ 4 ] ) : 0;
	dDelta = ( iArgs > 5 ) ? atof( aszArgs[ 5 ] ) : 0.01;
	dMax = ( iArgs > 6 ) ? atof( aszArgs[ 6 ] ) : 1;

	if( !DatAnswers.Open( szAnswers ) ) {
		cerr << "Couldn't open: " << szAnswers << endl;
		return 1; }
	if( !DatData.Open( szData ) ) {
		cerr << "Couldn't open: " << szData << endl;
		return 1; }

	return MainLLS( DatAnswers, DatData, fInvert, dMin, dDelta, dMax ); }

int MainPR( const CDat& DatAnswers, const CDat& DatData, bool fInvert, double dMin,
	double dDelta, double dMax ) {
	size_t					i, j, iOne, iTwo;
	vector<size_t>			veciGenes;
	vector<vector<size_t> >	vecveciResults;
	ETFPN					eTFPN;
	int						k, iMax;
	float					dAnswer, dValue;

	veciGenes.resize( DatAnswers.GetGenes( ) );
	for( i = 0; i < DatAnswers.GetGenes( ); ++i )
		veciGenes[ i ] = DatData.GetGene( DatAnswers.GetGene( i ) );
//	Balance( DatAnswers, DatData, veciGenes );

	vecveciResults.resize( (size_t)( ( dMax - dMin ) / dDelta ) );
	for( i = 0; i < vecveciResults.size( ); ++i )
		vecveciResults[ i ].resize( 4 );
	for( i = 0; i < DatAnswers.GetGenes( ); ++i ) {
		if( ( iOne = veciGenes[ i ] ) == -1 )
			continue;
		for( j = ( i + 1 ); j < DatAnswers.GetGenes( ); ++j ) {
			if( ( ( iTwo = veciGenes[ j ] ) == -1 ) ||
				CMeta::IsNaN( dAnswer = DatAnswers.Get( i, j ) ) ||
				CMeta::IsNaN( dValue = DatData.Get( iOne, iTwo ) ) )
				continue;
			if( fInvert )
				dValue = 1 - dValue;

			iMax = (int)( ( dValue - dMin ) / dDelta );
			if( iMax > (int)vecveciResults.size( ) )
				iMax = (int)vecveciResults.size( );
			eTFPN = (ETFPN)!dAnswer;
			for( k = 0; k < iMax; ++k )
				vecveciResults[ k ][ eTFPN ]++;
			eTFPN = (ETFPN)( 2 + !eTFPN );
			for( ; k < (int)vecveciResults.size( ); ++k )
				vecveciResults[ k ][ eTFPN ]++; } }

	cout << "Cut\tTP\tFP\tTN\tFN" << endl;
	for( i = 0; i < vecveciResults.size( ); ++i ) {
		cout << ( dMin + ( i * dDelta ) );
		for( j = 0; j < vecveciResults[ i ].size( ); ++j )
			cout << '\t' << (unsigned int)vecveciResults[ i ][ j ];
		cout << endl; }
	cout.flush( );

	return 0; }

void Balance( CDat& DatAnswers, const CDat& DatData, const vector<size_t>& veciGenes ) {
	size_t	i, j, iOne, iTwo;
	size_t	aiPN[] = { 0, 0 };
	float	dFrac, dAnswer;
	bool	fP;

	for( i = 0; i < veciGenes.size( ); ++i ) {
		if( ( iOne = veciGenes[ i ] ) == -1 )
			continue;
		for( j = ( i + 1 ); j < veciGenes.size( ); ++j ) {
			if( ( ( iTwo = veciGenes[ j ] ) == -1 ) ||
				CMeta::IsNaN( DatData.Get( iOne, iTwo ) ) ||
				CMeta::IsNaN( dAnswer = DatAnswers.Get( i, j ) ) )
				continue;
			aiPN[ (size_t)dAnswer ]++; } }

	fP = aiPN[ 1 ] > aiPN[ 0 ];
	dFrac = (float)aiPN[ !fP ] / aiPN[ fP ];
	for( i = 0; i < DatAnswers.GetGenes( ); ++i ) {
		if( ( iOne = veciGenes[ i ] ) == -1 )
			continue;
		for( j = ( i + 1 ); j < DatAnswers.GetGenes( ); ++j )
			if( !( ( iTwo = veciGenes[ j ] ) == -1 ) &&
				!CMeta::IsNaN( dAnswer = DatAnswers.Get( i, j ) ) && ( !!dAnswer == fP ) &&
				!CMeta::IsNaN( DatData.Get( iOne, iTwo ) ) &&
				( ( (float)rand( ) / RAND_MAX ) > dFrac ) )
				DatAnswers.Set( i, j, CMeta::GetNaN( ) ); } }

int MainLLS( const CDat& DatAnswers, const CDat& DatData, bool fInvert, double dMin,
	double dDelta, double dMax ) {
	size_t					i, j, iOne, iTwo;
	vector<size_t>			veciGenes, veciRec;
	vector<vector<bool> >	vecvecfGenes;
	vector<vector<size_t> >	vecveciResults;
	ETFPN					eTFPN;
	int						k, iMax;
	float					dAnswer, dValue;

	vecveciResults.resize( (size_t)( ( dMax - dMin ) / dDelta ) );
	veciGenes.resize( DatAnswers.GetGenes( ) );
	vecvecfGenes.resize( veciGenes.size( ) );
	for( i = 0; i < DatAnswers.GetGenes( ); ++i ) {
		veciGenes[ i ] = DatData.GetGene( DatAnswers.GetGene( i ) );
		vecvecfGenes[ i ].resize( vecveciResults.size( ) ); }

	for( i = 0; i < vecveciResults.size( ); ++i )
		vecveciResults[ i ].resize( 4 );
	for( i = 0; i < DatAnswers.GetGenes( ); ++i ) {
		if( ( iOne = veciGenes[ i ] ) == -1 )
			continue;
		for( j = ( i + 1 ); j < DatAnswers.GetGenes( ); ++j ) {
			if( ( ( iTwo = veciGenes[ j ] ) == -1 ) ||
				CMeta::IsNaN( dAnswer = DatAnswers.Get( i, j ) ) ||
				CMeta::IsNaN( dValue = DatData.Get( iOne, iTwo ) ) )
				continue;
			if( fInvert )
				dValue = 1 - dValue;

			iMax = (int)( ( dValue - dMin ) / dDelta );
			if( iMax > (int)vecveciResults.size( ) )
				iMax = (int)vecveciResults.size( );
			eTFPN = (ETFPN)!dAnswer;
			for( k = 0; k < iMax; ++k ) {
				vecveciResults[ k ][ eTFPN ]++;
				vecvecfGenes[ i ][ k ] = vecvecfGenes[ j ][ k ] = true; }
			eTFPN = (ETFPN)( 2 + !eTFPN );
			for( ; k < (int)vecveciResults.size( ); ++k )
				vecveciResults[ k ][ eTFPN ]++; } }

	veciRec.resize( vecveciResults.size( ) );
	for( i = 0; i < veciRec.size( ); ++i ) {
		veciRec[ i ] = 0;
		for( j = 0; j < vecvecfGenes.size( ); ++j )
			if( vecvecfGenes[ j ][ i ] )
				veciRec[ i ]++; }

	cout << "Cut\tGenes\tTP\tFP\tTN\tFN" << endl;
	for( i = 0; i < vecveciResults.size( ); ++i ) {
		cout << ( dMin + ( i * dDelta ) ) << '\t' << (unsigned int)veciRec[ i ];
		for( j = 0; j < vecveciResults[ i ].size( ); ++j )
			cout << '\t' << (unsigned int)vecveciResults[ i ][ j ];
		cout << endl; }
	cout.flush( );

	return 0; }

#include "stdafx.h"
#include "cmdline.h"

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
	CDat					DatAnswers, DatData;
	gengetopt_args_info		sArgs;
	size_t					i, j, iOne, iTwo, iGenes;
	vector<size_t>			veciGenes, veciRec;
	CFullMatrix<bool>		MatGenes;
	CFullMatrix<size_t>		MatResults;
	ETFPN					eTFPN;
	int						k, iMax;
	float					dAnswer, dValue;
	CGenome					Genome;
	CGenes					GenesIn( Genome );
	ifstream				ifsm;
	vector<bool>			vecfGenes;
	ofstream				ofsm;
	ostream*				postm;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( !DatAnswers.Open( sArgs.answers_arg ) ) {
		cerr << "Couldn't open: " << sArgs.answers_arg << endl;
		return 1; }
	if( sArgs.genes_arg ) {
		ifsm.open( sArgs.genes_arg );
		if( !GenesIn.Open( ifsm ) ) {
			cerr << "Couldn't open: " << sArgs.genes_arg << endl;
			return 1; }
		ifsm.close( );

		vecfGenes.resize( DatAnswers.GetGenes( ) );
		for( i = 0; i < vecfGenes.size( ); ++i )
			vecfGenes[ i ] = GenesIn.IsGene( DatAnswers.GetGene( i ) ); }
	if( !DatData.Open( sArgs.input_arg ) ) {
		cerr << "Couldn't open: " << sArgs.input_arg << endl;
		return 1; }

	veciGenes.resize( DatAnswers.GetGenes( ) );
	MatResults.Initialize( (size_t)( ( sArgs.max_arg - sArgs.min_arg ) / sArgs.delta_arg ),
		4 );
	MatGenes.Initialize( veciGenes.size( ), MatResults.GetRows( ) );
	for( i = 0; i < DatAnswers.GetGenes( ); ++i )
		veciGenes[ i ] = DatData.GetGene( DatAnswers.GetGene( i ) );

	for( iGenes = 0; !sArgs.inputs_num || ( iGenes < sArgs.inputs_num ); ++iGenes ) {
		MatResults.Clear( );
		MatGenes.Clear( );

		if( sArgs.inputs_num ) {
			ifsm.clear( );
			ifsm.open( sArgs.inputs[ iGenes ] );
			if( !GenesIn.Open( ifsm ) ) {
				cerr << "Couldn't open: " << sArgs.inputs[ iGenes ] << endl;
				return 1; }
			vecfGenes.resize( DatAnswers.GetGenes( ) );
			for( i = 0; i < vecfGenes.size( ); ++i )
				vecfGenes[ i ] = GenesIn.IsGene( DatAnswers.GetGene( i ) );
			cerr << "Processing " << sArgs.inputs[ iGenes ] << "..." << endl;
			ifsm.close( ); }

		for( i = 0; i < DatAnswers.GetGenes( ); ++i ) {
			if( ( iOne = veciGenes[ i ] ) == -1 )
				continue;
			for( j = ( i + 1 ); j < DatAnswers.GetGenes( ); ++j ) {
				if( ( ( iTwo = veciGenes[ j ] ) == -1 ) ||
					CMeta::IsNaN( dValue = DatData.Get( iOne, iTwo ) ) ||
					CMeta::IsNaN( dAnswer = DatAnswers.Get( i, j ) ) )
					continue;
				if( GenesIn.GetGenes( ) &&
					!( vecfGenes[ i ] || vecfGenes[ j ] ) )
//					( ( dAnswer && !( vecfGenes[ i ] && vecfGenes[ j ] ) ) ||
//					( !dAnswer && !( vecfGenes[ i ] || vecfGenes[ j ] ) ) ) )
					continue;
				if( sArgs.invert_flag )
					dValue = 1 - dValue;

				iMax = (int)( ( dValue - sArgs.min_arg ) / sArgs.delta_arg );
				if( iMax > (int)MatResults.GetRows( ) )
					iMax = (int)MatResults.GetRows( );
				eTFPN = (ETFPN)!dAnswer;
				for( k = 0; k < iMax; ++k ) {
					MatResults.Get( k, eTFPN )++;
					MatGenes.Set( i, k, true );
					MatGenes.Set( j, k, true ); }
				eTFPN = (ETFPN)( 2 + !eTFPN );
				for( ; k < (int)MatResults.GetRows( ); ++k )
					MatResults.Get( k, eTFPN )++; } }

		veciRec.resize( MatResults.GetRows( ) );
		for( i = 0; i < veciRec.size( ); ++i ) {
			veciRec[ i ] = 0;
			for( j = 0; j < MatGenes.GetRows( ); ++j )
				if( MatGenes.Get( j, i ) )
					veciRec[ i ]++; }

		if( sArgs.inputs_num ) {
			ofsm.open( ( (string)sArgs.directory_arg + '\\' +
				CMeta::Basename( sArgs.inputs[ iGenes ] ) + ".bins" ).c_str( ) );
			postm = &ofsm; }
		else
			postm = &cout;

		*postm << "Cut\tGenes\tTP\tFP\tTN\tFN" << endl;
		for( i = 0; i < MatResults.GetRows( ); ++i ) {
			*postm << ( sArgs.min_arg + ( i * sArgs.delta_arg ) ) << '\t' <<
				(unsigned int)veciRec[ i ];
			for( j = 0; j < MatResults.GetColumns( ); ++j )
				*postm << '\t' << (unsigned int)MatResults.Get( i, j );
			*postm << endl; }

		if( sArgs.inputs_num )
			ofsm.close( );
		else
			cout.flush( );

		if( !sArgs.inputs_num )
			break; }

	CMeta::Shutdown( );
	return 0; }

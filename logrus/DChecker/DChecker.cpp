#include "stdafx.h"
#include "cmdline.h"

enum ETFPN {
	ETFPN_TP	= 0,
	ETFPN_FP	= ETFPN_TP + 1,
	ETFPN_TN	= ETFPN_FP + 1,
	ETFPN_FN	= ETFPN_TN + 1
};

int main( int iArgs, char** aszArgs ) {
	CDat				Answers, Data;
	gengetopt_args_info	sArgs;
	size_t				i, j, iOne, iTwo, iGenes;
	vector<size_t>		veciGenes, veciRec;
	CFullMatrix<bool>	MatGenes;
	CFullMatrix<size_t>	MatResults;
	ETFPN				eTFPN;
	int					k, iMax;
	float				dAnswer, dValue;
	vector<bool>		vecfGenes;
	CBinaryMatrix		MatPairs;
	ofstream			ofsm;
	ostream*			postm;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( !Answers.Open( sArgs.answers_arg ) ) {
		cerr << "Couldn't open: " << sArgs.answers_arg << endl;
		return 1; }
	if( sArgs.genes_arg && !Answers.FilterGenes( sArgs.genes_arg, CDat::EFilterInclude ) ) {
		cerr << "Couldn't open: " << sArgs.genes_arg << endl;
		return 1; }
	if( sArgs.genet_arg && !Answers.FilterGenes( sArgs.genet_arg, CDat::EFilterTerm ) ) {
		cerr << "Couldn't open: " << sArgs.genet_arg << endl;
		return 1; }
	if( sArgs.genex_arg && !Answers.FilterGenes( sArgs.genex_arg, CDat::EFilterExclude ) ) {
		cerr << "Couldn't open: " << sArgs.genex_arg << endl;
		return 1; }
	if( !Data.Open( sArgs.input_arg ) ) {
		cerr << "Couldn't open: " << sArgs.input_arg << endl;
		return 1; }

	veciGenes.resize( Answers.GetGenes( ) );
	MatResults.Initialize( (size_t)( ( sArgs.max_arg - sArgs.min_arg ) / sArgs.delta_arg ) + 1,
		4 );
	MatGenes.Initialize( veciGenes.size( ), MatResults.GetRows( ) );
	for( i = 0; i < Answers.GetGenes( ); ++i )
		veciGenes[ i ] = Data.GetGene( Answers.GetGene( i ) );

	if( sArgs.inputs_num ) {
		MatPairs.Initialize( Answers.GetGenes( ) );
		for( i = 0; i < MatPairs.GetSize( ); ++i )
			for( j = ( i + 1 ); j < MatPairs.GetSize( ); ++j )
				MatPairs.Set( i, j, true );
		for( i = 0; i < sArgs.inputs_num; ++i ) {
			ifstream	ifsm( sArgs.inputs[ i ] );
			CGenome		Genome;
			CGenes		Genes( Genome );
			size_t		iGene;

			if( !Genes.Open( ifsm ) ) {
				cerr << "Couldn't open: " << sArgs.inputs[ i ] << endl;
				return 1; }
			for( j = 0; j < Genes.GetGenes( ); ++j ) {
				if( ( iGene = Answers.GetGene( Genes.GetGene( j ).GetName( ) ) ) == -1 )
					continue;
				for( k = 0; k < (int)MatPairs.GetSize( ); ++k )
					MatPairs.Set( k, iGene, false ); } } }

	for( iGenes = 0; !sArgs.inputs_num || ( iGenes < sArgs.inputs_num ); ++iGenes ) {
		MatResults.Clear( );
		MatGenes.Clear( );

		if( sArgs.inputs_num ) {
			CGenome		Genome;
			CGenes		Genes( Genome );
			ifstream	ifsm;

			ifsm.open( sArgs.inputs[ iGenes ] );
			if( !Genes.Open( ifsm ) ) {
				cerr << "Couldn't open: " << sArgs.inputs[ iGenes ] << endl;
				return 1; }
			vecfGenes.resize( Answers.GetGenes( ) );
			for( i = 0; i < vecfGenes.size( ); ++i )
				vecfGenes[ i ] = Genes.IsGene( Answers.GetGene( i ) );
			cerr << "Processing " << sArgs.inputs[ iGenes ] << "..." << endl;
			ifsm.close( ); }

		for( i = 0; i < Answers.GetGenes( ); ++i ) {
			if( ( iOne = veciGenes[ i ] ) == -1 )
				continue;
			for( j = ( i + 1 ); j < Answers.GetGenes( ); ++j ) {
				if( ( ( iTwo = veciGenes[ j ] ) == -1 ) ||
					CMeta::IsNaN( dValue = Data.Get( iOne, iTwo ) ) ||
					CMeta::IsNaN( dAnswer = Answers.Get( i, j ) ) )
					continue;
				if( !vecfGenes.empty( ) && !MatPairs.Get( i, j ) &&
					( ( dAnswer && !( vecfGenes[ i ] && vecfGenes[ j ] ) ) ||
					( !dAnswer && !( vecfGenes[ i ] || vecfGenes[ j ] ) ) ) )
					continue;
				if( sArgs.invert_flag )
					dValue = 1 - dValue;

				iMax = (int)ceil( ( dValue - sArgs.min_arg ) / sArgs.delta_arg );
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
			ofsm.open( ( (string)sArgs.directory_arg + '/' +
				CMeta::Basename( sArgs.inputs[ iGenes ] ) + ".bins" ).c_str( ) );
			postm = &ofsm; }
		else
			postm = &cout;

		*postm << "Cut\tGenes\tTP\tFP\tTN\tFN" << endl;
		for( i = 0; i < MatResults.GetRows( ); ++i ) {
			*postm << ( sArgs.min_arg + ( i * sArgs.delta_arg ) ) << '\t' << veciRec[ i ];
			for( j = 0; j < MatResults.GetColumns( ); ++j )
				*postm << '\t' << MatResults.Get( i, j );
			*postm << endl; }
		*postm << "#\tAUC\t" << CStatistics::WilcoxonRankSum( Data, Answers, vecfGenes, MatPairs,
			!!sArgs.invert_flag ) << endl;

		if( sArgs.inputs_num )
			ofsm.close( );
		else
			cout.flush( );

		if( !sArgs.inputs_num )
			break; }

	CMeta::Shutdown( );
	return 0; }

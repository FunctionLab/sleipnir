#include "stdafx.h"
#include "cmdline.h"

enum ETFPN {
	ETFPN_TP	= 0,
	ETFPN_FP	= ETFPN_TP + 1,
	ETFPN_TN	= ETFPN_FP + 1,
	ETFPN_FN	= ETFPN_TN + 1
};

struct SDatum {
	float	m_dValue;
	size_t	m_iOne;
	size_t	m_iTwo;

	SDatum( float dValue, size_t iOne, size_t iTwo ) : m_dValue(dValue), m_iOne(iOne), m_iTwo(iTwo) { }
};

struct SSorter {

	bool operator()( const SDatum& sOne, const SDatum& sTwo ) const {

		return ( sOne.m_dValue < sTwo.m_dValue ); }
};

int main( int iArgs, char** aszArgs ) {
	CDat				Answers, Data;
	gengetopt_args_info	sArgs;
	size_t				i, j, k, m, iOne, iTwo, iGenes, iPositives, iNegatives;
	vector<size_t>		veciGenes, veciRec;
	CFullMatrix<bool>	MatGenes;
	CFullMatrix<size_t>	MatResults;
	ETFPN				eTFPN;
	int					iMax;
	float				dAnswer, dValue;
	vector<bool>		vecfGenes;
	vector<float>		vecdScores;
	vector<size_t>		veciPositives, veciNegatives;
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
	MatResults.Initialize( sArgs.bins_arg ? ( sArgs.bins_arg + 1 ) :
		(size_t)( ( sArgs.max_arg - sArgs.min_arg ) / sArgs.delta_arg ) + 1, 4 );
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

		if( sArgs.bins_arg ) {
			vector<SDatum>	vecsData;
			size_t			iChunk;

			for( iPositives = iNegatives = i = 0; i < Answers.GetGenes( ); ++i ) {
				if( ( iOne = veciGenes[ i ] ) == -1 )
					continue;
				for( j = ( i + 1 ); j < Answers.GetGenes( ); ++j ) {
					if( ( ( iTwo = veciGenes[ j ] ) == -1 ) ||
						CMeta::IsNaN( dValue = Data.Get( iOne, iTwo ) ) ||
						CMeta::IsNaN( dAnswer = Answers.Get( i, j ) ) )
						continue;
					MatGenes.Set( i, 0, true );
					MatGenes.Set( j, 0, true );
					if( dAnswer )
						iPositives++;
					else
						iNegatives++;
					vecsData.push_back( SDatum( dValue, i, j ) ); } }
			sort( vecsData.begin( ), vecsData.end( ), SSorter( ) );
			veciPositives.resize( MatResults.GetRows( ) - 1 );
			veciNegatives.resize( veciPositives.size( ) );
			iChunk = (size_t)( 0.5 + ( (float)vecsData.size( ) / veciPositives.size( ) ) );
			for( i = j = 0; i < veciPositives.size( ); ++i,j += iChunk )
				for( k = 0; k < iChunk; ++k ) {
					if( ( j + k ) >= vecsData.size( ) )
						break;
					const SDatum&	sDatum	= vecsData[ j + k ];

					for( m = i; m > 0; --m ) {
						MatGenes.Set( sDatum.m_iOne, m, true );
						MatGenes.Set( sDatum.m_iTwo, m, true ); }
					if( Answers.Get( sDatum.m_iOne, sDatum.m_iTwo ) )
						veciPositives[ i ]++;
					else
						veciNegatives[ i ]++; }

			MatResults.Set( 0, ETFPN_TP, iPositives );
			MatResults.Set( 0, ETFPN_FP, iNegatives );
			MatResults.Set( 0, ETFPN_TN, 0 );
			MatResults.Set( 0, ETFPN_FN, 0 );
			for( i = 1; i < MatResults.GetRows( ); ++i ) {
				MatResults.Set( i, ETFPN_TP, MatResults.Get( i - 1, ETFPN_TP ) - veciPositives[ i - 1 ] );
				MatResults.Set( i, ETFPN_FP, MatResults.Get( i - 1, ETFPN_FP ) - veciNegatives[ i - 1 ] );
				MatResults.Set( i, ETFPN_TN, MatResults.Get( i - 1, ETFPN_TN ) + veciNegatives[ i - 1 ] );
				MatResults.Set( i, ETFPN_FN, MatResults.Get( i - 1, ETFPN_FN ) +
					veciPositives[ i - 1 ] ); } }
		else
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
					for( k = 0; (int)k < iMax; ++k ) {
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
			*postm << ( sArgs.bins_arg ? i : ( sArgs.min_arg + ( i * sArgs.delta_arg ) ) ) << '\t' <<
				veciRec[ i ];
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

#include "stdafx.h"
#include "cmdline.h"

struct SDatum {
	float	m_dDiff;
	size_t	m_iOne;
	size_t	m_iTwo;

	SDatum( float dDiff, size_t iOne, size_t iTwo ) : m_dDiff(dDiff), m_iOne(iOne), m_iTwo(iTwo) { }
};

struct SSorter {

	bool operator()( const SDatum& sOne, const SDatum& sTwo ) const {

		return ( sOne.m_dDiff < sTwo.m_dDiff ); }
};

int main( int iArgs, char** aszArgs ) {
	CDat				Answers, Data;
	gengetopt_args_info	sArgs;
	vector<size_t>		veciGenes;
	size_t				i, j, iOne, iTwo, iNumber;
	float				dValue, dAnswer;
	vector<SDatum>		vecsData;
	ifstream			ifsm;
	CGenome				Genome;
	CGene*				pOne;
	CGene*				pTwo;

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
	Data.Normalize( );
	if( sArgs.invert_flag )
		Data.Invert( );

	if( sArgs.features_arg ) {
		ifsm.open( sArgs.features_arg );
		if( !Genome.Open( ifsm ) ) {
			cerr << "Could not open: " << sArgs.features_arg << endl;
			return 1; }
		ifsm.close( ); }

	veciGenes.resize( Answers.GetGenes( ) );
	for( i = 0; i < Answers.GetGenes( ); ++i )
		veciGenes[ i ] = Data.GetGene( Answers.GetGene( i ) );
	for( i = 0; i < Answers.GetGenes( ); ++i ) {
		if( ( iOne = veciGenes[ i ] ) == -1 )
			continue;
		for( j = ( i + 1 ); j < Answers.GetGenes( ); ++j ) {
			if( ( ( iTwo = veciGenes[ j ] ) == -1 ) ||
				CMeta::IsNaN( dValue = Data.Get( iOne, iTwo ) ) ||
				CMeta::IsNaN( dAnswer = Answers.Get( i, j ) ) )
				continue;
			if( sArgs.positives_flag && !dAnswer )
				continue;
			vecsData.push_back( SDatum( fabs( dValue - dAnswer ), i, j ) ); } }
	sort( vecsData.begin( ), vecsData.end( ), SSorter( ) );

	if( ( ( iNumber = sArgs.number_arg ) < 0 ) || ( iNumber >= vecsData.size( ) ) )
		iNumber = vecsData.size( );
	for( i = 0; i < iNumber; ++i ) {
		const SDatum&	Datum	= vecsData[ i ];

		cout << Answers.GetGene( Datum.m_iOne ) << '\t' << Answers.GetGene( Datum.m_iTwo ) << '\t' <<
			Answers.Get( Datum.m_iOne, Datum.m_iTwo ) << '\t' << Data.Get( veciGenes[ Datum.m_iOne ],
			veciGenes[ Datum.m_iTwo ] ) << endl;
		if( Genome.GetGenes( ) ) {
			iOne = Genome.GetGene( Answers.GetGene( Datum.m_iOne ) );
			pOne = ( iOne == -1 ) ? NULL : &Genome.GetGene( iOne );
			iTwo = Genome.GetGene( Answers.GetGene( Datum.m_iTwo ) );
			pTwo = ( iTwo == -1 ) ? NULL : &Genome.GetGene( iTwo );

			cout << '\t';
			if( pOne && pOne->GetSynonyms( ) )
				cout << pOne->GetSynonym( 0 );
			cout << '\t';
			if( pOne )
				cout << pOne->GetGloss( );
			cout << endl << '\t';
			if( pTwo && pTwo->GetSynonyms( ) )
				cout << pTwo->GetSynonym( 0 );
			cout << '\t';
			if( pTwo )
				cout << pTwo->GetGloss( );
			cout << endl; } }

	CMeta::Shutdown( );
	return 0; }

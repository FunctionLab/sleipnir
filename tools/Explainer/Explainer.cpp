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

static const char	c_szExclude[]	= "exclude";
static const char	c_szInclude[]	= "include";
static const char	c_szOnly[]		= "only";
static const char	c_szUnknown[]	= "GO:0008150";

struct SDatum {
	float	m_dDiff;
	size_t	m_iOne;
	size_t	m_iTwo;

	SDatum( float dDiff, size_t iOne, size_t iTwo ) : m_dDiff(dDiff), m_iOne(iOne), m_iTwo(iTwo) { }
};

struct SSorter {
	bool	m_fReverse;

	SSorter( bool fReverse ) : m_fReverse(fReverse) { }

	bool operator()( const SDatum& sOne, const SDatum& sTwo ) const {
		bool	fRet;

		fRet = sOne.m_dDiff > sTwo.m_dDiff;
		return ( m_fReverse ? !fRet : fRet ); }
};

int main( int iArgs, char** aszArgs ) {
	CDat				Answers, Data;
	gengetopt_args_info	sArgs;
	vector<size_t>		veciGenes;
	size_t				i, j, k, iOne, iTwo, iNumber;
	float				dValue, dAnswer;
	vector<SDatum>		vecsData;
	ifstream			ifsm;
	COntologyGO			GOBP;
	CGenome				Genome;
	CGene*				pOne;
	CGene*				pTwo;
	bool				fOne, fTwo;
	string				strOne, strTwo;
	int					iRet;

	iRet = cmdline_parser2( iArgs, aszArgs, &sArgs, 0, 1, 0 );
	if( sArgs.config_arg )
		iRet = cmdline_parser_configfile( sArgs.config_arg, &sArgs, 0, 0, 1 ) && iRet;
	if( iRet ) {
		cmdline_parser_print_help( );
		return iRet; }
	CMeta Meta( sArgs.verbosity_arg );

	if( !strcmp( c_szInclude, sArgs.unknowns_arg ) || !strcmp( c_szOnly, sArgs.unknowns_arg ) )
		sArgs.everything_flag = true;

	if( !Answers.Open( sArgs.answers_arg, sArgs.memmap_flag && !sArgs.genes_arg && !sArgs.genet_arg && !sArgs.genex_arg ) ) {
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
	if( !Data.Open( sArgs.input_arg, sArgs.memmap_flag && !sArgs.normalize_flag && !sArgs.invert_flag ) ) {
		cerr << "Couldn't open: " << sArgs.input_arg << endl;
		return 1; }
	if( sArgs.normalize_flag )
		Data.Normalize( CDat::ENormalizeMinMax );
	if( sArgs.invert_flag )
		Data.Invert( );

	if( sArgs.features_arg ) {
		ifsm.open( sArgs.features_arg );
		if( !Genome.Open( ifsm ) ) {
			cerr << "Could not open: " << sArgs.features_arg << endl;
			return 1; }
		ifsm.close( ); }

	if( sArgs.go_onto_arg ) {
		ifstream	ifsmGenes;

		ifsm.clear( );
		ifsm.open( sArgs.go_onto_arg );
		if( sArgs.go_anno_arg )
			ifsmGenes.open( sArgs.go_anno_arg );
		if( !GOBP.Open( ifsm, ifsmGenes, Genome, COntologyGO::c_szBiologicalProcess ) ) {
			cerr << "Could not open: " << sArgs.go_onto_arg << ", " << sArgs.go_anno_arg << endl;
			return 1; }
		ifsm.close( ); }

	veciGenes.resize( Data.GetGenes( ) );
	for( i = 0; i < Data.GetGenes( ); ++i )
		veciGenes[ i ] = Answers.GetGene( Data.GetGene( i ) );
	for( i = 0; i < Data.GetGenes( ); ++i ) {
		if( !( i % 1000 ) )
			cerr << "Gene " << i << '/' << Data.GetGenes( ) << endl;
		if( !sArgs.everything_flag && ( ( iOne = veciGenes[ i ] ) == -1 ) )
			continue;
		for( j = ( i + 1 ); j < Data.GetGenes( ); ++j ) {
			if( CMeta::IsNaN( dValue = Data.Get( i, j ) ) )
				continue;
			if( !sArgs.everything_flag && ( ( ( iTwo = veciGenes[ j ] ) == -1 ) ||
				CMeta::IsNaN( dAnswer = Answers.Get( iOne, iTwo ) ) ||
				( sArgs.positives_flag && ( dAnswer <= 0 ) ) ) )
				continue;
			if( ( (float)rand( ) / RAND_MAX ) > sArgs.fraction_arg )
				continue;
			if( sArgs.everything_flag )
				dAnswer = dValue ? ( dValue - ( 1 / dValue ) ) : -FLT_MAX;
			vecsData.push_back( SDatum( fabs( dValue - dAnswer ), i, j ) ); } }
	sort( vecsData.begin( ), vecsData.end( ), SSorter( !!sArgs.reverse_flag ) );

	if( ( ( iNumber = sArgs.count_arg ) < 0 ) || ( iNumber >= vecsData.size( ) ) )
		iNumber = vecsData.size( );
	for( i = 0; i < iNumber; ++i ) {
		const SDatum&	Datum	= vecsData[ i ];

		pOne = pTwo = NULL;
		strOne = Data.GetGene( Datum.m_iOne );
		strTwo = Data.GetGene( Datum.m_iTwo );
		if( Genome.GetGenes( ) ) {
			iOne = Genome.GetGene( strOne );
			pOne = ( iOne == -1 ) ? NULL : &Genome.GetGene( iOne );
			iTwo = Genome.GetGene( strTwo );
			pTwo = ( iTwo == -1 ) ? NULL : &Genome.GetGene( iTwo );
			fOne = fTwo = true;
			if( pOne )
				for( j = 0; j < pOne->GetOntologies( ); ++j )
					if( pOne->GetAnnotations( j ) && ( ( pOne->GetAnnotations( j ) > 1 ) ||
						( pOne->GetOntology( j )->GetID( pOne->GetAnnotation( j, 0 ) ) != c_szUnknown ) ) ) {
						fOne = false;
						break; }
			if( pTwo )
				for( j = 0; j < pTwo->GetOntologies( ); ++j )
					if( pTwo->GetAnnotations( j ) && ( ( pTwo->GetAnnotations( j ) > 1 ) ||
						( pTwo->GetOntology( j )->GetID( pTwo->GetAnnotation( j, 0 ) ) != c_szUnknown ) ) ) {
						fTwo = false;
						break; }
			if( fOne || fTwo ) {
				if( !strcmp( c_szExclude, sArgs.unknowns_arg ) )
					continue; }
			else if( !strcmp( c_szOnly, sArgs.unknowns_arg ) )
				continue; }

		cout << strOne << '\t' << strTwo << '\t' << Data.Get( Datum.m_iOne, Datum.m_iTwo ) << '\t';
		dAnswer = ( ( ( j = veciGenes[ Datum.m_iOne ] ) == -1 ) || ( ( k = veciGenes[ Datum.m_iTwo ] ) == -1 ) ) ?
			CMeta::GetNaN( ) : Answers.Get( j, k );
		if( CMeta::IsNaN( dAnswer ) )
			cout << '-';
		else
			cout << dAnswer;
		cout << endl;
		if( pOne ) {
			cout << '\t';
			if( pOne->GetSynonyms( ) )
				cout << pOne->GetSynonym( 0 );
			cout << '\t' << pOne->GetGloss( ) << endl; }
		if( pTwo ) {
			cout << '\t';
			if( pTwo->GetSynonyms( ) )
				cout << pTwo->GetSynonym( 0 );
			cout << '\t' << pTwo->GetGloss( ) << endl; } }

	return 0; }

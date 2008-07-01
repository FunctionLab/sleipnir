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
	float	m_dAnswer;

	SDatum( float dValue, size_t iOne, size_t iTwo, float dAnswer ) : m_dValue(dValue), m_iOne(iOne), m_iTwo(iTwo),
		m_dAnswer(dAnswer) { }
};

struct SSorter {
	bool	m_fInvert;

	SSorter( bool fInvert ) : m_fInvert(fInvert) { }

	bool operator()( const SDatum& sOne, const SDatum& sTwo ) const {

		return ( m_fInvert ? ( sOne.m_dValue > sTwo.m_dValue ) : ( sOne.m_dValue < sTwo.m_dValue ) ); }
};

int main( int iArgs, char** aszArgs ) {
	CDat				Answers, Data;
	gengetopt_args_info	sArgs;
	size_t				i, j, k, m, iOne, iTwo, iGenes, iPositives, iNegatives, iBins;
	vector<size_t>		veciGenes, veciRec, veciRecTerm;
	CFullMatrix<bool>	MatGenes;
	CFullMatrix<size_t>	MatResults;
	ETFPN				eTFPN;
	int					iMax;
	float				dAnswer, dValue;
	vector<bool>		vecfHere;
	vector<float>		vecdScores, vecdSSE;
	vector<size_t>		veciPositives, veciNegatives, veciGenesTerm;
	ofstream			ofsm;
	ostream*			postm;
	map<float,size_t>	mapValues;
	bool				fMapAnswers;
	CGenome				Genome;
	CGenes				GenesTm( Genome );

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta = CMeta( sArgs.verbosity_arg );

	fMapAnswers = !!sArgs.memmap_flag && !( sArgs.genes_arg || sArgs.genet_arg || sArgs.genex_arg || sArgs.genee_arg );
	if( !Answers.Open( sArgs.answers_arg, fMapAnswers ) ) {
		cerr << "Couldn't open: " << sArgs.answers_arg << endl;
		return 1; }
	if( sArgs.genes_arg && !Answers.FilterGenes( sArgs.genes_arg, CDat::EFilterInclude ) ) {
		cerr << "Couldn't open: " << sArgs.genes_arg << endl;
		return 1; }
	if( sArgs.genee_arg && !Answers.FilterGenes( sArgs.genee_arg, CDat::EFilterEdge ) ) {
		cerr << "Couldn't open: " << sArgs.genee_arg << endl;
		return 1; }
	if( sArgs.genet_arg ) {
		if( !( Answers.FilterGenes( sArgs.genet_arg, CDat::EFilterTerm ) &&
			GenesTm.Open( sArgs.genet_arg ) ) ) {
			cerr << "Couldn't open: " << sArgs.genet_arg << endl;
			return 1; }
		veciGenesTerm.reserve( GenesTm.GetGenes( ) );
		for( i = 0; i < GenesTm.GetGenes( ); ++i )
			if( ( j = Answers.GetGene( GenesTm.GetGene( i ).GetName( ) ) ) != -1 )
				veciGenesTerm.push_back( j ); }
	if( sArgs.genex_arg && !Answers.FilterGenes( sArgs.genex_arg, CDat::EFilterExclude ) ) {
		cerr << "Couldn't open: " << sArgs.genex_arg << endl;
		return 1; }
	if( !Data.Open( sArgs.input_arg, !!sArgs.memmap_flag ) ) {
		cerr << "Couldn't open: " << sArgs.input_arg << endl;
		return 1; }
	if( sArgs.normalize_flag )
		Data.Normalize( CDat::ENormalizeMinMax );

	veciGenes.resize( Answers.GetGenes( ) );
	for( i = 0; i < Answers.GetGenes( ); ++i )
		veciGenes[ i ] = Data.GetGene( Answers.GetGene( i ) );
	if( sArgs.finite_flag ) {
		vector<float>	vecdValues;
		{
			set<float>		setdValues;

			for( i = 0; i < Answers.GetGenes( ); ++i ) {
				if( ( iOne = veciGenes[ i ] ) == -1 )
					continue;
				for( j = ( i + 1 ); j < Answers.GetGenes( ); ++j ) {
					if( ( ( iTwo = veciGenes[ j ] ) == -1 ) ||
						CMeta::IsNaN( dValue = Data.Get( iOne, iTwo ) ) ||
						CMeta::IsNaN( Answers.Get( i, j ) ) )
						continue;
					if( sArgs.invert_flag )
						dValue = 1 - dValue;
					setdValues.insert( dValue ); } }
			vecdValues.resize( setdValues.size( ) );
			copy( setdValues.begin( ), setdValues.end( ), vecdValues.begin( ) );
		}
		sort( vecdValues.begin( ), vecdValues.end( ) );
		for( i = 0; i < vecdValues.size( ); ++i )
			mapValues[ vecdValues[ i ] ] = i;
		iBins = mapValues.size( ); }
	else
		iBins = sArgs.bins_arg;
	MatResults.Initialize( iBins ? ( iBins + 1 ) :
		(size_t)( ( sArgs.max_arg - sArgs.min_arg ) / sArgs.delta_arg ) + 1, 4 );
	MatGenes.Initialize( veciGenes.size( ), MatResults.GetRows( ) );

	for( iGenes = 0; !sArgs.inputs_num || ( iGenes < sArgs.inputs_num ); ++iGenes ) {
		MatResults.Clear( );
		MatGenes.Clear( );

		if( sArgs.inputs_num ) {
			CGenes		Genes( Genome );
			ifstream	ifsm;

			ifsm.open( sArgs.inputs[ iGenes ] );
			if( !Genes.Open( ifsm ) ) {
				cerr << "Couldn't open: " << sArgs.inputs[ iGenes ] << endl;
				return 1; }
			vecfHere.resize( Answers.GetGenes( ) );
			for( i = 0; i < vecfHere.size( ); ++i )
				vecfHere[ i ] = Genes.IsGene( Answers.GetGene( i ) );
			cerr << "Processing " << sArgs.inputs[ iGenes ] << "..." << endl;
			ifsm.close( ); }

		if( mapValues.size( ) ) {
			for( i = 0; i < Answers.GetGenes( ); ++i ) {
				if( ( iOne = veciGenes[ i ] ) == -1 )
					continue;
				for( j = ( i + 1 ); j < Answers.GetGenes( ); ++j ) {
					if( ( ( iTwo = veciGenes[ j ] ) == -1 ) ||
						CMeta::IsNaN( dValue = Data.Get( iOne, iTwo ) ) ||
						CMeta::IsNaN( dAnswer = Answers.Get( i, j ) ) )
						continue;
					if( !( vecfHere.empty( ) ||
						( dAnswer && vecfHere[ i ] && vecfHere[ j ] ) ||
						( !dAnswer && ( vecfHere[ i ] || vecfHere[ j ] ) ) ) )
						continue;
					if( sArgs.invert_flag )
						dValue = 1 - dValue;
					for( k = 0; k <= mapValues[ dValue ]; ++k ) {
						MatGenes.Set( i, k, true );
						MatGenes.Set( j, k, true );
						MatResults.Get( k, dAnswer ? ETFPN_TP : ETFPN_FP )++; }
					for( ; k < MatResults.GetRows( ); ++k )
						MatResults.Get( k, dAnswer ? ETFPN_FN : ETFPN_TN )++; } } }
		else if( iBins ) {
			vector<SDatum>	vecsData;
			size_t			iChunk;

			for( iPositives = iNegatives = i = 0; i < Answers.GetGenes( ); ++i ) {
				if( ( iOne = veciGenes[ i ] ) == -1 )
					continue;
				for( j = ( i + 1 ); j < Answers.GetGenes( ); ++j ) {
					if( ( ( iTwo = veciGenes[ j ] ) == -1 ) ||
						CMeta::IsNaN( dAnswer = Answers.Get( i, j ) ) ||
						CMeta::IsNaN( dValue = Data.Get( iOne, iTwo ) ) )
						continue;
					if( !( vecfHere.empty( ) ||
						( dAnswer && vecfHere[ i ] && vecfHere[ j ] ) ||
						( !dAnswer && ( vecfHere[ i ] || vecfHere[ j ] ) ) ) )
						continue;

					MatGenes.Set( i, 0, true );
					MatGenes.Set( j, 0, true );
					if( dAnswer )
						iPositives++;
					else
						iNegatives++;
					vecsData.push_back( SDatum( dValue, i, j, dAnswer ) ); } }
			sort( vecsData.begin( ), vecsData.end( ), SSorter( !!sArgs.invert_flag ) );
			iChunk = (size_t)( 0.5 + ( (float)vecsData.size( ) / ( MatResults.GetRows( ) - 1 ) ) );
			if( sArgs.sse_flag ) {
				vecdSSE.resize( MatResults.GetRows( ) );
				veciPositives.resize( vecdSSE.size( ) );
				for( i = 1,j = 0; i < vecdSSE.size( ); ++i,j += iChunk ) {
					veciPositives[ veciPositives.size( ) - i - 1 ] = veciPositives[ veciPositives.size( ) - i ];
					vecdSSE[ vecdSSE.size( ) - i - 1 ] = vecdSSE[ vecdSSE.size( ) - i ];
					for( k = 0; k < iChunk; ++k ) {
						if( ( j + k ) >= vecsData.size( ) )
							break;
						const SDatum&	sDatum	= vecsData[ vecsData.size( ) - ( j + k ) - 1 ];

						for( m = 0; m < ( vecdSSE.size( ) - i ); ++m ) {
							MatGenes.Set( sDatum.m_iOne, m, true );
							MatGenes.Set( sDatum.m_iTwo, m, true ); }
						dValue = sDatum.m_dValue - sDatum.m_dAnswer;
						veciPositives[ veciPositives.size( ) - i - 1 ]++;
						vecdSSE[ vecdSSE.size( ) - i - 1 ] += dValue * dValue; } } }
			else {
				veciPositives.resize( MatResults.GetRows( ) - 1 );
				veciNegatives.resize( veciPositives.size( ) );
				for( i = 0; i < veciNegatives.size( ); ++i )
					veciNegatives[ i ] = veciPositives[ i ] = 0;
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
						veciPositives[ i - 1 ] ); } } }
		else
			for( i = 0; i < Answers.GetGenes( ); ++i ) {
				if( !( i % 1000 ) )
					cerr << "Processing gene " << i << '/' << Answers.GetGenes( ) << endl;
				if( ( iOne = veciGenes[ i ] ) == -1 )
					continue;
				for( j = ( i + 1 ); j < Answers.GetGenes( ); ++j ) {
					if( ( ( iTwo = veciGenes[ j ] ) == -1 ) ||
						CMeta::IsNaN( dAnswer = Answers.Get( i, j ) ) ||
						CMeta::IsNaN( dValue = Data.Get( iOne, iTwo ) ) )
						continue;
					if( !( vecfHere.empty( ) ||
						( dAnswer && vecfHere[ i ] && vecfHere[ j ] ) ||
						( !dAnswer && ( vecfHere[ i ] || vecfHere[ j ] ) ) ) )
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
		for( iPositives = iNegatives = i = 0; i < Answers.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < Answers.GetGenes( ); ++j ) {
				if( CMeta::IsNaN( dAnswer = Answers.Get( i, j ) ) ||
					!( vecfHere.empty( ) ||
					( dAnswer && vecfHere[ i ] && vecfHere[ j ] ) ||
					( !dAnswer && ( vecfHere[ i ] || vecfHere[ j ] ) ) ) )
					continue;
				if( dAnswer )
					iPositives++;
				else
					iNegatives++; }

		veciRec.resize( MatResults.GetRows( ) );
		veciRecTerm.resize( MatResults.GetRows( ) );
		for( i = 0; i < veciRec.size( ); ++i ) {
			veciRec[ i ] = veciRecTerm[ i ] = 0;
			for( j = 0; j < MatGenes.GetRows( ); ++j )
				if( MatGenes.Get( j, i ) ) {
					veciRec[ i ]++;
					if( vecfHere.size( ) && vecfHere[ j ] )
						veciRecTerm[ i ]++; }
			for( j = 0; j < veciGenesTerm.size( ); ++j )
				if( MatGenes.Get( veciGenesTerm[ j ], i ) &&
					( vecfHere.empty( ) || !vecfHere[ veciGenesTerm[ j ] ] ) )
					veciRecTerm[ i ]++; }

		if( sArgs.inputs_num ) {
			ofsm.open( ( (string)sArgs.directory_arg + '/' +
				CMeta::Basename( sArgs.inputs[ iGenes ] ) + ".bins" ).c_str( ) );
			postm = &ofsm; }
		else
			postm = &cout;

		if( !sArgs.sse_flag ) {
			*postm << "#	P	" << iPositives << endl;
			*postm << "#	N	" << iNegatives << endl; }
		*postm << "Cut	Genes	" << ( sArgs.sse_flag ? "Pairs	SSE" : "TP	FP	TN	FN" ) << endl;
		for( i = 0; i < MatResults.GetRows( ); ++i ) {
			*postm << ( iBins ? i : ( sArgs.min_arg + ( i * sArgs.delta_arg ) ) ) << '\t' <<
				veciRec[ i ];
			if( sArgs.sse_flag )
				*postm << '\t' << veciPositives[ i ] << '\t' << vecdSSE[ i ];
			else
				for( j = 0; j < MatResults.GetColumns( ); ++j )
					*postm << '\t' << MatResults.Get( i, j );
			if( veciGenesTerm.size( ) || vecfHere.size( ) )
				*postm << '\t' << veciRecTerm[ i ];
			*postm << endl; }
		if( !sArgs.sse_flag )
			*postm << "#	AUC	" << CStatistics::WilcoxonRankSum( Data, Answers, vecfHere,
				!!sArgs.invert_flag ) << endl;

		if( sArgs.inputs_num )
			ofsm.close( );
		else
			cout.flush( );

		if( !sArgs.inputs_num )
			break; }

	return 0; }

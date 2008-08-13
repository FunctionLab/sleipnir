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

struct SSorter {
	const vector<size_t>&	m_veciSizes;

	SSorter( const vector<size_t>& veciSizes ) : m_veciSizes(veciSizes) { }

	bool operator()( size_t iOne, size_t iTwo ) const {

		return ( m_veciSizes[ iOne ] > m_veciSizes[ iTwo ] ); }
};

size_t find_value( size_t iOne, size_t iTwo, const CDataPair& Dat, size_t iDefault, bool fRandom ) {
	float	d;

	d = ( ( iOne == -1 ) || ( iTwo == -1 ) ) ? CMeta::GetNaN( ) : Dat.Get( iOne, iTwo );
	return ( CMeta::IsNaN( d ) ? ( fRandom ? ( rand( ) % Dat.GetValues( ) ) : iDefault ) :
		Dat.Quantize( d ) ); }

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info			sArgs;
	CDataPair					DatOne, DatTwo;
	size_t						i, j, iDatOne, iDatTwo, iGeneOne, iGeneTwo, iCountOne, iCountTwo;
	size_t						iCountJoint, iJoint, iValueOne, iValueTwo;
	vector<size_t>				veciGenesOne, veciGenesTwo, veciOne, veciTwo, veciDefaults, veciSizes;
	vector<vector<size_t> >		vecveciJoint;
	vector<string>				vecstrInputs, vecstrGenes;
	float						dOne, dJoint, dMI;
	map<string, size_t>			mapZeros;
	vector<float>				vecdOne, vecdTwo;
	float*						adOne;
	float*						adTwo;
	IMeasure*					pMeasure;
	CMeasurePearson				Pearson;
	CMeasureEuclidean			Euclidean;
	CMeasureKendallsTau			KendallsTau;
	CMeasureKolmogorovSmirnov	KolmSmir;
	CMeasureHypergeometric		Hypergeom;
	CMeasureQuickPearson		PearQuick;
	CMeasureInnerProduct		InnerProd;
	CMeasureBinaryInnerProduct	BinInnerProd;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg, sArgs.random_arg );

	CMeasureSigmoid				EuclideanSig( &Euclidean, false, 1.0f / sArgs.inputs_num );
	IMeasure*					apMeasures[]	= { &Pearson, &EuclideanSig, &KendallsTau,
		&KolmSmir, &Hypergeom, &PearQuick, &InnerProd, &BinInnerProd, NULL };

	vecstrInputs.resize( sArgs.inputs_num );
	copy( sArgs.inputs, sArgs.inputs + sArgs.inputs_num, vecstrInputs.begin( ) );
	if( sArgs.only_arg == -1 ) {
		vector<size_t>	veciIndices, veciSizes;

		veciSizes.resize( vecstrInputs.size( ) );
		veciIndices.resize( vecstrInputs.size( ) );
		for( i = 0; i < vecstrInputs.size( ); ++i ) {
			CDat		Dat;
			ifstream	ifsm;

			veciIndices[ i ] = i;
			ifsm.open( vecstrInputs[ i ].c_str( ) );
			Dat.OpenGenes( ifsm, true );
			veciSizes[ i ] = Dat.GetGenes( ); }
		sort( veciIndices.begin( ), veciIndices.end( ), SSorter( veciSizes ) );
		CMeta::Permute( vecstrInputs, veciIndices ); }
	{
		CDataset	Data;

		Data.OpenGenes( vecstrInputs );
		vecstrGenes.resize( Data.GetGenes( ) );
		copy( Data.GetGeneNames( ).begin( ), Data.GetGeneNames( ).end( ), vecstrGenes.begin( ) );
	}
	veciSizes.resize( vecstrInputs.size( ) );
	for( i = 0; i < veciSizes.size( ); ++i ) {
		DatOne.OpenQuants( vecstrInputs[ i ].c_str( ) );
		veciSizes[ i ] = DatOne.GetValues( ); }

	pMeasure = NULL;
	if( sArgs.distance_arg )
		for( i = 0; apMeasures[ i ]; ++i )
			if( !strcmp( apMeasures[ i ]->GetName( ), sArgs.distance_arg ) ) {
				pMeasure = apMeasures[ i ];
				break; }

	if( sArgs.zeros_arg ) {
		ifstream		ifsm;
		vector<string>	vecstrZeros;
		char			acLine[ 1024 ];

		ifsm.open( sArgs.zeros_arg );
		if( !ifsm.is_open( ) ) {
			cerr << "Could not open: " << sArgs.zeros_arg << endl;
			return 1; }
		while( !ifsm.eof( ) ) {
			ifsm.getline( acLine, ARRAYSIZE(acLine) - 1 );
			acLine[ ARRAYSIZE(acLine) - 1 ] = 0;
			vecstrZeros.clear( );
			CMeta::Tokenize( acLine, vecstrZeros );
			if( vecstrZeros.empty( ) )
				continue;
			mapZeros[ vecstrZeros[ 0 ] ] = atoi( vecstrZeros[ 1 ].c_str( ) ); } }
	veciDefaults.resize( vecstrInputs.size( ) );
	for( i = 0; i < veciDefaults.size( ); ++i ) {
		map<string, size_t>::const_iterator	iterZero;

		if( ( iterZero = mapZeros.find( CMeta::Deextension( CMeta::Basename(
			vecstrInputs[ i ].c_str( ) ) ) ) ) != mapZeros.end( ) )
			veciDefaults[ i ] = iterZero->second;
		else
			veciDefaults[ i ] = sArgs.randomize_flag ? -1 :
				( sArgs.zero_flag ? 0 : veciSizes[ i ]++ ); }

	if( sArgs.table_flag ) {
		for( i = 0; i < vecstrInputs.size( ); ++i )
			cout << '\t' << vecstrInputs[ i ];
		cout << endl; }
	for( iDatOne = ( ( sArgs.only_arg == -1 ) ? 0 : sArgs.only_arg );
		iDatOne < ( ( sArgs.only_arg == -1 ) ? vecstrInputs.size( ) : ( sArgs.only_arg + 1 ) );
		++iDatOne ) {
		cerr << "Processing " << iDatOne << '/' << vecstrInputs.size( ) << endl;
		if( sArgs.table_flag ) {
			cout << vecstrInputs[ iDatOne ];
			for( i = 0; i < iDatOne; ++i )
				cout << '\t'; }

		DatOne.Open( vecstrInputs[ iDatOne ].c_str( ), false, !!sArgs.memmap_flag );
		veciGenesOne.resize( vecstrGenes.size( ) );
		for( i = 0; i < vecstrGenes.size( ); ++i )
			veciGenesOne[ i ] = DatOne.GetGene( vecstrGenes[ i ] );
		veciOne.resize( veciSizes[ iDatOne ] );
		fill( veciOne.begin( ), veciOne.end( ), 0 );
		for( i = iCountOne = 0; i < vecstrGenes.size( ); ++i ) {
			iGeneOne = veciGenesOne[ i ];
			for( j = ( i + 1 ); j < vecstrGenes.size( ); ++j )
				if( ( iValueOne = find_value( iGeneOne, veciGenesOne[ j ], DatOne,
					veciDefaults[ iDatOne ], !!sArgs.randomize_flag ) ) != -1 ) {
					iCountOne++;
					veciOne[ iValueOne ]++; } }

		if( !pMeasure ) {
/*
cout << iCountOne<< ':';
for( i = 0; i < veciOne.size( ); ++i )
cout << ' ' << veciOne[ i ];
cout << endl;
//*/
			for( dMI = 0,i = 0; i < veciOne.size( ); ++i )
				if( dOne = (float)veciOne[ i ] / iCountOne )
					dMI += dOne * log( 1 / dOne );
			dMI /= log( 2.0f );
			if( sArgs.table_flag )
				cout << '\t' << dMI;
			else
				cout << vecstrInputs[ iDatOne ] << '\t' << vecstrInputs[ iDatOne ] << '\t' << dMI <<
					endl; }
		if( ( iDatOne + 1 ) == vecstrInputs.size( ) )
			break;

		for( iDatTwo = ( iDatOne + 1 ); iDatTwo < vecstrInputs.size( ); ++iDatTwo ) {
			DatTwo.Open( vecstrInputs[ iDatTwo ].c_str( ), false, !!sArgs.memmap_flag );
			veciGenesTwo.resize( vecstrGenes.size( ) );
			for( i = 0; i < veciGenesTwo.size( ); ++i )
				veciGenesTwo[ i ] = DatTwo.GetGene( vecstrGenes[ i ] );
			veciTwo.resize( veciSizes[ iDatTwo ] );
			fill( veciTwo.begin( ), veciTwo.end( ), 0 );

			vecveciJoint.resize( veciOne.size( ) );
			for( i = 0; i < vecveciJoint.size( ); ++i ) {
				vecveciJoint[ i ].resize( veciTwo.size( ) );
				fill( vecveciJoint[ i ].begin( ), vecveciJoint[ i ].end( ), 0 ); }

			if( pMeasure ) {
				vecdOne.clear( );
				vecdTwo.clear( ); }
			for( i = iCountTwo = iCountJoint = 0; i < vecstrGenes.size( ); ++i ) {
				iGeneOne = veciGenesOne[ i ];
				iGeneTwo = veciGenesTwo[ i ];
				for( j = ( i + 1 ); j < vecstrGenes.size( ); ++j )
					if( ( iValueTwo = find_value( iGeneTwo, veciGenesTwo[ j ], DatTwo,
						veciDefaults[ iDatTwo ], !!sArgs.randomize_flag ) ) != -1 ) {
						iCountTwo++;
						veciTwo[ iValueTwo ]++;
						if( ( iValueOne = find_value( iGeneOne, veciGenesOne[ j ], DatOne,
							veciDefaults[ iDatOne ], !!sArgs.randomize_flag ) ) != -1 ) {
							if( pMeasure ) {
								vecdOne.push_back( (float)iValueOne );
								vecdTwo.push_back( (float)iValueTwo ); }
							iCountJoint++;
							vecveciJoint[ iValueOne ][ iValueTwo ]++; } } }

			if( pMeasure ) {
				adOne = new float[ vecdOne.size( ) ];
				adTwo = new float[ vecdTwo.size( ) ];
				for( i = 0; i < vecdOne.size( ); ++i ) {
					adOne[ i ] = vecdOne[ i ];
					adTwo[ i ] = vecdTwo[ i ]; }
				dMI = (float)pMeasure->Measure( adOne, vecdOne.size( ), adTwo, vecdTwo.size( ) );
				delete[] adTwo;
				delete[] adOne; }
			else {
/*
cout << iCountOne << ':';
for( i = 0; i < veciOne.size( ); ++i )
cout << ' ' << veciOne[ i ];
cout << endl;
cout << iCountTwo << ':';
for( i = 0; i < veciTwo.size( ); ++i )
cout << ' ' << veciTwo[ i ];
cout << endl << endl << iCountJoint << ':' << endl;
for( i = 0; i < vecveciJoint.size( ); ++i ) {
for( j = 0; j < vecveciJoint[ i ].size( ); ++j )
cout << vecveciJoint[ i ][ j ] << ' ';
cout << endl; }
cout << endl;
//*/
				for( dMI = 0,i = 0; i < veciOne.size( ); ++i ) {
					dOne = (float)veciOne[ i ] / iCountOne;
					for( j = 0; j < veciTwo.size( ); ++j )
						if( iJoint = vecveciJoint[ i ][ j ] ) {
							dJoint = (float)iJoint / iCountJoint;
							dMI += dJoint * log( dJoint * iCountTwo / dOne / veciTwo[ j ] ); } }
				dMI -= ( veciOne.size( ) - 1 ) * ( veciTwo.size( ) - 1 ) / ( 2.0f * ( iCountOne +
					iCountTwo ) );
				dMI = ( dMI < 0 ) ? 0 : ( dMI / log( 2.0f ) ); }
			if( sArgs.table_flag )
				cout << '\t' << dMI;
			else
				cout << vecstrInputs[ iDatOne ] << '\t' << vecstrInputs[ iDatTwo ] << '\t' << dMI <<
					endl; }
		if( sArgs.table_flag )
			cout << endl; }

	return 0; }

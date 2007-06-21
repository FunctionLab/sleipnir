#include "stdafx.h"
#include "cmdline.h"

size_t find_value( size_t iOne, size_t iTwo, const CDataPair& Dat, size_t iDefault ) {
	float	d;

	d = ( ( iOne == -1 ) || ( iTwo == -1 ) ) ? CMeta::GetNaN( ) : Dat.Get( iOne, iTwo );
	return ( CMeta::IsNaN( d ) ? iDefault : Dat.Quantize( d ) ); }

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info		sArgs;
	CDataPair				DatOne, DatTwo;
	size_t					i, j, iDatOne, iDatTwo, iGeneOne, iGeneTwo, iCountOne, iCountTwo;
	size_t					iCountJoint, iJoint, iValueOne, iValueTwo, iCountOneBase;
	vector<size_t>			veciGenesOne, veciGenesTwo, veciOne, veciTwo, veciOneBase, veciDefaults;
	vector<vector<size_t> >	vecveciJoint;
	float					dOne, dJoint, dMI;
	map<string, size_t>		mapZeros;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

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
	veciDefaults.resize( sArgs.inputs_num );
	for( i = 0; i < veciDefaults.size( ); ++i ) {
		map<string, size_t>::const_iterator	iterZero;

		if( ( iterZero = mapZeros.find( CMeta::Deextension( CMeta::Basename(
			sArgs.inputs[ i ] ) ) ) ) != mapZeros.end( ) )
			veciDefaults[ i ] = iterZero->second;
		else
			veciDefaults[ i ] = sArgs.zero_flag ? 0 : -1; }

	if( sArgs.table_flag ) {
		for( i = 0; i < sArgs.inputs_num; ++i )
			cout << '\t' << sArgs.inputs[ i ];
		cout << endl; }
	for( iDatOne = ( ( sArgs.only_arg == -1 ) ? 0 : sArgs.only_arg );
		iDatOne < ( ( sArgs.only_arg == -1 ) ? sArgs.inputs_num : ( sArgs.only_arg + 1 ) );
		++iDatOne ) {
		if( sArgs.table_flag ) {
			cout << sArgs.inputs[ iDatOne ];
			for( i = 0; i <= iDatOne; ++i )
				cout << '\t'; }
		if( ( iDatOne + 1 ) == sArgs.inputs_num )
			break;
		DatOne.Open( sArgs.inputs[ iDatOne ], false, !!sArgs.memmap_flag );
		veciGenesOne.resize( DatOne.GetGenes( ) );
		veciOneBase.resize( DatOne.GetValues( ) );
		fill( veciOneBase.begin( ), veciOneBase.end( ), 0 );
		for( i = iCountOneBase = 0; i < DatOne.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < DatOne.GetGenes( ); ++j )
				if( ( iValueOne = find_value( i, j, DatOne, veciDefaults[ iDatOne ] ) ) != -1 ) {
					iCountOneBase++;
					veciOneBase[ iValueOne ]++; }
		for( iDatTwo = ( iDatOne + 1 ); iDatTwo < sArgs.inputs_num; ++iDatTwo ) {
			DatTwo.Open( sArgs.inputs[ iDatTwo ], false, !!sArgs.memmap_flag );
			for( i = 0; i < veciGenesOne.size( ); ++i )
				veciGenesOne[ i ] = DatTwo.GetGene( DatOne.GetGene( i ) );
			veciGenesTwo.resize( DatTwo.GetGenes( ) );
			for( i = 0; i < veciGenesTwo.size( ); ++i )
				veciGenesTwo[ i ] = DatOne.GetGene( DatTwo.GetGene( i ) );
			veciTwo.resize( DatTwo.GetValues( ) );
			fill( veciTwo.begin( ), veciTwo.end( ), 0 );

			iCountOne = iCountOneBase;
			veciOne.resize( veciOneBase.size( ) );
			copy( veciOneBase.begin( ), veciOneBase.end( ), veciOne.begin( ) );
			vecveciJoint.resize( veciOne.size( ) );
			for( i = 0; i < vecveciJoint.size( ); ++i ) {
				vecveciJoint[ i ].resize( veciTwo.size( ) );
				fill( vecveciJoint[ i ].begin( ), vecveciJoint[ i ].end( ), 0 ); }

			for( i = iCountTwo = iCountJoint = 0; i < DatTwo.GetGenes( ); ++i ) {
				iGeneOne = veciGenesTwo[ i ];
				for( j = ( i + 1 ); j < DatTwo.GetGenes( ); ++j ) {
					iGeneTwo = veciGenesTwo[ j ];
					if( ( iValueTwo = find_value( i, j, DatTwo, veciDefaults[ iDatTwo ] ) ) != -1 ) {
						iCountTwo++;
						veciTwo[ iValueTwo ]++;
						if( ( iValueOne = find_value( iGeneOne, iGeneTwo, DatOne,
							veciDefaults[ iDatOne ] ) ) != -1 ) {
							if( ( iGeneOne == -1 ) || ( iGeneTwo == -1 ) ) {
								iCountOne++;
								veciOne[ iValueOne ]++; }
							iCountJoint++;
							vecveciJoint[ iValueOne ][ iValueTwo ]++; } } } }
			for( i = 0; i < veciGenesOne.size( ); ++i ) {
					iGeneOne = veciGenesOne[ i ];
					for( j = ( i + 1 ); j < veciGenesOne.size( ); ++j ) {
						iGeneTwo = veciGenesOne[ j ];
						if( !( ( iGeneOne == -1 ) || ( iGeneTwo == -1 ) ) )
							continue;
						if( ( iValueTwo = find_value( -1, -1, DatTwo, veciDefaults[ iDatTwo ] ) ) != -1 ) {
							iCountTwo++;
							veciTwo[ iValueTwo ]++;
							if( ( iValueOne = find_value( i, j, DatOne, veciDefaults[ iDatOne ] ) ) != -1 ) {
								iCountJoint++;
								vecveciJoint[ iValueOne ][ iValueTwo ]++; } } } }

			for( dMI = 0,i = 0; i < veciOne.size( ); ++i ) {
				dOne = (float)veciOne[ i ] / iCountOne;
				for( j = 0; j < veciTwo.size( ); ++j )
					if( iJoint = vecveciJoint[ i ][ j ] ) {
						dJoint = (float)iJoint / iCountJoint;
						dMI += dJoint * log( dJoint * iCountTwo / dOne / veciTwo[ j ] ); } }
			dMI -= ( veciOne.size( ) - 1 ) * ( veciTwo.size( ) - 1 ) / ( 2.0f * ( iCountOne + iCountTwo ) );
			dMI = ( dMI < 0 ) ? 0 : ( dMI / log( 2.0f ) );
			if( sArgs.table_flag )
				cout << '\t' << dMI;
			else
				cout << sArgs.inputs[ iDatOne ] << '\t' << sArgs.inputs[ iDatTwo ] << '\t' << dMI <<
					endl; }
		if( sArgs.table_flag )
			cout << endl; }

	CMeta::Shutdown( );
	return 0; }

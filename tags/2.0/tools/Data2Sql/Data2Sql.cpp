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

int main( int iArgs, char** aszArgs ) {
	static const size_t					c_iBuffer	= 1024;
	gengetopt_args_info					sArgs;
	ifstream							ifsm;
	istream*							pistm;
	size_t								iFile, i, j, iOne, iTwo, iFirst, iSecond, iCount;
	float								d;
	map<string, size_t>					mapstriGenes;
	map<string, size_t>::const_iterator	iterGene;
	vector<string>						vecstrLine;
	char								acBuffer[ c_iBuffer ];
	vector<size_t>						veciGenes;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg );

	if( sArgs.input_arg ) {
		ifsm.open( sArgs.input_arg );
		pistm = &ifsm; }
	else
		pistm = &cin;
	while( !pistm->eof( ) ) {
		pistm->getline( acBuffer, c_iBuffer - 1 );
		acBuffer[ c_iBuffer - 1 ] = 0;
		vecstrLine.clear( );
		CMeta::Tokenize( acBuffer, vecstrLine );
		if( vecstrLine.size( ) != 2 ) {
			cerr << "Ignoring line: " << acBuffer << endl;
			continue; }
		mapstriGenes[ vecstrLine[ 1 ] ] = atoi( vecstrLine[ 0 ].c_str( ) ); }
	if( sArgs.input_arg )
		ifsm.close( );

	for( iCount = iFile = 0; iFile < sArgs.inputs_num; ++iFile ) {
		CDataPair	Dat;

		if( sArgs.datasets_flag ) {
			cout << ( iFile + 1 ) << '\t' << CMeta::Deextension( CMeta::Basename( sArgs.inputs[ iFile ] ) ) <<
				endl;
			continue; }
		if( !Dat.Open( sArgs.inputs[ iFile ], false, !!sArgs.memmap_flag ) ) {
			cerr << "Could not open: " << sArgs.inputs[ iFile ] << endl;
			return 1; }
		veciGenes.resize( Dat.GetGenes( ) );
		for( i = 0; i < veciGenes.size( ); ++i )
#ifdef _MSC_VER
			(size_t)
#endif // _MSC_VER
			veciGenes[ i ] = ( ( iterGene = mapstriGenes.find( Dat.GetGene( i ) ) ) ==
				mapstriGenes.end( ) ) ? -1 : iterGene->second;
		for( i = 0; i < veciGenes.size( ); ++i ) {
			if( !( i % 100 ) )
				cerr << i << '/' << veciGenes.size( ) << endl;
			if( ( iOne = veciGenes[ i ] ) == -1 )
				continue;
			for( j = ( i + 1 ); j < veciGenes.size( ); ++j )
				if( ( ( iTwo = veciGenes[ j ] ) != -1 ) && !CMeta::IsNaN( d = Dat.Get( i, j ) ) ) {
					if( iOne < iTwo ) {
						iFirst = iOne;
						iSecond = iTwo; }
					else {
						iFirst = iTwo;
						iSecond = iOne; }
					if( iCount % sArgs.block_arg )
						cout << ',';
					else
						cout << "INSERT INTO " << sArgs.table_arg << " VALUES " << endl;
					cout << '(' << ( iFile + 1 ) << ',' << iFirst << ',' << iSecond << ',' <<
						Dat.Quantize( d ) << ')';
					if( !( ++iCount % sArgs.block_arg ) )
						cout << ';' << endl; } } }
	if( iCount )
		cout << ';' << endl;

	return 0; }

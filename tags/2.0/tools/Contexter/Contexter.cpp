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

static const char	c_szXDSL[]	= ".xdsl";
static const char	c_szDSL[]	= ".dsl";

int main_database( const gengetopt_args_info& );
int main_datfile( const gengetopt_args_info& );
size_t count_contexts( const char* );
bool open_contexts( const char*, size_t, size_t, CCompactFullMatrix& );

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	int					iRet;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg );

	iRet = sArgs.dat_flag ? main_datfile( sArgs ) : main_database( sArgs );

	return iRet; }

size_t count_contexts( const char* szContexts ) {
	static const size_t	c_iBuffer	= 1024;
	char				acBuffer[ c_iBuffer ];
	ifstream			ifsm;
	vector<string>		vecstrLine;
	size_t				iRet, iCur;

	iRet = -1;
	ifsm.open( szContexts );
	if( !ifsm.is_open( ) ) {
		cerr << "Could not open: " << szContexts << endl;
		return -1; }
	while( !ifsm.eof( ) ) {
		ifsm.getline( acBuffer, c_iBuffer - 1 );
		acBuffer[ c_iBuffer - 1 ] = 0;
		if( !acBuffer[ 0 ] )
			continue;
		vecstrLine.clear( );
		CMeta::Tokenize( acBuffer, vecstrLine );
		if( vecstrLine.empty( ) ) {
			cerr << "Invalid line in " << szContexts << ": " << acBuffer << endl;
			return -1; }
		iCur = atoi( vecstrLine[ 0 ].c_str( ) );
		if( ( iRet == -1 ) || ( iCur > iRet ) )
			iRet = iCur; }
	ifsm.close( );

	return iRet; }

bool open_contexts( const char* szContexts, size_t iContexts, size_t iGenes,
	CCompactFullMatrix& MatContexts ) {
	static const size_t	c_iBuffer	= 1024;
	char				acBuffer[ c_iBuffer ];
	ifstream			ifsm;
	vector<string>		vecstrLine;

	if( ( iContexts == -1 ) && ( ( iContexts = count_contexts( szContexts ) ) == -1 ) )
		return false;

	ifsm.open( szContexts );
	if( !ifsm.is_open( ) ) {
		cerr << "Could not open: " << szContexts << endl;
		return false; }
	MatContexts.Initialize( iContexts, iGenes, 2, true );
	while( !ifsm.eof( ) ) {
		size_t	iContext, iGene;

		ifsm.getline( acBuffer, c_iBuffer - 1 );
		acBuffer[ c_iBuffer - 1 ] = 0;
		if( !acBuffer[ 0 ] )
			continue;
		vecstrLine.clear( );
		CMeta::Tokenize( acBuffer, vecstrLine );
		if( vecstrLine.size( ) != 2 ) {
			cerr << "Invalid line in " << szContexts << ": " << acBuffer << endl;
			return false; }
		iContext = atoi( vecstrLine[ 0 ].c_str( ) ) - 1;
		iGene = atoi( vecstrLine[ 1 ].c_str( ) ) - 1;
		if( iContext >= MatContexts.GetRows( ) ) {
			cerr << "Invalid context on line: " << acBuffer << endl;
			return false; }
		if( iGene >= MatContexts.GetColumns( ) ) {
			cerr << "Invalid gene on line: " << acBuffer << endl;
			return false; }
		MatContexts.Set( iContext, iGene, 1 ); }
	ifsm.close( );

	return true; }

int main_database( const gengetopt_args_info& sArgs ) {
	static const size_t					c_iBuffer	= 1024;
	char								acBuffer[ c_iBuffer ];
	CBayesNetSmile						BNSmile;
	CBayesNetMinimal					BNDefault;
	vector<CBayesNetMinimal>			vecBNs;
	ifstream							ifsm;
	istream*							pistm;
	map<size_t, string>					mapistrBNs;
	map<size_t, string>::const_iterator	iterBN;
	size_t								i, iMax, iGene;
	CDatabase							Database;
	uint32_t							iSize;
	float*								adGenes;
	CDat								Dat;
	CCompactFullMatrix					MatContexts;
	vector<string>						vecstrLine;

	if( !Database.Open( sArgs.database_arg ) ) {
		cerr << "Could not open: " << sArgs.database_arg << endl;
		return 1; }
	if( sArgs.minimal_in_flag ) {
		ifsm.open( sArgs.networks_arg, ios_base::binary );
		if( !BNDefault.Open( ifsm ) ) {
			cerr << "Could not read: " << sArgs.networks_arg << endl;
			return 1; }
		ifsm.read( (char*)&iSize, sizeof(iSize) );
		vecBNs.resize( iSize );
		for( i = 0; i < vecBNs.size( ); ++i )
			if( !vecBNs[ i ].Open( ifsm ) ) {
				cerr << "Could not read: " << sArgs.networks_arg << " (" << i << ")" << endl;
				return 1; }
		ifsm.close( ); }
	else {
		if( sArgs.default_arg && !( BNSmile.Open( sArgs.default_arg ) && BNDefault.Open( BNSmile ) ) ) {
			cerr << "Could not open: " << sArgs.default_arg << endl;
			return 1; }
		BNDefault.SetID( sArgs.default_arg );
		if( sArgs.input_arg ) {
			ifsm.open( sArgs.input_arg );
			pistm = &ifsm; }
		else
			pistm = &cin;
		iMax = 0;
		while( !pistm->eof( ) ) {
			pistm->getline( acBuffer, c_iBuffer - 1 );
			acBuffer[ c_iBuffer - 1 ] = 0;
			if( !acBuffer[ 0 ] )
				continue;
			vecstrLine.clear( );
			CMeta::Tokenize( acBuffer, vecstrLine );
			if( vecstrLine.size( ) < 2 ) {
				cerr << "Ignoring line: " << acBuffer << endl;
				continue; }
			if( ( i = atoi( vecstrLine[ 0 ].c_str( ) ) ) > iMax )
				iMax = i;
			mapistrBNs[ i ] = vecstrLine[ 1 ]; }
		if( sArgs.input_arg )
			ifsm.close( );
		vecBNs.resize( iMax );
		for( iterBN = mapistrBNs.begin( ); iterBN != mapistrBNs.end( ); ++iterBN ) {
			if( !( BNSmile.Open( ( (string)sArgs.networks_arg + '/' + CMeta::Filename( iterBN->second ) +
				( sArgs.xdsl_flag ? c_szXDSL : c_szDSL ) ).c_str( ) ) &&
				vecBNs[ iterBN->first - 1 ].Open( BNSmile ) ) ) {
				cerr << "Could not open: " << iterBN->second << endl;
				return 1; }
			vecBNs[ iterBN->first - 1 ].SetID( iterBN->second ); } }

	if( sArgs.minimal_out_arg ) {
		ofstream	ofsm;

		ofsm.open( sArgs.minimal_out_arg, ios_base::binary );
		BNDefault.Save( ofsm );
		iSize = (uint32_t)vecBNs.size( );
		ofsm.write( (const char*)&iSize, sizeof(iSize) );
		for( i = 0; i < vecBNs.size( ); ++i )
			vecBNs[ i ].Save( ofsm ); }

	if( !open_contexts( sArgs.contexts_arg, vecBNs.size( ), Database.GetGenes( ), MatContexts ) )
		return 1;

	vecstrLine.resize( Database.GetGenes( ) );
	for( i = 0; i < vecstrLine.size( ); ++i )
		vecstrLine[ i ] = Database.GetGene( i );
	Dat.Open( vecstrLine );
	adGenes = new float[ Database.GetGenes( ) ];
	for( iGene = 0; ( iGene + 1 ) < Database.GetGenes( ); ++iGene ) {
		vector<unsigned char>	vecbData;

		if( !Database.Get( iGene, vecbData ) ) {
			cerr << "Could not retrieve gene: " << iGene << " (" << Database.GetGene( iGene ) << ")" << endl;
			return 1; }
		if( !vecBNs[ sArgs.context_arg - 1 ].Evaluate( vecbData, adGenes, Database.GetGenes( ), iGene + 1 ) ) {
			cerr << "Inference failed on gene: " << iGene << " (" << Database.GetGene( iGene ) << ")" << endl;
			return 1; }

		for( i = ( iGene + 1 ); i < Database.GetGenes( ); ++i )
			Dat.Set( iGene, i, adGenes[ i ] ); }
	delete[] adGenes;

	return 0; }

int main_datfile( const gengetopt_args_info& sArgs ) {
	static const size_t					c_iBuffer	= 1024;
	char								acBuffer[ c_iBuffer ];
	CDat								Dat;
	size_t								i, j, iGene, iIn, iOut;
	float								d, dIn, dOut;
	CCompactFullMatrix					MatContexts;
	vector<size_t>						veciGenes;
	map<string, size_t>					mapstriGenes;
	map<string, size_t>::const_iterator	iterGene;
	ifstream							ifsm;
	vector<string>						vecstrLine, vecstrGenes;

	if( sArgs.input_arg ) {
		if( !Dat.Open( sArgs.input_arg, !!sArgs.memmap_flag ) ) {
			cerr << "Could not open: " << sArgs.input_arg << endl;
			return 1; } }
	else {
		if( !Dat.Open( cin, CDat::EFormatText ) ) {
			cerr << "Could not open DAT" << endl;
			return 1; } }

	ifsm.open( sArgs.genes_arg );
	if( !ifsm.is_open( ) ) {
		cerr << "Could not open: " << sArgs.genes_arg << endl;
		return -1; }
	while( !ifsm.eof( ) ) {
		ifsm.getline( acBuffer, c_iBuffer - 1 );
		acBuffer[ c_iBuffer - 1 ] = 0;
		if( !acBuffer[ 0 ] )
			continue;
		vecstrLine.clear( );
		CMeta::Tokenize( acBuffer, vecstrLine );
		if( vecstrLine.size( ) != 2 ) {
			cerr << "Invalid line in " << sArgs.genes_arg << ": " << acBuffer << endl;
			return -1; }
		vecstrGenes.push_back( vecstrLine[ 1 ] );
		mapstriGenes[ vecstrLine[ 1 ] ] = atoi( vecstrLine[ 0 ].c_str( ) ); }
	ifsm.close( );

	if( !open_contexts( sArgs.contexts_arg, -1, Dat.GetGenes( ), MatContexts ) )
		return false;

	if( sArgs.lookup_arg ) {
		veciGenes.resize( Dat.GetGenes( ) );
		for( i = 0; i < Dat.GetGenes( ); ++i )
			veciGenes[ i ] = ( ( iterGene = mapstriGenes.find( Dat.GetGene( i ) ) ) == mapstriGenes.end( ) ) ?
				-1 : ( iterGene->second - 1 );
		if( ( iGene = Dat.GetGene( sArgs.lookup_arg ) ) == -1 ) {
			cerr << "Unknown gene: " << sArgs.lookup_arg << endl;
			return 1; }
		dIn = dOut = 0;
		for( iIn = iOut = i = 0; i < Dat.GetGenes( ); ++i ) {
			if( ( i == iGene ) || CMeta::IsNaN( d = Dat.Get( iGene, i ) ) || ( veciGenes[ i ] == -1 ) )
				continue;
			if( MatContexts.Get( sArgs.context_arg - 1, veciGenes[ i ] ) ) {
				iIn++;
				dIn += d; }
			iOut++;
			dOut += d; }
		cout << Dat.GetGene( iGene ) << '\t' << ( dIn / iIn ) << '\t' << ( dOut / iOut ) << endl; }
	else {
		veciGenes.resize( vecstrGenes.size( ) );
		for( i = 0; i < vecstrGenes.size( ); ++i )
			veciGenes[ i ] = Dat.GetGene( vecstrGenes[ i ] );
		for( i = 0; i < vecstrGenes.size( ); ++i ) {
			if( !( i % 100 ) )
				cerr << "Gene " << i << '/' << vecstrGenes.size( ) << endl;
			if( ( iGene = veciGenes[ i ] ) == -1 ) {
				cout << vecstrGenes[ i ] << '\t' << 0 << '\t' << 0 << endl;
				continue; }
			dIn = dOut = 0;
			for( iIn = iOut = j = 0; j < vecstrGenes.size( ); ++j ) {
				if( ( i == j ) || ( veciGenes[ j ] == -1 ) ||
					CMeta::IsNaN( d = Dat.Get( iGene, veciGenes[ j ] ) ) )
					continue;
				if( MatContexts.Get( sArgs.context_arg - 1, j ) ) {
					iIn++;
					dIn += d; }
				iOut++;
				dOut += d; }
			cout << vecstrGenes[ i ] << '\t' << ( dIn / iIn ) << '\t' << ( dOut / iOut ) << endl; } }

	return 0; }

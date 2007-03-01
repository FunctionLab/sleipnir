#include "stdafx.h"
#include "cmdline.h"

const char	c_szTxt[]	= "txt";
const char	c_szBin[]	= "bin";
const char	c_szDat[]	= "dat";

bool OpenBin( istream&, ostream& );
bool OpenText( istream&, ostream& );
bool OpenDats( const CDataPair&, const vector<string>&, ostream&, const CGenes&,
	const CGenes& );

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	ifstream			ifsm;
	vector<string>		vecstrDats;
	size_t				i;
	ofstream			ofsm;
	CGenome				Genome;
	CGenes				GenesEx( Genome ), GenesIn( Genome );
	CDataPair			Answers;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( sArgs.genex_arg ) {
		ifsm.open( sArgs.genex_arg );
		GenesEx.Open( ifsm );
		ifsm.close( ); }
	if( sArgs.genes_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.genes_arg );
		GenesIn.Open( ifsm );
		ifsm.close( ); }

	ifsm.clear( );
	if( !strcmp( sArgs.from_arg, c_szBin ) ) {
		if( sArgs.output_arg )
			ofsm.open( sArgs.output_arg );
		ifsm.open( sArgs.input_arg, ios_base::binary );
		if( !OpenBin( ifsm, sArgs.output_arg ? (ostream&)ofsm : cout ) ) {
			cerr << "Couldn't open: " << sArgs.input_arg << endl;
			return 1; }
		ifsm.close( ); }
	else if( !strcmp( sArgs.from_arg, c_szDat ) ) {
		ofsm.open( sArgs.output_arg, ios_base::binary );
		if( !Answers.Open( sArgs.answers_arg, true ) ) {
			cerr << "Couldn't open: " << sArgs.answers_arg << endl;
			return 1; }
		vecstrDats.resize( sArgs.inputs_num );
		for( i = 0; i < vecstrDats.size( ); ++i )
			vecstrDats[ i ] = sArgs.inputs[ i ];
		if( !OpenDats( Answers, vecstrDats, ofsm, GenesIn, GenesEx ) ) {
			cerr << "Couldn't open DAT files" << endl;
			return 1; } }
	else {
		ofsm.open( sArgs.output_arg, ios_base::binary );
		if( sArgs.input_arg )
			ifsm.open( sArgs.input_arg );
		if( !OpenText( sArgs.input_arg ? (istream&)ifsm : cin, ofsm ) ) {
			cerr << "Couldn't open: " << ( sArgs.input_arg ? sArgs.input_arg : "stdin" ) <<
				endl;
			return 1; }
		if( sArgs.input_arg )
			ifsm.close( ); }

	if( sArgs.output_arg )
		ofsm.close( );
	else
		cout.flush( );

	CMeta::Shutdown( );
	return 0; }

bool OpenDats( const CDataPair& Answers, const vector<string>& vecstrDATs, ostream& ostm,
	const CGenes& GenesIn, const CGenes& GenesEx ) {
	uint32_t		i, j, k, iPairs;
	CBinaryMatrix	Pairs;
	vector<size_t>	veciGenes;
	float			d;

	Pairs.Initialize( Answers.GetGenes( ) );
	for( i = 0; i < Pairs.GetSize( ); ++i )
		for( j = ( i + 1 ); j < Pairs.GetSize( ); ++j )
			Pairs.Set( i, j, false );

	veciGenes.resize( Answers.GetGenes( ) );
	for( i = 0; i < vecstrDATs.size( ); ++i ) {
		CDat	Dat;

		if( !Dat.Open( vecstrDATs[ i ].c_str( ) ) ) {
			cerr << "Couldn't open: " << vecstrDATs[ i ] << endl;
			return false; }
		cerr << "OpenDats( ) testing " << vecstrDATs[ i ] << endl;
		for( j = 0; j < veciGenes.size( ); ++j )
			veciGenes[ j ] = Dat.GetGene( Answers.GetGene( j ) );
		for( j = 0; j < Pairs.GetSize( ); ++j )
			if( veciGenes[ j ] != -1 )
				for( k = ( j + 1 ); k < Pairs.GetSize( ); ++k )
					if( ( veciGenes[ k ] != -1 ) && !Pairs.Get( j, k ) &&
						!CMeta::IsNaN( Dat.Get( veciGenes[ j ], veciGenes[ k ] ) ) )
						Pairs.Set( j, k, true ); }

	if( GenesEx.GetGenes( ) ) {
		for( i = 0; i < veciGenes.size( ); ++i )
			veciGenes[ i ] = GenesEx.IsGene( Answers.GetGene( i ) );
		for( i = 0; i < Pairs.GetSize( ); ++i ) {
			if( veciGenes[ i ] ) {
				for( j = ( i + 1 ); j < Pairs.GetSize( ); ++j )
					Pairs.Set( i, j, false );
				continue; }
			for( j = ( i + 1 ); j < Pairs.GetSize( ); ++j )
				if( veciGenes[ j ] )
					Pairs.Set( i, j, false ); } }
	if( GenesIn.GetGenes( ) ) {
		for( i = 0; i < veciGenes.size( ); ++i )
			veciGenes[ i ] = GenesIn.IsGene( Answers.GetGene( i ) );
		for( i = 0; i < Pairs.GetSize( ); ++i ) {
			if( !veciGenes[ i ] )
				for( j = ( i + 1 ); j < Pairs.GetSize( ); ++j )
					if( !veciGenes[ j ] )
						Pairs.Set( i, j, false ); } }

	cerr << "OpenDats( ) storing answers" << endl;
	k = 2 * sizeof(iPairs);
	ostm.seekp( k, ios_base::beg );
	for( iPairs = i = 0; i < Pairs.GetSize( ); ++i )
		for( j = ( i + 1 ); j < Pairs.GetSize( ); ++j )
			if( CMeta::IsNaN( d = Answers.Get( i, j ) ) )
				Pairs.Set( i, j, false );
			else if( Pairs.Get( i, j ) ) {
				d = d ? 1 : -1.0f;
				ostm.write( (char*)&d, sizeof(d) );
				ostm.seekp( (ostream::off_type)( vecstrDATs.size( ) * sizeof(float) ),
					ios_base::cur );
				ostm.write( (char*)&k, sizeof(k) );
				ostm.write( (char*)&i, sizeof(i) );
				ostm.write( (char*)&j, sizeof(j) );
				iPairs++; }

	ostm.seekp( 0, ios_base::beg );
	i = (uint32_t)vecstrDATs.size( );
	ostm.write( (char*)&i, sizeof(i) );
	ostm.write( (char*)&iPairs, sizeof(iPairs) );
	for( i = 0; i < vecstrDATs.size( ); ++i ) {
		CDat	Dat;

		if( !Dat.Open( vecstrDATs[ i ].c_str( ) ) ) {
			cerr << "Couldn't open: " << vecstrDATs[ i ] << endl;
			return false; }
		cerr << "OpenDats( ) storing " << vecstrDATs[ i ] << endl;
		for( j = 0; j < veciGenes.size( ); ++j )
			veciGenes[ j ] = Dat.GetGene( Answers.GetGene( j ) );
		ostm.seekp( sizeof(float) + ( 2 * sizeof(iPairs) ) + ( i * sizeof(float) ),
			ios_base::beg );
		for( j = 0; j < Pairs.GetSize( ); ++j )
			for( k = ( j + 1 ); k < Pairs.GetSize( ); ++k )
				if( Pairs.Get( j, k ) ) {
					d = ( ( veciGenes[ j ] == -1 ) || ( veciGenes[ k ] == -1 ) ||
						CMeta::IsNaN( d = Dat.Get( veciGenes[ j ], veciGenes[ k ] ) ) ) ?
						0 : ( 2 * d ) - 1;
					ostm.write( (char*)&d, sizeof(d) );
					ostm.seekp( (ostream::off_type)( ( 3 * sizeof(iPairs) ) +
						( vecstrDATs.size( ) * sizeof(float) ) ), ios_base::cur ); } }
	ostm.seekp( 0, ios_base::end );
	i = (uint32_t)Answers.GetGenes( );
	ostm.write( (char*)&i, sizeof(i) );
	for( j = i = 0; i < Answers.GetGenes( ); ++i ) {
		const string&	strGene	= Answers.GetGene( i );

		ostm.write( strGene.c_str( ), (streamsize)strGene.length( ) );
		ostm.write( (char*)&j, 1 ); }

	return true; }

bool OpenText( istream& istm, ostream& ostm ) {

	return false; }

bool OpenBin( istream& istm, ostream& ostm ) {
	static const size_t	c_iSize	= 512;
	char			sz[ c_iSize ];
	char*			pc;
	uint32_t		i, j, k, iWords, iDocs, iGenes;
	float*			ad;
	vector<string>	vecstrGenes;

	istm.read( (char*)&iWords, sizeof(iWords) );
	istm.read( (char*)&iDocs, sizeof(iDocs) );

	istm.seekg( (istream::off_type)( iDocs * ( ( ( iWords + 1 ) * sizeof(float) ) +
		( 3 * sizeof(iWords) ) ) ), ios_base::cur );
	istm.read( (char*)&iGenes, sizeof(iGenes) );
	vecstrGenes.resize( iGenes );
	for( i = 0; i < vecstrGenes.size( ); ++i ) {
		for( pc = sz; ; ++pc ) {
			istm.read( pc, 1 );
			if( !*pc )
				break; }
		vecstrGenes[ i ] = sz; }
	istm.seekg( 2 * sizeof(iWords), ios_base::beg );

	ad = new float[ iWords + 1 ];
	for( i = 0; i < iDocs; ++i ) {
		istm.read( (char*)ad, (streamsize)( iWords + 1 ) * sizeof(*ad) );
		cout << ad[ 0 ];
		for( j = 1; j <= iWords; ++j )
			cout << '\t' << (unsigned int)j << ':' << ad[ j ];
		istm.read( (char*)&j, sizeof(j) );
		if( j == ( 2 * sizeof(iWords) ) ) {
			istm.read( (char*)&j, sizeof(j) );
			istm.read( (char*)&k, sizeof(k) );
			cout << " # " << vecstrGenes[ j ] << '\t' << vecstrGenes[ k ]; }
		else if( j ) {
			istm.read( sz, (streamsize)j );
			sz[ j ] = 0;
			cout << " # " << sz; }
		cout << endl; }

	return true; }

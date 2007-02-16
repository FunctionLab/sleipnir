#include "stdafx.h"
#include "cmdline.h"

struct SFeature {
	string	m_strName;
	size_t	m_iDefault;

	SFeature( const string& strName, size_t iDefault ) : m_strName(strName), m_iDefault(iDefault) { }
};

struct SDatum {
	string				m_strName;
	map<string,size_t>	m_mapstriFeatures;
};

int main( int iArgs, char** aszArgs ) {
	static const size_t	c_iBuf	= 1024;
	gengetopt_args_info	sArgs;
	size_t				i, j;
	CPCL				DataOut;
	vector<SFeature>	vecsFeatures;
	vector<string>		vecstrLine, vecstrToken;
	ifstream			ifsm;
	char				szBuf[ c_iBuf ];
	map<string,SDatum>	mapValues;
	CGenome				Genome;
	CGenes				Genes( Genome );

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	ifsm.open( sArgs.features_arg );
	if( !ifsm.is_open( ) ) {
		cerr << "Could not open: " << sArgs.features_arg << endl;
		return 1; }
	while( ifsm.peek( ) != EOF ) {
		ifsm.getline( szBuf, c_iBuf - 1 );
		vecstrLine.clear( );
		CMeta::Tokenize( szBuf, vecstrLine );
		if( vecstrLine.size( ) != 0 )
			vecsFeatures.push_back( SFeature( vecstrLine[ 0 ], ( vecstrLine.size( ) == 1 ) ? -1 :
				atoi( vecstrLine[ 1 ].c_str( ) ) ) ); }

	ifsm.clear( );
	ifsm.open( sArgs.data_arg );
	if( !ifsm.is_open( ) ) {
		cerr << "Could not open: " << sArgs.data_arg << endl;
		return 1; }
	while( ifsm.peek( ) != EOF ) {
		ifsm.getline( szBuf, c_iBuf - 1 );
		vecstrLine.clear( );
		CMeta::Tokenize( szBuf, vecstrLine );
		if( vecstrLine.size( ) == 0 )
			continue;
		{
			SDatum&	sCur	= mapValues[ vecstrLine[ 0 ] ];

			sCur.m_strName = vecstrLine[ 0 ];
			for( i = 1; i < vecstrLine.size( ); ++i ) {
				vecstrToken.clear( );
				CMeta::Tokenize( vecstrLine[ i ].c_str( ), vecstrToken, "|" );
				if( vecstrToken.size( ) != 2 ) {
					cerr << "Illegal token in " << sArgs.data_arg << ": " << szBuf << endl;
					return 1; }
				sCur.m_mapstriFeatures[ vecstrToken[ 0 ] ] = atoi( vecstrToken[ 1 ].c_str( ) ); } } }

	if( sArgs.input_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.input_arg ); }
	if( !Genes.Open( sArgs.input_arg ? ifsm : cin ) ) {
		cerr << "Could not open: " << ( sArgs.input_arg ? sArgs.input_arg : "input genes" ) << endl;
		return 1; }
	if( sArgs.input_arg )
		ifsm.close( );

	cout << '#';
	for( i = 0; i < vecsFeatures.size( ); ++i )
		cout << '\t' << vecsFeatures[ i ].m_strName;
	cout << endl;
	for( i = 0; i < sArgs.inputs_num; ++i ) {
		const SDatum&	sDatum	= mapValues[ CMeta::Basename( sArgs.inputs[ i ] ) ];

		ifsm.clear( );
		ifsm.open( sArgs.inputs[ i ] );
		if( !ifsm.is_open( ) ) {
			cerr << "Could not open: " << sArgs.inputs[ i ] << endl;
			return 1; }
		while( ifsm.peek( ) != EOF ) {
			map<string,size_t>	mapstriCur;

			ifsm.getline( szBuf, c_iBuf - 1 );
			vecstrLine.clear( );
			CMeta::Tokenize( szBuf, vecstrLine );
			if( vecstrLine.size( ) == 0 )
				continue;
			if( vecstrLine.size( ) < 2 ) {
				cerr << "Illegal line in " << sArgs.inputs[ i ] << ": " << szBuf << endl;
				return 1; }

			for( j = 2; j < vecstrLine.size( ); ++j ) {
				vecstrToken.clear( );
				CMeta::Tokenize( vecstrLine[ j ].c_str( ), vecstrToken, "|" );
				if( vecstrToken.size( ) != 2 ) {
					cerr << "Illegal token in " << sArgs.data_arg << ": " << szBuf << endl;
					return 1; }
				mapstriCur[ vecstrToken[ 0 ] ] = atoi( vecstrToken[ 1 ].c_str( ) ); }

			cout << ( Genes.IsGene( vecstrLine[ 0 ] ) ? 1 : -1 );
			for( j = 0; j < vecsFeatures.size( ); ++j ) {
				const SFeature&						sFeature	= vecsFeatures[ j ];
				map<string,size_t>::const_iterator	iterCur;
				size_t								iCur;

				if( ( iterCur = mapstriCur.find( sFeature.m_strName ) ) != mapstriCur.end( ) )
					iCur = iterCur->second;
				else if( ( iterCur = sDatum.m_mapstriFeatures.find( sFeature.m_strName ) ) !=
					sDatum.m_mapstriFeatures.end( ) )
					iCur = iterCur->second;
				else
					iCur = sFeature.m_iDefault;
				cout << '\t';
				if( iCur != -1 )
					cout << iCur; }
			cout << "	#	" << vecstrLine[ 0 ] << '\t' << sArgs.inputs[ i ];
			cout << endl; } }

	CMeta::Shutdown( );
	return 0; }

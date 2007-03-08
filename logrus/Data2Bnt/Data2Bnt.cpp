#include "stdafx.h"
#include "cmdline.h"

static const char	c_szERROR[]	= "ERROR";
static const char	c_cComment	= '#';
static const char	c_cDot		= '.';

struct SFeature {
	string			m_strName;
	vector<string>	m_vecstrValues;
	size_t			m_iDefault;

	SFeature( const string& strName ) : m_strName(strName), m_iDefault(-1) { }

	size_t quantize( const string& strValue ) const {
		size_t	i;

		for( i = 0; i < m_vecstrValues.size( ); ++i )
			if( strValue == m_vecstrValues[ i ] )
				return i;

		return -1; }
};

struct SDatum {
	string				m_strName;
	map<size_t,size_t>	m_mapiiFeatures;
};

int main( int iArgs, char** aszArgs ) {
	static const size_t	c_iBuf	= 1024;
	gengetopt_args_info	sArgs;
	size_t				i, j, iFeature;
	CPCL				DataOut;
	vector<SFeature>	vecsFeatures;
	vector<string>		vecstrLine, vecstrToken;
	ifstream			ifsm;
	char				szBuf[ c_iBuf ];
	map<string,SDatum>	mapValues;
	CGenome				Genome;
	CGenes				Genes( Genome );
	string				strName;
	vector<float>		vecdQuants;
	float				d;

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
		if( ( vecstrLine.size( ) != 0 ) && ( vecstrLine[ 0 ][ 0 ] != c_cComment ) ) {
			vecsFeatures.push_back( SFeature( vecstrLine[ 0 ] ) );
			{
				SFeature&	sFeature	= vecsFeatures[ vecsFeatures.size( ) - 1 ];

				vecstrToken.clear( );
				CMeta::Tokenize( vecstrLine[ 1 ].c_str( ), vecstrToken, "|" );
				sFeature.m_vecstrValues.resize( vecstrToken.size( ) );
				copy( vecstrToken.begin( ), vecstrToken.end( ), sFeature.m_vecstrValues.begin( ) );

				if( vecstrLine.size( ) > 2 )
					sFeature.m_iDefault = sFeature.quantize( vecstrLine[ 2 ] );
			} } }
	ifsm.close( );

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
				for( j = 0; j < vecsFeatures.size( ); ++j )
					if( vecstrToken[ 0 ] == vecsFeatures[ j ].m_strName )
						break;
				if( j >= vecsFeatures.size( ) ) {
					cerr << "Unknown feature: " << vecstrLine[ i ] << endl;
					return 1; }
				sCur.m_mapiiFeatures[ j ] = vecsFeatures[ j ].quantize( vecstrToken[ 1 ] ); }
		} }
	ifsm.close( );

	ifsm.clear( );
	ifsm.open( sArgs.quants_arg );
	if( !ifsm.is_open( ) ) {
		cerr << "Could not open: " << sArgs.quants_arg << endl;
		return 1; }
	ifsm.getline( szBuf, c_iBuf - 1 );
	vecstrLine.clear( );
	CMeta::Tokenize( szBuf, vecstrLine );
	vecdQuants.resize( vecstrLine.size( ) );
	for( i = 0; i < vecdQuants.size( ); ++i )
		vecdQuants[ i ] = (float)atof( vecstrLine[ i ].c_str( ) );
	ifsm.close( );

	if( sArgs.input_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.input_arg ); }
	if( !Genes.Open( sArgs.input_arg ? ifsm : cin ) ) {
		cerr << "Could not open: " << ( sArgs.input_arg ? sArgs.input_arg : "input genes" ) << endl;
		return 1; }
	if( sArgs.input_arg )
		ifsm.close( );

	if( sArgs.comments_flag ) {
		cout << '#';
		for( i = 1; i < vecsFeatures.size( ); ++i )
			cout << '\t' << vecsFeatures[ i ].m_strName;
		cout << endl; }
	for( i = 0; i < sArgs.inputs_num; ++i ) {
		strName = CMeta::Basename( sArgs.inputs[ i ] );
		if( ( j = strName.rfind( c_cDot ) ) != string::npos )
			strName = strName.substr( 0, j );
		const SDatum&	sDatum	= mapValues[ strName ];

		ifsm.clear( );
		ifsm.open( sArgs.inputs[ i ] );
		if( !ifsm.is_open( ) ) {
			cerr << "Could not open: " << sArgs.inputs[ i ] << endl;
			return 1; }
		while( ifsm.peek( ) != EOF ) {
			map<size_t,size_t>	mapiiCur;

			ifsm.getline( szBuf, c_iBuf - 1 );
			vecstrLine.clear( );
			CMeta::Tokenize( szBuf, vecstrLine );
			if( ( vecstrLine.size( ) == 0 ) || ( vecstrLine[ 0 ].find( c_szERROR ) == 0 ) )
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
				for( iFeature = 0; iFeature < vecsFeatures.size( ); ++iFeature )
					if( vecstrToken[ 0 ] == vecsFeatures[ iFeature ].m_strName )
						break;
				if( iFeature >= vecsFeatures.size( ) ) {
					cerr << "Unknown feature: " << vecstrLine[ j ] << endl;
					return 1; }
				mapiiCur[ iFeature ] = vecsFeatures[ iFeature ].quantize( vecstrToken[ 1 ] ); }

			cout << ( Genes.IsGene( vecstrLine[ 0 ] ) ? 2 : 1 ) << '\t';
			d = (float)atof( vecstrLine[ 1 ].c_str( ) );
			for( j = 0; j < vecdQuants.size( ); ++j )
				if( d <= vecdQuants[ j ] )
					break;
			if( j == vecdQuants.size( ) )
				j--;
			cout << ( j + 1 );
			for( j = 2; j < vecsFeatures.size( ); ++j ) {
				const SFeature&						sFeature	= vecsFeatures[ j ];
				map<size_t,size_t>::const_iterator	iterCur;
				size_t								iCur;

				if( ( iterCur = mapiiCur.find( j ) ) != mapiiCur.end( ) )
					iCur = iterCur->second;
				else if( ( iterCur = sDatum.m_mapiiFeatures.find( j ) ) != sDatum.m_mapiiFeatures.end( ) )
					iCur = iterCur->second;
				else
					iCur = sFeature.m_iDefault;
				cout << '\t';
				if( iCur != -1 )
					cout << ( iCur + 1 ); }
			if( sArgs.comments_flag )
				cout << "	#	" << vecstrLine[ 0 ] << '\t' << sArgs.inputs[ i ];
			cout << endl; }
		ifsm.close( ); }

	CMeta::Shutdown( );
	return 0; }

#include "stdafx.h"
#include "cmdline.h"

static const char	c_szERROR[]		= "ERROR";
static const char	c_szAnnotated[]	= "annotated";
static const char	c_szValue[]		= "value";
static const char	c_cComment		= '#';
static const char	c_cDot			= '.';

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

void output_row( const string&, const CGenes&, size_t, const SDatum&, const vector<SFeature>&,
	const map<size_t,size_t>&, bool, bool, bool );
size_t get_feature( const string&, const vector<SFeature>& );

int main( int iArgs, char** aszArgs ) {
	static const size_t	c_iBuf	= 1024;
	gengetopt_args_info	sArgs;
	size_t				i, j, k, iFeature;
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
	map<size_t,size_t>	mapiiCur;

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
			char*	pc;

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
				k = strtol( vecstrToken[ 1 ].c_str( ), &pc, 10 );
				sCur.m_mapiiFeatures[ j ] = ( pc != vecstrToken[ 1 ].c_str( )) ? k :
					vecsFeatures[ j ].quantize( vecstrToken[ 1 ] ); }
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

	if( sArgs.genome_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.genome_arg );
		if( !Genome.Open( ifsm ) ) {
			cerr << "Could not open: " << sArgs.genome_arg << endl;
			return 1; }
		ifsm.close( ); }

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
		CDat	Dat;

		cerr << "Processing: " << sArgs.inputs[ i ] << endl;
		strName = CMeta::Basename( sArgs.inputs[ i ] );
		if( ( j = strName.rfind( c_cDot ) ) != string::npos )
			strName = strName.substr( 0, j );
		const SDatum&	sDatum	= mapValues[ strName ];

		if( Dat.Open( sArgs.inputs[ i ] ) ) {
			size_t								iValue, iDefault, iJ;
			map<size_t,size_t>::const_iterator	iterDefault;

			Dat.Normalize( false );
			for( j = 0; j < Dat.GetGenes( ); ++j ) {
				iJ = get_feature( Dat.GetGene( j ), vecsFeatures );
				for( k = ( j + 1 ); k < Dat.GetGenes( ); ++k ) {
					if( CMeta::IsNaN( d = Dat.Get( j, k ) ) )
						continue;
					iValue = CMeta::Quantize( d, vecdQuants );
					mapiiCur.clear( );
					mapiiCur[ get_feature( Dat.GetGene( k ), vecsFeatures ) ] = 1;
					output_row( Dat.GetGene( j ), Genes, iValue, sDatum, vecsFeatures, mapiiCur,
						!!sArgs.comments_flag, !!sArgs.sparse_flag, false );
					mapiiCur.clear( );
					mapiiCur[ iJ ] = 1;
					output_row( Dat.GetGene( k ), Genes, iValue, sDatum, vecsFeatures, mapiiCur,
						!!sArgs.comments_flag, !!sArgs.sparse_flag, false ); } }
			for( iDefault = -1,j = 0; j < vecsFeatures.size( ); ++j )
				if( vecsFeatures[ j ].m_strName == c_szValue ) {
					iDefault = ( ( iterDefault = sDatum.m_mapiiFeatures.find( j ) ) ==
						sDatum.m_mapiiFeatures.end( ) ) ? vecsFeatures[ j ].m_iDefault : iterDefault->second;
					break; }
			if( iDefault != -1 ) {
				vector<size_t>	veciGenes;
				size_t			iOne, iTwo;

				veciGenes.resize( sArgs.genome_arg ? Genome.GetGenes( ) : Dat.GetGenes( ) );
				for( j = 0; j < veciGenes.size( ); ++j )
					veciGenes[ j ] = sArgs.genome_arg ? Dat.GetGene( Genome.GetGene( j ).GetName( ) ) : j;
				for( j = 0; j < veciGenes.size( ); ++j ) {
					if( ( iOne = veciGenes[ j ] ) == -1 )
						continue;
					iJ = get_feature( Dat.GetGene( iOne ), vecsFeatures );
					for( k = ( j + 1 ); k < veciGenes.size( ); ++k ) {
						if( ( ( iTwo = veciGenes[ k ] ) == -1 ) || !CMeta::IsNaN( Dat.Get( iOne, iTwo ) ) ||
							( ( (float)rand( ) / RAND_MAX ) > sArgs.fraction_arg ) )
							continue;
						mapiiCur.clear( );
						mapiiCur[ get_feature( Dat.GetGene( iTwo ), vecsFeatures ) ] = 1;
						output_row( Dat.GetGene( iOne ), Genes, iDefault, sDatum, vecsFeatures, mapiiCur,
							!!sArgs.comments_flag, !!sArgs.sparse_flag, true );
						mapiiCur.clear( );
						mapiiCur[ iJ ] = 1;
						output_row( Dat.GetGene( iTwo ), Genes, iDefault, sDatum, vecsFeatures, mapiiCur,
							!!sArgs.comments_flag, !!sArgs.sparse_flag, true ); } } } }
		else {
			ifsm.clear( );
			ifsm.open( sArgs.inputs[ i ] );
			if( !ifsm.is_open( ) ) {
				cerr << "Could not open: " << sArgs.inputs[ i ] << endl;
				return 1; }
			while( ifsm.peek( ) != EOF ) {
				mapiiCur.clear( );
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

				output_row( vecstrLine[ 0 ], Genes, CMeta::Quantize( (float)atof( vecstrLine[ 1 ].c_str( ) ),
					vecdQuants ), sDatum, vecsFeatures, mapiiCur, !!sArgs.comments_flag, !!sArgs.sparse_flag,
					false ); }
			ifsm.close( ); } }

	CMeta::Shutdown( );
	return 0; }

void output_row( const string& strGene, const CGenes& Genes, size_t iValue, const SDatum& sDatum,
	const vector<SFeature>& vecsFeatures, const map<size_t,size_t>& mapiiCur, bool fComments,
	bool fSparse, bool fDefault ) {
	size_t	i;

	if( fSparse ) {
		if( Genes.IsGene( strGene ) )
			cout << get_feature( c_szAnnotated, vecsFeatures ) << "|1	";
		if( iValue )
			cout << get_feature( c_szValue, vecsFeatures ) << '|' << iValue << '\t'; }
	else
		cout << ( Genes.IsGene( strGene ) ? 2 : 1 ) << '\t' << ( 1 + iValue );
	for( i = 2; i < vecsFeatures.size( ); ++i ) {
		const SFeature&						sFeature	= vecsFeatures[ i ];
		map<size_t,size_t>::const_iterator	iterCur;
		size_t								iCur;

		if( ( iterCur = mapiiCur.find( i ) ) != mapiiCur.end( ) )
			iCur = iterCur->second;
		else if( ( iterCur = sDatum.m_mapiiFeatures.find( i ) ) != sDatum.m_mapiiFeatures.end( ) )
			iCur = iterCur->second;
		else
			iCur = sFeature.m_iDefault;
		if( fSparse ) {
			if( ( iCur != -1 ) && ( iCur > 0 ) )
				cout << i << '|' << iCur << '\t'; }
		else {
			cout << '\t';
			if( iCur != -1 )
				cout << ( iCur + 1 ); } }
	if( fComments ) {
		cout << "	#	" << strGene << '\t' << sDatum.m_strName;
		if( fDefault )
			cout << "	default"; }
	cout << endl; }

size_t get_feature( const string& strName, const vector<SFeature>& vecsFeatures ) {
	size_t	i;

	for( i = 0; i < vecsFeatures.size( ); ++i )
		if( vecsFeatures[ i ].m_strName == strName )
			return i;

	return -1; }

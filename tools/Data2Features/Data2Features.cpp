#include "stdafx.h"
#include "cmdline.h"

static const char	c_szERROR[]		= "ERROR";
static const char	c_szAnnotated[]	= "annotated";
static const char	c_szValue[]		= "value";
static const char	c_cComment		= '#';
static const char	c_cDot			= '.';
static const size_t	c_iBuf			= 1024;

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
	gengetopt_args_info	sArgs;
	size_t				i, j, k;
	vector<SFeature>	vecsFeatures;
	vector<string>		vecstrLine, vecstrToken;
	ifstream			ifsm;
	char				szBuf[ c_iBuf ];
	map<string,SDatum>	mapValues;
	CGenome				Genome;
	CGenes				GenesPos( Genome ), GenesNeg( Genome );
//	vector<float>		vecdQuants;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	ifsm.open( sArgs.environment_arg );
	if( !ifsm.is_open( ) ) {
		cerr << "Could not open: " << sArgs.environment_arg << endl;
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

/*
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
*/

	if( sArgs.genome_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.genome_arg );
		if( !Genome.Open( ifsm ) ) {
			cerr << "Could not open: " << sArgs.genome_arg << endl;
			return 1; }
		ifsm.close( ); }

	if( sArgs.positives_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.positives_arg ); }
	if( !GenesPos.Open( sArgs.positives_arg ? ifsm : cin ) ) {
		cerr << "Could not open: " << ( sArgs.positives_arg ? sArgs.positives_arg : "input genes" ) << endl;
		return 1; }
	if( sArgs.positives_arg )
		ifsm.close( );
	ifsm.clear( );
	ifsm.open( sArgs.negatives_arg );
	if( !GenesNeg.Open( ifsm ) ) {
		cerr << "Could not open: " << sArgs.negatives_arg << endl;
		return 1; }
	ifsm.close( );

	cout << "<?xml version='1.0' encoding='utf-8'?>" << endl;
	cout << "<dataset name='datasets'>" << endl;
	cout << "  <header>" << endl;
	cout << "    <attributes>" << endl;
/*
	cout << "      <attribute name='annotated' type='nominal'>" << endl;
	cout << "        <labels>" << endl;
	cout << "          <label>1</label>" << endl;
	cout << "          <label>2</label>" << endl;
	cout << "        </labels>" << endl;
	cout << "      </attribute>" << endl;
*/
	cout << "      <attribute name='value' type='numeric'/>" << endl;
	for( i = 0; i < vecsFeatures.size( ); ++i ) {
		cout << "      <attribute name='" << vecsFeatures[ i ].m_strName << "' type='nominal'>" << endl;
		cout << "        <labels>" << endl;
		for( j = 0; j < vecsFeatures[ i ].m_vecstrValues.size( ); ++j )
			cout << "          <label>" << ( j + 1 ) << "</label>" << endl;
		cout << "        </labels>" << endl;
		cout << "      </attribute>" << endl; }
/*
	for( i = 0; i < Genome.GetGenes( ); ++i ) {
		cout << "      <attribute name='" << Genome.GetGene( i ).GetName( ) << "' type='nominal'>" << endl;
		cout << "        <labels>" << endl;
		cout << "          <label>1</label>" << endl;
		cout << "          <label>2</label>" << endl;
		cout << "        </labels>" << endl;
		cout << "      </attribute>" << endl; }
*/
	cout << "    </attributes>" << endl;
	cout << "  </header>" << endl;

	cout << "  <body>" << endl;
	cout << "    <instances>" << endl;
	for( i = 0; i < sArgs.inputs_num; ++i ) {
		CPCL				PCL;
		CDat				Dat;
		size_t				iCount;
		vector<size_t>		veciGenes;
		float				d, dAverage;
		map<size_t,size_t>	mapiiCur;
		int					iRet;
		string				strName;

		cerr << "Processing: " << sArgs.inputs[ i ] << endl;

		strName = CMeta::Basename( sArgs.inputs[ i ] );
		while( ( j = strName.rfind( c_cDot ) ) != string::npos )
			strName = strName.substr( 0, j );
		const SDatum&	sDatum	= mapValues[ strName ];

		if( !Dat.Open( sArgs.inputs[ i ], !!sArgs.memmap_flag ) && ( iRet = CPCL::Distance( sArgs.inputs[ i ],
			sArgs.skip_arg, sArgs.distance_arg, !!sArgs.normalize_flag, !!sArgs.zscore_flag, false,
			sArgs.genome_arg, CMeta::GetNaN( ), -1, PCL, Dat ) ) ) {
			cerr << "Could not open: " << sArgs.inputs[ i ] << endl;
			cmdline_parser_print_help( );
			return iRet; }

		veciGenes.resize( GenesPos.GetGenes( ) );
		for( j = 0; j < veciGenes.size( ); ++j )
			veciGenes[ j ] = Dat.GetGene( GenesPos.GetGene( j ).GetName( ) );
		dAverage = 0;
		for( iCount = j = 0; j < veciGenes.size( ); ++j ) {
			if( veciGenes[ j ] == -1 )
				continue;
			for( k = ( j + 1 ); k < veciGenes.size( ); ++k )
				if( ( veciGenes[ k ] != -1 ) &&
					!CMeta::IsNaN( d = Dat.Get( veciGenes[ j ], veciGenes[ k ] ) ) ) {
					iCount++;
					dAverage += d; } }
		if( iCount )
			dAverage /= iCount;
/*
		adCentroid = new float[ PCL.GetExperiments( ) ];
		for( iCount = i = 0; i < GenesPos.GetGenes( ); ++i ) {
			if( ( j = PCL.GetGene( GenesPos.GetGene( i ).GetName( ) ) ) == -1 )
				continue;
			iCount++;
			for( k = 0; k < PCL.GetExperiments( ); ++k )
				adCentroid[ k ] += PCL.Get( j, k ); }
		for( i = 0; i < PCL.GetExperiments( ); ++i )
			adCentroid[ i ] /= iCount;
*/

		cout << "      <instance>" << endl;
		cout << "        <value>" << dAverage << "</value>" << endl;
		for( j = 0; j < vecsFeatures.size( ); ++j ) {
			const SFeature&						sFeature	= vecsFeatures[ j ];
			map<size_t,size_t>::const_iterator	iterCur;
			size_t								iCur;

			if( ( iterCur = sDatum.m_mapiiFeatures.find( j ) ) != sDatum.m_mapiiFeatures.end( ) )
				iCur = iterCur->second;
			else
				iCur = sFeature.m_iDefault;
			cout << "        <value>" << ( iCur + 1 ) << "</value>" << endl; }
		cout << "      </instance>" << endl; }
	cout << "    </instances>" << endl;
	cout << "  </body>" << endl;
	cout << "</dataset>" << endl;

	CMeta::Shutdown( );
	return 0; }

#include "stdafx.h"
#include "cmdline.h"
#include "BNServer.h"
#include "dot.h"

static const char				c_szXDSL[]						= ".xdsl";
static const char				c_szDSL[]						= ".dsl";
static const char				c_szDOT[]						= ".dot";
static const char				c_szSVG[]						= ".svg";
const float						CBNServer::c_dCutoff			= 0.1f;
const float						CBNServer::c_adColorMin[]		= {0, 1, 0};
const float						CBNServer::c_adColorMax[]		= {1, 0, 0};
const CBNServer::TPFNProcessor	CBNServer::c_apfnProcessors[]	=
	{&CBNServer::ProcessInference, &CBNServer::ProcessData, &CBNServer::ProcessGraph,
	&CBNServer::ProcessContexts};

struct SPixie {
	size_t	m_iNode;
	float	m_dScore;

	SPixie( size_t iNode, float dScore ) : m_iNode(iNode), m_dScore(dScore) { }

	bool operator<( const SPixie& sPixie ) const {

		return ( m_dScore < sPixie.m_dScore ); }
};

int main( int iArgs, char** aszArgs ) {
	static const size_t					c_iBuffer	= 1024;
	gengetopt_args_info					sArgs;
	CServer								Server;
	CBayesNetSmile						BNSmile;
	CBayesNetMinimal					BNDefault;
	vector<CBayesNetMinimal>			vecBNs;
	ifstream							ifsm;
	istream*							pistm;
	char								acBuffer[ c_iBuffer ];
	vector<string>						vecstrLine;
	map<size_t, string>					mapistrBNs;
	map<size_t, string>::const_iterator	iterBN;
	size_t								i, iMax;
	CDatabase							Database;
	uint32_t							iSize;
	CCompactFullMatrix					MatContexts;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );
	EnableXdslFormat( );

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

	ifsm.clear( );
	ifsm.open( sArgs.contexts_arg );
	if( !ifsm.is_open( ) ) {
		cerr << "Could not open: " << sArgs.contexts_arg << endl;
		return 1; }
	MatContexts.Initialize( vecBNs.size( ), Database.GetGenes( ), 2, true );
	while( !ifsm.eof( ) ) {
		size_t	iContext, iGene;

		ifsm.getline( acBuffer, c_iBuffer - 1 );
		acBuffer[ c_iBuffer - 1 ] = 0;
		if( !acBuffer[ 0 ] )
			continue;
		vecstrLine.clear( );
		CMeta::Tokenize( acBuffer, vecstrLine );
		if( vecstrLine.size( ) != 2 ) {
			cerr << "Invalid line in " << sArgs.contexts_arg << ": " << acBuffer << endl;
			return 1; }
		iContext = atoi( vecstrLine[ 0 ].c_str( ) ) - 1;
		iGene = atoi( vecstrLine[ 1 ].c_str( ) ) - 1;
		if( iContext >= MatContexts.GetRows( ) ) {
			cerr << "Invalid context on line: " << acBuffer << endl;
			return 1; }
		if( iGene >= MatContexts.GetColumns( ) ) {
			cerr << "Invalid gene on line: " << acBuffer << endl;
			return 1; }
		MatContexts.Set( iContext, iGene, 1 ); }
	ifsm.close( );

	CBNServer	BNServer( BNDefault, vecBNs, MatContexts, 0, Database, "", sArgs.files_arg,
		sArgs.graphviz_arg );

	Server.Initialize( sArgs.port_arg, sArgs.timeout_arg, &BNServer );
#ifdef WIN32
	pthread_win32_process_attach_np( );
#endif // WIN32
	Server.Start( );
#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32

	CMeta::Shutdown( );
	return 0; }

CBNServer::CBNServer( const CBayesNetMinimal& BNDefault, const vector<CBayesNetMinimal>& vecBNs,
	const CCompactFullMatrix& MatContexts, SOCKET iSocket, const CDatabase& Database,
	const string& strConnection, const char* szFiles, const char* szGraphviz ) : m_BNDefault(BNDefault),
	m_vecBNs(vecBNs), m_MatContexts(MatContexts), m_iSocket(iSocket), m_Database(Database),
	m_iGenes(Database.GetGenes( )), m_strConnection(strConnection), m_adGenes(NULL), m_adContexts(NULL),
	m_strGraphviz(szGraphviz), m_strFiles(szFiles) {

	if( m_strConnection.length( ) > 0 )
		cerr << "New connection from: " << m_strConnection << endl; }

CBNServer::~CBNServer( ) {

	if( m_adGenes )
		delete[] m_adGenes;
	if( m_adContexts )
		delete[] m_adContexts; }

IServerClient* CBNServer::NewInstance( SOCKET iSocket, uint32_t iHost, uint16_t sPort,
	const CPropertyFile* pConfig ) {
	string	strConnection;
	char	acBuffer[ 16 ];
	in_addr	sAddr;

#pragma warning(disable : 4996)
	sprintf( acBuffer, "%hu", sPort );
#pragma warning(default : 4996)
	sAddr.s_addr = htonl( iHost );
	strConnection = (string)inet_ntoa( sAddr ) + ":" + acBuffer;
	return new CBNServer( m_BNDefault, m_vecBNs, m_MatContexts, iSocket, m_Database, strConnection,
		m_strFiles.c_str( ), m_strGraphviz.c_str( ) ); }

void CBNServer::Destroy( ) {

	cerr << "Disconnected: " << m_strConnection << endl;

	delete this; }

bool CBNServer::ProcessMessage( const vector<unsigned char>& vecbMessage ) {
	size_t	iProcessed, iOffset;

	for( iOffset = 0; iOffset < vecbMessage.size( ); iOffset += ( iProcessed + 1 ) ) {
		if( vecbMessage[ iOffset ] >= ARRAYSIZE(c_apfnProcessors) ) {
			cerr << m_strConnection << " unknown opcode: " << (int)vecbMessage[ iOffset ] << endl;
			return false; }
		if( ( iProcessed = (this->*c_apfnProcessors[ vecbMessage[ iOffset ] ])( vecbMessage,
			iOffset + 1 ) ) == -1 )
			return false; }

	return true; }
	
size_t CBNServer::ProcessInference( const vector<unsigned char>& vecbMessage, size_t iOffset ) {
	size_t		iStart;
	uint32_t	iGene, iContext;

	if( ( iOffset + sizeof(iContext) ) > vecbMessage.size( ) )
		return -1;
	iStart = iOffset;
	iContext = *(uint32_t*)&vecbMessage[ iOffset ];
	for( iOffset += sizeof(iContext); ( iOffset + sizeof(iGene) ) <= vecbMessage.size( );
		iOffset += sizeof(iGene) ) {
		iGene = *(uint32_t*)&vecbMessage[ iOffset ];
		if( !Get( iGene, iContext ) )
			return -1; }

	return ( iOffset - iStart ); }

bool CBNServer::Get( size_t iGene, size_t iContext, float* adValues ) {
	const CBayesNetMinimal&	BNet		= ( iContext && m_vecBNs.size( ) ) ?
											m_vecBNs[ ( iContext - 1 ) % m_vecBNs.size( ) ] : m_BNDefault;
	vector<unsigned char>	vecbData;
	uint32_t				iSize;
	float*					adTarget;

	cerr << m_strConnection << " inferring " << iGene  << " (" << m_Database.GetGene( iGene - 1 ) << ") in " <<
		iContext << " (" << BNet.GetID( ) << ")" << endl;
	if( ( iGene < 1 ) || !m_Database.Get( iGene - 1, vecbData ) )
		return false;
	if( !( adTarget = adValues ) ) {
		if( !m_adGenes )
			m_adGenes = new float[ m_iGenes ];
		adTarget = m_adGenes; }
	if( !BNet.Evaluate( vecbData, adTarget, m_iGenes ) )
		return false;
	adTarget[ iGene - 1 ] = CMeta::GetNaN( );

	if( !adValues ) {
		iSize = (uint32_t)( m_iGenes * sizeof(*m_adGenes) );
		send( m_iSocket, (char*)&iSize, sizeof(iSize), 0 );
		send( m_iSocket, (char*)m_adGenes, iSize, 0 ); }

	return true; }

bool CBNServer::Get( size_t iGene, const vector<size_t>& veciGenes, size_t iContext, float* adValues ) {
	const CBayesNetMinimal&	BNet		= ( iContext && m_vecBNs.size( ) ) ?
											m_vecBNs[ iContext % m_vecBNs.size( ) ] : m_BNDefault;
	vector<unsigned char>	vecbData;

	return ( ( iGene >= 1 ) && m_Database.Get( iGene - 1, veciGenes, vecbData ) &&
		BNet.Evaluate( vecbData, adValues, veciGenes.size( ) ) ); }

size_t CBNServer::ProcessData( const vector<unsigned char>& vecbMessage, size_t iOffset ) {
	uint32_t				iSize, iOne, iTwo;
	vector<unsigned char>	vecbData, vecbSend;
	size_t					i;
	unsigned char			b, bValue;

	iSize = sizeof(iOne) + sizeof(iTwo);
	if( ( iOffset + iSize ) > vecbMessage.size( ) )
		return -1;
	iOne = *(uint32_t*)&vecbMessage[ iOffset ];
	iTwo = *(uint32_t*)&vecbMessage[ iOffset + sizeof(iOne) ];
	cerr << m_strConnection << " requested " << iOne << ',' << iTwo << endl;
	if( !m_Database.Get( iOne - 1, iTwo - 1, vecbData ) )
		return -1;

	vecbSend.resize( vecbData.size( ) * 2 );
	for( i = 0; i < vecbData.size( ); ++i ) {
		b = vecbData[ i ];
		if( ( bValue = ( b & 0xF ) ) == 0xF )
			bValue = -1;
		vecbSend[ 2 * i ] = bValue;
		if( ( bValue = ( ( b >> 4 ) & 0xF ) ) == 0xF )
			bValue = -1;
		vecbSend[ ( 2 * i ) + 1 ] = bValue; }
	iOne = (uint32_t)vecbSend.size( );
	send( m_iSocket, (char*)&iOne, sizeof(iOne), 0 );
	send( m_iSocket, (char*)&vecbSend[ 0 ], (int)vecbSend.size( ), 0 );

	return iSize; }

size_t CBNServer::ProcessGraph( const vector<unsigned char>& vecbMessage, size_t iOffset ) {
	size_t						iRet, i, iSize;
	uint32_t					iGene, iContext, iLimit;
	vector<size_t>				veciQuery;
	set<size_t>					setiQuery;
	CDat						DatGraph;
	vector<bool>				vecfQuery;
	set<size_t>::const_iterator	iterQuery;

	iSize = sizeof(iContext) + sizeof(iLimit);
	if( ( iOffset + iSize ) > vecbMessage.size( ) )
		return -1;
	iContext = *(uint32_t*)&vecbMessage[ iOffset ];
	iLimit = *(uint32_t*)&vecbMessage[ iOffset + sizeof(iContext) ];
	for( i = iOffset + iSize; ( i + sizeof(iGene) ) <= vecbMessage.size( ); i += sizeof(iGene) ) {
		iGene = *(uint32_t*)&vecbMessage[ i ];
		setiQuery.insert( iGene ); }
	iRet = i - iOffset;

	veciQuery.reserve( setiQuery.size( ) );
	for( iterQuery = setiQuery.begin( ); iterQuery != setiQuery.end( ); ++iterQuery )
		veciQuery.push_back( *iterQuery );
	if( !( PixieCreate( veciQuery, iContext, iLimit, vecfQuery, DatGraph ) &&
		PixieGraph( DatGraph, vecfQuery ) ) )
		return -1;

	return iRet; }

bool CBNServer::PixieCreate( const vector<size_t>& veciQuery, size_t iContext, size_t iLimit,
	vector<bool>& vecfQuery, CDat& DatGraph ) const {
	vector<float>			vecdNeighbors;
	priority_queue<SPixie>	pqueNeighbors;
	float					d, dAve, dStd, dMin, dCutoff;
	size_t					i, j, iN, iMinOne, iMinTwo;
	vector<size_t>			veciDegree, veciGenes, veciNeighbors;
	bool					fDone;
	CDataMatrix				MatQuery, MatNeighbors;
	vector<string>			vecstrGenes;
	set<size_t>				setiQuery;

	cerr << m_strConnection << " PIXIE query " << iContext << ':';
	for( i = 0; i < veciQuery.size( ); ++i )
		cerr << ' ' << veciQuery[ i ];
	cerr << endl;

	for( i = 0; i < veciQuery.size( ); ++i )
		setiQuery.insert( veciQuery[ i ] );
	vecdNeighbors.resize( m_Database.GetGenes( ) );
	fill( vecdNeighbors.begin( ), vecdNeighbors.end( ), 0.0f );
	MatQuery.Initialize( veciQuery.size( ), vecdNeighbors.size( ) );
	dAve = dStd = 0;
	for( iN = i = 0; i < veciQuery.size( ); ++i ) {
		if( !((CBNServer*)this)->Get( veciQuery[ i ], iContext, MatQuery.Get( i ) ) )
			return false;
		for( j = 0; j < m_Database.GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = MatQuery.Get( i, j ) ) ) {
				iN++;
				dAve += d;
				dStd += d * d;
				if( setiQuery.find( j + 1 ) == setiQuery.end( ) )
					vecdNeighbors[ j ] += d; } }
	for( i = 0; i < vecdNeighbors.size( ); ++i )
		if( ( d = vecdNeighbors[ i ] ) > 0 )
			pqueNeighbors.push( SPixie( i, d ) );
	while( !pqueNeighbors.empty( ) && ( ( veciQuery.size( ) + veciNeighbors.size( ) ) < iLimit ) ) {
		veciNeighbors.push_back( pqueNeighbors.top( ).m_iNode );
		pqueNeighbors.pop( ); }

	vecstrGenes.resize( veciQuery.size( ) + veciNeighbors.size( ) );
	vecfQuery.resize( vecstrGenes.size( ) );
	for( i = 0; i < veciQuery.size( ); ++i ) {
		vecstrGenes[ i ] = m_Database.GetGene( veciQuery[ i ] - 1 );
		vecfQuery[ i ] = true; }
	for( i = 0; i < veciNeighbors.size( ); ++i ) {
		vecstrGenes[ veciQuery.size( ) + i ] = m_Database.GetGene( veciNeighbors[ i ] );
		vecfQuery[ veciQuery.size( ) + i ] = false; }
	MatNeighbors.Initialize( veciNeighbors.size( ), veciNeighbors.size( ) );
	for( i = 0; i < veciNeighbors.size( ); ++i )
		if( !((CBNServer*)this)->Get( veciNeighbors[ i ] + 1, veciNeighbors, iContext,
			MatNeighbors.Get( i ) ) )
			return false;
	DatGraph.Open( vecstrGenes );
	for( i = 0; i < veciQuery.size( ); ++i ) {
		for( j = ( i + 1 ); j < veciQuery.size( ); ++j )
			DatGraph.Set( i, j, MatQuery.Get( i, veciQuery[ j ] - 1 ) );
		for( j = 0; j < veciNeighbors.size( ); ++j )
			DatGraph.Set( i, veciQuery.size( ) + j, MatQuery.Get( i, veciNeighbors[ j ] ) ); }
	for( i = 0; i < veciNeighbors.size( ); ++i ) {
		const float*	adOne	= MatNeighbors.Get( i );
		size_t			iOne	= veciQuery.size( ) + i;

		for( j = ( i + 1 ); j < veciNeighbors.size( ); ++j )
			if( !CMeta::IsNaN( d = adOne[ j ] ) ) {
				iN++;
				dAve += d;
				dStd += d * d;
				DatGraph.Set( iOne, veciQuery.size( ) + j, adOne[ j ] ); } }

	dAve /= iN;
	dStd = sqrt( ( dStd / iN ) - ( dAve * dAve ) );
	dCutoff = max( dAve + dStd, c_dCutoff );
	veciDegree.resize( DatGraph.GetGenes( ) );
	fill( veciDegree.begin( ), veciDegree.end( ), veciDegree.size( ) - 1 );
	for( fDone = false; !fDone; ) {
		fDone = true;
		dMin = FLT_MAX;
		for( i = 0; i < DatGraph.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < DatGraph.GetGenes( ); ++j )
				if( !CMeta::IsNaN( d = DatGraph.Get( i, j ) ) && ( d < dCutoff ) && ( d < dMin ) &&
					( veciDegree[ i ] > c_iDegree ) && ( veciDegree[ j ] > c_iDegree ) ) {
					fDone = false;
					dMin = d;
					iMinOne = i;
					iMinTwo = j; }
		if( !fDone ) {
			veciDegree[ iMinOne ]--;
			veciDegree[ iMinTwo ]--;
			DatGraph.Set( iMinOne, iMinTwo, CMeta::GetNaN( ) ); } }

	return true; }

bool CBNServer::PixieGraph( const CDat& DatGraph, const vector<bool>& vecfQuery ) const {
	static const size_t	c_iBuffer	= 1024;
	char		acBuffer[ c_iBuffer ];
	string		strCmd, strDotIn, strDotOut, strSvg;
	ofstream	ofsm;
	CDot		DotOut( DatGraph );
	CGenome		Genome;
	ostrstream	ossm;
	uint32_t	iSize;

	sprintf_s( acBuffer, ( m_strFiles + "/inXXXXXX" ).c_str( ) );
	if( _mktemp_s( acBuffer ) )
		return false;
	strDotIn = acBuffer;
	strDotIn += c_szDOT;
	ofsm.open( strDotIn.c_str( ) );
	if( !ofsm.is_open( ) )
		return false;
	Genome.Open( DatGraph.GetGeneNames( ) );
	DatGraph.SaveDOT( ofsm, CMeta::GetNaN( ), &Genome );
	ofsm.close( );

	sprintf_s( acBuffer, ( m_strFiles + "/outXXXXXX" ).c_str( ) );
	if( _mktemp_s( acBuffer ) )
		return false;
	strDotOut = acBuffer;
	strDotOut += c_szDOT;
	strCmd = m_strGraphviz + " -Tdot -o" + strDotOut + ' ' + strDotIn;
	system( strCmd.c_str( ) );

	if( !( DotOut.Open( strDotOut.c_str( ) ) && DotOut.Save( ossm, vecfQuery ) ) )
		return false;
	iSize = ossm.pcount( );
	send( m_iSocket, (const char*)&iSize, sizeof(iSize), 0 );
	send( m_iSocket, ossm.str( ), ossm.pcount( ), 0 );

	return true; }

size_t CBNServer::ProcessContexts( const vector<unsigned char>& vecbMessage, size_t iOffset ) {
	uint32_t				iGene;
	size_t					iPlace, iContext;
	vector<unsigned char>	vecbData;

	if( !m_adContexts )
		m_adContexts = new float[ 2 * m_vecBNs.size( ) ];
	for( iPlace = iOffset; ( iPlace + sizeof(iGene) ) <= vecbMessage.size( ); iPlace += sizeof(iGene) ) {
		iGene = *(uint32_t*)&vecbMessage[ iPlace ];
		if( ( iGene < 1 ) || !m_Database.Get( iGene - 1, vecbData ) )
			return -1;
		cerr << m_strConnection << " contexting " << iGene  << " (" << m_Database.GetGene( iGene - 1 ) <<
			")" << endl;
		for( iContext = 0; iContext < m_vecBNs.size( ); ++iContext ) {
			if( !( iContext % 10 ) )
				cerr << m_strConnection << " context " << iContext << '/' << m_vecBNs.size( ) << endl;
			if( !Get( vecbData, iContext ) )
				return -1; }
		iGene = (uint32_t)( m_vecBNs.size( ) * 2 * sizeof(*m_adContexts) );
		send( m_iSocket, (const char*)&iGene, sizeof(iGene), 0 );
		send( m_iSocket, (const char*)m_adContexts, iGene, 0 ); }

	return ( iPlace - iOffset ); }

bool CBNServer::Get( const vector<unsigned char>& vecbData, size_t iContext ) {
	size_t	i, iIn, iOut;
	float	d, dIn, dOut;

	if( !m_adGenes )
		m_adGenes = new float[ m_iGenes ];
	if( !m_vecBNs[ iContext ].Evaluate( vecbData, m_adGenes, m_iGenes ) )
		return false;

// Score is: (sum(w(g, h), h in C)/|C|)/(sum(w(g, h), h in G)/|G|)
//  for gene of interest g, weight w, context C, and genome G
// This is equivalent to: |G|/|C|/sum(w(g, h), h notin C)
	dIn = dOut = 0;
	for( iIn = iOut = i = 0; i < m_iGenes; ++i ) {
		if( CMeta::IsNaN( d = m_adGenes[ i ] ) )
			continue;
		if( m_MatContexts.Get( iContext, i ) ) {
			iIn++;
			dIn += d; }
		iOut++;
		dOut += d; }
	m_adContexts[ iContext << 1 ] = iIn ? ( dIn / iIn ) : 0;
	m_adContexts[ ( iContext << 1 ) | 1 ] = iOut ? ( dOut / iOut ) : 0;

	return true; }

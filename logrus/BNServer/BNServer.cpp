#include "stdafx.h"
#include "cmdline.h"

static const char	c_szXDSL[]	= ".xdsl";
static const char	c_szDSL[]	= ".dsl";

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

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );
	EnableXdslFormat( );

	if( sArgs.default_arg && !( BNSmile.Open( sArgs.default_arg ) && BNDefault.Open( BNSmile ) ) ) {
		cerr << "Could not open: " << sArgs.default_arg << endl;
		return 1; }
	if( sArgs.input_arg ) {
		ifsm.open( sArgs.input_arg );
		pistm = &ifsm; }
	else
		pistm = &cin;
	iMax = 0;
	while( !pistm->eof( ) ) {
		pistm->getline( acBuffer, c_iBuffer - 1 );
		acBuffer[ c_iBuffer - 1 ] = 0;
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
	for( iterBN = mapistrBNs.begin( ); iterBN != mapistrBNs.end( ); ++iterBN )
		cerr << ( ( BNSmile.Open( ( (string)sArgs.directory_arg + '/' + CMeta::Filename( iterBN->second ) +
			( sArgs.xdsl_flag ? c_szXDSL : c_szDSL ) ).c_str( ) ) &&
			vecBNs[ iterBN->first - 1 ].Open( BNSmile ) ) ? "Opened" : "Could not open" ) << ": " <<
			iterBN->second << endl;

	CBNServer	BNServer( BNDefault, vecBNs, 0, sArgs.host_arg, sArgs.dbport_arg,
		sArgs.database_arg ? sArgs.database_arg : "", sArgs.username_arg,
		sArgs.password_arg ? sArgs.password_arg : "", "" );

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
	SOCKET iSocket, const string& strHost, size_t iPort, const string& strDatabase, const string& strUsername,
	const string& strPassword, const string& strConnection ) : m_BNDefault(BNDefault), m_vecBNs(vecBNs),
	m_iSocket(iSocket), m_strHost(strHost), m_iPort(iPort), m_strDatabase(strDatabase),
	m_strUsername(strUsername), m_strPassword(strPassword), m_strConnection(strConnection), m_adValues(NULL) {

	cerr << "New connection from: " << m_strConnection << endl; }

CBNServer::~CBNServer( ) {

	if( m_adValues )
		delete[] m_adValues; }

IServerClient* CBNServer::NewInstance( SOCKET iSocket, uint32_t iHost, uint16_t sPort,
	const CPropertyFile* pConfig ) {
	string	strConnection;
	char	acBuffer[ 16 ];
	in_addr	sAddr;

	sprintf( acBuffer, "%hu", sPort );
	sAddr.s_addr = htonl( iHost );
	strConnection = (string)inet_ntoa( sAddr ) + ":" + acBuffer;
	return new CBNServer( m_BNDefault, m_vecBNs, iSocket, m_strHost, m_iPort, m_strDatabase, m_strUsername,
		m_strPassword, strConnection ); }

void CBNServer::Destroy( ) {

	cerr << "Disconnected: " << m_strConnection << endl;

	delete this; }

bool CBNServer::ProcessMessage( const vector<unsigned char>& vecbMessage ) {
	size_t		iOffset, iStep;
	uint32_t	iGene, iContext;

	if( m_MSQLConnection.ping( ) && !Connect( ) )
		return false;

	iStep = sizeof(iGene) + sizeof(iContext);
	for( iOffset = 0; iOffset < vecbMessage.size( ); iOffset += iStep ) {
		if( ( vecbMessage.size( ) - iOffset ) < iStep )
			return false;
		iGene = *(uint32_t*)&vecbMessage[ iOffset ];
		iContext = *(uint32_t*)&vecbMessage[ iOffset + sizeof(iGene) ];
		cerr << m_strConnection << " requested " << iGene << ':' << iContext << endl;
		if( !Get( iGene, iContext ) )
			return false; }

	return true; }

bool CBNServer::Connect( ) {
	bool	fRet;

	try {
		if( fRet = m_MSQLConnection.connect( m_strDatabase.c_str( ), m_strHost.c_str( ), m_strUsername.c_str( ),
			m_strPassword.c_str( ), (uint)m_iPort ) ) {
			Query	MSQLQuery	= m_MSQLConnection.query( );
			Result	MSQLResult;

			MSQLQuery << "SELECT max(genetwo_id) FROM datapairs";
			if( MSQLResult = MSQLQuery.store( ) )
				m_iGenes = MSQLResult.at( 0 )[ (size_t)0 ];
			else
				fRet = false; } }
	catch( ConnectionFailed ) {
		cerr << "Connection failed: " << m_strConnection << " to " << m_strHost << ':' << m_iPort << ' ' <<
			m_strUsername << '@' << m_strDatabase << endl;
		return false; }

	return fRet; }

bool CBNServer::Get( size_t iGene, size_t iContext ) {
	Query					MSQLQuery	= m_MSQLConnection.query( );
	const CBayesNetMinimal&	BNet		= iContext ? m_vecBNs[ iContext ] : m_BNDefault;
	Result					MSQLResult;
	Row						MSQLRow;
	size_t					i, iPrev, iDataset, iOne, iTwo, iOther;
	unsigned char			bValue;
	vector<unsigned char>	vecbDatum;
	uint32_t				iSize;
	float					d;

	MSQLQuery << "SELECT * FROM datapairs WHERE geneone_id = " << iGene << " OR genetwo_id = " << iGene <<
		" ORDER BY geneone_id, genetwo_id";
	if( !( MSQLResult = MSQLQuery.store( ) ) )
		return false;

	if( !m_adValues )
		m_adValues = new float[ m_iGenes ];
	d = CMeta::GetNaN( );
	memset( m_adValues, *(int*)&d, m_iGenes * sizeof(*m_adValues) );
	vecbDatum.resize( m_BNDefault.GetNodes( ) );
	iPrev = -1;
	for( i = 0; i < MSQLResult.rows( ); ++i ) {
		MSQLRow = MSQLResult.at( (uint)i );
		iDataset = MSQLRow[ (size_t)0 ];
		iOne = MSQLRow[ 1 ];
		iTwo = MSQLRow[ 2 ];
		bValue = MSQLRow[ 3 ];
		iOther = ( iOne < iGene ) ? iOne : iTwo;

		if( iOther != iPrev ) {
			if( iPrev != -1 )
				m_adValues[ iPrev - 1 ] = BNet.Evaluate( vecbDatum );
			iPrev = iOther;
			fill( vecbDatum.begin( ), vecbDatum.end( ), -1 ); }
		vecbDatum[ iDataset
// HACK HACK REMOVE ME
+ 1
// HACK HACK REMOVE ME
 ] = bValue; }
	if( iPrev != -1 )
		m_adValues[ iPrev - 1 ] = BNet.Evaluate( vecbDatum );

	iSize = (uint32_t)( m_iGenes * sizeof(*m_adValues) );
	send( m_iSocket, (char*)&iSize, sizeof(iSize), 0 );
	send( m_iSocket, (char*)m_adValues, iSize, 0 );

	return true; }

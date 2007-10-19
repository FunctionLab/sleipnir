#include "stdafx.h"
#include "cmdline.h"

static const char	c_szXDSL[]	= ".xdsl";
static const char	c_szDSL[]	= ".dsl";

int main( int iArgs, char** aszArgs ) {
	static const size_t					c_iBuffer	= 1024;
	gengetopt_args_info					sArgs;
	CServer								Server;
	CBNServer							BNServer( 0 );
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

CBNServer::CBNServer( SOCKET iSocket ) : m_iSocket(iSocket) { }

IServerClient* CBNServer::NewInstance( SOCKET iSocket, const CPropertyFile* pConfig ) {

	return new CBNServer( iSocket ); }

void CBNServer::Destroy( ) {

	delete this; }

bool CBNServer::ProcessMessage( const vector<unsigned char>& vecbMessage ) {
	size_t	i;

	for( i = 0; i < vecbMessage.size( ); ++i )
		cerr << vecbMessage[ i ];
	cerr << endl;

	return true; }

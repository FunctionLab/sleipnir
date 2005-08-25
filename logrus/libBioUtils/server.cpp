#include "stdafx.h"
#include "server.h"
#include "serverclient.h"

namespace libBioUtils {

const char*	CServerImpl::c_szPort		= "port";
const char*	CServerImpl::c_szTimeout	= "timeout";
CServer*	CServerImpl::s_pServer		= NULL;

bool CServer::Initialize( const char* szConfig, IServerClient* pClient ) {
	SVariant			Var;
#ifndef WIN32
	struct sigaction	Sigact;
#endif // WIN32

	m_pClient = pClient;
	m_Config.Open( szConfig );
	
	Var = m_Config.Get( c_szPort );
	if( Var.m_eType != SVariant::eInt )
		return false;
	m_iPort = Var.m_i;

	Var = m_Config.Get( c_szTimeout );
	if( Var.m_eType != SVariant::eInt )
		return false;
	m_iTimeout = Var.m_i;

#ifndef WIN32
	Sigact.sa_handler = Alarm;
	memset( &Sigact.sa_mask, 0, sizeof(Sigact.sa_mask) );
	Sigact.sa_flags = 0;
	sigaction( SIGALRM, &Sigact, NULL );
#endif // WIN32

	return true; }

#ifndef WIN32
void CServerImpl::Alarm( int iSig ) { }
#endif // WIN32

bool CServer::Start( ) {
	sockaddr_in	Addr;
	char		cOn;

	m_fStop = false;
	m_iSocket = socket( PF_INET, SOCK_STREAM, 0 );
	cOn = true;
	setsockopt( m_iSocket, SOL_SOCKET, SO_REUSEADDR, &cOn, sizeof(cOn) );

	Addr.sin_family = AF_INET;
	Addr.sin_port = htons( m_iPort );
	Addr.sin_addr.s_addr = INADDR_ANY;
	s_pServer = this;
	if( bind( m_iSocket, (const sockaddr*)&Addr, sizeof(Addr) ) ||
		listen( m_iSocket, INT_MAX ) ) {
		g_CatBioUtils.error( "CServer::Start( ) bind failed: %s", strerror( errno ) );
		return false; }
	g_CatBioUtils.notice( "CServer::Start( ) bound to port %d", m_iPort );

	return true; }

void CServer::Listen( ) {
	SOCKET			iClient;
	socklen_t		iSize;
	IServerClient*	pClient;
	CThread			ThreadClient;
	sockaddr_in		Addr;
#ifndef WIN32
	itimerval		Time;

	memset( &Time, 0, sizeof(Time) );
#endif // WIN32
	while( !m_fStop ) {
#ifndef WIN32
		Time.it_value.tv_usec = 1000 * m_iTimeout;
		setitimer( ITIMER_REAL, &Time, NULL );
		Time.it_value.tv_usec = 0;
#endif // WIN32
		iSize = sizeof(Addr);
		if( ( iClient = accept( m_iSocket, (sockaddr*)&Addr, &iSize ) ) == -1 )
			continue;
#ifndef WIN32
		setitimer( ITIMER_REAL, &Time, NULL );
#endif // WIN32

		pClient = m_pClient->NewInstance( );
		pClient->SetParent( this );
		pClient->SetSocket( iClient );
		pClient->SetConfig( &m_Config );

		iSize = ntohl( Addr.sin_addr.s_addr );
		g_CatBioUtils.info( "CServer::Listen( ) client 0x%08x connected from %d.%d.%d.%d",
			pClient, ( iSize >> 24 ) & 0xFF, ( iSize >> 16 ) & 0xFF, ( iSize >> 8 ) & 0xFF,
			iSize & 0xFF );
		ThreadClient.Run( pClient->GetStartRoutine( ), pClient ); }

	g_CatBioUtils.info( "CServer::Listen( ) preparing to shutdown..." ); }

void CServer::Stop( ) {

	m_fStop = true; }

}

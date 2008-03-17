#include "stdafx.h"
#include "serverclient.h"
#include "server.h"
#include "xmlvacuum.h"

namespace Sleipnir {

void* CServerClientImpl::StartRoutine( void* pArg ) {
	CServerClientImpl*	pClient;

	pClient = (CServerClientImpl*)pArg;
	pClient->StartRoutine( );
	delete pClient;

	return NULL; }

void CServerClientImpl::StartRoutine( ) {

	do
		ReadMessage( );
	while( m_vecbMessage.size( ) && m_pClient->ProcessMessage( m_vecbMessage ) );

	g_CatSleipnir.info( "Client 0x%08x disconnected", m_pClient ); }

CServerClientImpl::CServerClientImpl( SOCKET iSocket, IServerClient* pClient ) : m_iSocket(iSocket),
	m_pClient(pClient) { }

CServerClientImpl::~CServerClientImpl( ) {

	closesocket( m_iSocket );
	m_pClient->Destroy( ); }

void CServerClientImpl::ReadMessage( ) {
	unsigned char	acBuffer[ c_iBuffer ];
	size_t			iRead, iTotal, iPrev;

	m_vecbMessage.clear( );
	if( !( iRead = recv( m_iSocket, (char*)acBuffer, c_iBuffer, 0 ) ) || ( iRead == -1 ) ||
		( iRead < sizeof(uint32_t) ) )
		return;
	iTotal = *(uint32_t*)acBuffer;
	m_vecbMessage.resize( max( iTotal, iPrev = iRead - sizeof(uint32_t) ) );
	copy( acBuffer + sizeof(uint32_t), acBuffer + iRead, m_vecbMessage.begin( ) );
	while( ( iPrev < iTotal ) && ( iRead = recv( m_iSocket, (char*)acBuffer, c_iBuffer, 0 ) ) ) {
		if( ( iPrev + iRead ) >= m_vecbMessage.size( ) )
			m_vecbMessage.resize( iPrev + iRead );
		copy( acBuffer, acBuffer + iRead, m_vecbMessage.begin( ) + iPrev );
		iPrev += iRead; } }

}

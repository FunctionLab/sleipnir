#include "stdafx.h"
#include "serverclient.h"
#include "server.h"
#include "xmlvacuum.h"

namespace libBioUtils {

#ifdef XML_ENABLED
const XMLCh	CServerClientImpl::c_pxszLS[]	= { chLatin_L, chLatin_S, 0 };
#endif // XML_ENABLED

void* CServerClientImpl::StartRoutine( void* pArg ) {
	CServerClientImpl*	pClient;

	pClient = (CServerClientImpl*)pArg;
	pClient->StartRoutine( );
	delete pClient;

	return NULL; }

void CServerClientImpl::StartRoutine( ) {

	do
		ReadMessage( );
	while(
#ifdef XML_ENABLED
			( m_pxdocMessage && m_pxdocMessage->getDocumentElement( ) )
#else // XML_ENABLED
			m_vecbMessage.size( )
#endif // XML_ENABLED
			&& m_pClient->ProcessMessage(
#ifdef XML_ENABLED
			m_pxdocMessage
#else // XML_ENABLED
			m_vecbMessage
#endif // XML_ENABLED
			) );

	g_CatBioUtils.info( "Client 0x%08x disconnected", m_pClient ); }

CServerClientImpl::CServerClientImpl( SOCKET iSocket, IServerClient* pClient ) : m_iSocket(iSocket),
	m_pClient(pClient) {
#ifdef XML_ENABLED
	DOMImplementation*	pximp;

	m_pxdocMessage = NULL;
	pximp = DOMImplementationRegistry::getDOMImplementation( c_pxszLS );
	m_pxbld = ((DOMImplementationLS*)pximp)->createDOMBuilder(
		DOMImplementationLS::MODE_SYNCHRONOUS, NULL );
	m_pxbld->setFeature( XMLUni::fgDOMComments, false );
	m_pxbld->setFeature( XMLUni::fgDOMWhitespaceInElementContent, false );
#endif // XML_ENABLED
}

CServerClientImpl::~CServerClientImpl( ) {

	closesocket( m_iSocket );
	m_pClient->Destroy( );
#ifdef XML_ENABLED
	if( m_pxbld )
		m_pxbld->release( );
#endif // XML_ENABLED
}

void CServerClientImpl::ReadMessage( ) {
#ifdef XML_ENABLED
	string				strXml;
	CXmlVacuum			Vacuum( (int)m_iSocket );
	CStringInputSource	In( strXml );
	
	strXml = Vacuum.Suck( );
	m_pxdocMessage = m_pxbld->parse( In );
#else // XML_ENABLED
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
		iPrev += iRead; }
#endif // XML_ENABLED
}

#ifdef XML_ENABLED
CStringInputSource::CStringInputSource( const string& str ) : m_str( str ) { }

const XMLCh* CStringInputSource::getEncoding( ) const {

	return NULL; }

const XMLCh* CStringInputSource::getPublicId( ) const {

	return NULL; }

const XMLCh* CStringInputSource::getSystemId( ) const {

	return NULL; }

const XMLCh* CStringInputSource::getBaseURI( ) const {

	return NULL; }

void CStringInputSource::setEncoding( const XMLCh* ) { }

void CStringInputSource::setPublicId( const XMLCh* ) { }

void CStringInputSource::setSystemId( const XMLCh* ) { }

void CStringInputSource::setBaseURI( const XMLCh* ) { }

void CStringInputSource::setIssueFatalErrorIfNotFound( bool ) { }

bool CStringInputSource::getIssueFatalErrorIfNotFound( ) const {
	
	return false; }

BinInputStream* CStringInputSource::makeStream( ) const {
	
	return new CBinStringInputStream( m_str ); }

void CStringInputSource::release( ) {
	
	delete this; }

CBinStringInputStream::CBinStringInputStream( const string& str ) : m_str( str ), m_iPos(-2),
	m_fHalf(false) { }

unsigned int CBinStringInputStream::curPos( ) const {
	
	return ( m_iPos + 2 ); }

unsigned int CBinStringInputStream::readBytes( XMLByte* const pxbytFill, unsigned int iCount ) {
	size_t	iByte;

	iByte = 0;
	if( ( m_iPos < 0 ) && iCount ) {
		pxbytFill[ iByte++ ] = 0xFE;
		m_iPos++;
		if( iCount > 1 ) {
			pxbytFill[ iByte++ ] = 0xFF;
			m_iPos++; } }
	for( ; ( m_iPos < (int)m_str.size( ) ) && ( iByte < iCount ); ++iByte ) {
		pxbytFill[ iByte ] = ( m_fHalf ? ( m_str[ m_iPos++ ] & 0xFF ) : 0 );
		m_fHalf = !m_fHalf; }

	return iByte; }
#endif // XML_ENABLED

}

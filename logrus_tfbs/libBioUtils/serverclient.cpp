#include "stdafx.h"
#include "serverclient.h"
#include "xmlvacuum.h"

namespace libBioUtils {

const XMLCh	CServerClientImpl::c_pxszLS[]	= { chLatin_L, chLatin_S, 0 };

void* CServerClientImpl::StartRoutine( void* pArg ) {
	CServerClient*	pClient;
	
	pClient = (CServerClient*)pArg;
	while( pClient->m_fOpen ) {
		pClient->ReadMessage( );
		if( !( pClient->m_pxdocMessage &&
			pClient->m_pxdocMessage->getDocumentElement( ) ) )
			break;
		pClient->ProcessMessage( ); }

	g_CatBioUtils.info( "Client 0x%08x disconnected", pClient );

	closesocket( pClient->m_iSocket );
	pClient->Destroy( );

	return pArg; }

CServerClientImpl::CServerClientImpl( ) {
	DOMImplementation*	pximp;

	m_pParent = NULL;

	m_pxdocMessage = NULL;
	pximp = DOMImplementationRegistry::getDOMImplementation( c_pxszLS );
	m_pxbld = ((DOMImplementationLS*)pximp)->createDOMBuilder(
		DOMImplementationLS::MODE_SYNCHRONOUS, NULL );
	m_pxbld->setFeature( XMLUni::fgDOMComments, false );
	m_pxbld->setFeature( XMLUni::fgDOMWhitespaceInElementContent, false ); }

CServerClientImpl::~CServerClientImpl( ) {
	
	if( m_pxbld )
		m_pxbld->release( ); }

void CServerClient::SetParent( CServer* pParent ) {
	
	m_pParent = pParent; }

void CServerClient::SetSocket( SOCKET iSocket ) {
	
	m_iSocket = iSocket;
	m_fOpen = true; }

void CServerClient::SetConfig( const CPropertyFile* pConfig ) {
	
	m_pConfig = pConfig; }

IServerClient::TPFnStart CServerClient::GetStartRoutine( ) {
	
	return StartRoutine; }

void CServerClientImpl::ReadMessage( ) {
	string				strXml;
	CXmlVacuum			Vacuum( (int)m_iSocket );
	CStringInputSource	In( strXml );
	
	strXml = Vacuum.Suck( );
	m_pxdocMessage = m_pxbld->parse( In ); }

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

}

#ifndef SERVERCLIENT_H
#define SERVERCLIENT_H

#include "typesi.h"

namespace Sleipnir {

class CServer;

/*!
 * \brief
 * Provide a simple interface for objects handling network requests on a server thread.
 * 
 * When a CServer detects an incoming connection, it creates a new thread and passes control of that
 * thread to an IServerClient to communicate with the incoming client.  By implementing IServerClient,
 * a Sleipnir user can quickly and easily create multithreaded network servers.
 * 
 * Sleipnir servers always assume that TCP/IP messages are passed preceded by a four-byte unsigned
 * integer indicating the size of the message.  For example, to send the string, "Hello, world!" to a
 * CServer object, a client should send the integer 13 to the server, followed by the 13 bytes of the
 * message (assuming ASCII encoding, of course).  It is not necessary for Sleipnir servers to respond to
 * their clients in the same manner - this is left up to the implementation of ProcessMessage - but it is
 * good practice and can greatly simplify TCP/IP communication.
 * 
 * For example, to create an echo server client, one might implement:
 * \code
 * class CEchoServerClient : public IServerClient {
 * public:
 *   SOCKET  m_iSocket;
 * 
 *   CEchoServerClient( SOCKET iSocket ) : m_iSocket(iSocket) { }
 * 
 *   IServerClient* NewInstance( SOCKET iSocket, uint32_t iHost, uint16_t iPort ) {
 *     return new CEchoServerClient( iSocket ); }
 * 
 *   void Destroy( ) {
 *     delete this; }
 * 
 *   bool ProcessMessage( const vector<unsigned char>& vecbMessage ) {
 *     size_t    i;
 *     uint32_t  iSize;
 * 
 *     iSize = vecbMessage.size( );
 *     send( m_iSocket, &iSize, sizeof(iSize), 0 );
 *     for( i = 0; i < vecbMessage.size( ); ++i )
 *       if( vecbMessage[ i ] == EOF )
 *         return false;
 *       else
 *         send( m_iSocket, &vecbMessage[ i ], 1, 0 );
 * 
 *     return true; }
 * }
 * \endcode
 * See CServer for an example server setup using this client object.
 * 
 * \remarks
 * CServer will call IServerClient::NewInstance for each new incoming request, creating a new client object
 * to handle that thread.  When the connection is closed by IServerClient::ProcessMessage returning false,
 * the server will call IServerClient::Destroy to clean up that object.  The original server client object
 * (ESC in the example above) is only used to create additional new instances.
 */
class IServerClient {
public:
	/*!
	 * \brief
	 * Requests a clone of the current server client object to handle an incoming connection on the given
	 * socket.
	 * 
	 * \param iSocket
	 * Socket to be used for communication with the incoming request.
	 * 
	 * \param iHost
	 * Host byte-ordered encoding of the incoming client's IP address.
	 * 
	 * \param iPort
	 * Host byte-ordered encoding of the incoming client's IP port.
	 * 
	 * \returns
	 * A clone of the current object configured to handle the given incoming connection request.
	 * 
	 * CServer calls NewInstance whenever an incoming connection request is detected.  The returned new
	 * server client object should be configured to store iSocket and, optionally, the host and port
	 * information so as to communicate with the new client during ProcessMessage.
	 * 
	 * \remarks
	 * iHost and iPort are ntoh decoded versions of sin_addr.s_addr and sin_port, respectively.
	 */
	virtual IServerClient* NewInstance( SOCKET iSocket, uint32_t iHost, uint16_t iPort ) = 0;
	/*!
	 * \brief
	 * Called when a new message has been received from the client and should be processed by the server.
	 * 
	 * \param vecbMessage
	 * Contents of the message sent by the clients (without the preceding byte size indicator).
	 * 
	 * \returns
	 * True if the message was handled and the server should continue the connection, false to close the
	 * connection.
	 * 
	 * \remarks
	 * Sleipnir puts no constraints on the contents of messages passed to and from network servers - the bytes
	 * in vecbMessage can be ASCII text or any other arbitrary data.
	 */
	virtual bool ProcessMessage( const std::vector<unsigned char>& vecbMessage ) = 0;
	/*!
	 * \brief
	 * Indicates that the server client object should delete itself.
	 * 
	 * \remarks
	 * Provided so that delete doesn't need to be called directly on an opaque interface.
	 */
	virtual void Destroy( ) = 0;
};

}

#endif // SERVERCLIENT_H

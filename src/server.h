#ifndef SERVER_H
#define SERVER_H

#include "serveri.h"

namespace Sleipnir {

/*!
 * \brief
 * Provide a basic multithreaded TCP/IP server for simple communication tasks.
 * 
 * A server object listens for incoming connections on a specified port.  When a connection is made, a new
 * thread is spawned and control is given to an IServerClient object, a handler that knows how to respond
 * to incoming messages.  This means that by implementing the simple IServerClient interface, a server
 * program can be created quickly and easily to service network requests in a robust, multithreaded manner.
 * 
 * Sleipnir servers always assume that TCP/IP messages are passed preceded by a four-byte unsigned
 * integer indicating the size of the message.  For example, to send the string, "Hello, world!" to a
 * CServer object, a client should send the integer 13 to the server, followed by the 13 bytes of the
 * message (assuming ASCII encoding, of course).  It is not necessary for Sleipnir servers to respond to
 * their clients in the same manner - this is left up to the implementation of IServerClient::ProcessMessage -
 * but it is good practice and can greatly simplify TCP/IP communication.
 * 
 * Example usage for a server on port 1234 might resemble:
 * \code
 * CServer            Server;
 * CEchoServerClient  ESC;
 * 
 * Server.Initialize( 1234, 100, &ESC );
 * Server.Start( );
 * \endcode
 * See IServerClient for an example implementation of an echo server client.
 * 
 * \remarks
 * CServer will call IServerClient::NewInstance for each new incoming request, creating a new client object
 * to handle that thread.  When the connection is closed by IServerClient::ProcessMessage returning false,
 * the server will call IServerClient::Destroy to clean up that object.  The original server client object
 * (ESC in the example above) is only used to create additional new instances.
 */
class CServer : CServerImpl {
public:
	bool Initialize( size_t iPort, size_t iTimeout, IServerClient* pServerClient );
	bool Start( );

	/*!
	 * \brief
	 * Signals that the server should stop listening for incoming connections.
	 * 
	 * \remarks
	 * Results in an exit from Start at the next timeout interval.
	 */
	void Stop( ) {

		m_fStop = true; }
};

}

#endif // SERVER_H

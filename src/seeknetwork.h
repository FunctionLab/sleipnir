#ifndef SEEKNETWORK_H
#define SEEKNETWORK_H

#include "seekbasic.h"

//additional network include files
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>
#include <string.h>

namespace Sleipnir {

/*!
 * \brief Utilities for sending and receiving data over the network
 *
 * This class provides static utility functions to facilitate the exchange of messages between the Seek
 * client and the Seek server. In order to allow this exchange to occur, all messages
 * must conform to a uniform standard.
 *
 * \section sec_out Outgoing messages
 *
 * On the sending end, all outgoing messages must first begin with a message header that specifies the
 * length and the type of the message. Then the body of the message follows.
 *
 * The supported outgoing messages are: an array of \c chars (such as a \c string), an array of \c floats.
 * The outgoing message is structured as follows:
 * \li Byte #1-4: An \c unsigned \c integer that specifies the size of one element (\a S). (1 for a \c char, 4 for a \c float)
 * \li Byte #5-8: An \c unsigned \c integer that specifies the total number of elements to be sent (\a N). (1 for a single-value,
 * otherwise the size of the array)
 * \li Byte #9 and onward: \a S times \a N bytes specifying the array content
 *
 * IMPORTANT:
 * <b>Outgoing messages are always encoded using bytes in the Little Endian order.</b>
 *
 *
 * \section sec_in Incoming messages
 *
 * On the receiving end, CSeekNetwork also supports the receiving of a \c char array (or a \c string) or a \c float array.
 *
 * In order to be properly recognized, the incoming message should be structured as follows:
 *
 * For a \c char array:
 * \li Byte #1-4: A \c signed \c integer that specifies the length of the \c char array to receive (\a NC)
 * \li Byte #5 and onward: \a NC bytes specifying the \c char array.
 *
 * For a \c float array:
 * \li Byte #1-4: A \c signed \c integer that specifies the length of the \c float array to receive (\a NF)
 * \li Byte #5 and onward: \a NF times 4 bytes specifying the \c float array.
 *
 * IMPORTANT:
 * <b>For an incoming message to be properly recognized, the incoming message should also be encoded with bytes in the Little Endian order.
 * </b>
 */
class CSeekNetwork{
public:
	/*!
	 * \brief Send a string
	 *
	 * Encodes an outgoing message and sends it to the client
	 *
	 * \param new_fd The client socket
	 * \param str The string to be sent to the client
	 *
	 * \remarks Assumes that the client connection has been established.
	 */
	static int Send(int, const string&);

	/*!
	 * \brief Send a float array
	 *
	 * Encodes an outgoing message and sends it to the client
	 *
	 * \param new_fd The client socket
	 * \param str The array of floats to be sent to the client
	 *
	 * \remarks Assumes that the client connection has been established.
	 */
	static int Send(int, const vector<float>&);

	/*!
	 * \brief Low-level send function
	 *
	 * \param new_fd The client socket
	 * \param c The message
	 * \param size The message length
	 * \return -1 if an error occurs or \c size if the sending is successful
	 *
	 * \remarks Assumes that the client connection has been established.
	 */
	static int Send(int, char*, int);

	/*!
	 * \brief Clear a char array
	 *
	 * Clears a char array by zeroing all bytes
	 *
	 * \param b The char array
	 * \param size The size of the char array
	 */
	static void Clear(char*, int);

	/*!
	 * \brief Copy a char array
	 *
	 * Copies the entire source array (0...N) to the destination array beginning at the index \c beg
	 *
	 * \param d The destination
	 * \param s The source
	 * \param beg The position on the destination array where the pasting starts
	 * \param num The size of the source array
	 * \return \c beg + \c num
	 */
	static int Copy(char*, char*, int, int);

	/*!
	 * \brief Receive a string
	 *
	 * Receive a string from the client
	 *
	 * \param new_fd The client socket
	 * \param s The string where the message will be received to
	 *
	 * \remarks Assumes that the client connection has been established.
	 */
	static int Receive(int, string&);

	/*!
	 * \brief Receive a float array
	 *
	 * Receive a float array from the client
	 *
	 * \param new_fd The client socket
	 * \param f The float array where the message will be received to
	 *
	 * \remarks Assumes that the client connection has been established.
	 */
	static int Receive(int, vector<float>&);
};

}	
#endif

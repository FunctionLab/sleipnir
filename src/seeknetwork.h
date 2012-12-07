#ifndef SEEKNETWORK_H
#define SEEKNETWORK_H

#include "seekbasic.h"

//additinal network include files
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>
#include <string.h>

namespace Sleipnir {

class CSeekNetwork{
public:
	static int Send(int, const string&);
	static int Send(int, const vector<float>&);
	static int Send(int, char*, int);
	static void Clear(char*, int);
	static int Copy(char*, char*, int, int);
	static int Receive(int, string&);
};

}	
#endif

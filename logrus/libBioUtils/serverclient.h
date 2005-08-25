#ifndef SERVERCLIENT_H
#define SERVERCLIENT_H

namespace libBioUtils {

class CPropertyFile;
class CServer;

class IServerClient {
public:
	typedef void* (*TPFnStart)( void* );
	
	virtual IServerClient* NewInstance( ) = 0;
	virtual void SetParent( CServer* ) = 0;
	virtual void SetSocket( SOCKET ) = 0;
	virtual void SetConfig( const CPropertyFile* ) = 0;
	virtual TPFnStart GetStartRoutine( ) = 0;
	virtual void Destroy( ) = 0;
};

}

#endif // SERVERCLIENT_H

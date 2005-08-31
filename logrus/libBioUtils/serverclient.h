#ifndef SERVERCLIENT_H
#define SERVERCLIENT_H

#include <xercesc/dom/DOMInputSource.hpp>
#include <xercesc/util/BinInputStream.hpp>

#include "serverclienti.h"

namespace libBioUtils {

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

class CServerClient : public IServerClient, public CServerClientImpl {
public:
	void SetParent( CServer* );
	void SetSocket( SOCKET );
	void SetConfig( const CPropertyFile* );
	TPFnStart GetStartRoutine( );
};

class CStringInputSource : public xercesc::DOMInputSource {
public:
	CStringInputSource( const std::string& );
	
	const XMLCh* getEncoding( ) const;
	const XMLCh* getPublicId( ) const;
	const XMLCh* getSystemId( ) const;
	const XMLCh* getBaseURI( ) const;
	void setEncoding( const XMLCh* const );
	void setPublicId( const XMLCh* const );
	void setSystemId( const XMLCh* const );
	void setBaseURI( const XMLCh* const );
	void setIssueFatalErrorIfNotFound( const bool );
	bool getIssueFatalErrorIfNotFound( ) const;
	xercesc::BinInputStream* makeStream( ) const;
	void release( );

private:
	const std::string&	m_str;
};

class CBinStringInputStream : public xercesc::BinInputStream {
public:
	CBinStringInputStream( const std::string& );
	
	unsigned int curPos( ) const;
	unsigned int readBytes( XMLByte* const, const unsigned int );

private:
	const std::string&	m_str;
	int					m_iPos;
	bool				m_fHalf;
};

}

#endif // SERVERCLIENT_H

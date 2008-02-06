#ifndef SERVERCLIENTI_H
#define SERVERCLIENTI_H

#include "xml.h"

#include <xercesc/dom/DOMBuilder.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/util/XMLChar.hpp>

namespace libBioUtils {

class CPropertyFile;
class CServer;

class CServerClientImpl {
protected:
	static const XMLCh	c_pxszLS[];

	static void* StartRoutine( void* );

	CServerClientImpl( );
	~CServerClientImpl( );

	virtual void ProcessMessage( ) = 0;

	SOCKET					m_iSocket;
	CServer*				m_pParent;
	xercesc::DOMDocument*	m_pxdocMessage;
	xercesc::DOMBuilder*	m_pxbld;
	const CPropertyFile*	m_pConfig;
	bool					m_fOpen;

private:
	void ReadMessage( );
};

}

#endif // SERVERCLIENTI_H

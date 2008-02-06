#ifndef STDAFX_H
#define STDAFX_H

#include <math.h>
#include <pthread.h>
#ifdef WIN32
#include <winsock2.h>
#else // WIN32
#include <sys/socket.h>

#define _strdup		strdup
#define SOCKET		int
#define sprintf_s	sprintf
#endif // WIN32

#include <algorithm>
#include <fstream>
using namespace std;

#pragma warning(disable:4275)

#ifdef USE_XML
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/util/XMLString.hpp>
XERCES_CPP_NAMESPACE_USE

#include <xercesc/util/PlatformUtils.hpp>
#include <xalanc/XPath/XPathEvaluator.hpp>
XALAN_USING_XALAN(XPathEvaluator)
XALAN_USING_XERCES(XMLPlatformUtils)
#endif // USE_XML

#include <history.h>
#include <readline.h>

#include "annotation.h"
#include "genome.h"
#include "meta.h"
#include "server.h"
#include "serverclient.h"
#include "xmlwriter.h"
using namespace libBioUtils;

#endif // STDAFX_H

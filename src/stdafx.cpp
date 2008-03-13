#include "stdafx.h"

namespace Sleipnir {

const char	c_szSleipnir[]	= "Sleipnir";
#ifdef USE_LOG4CPP_STUB
Category	g_CatSleipnir;
#else // USE_LOG4CPP_STUB
Category&	g_CatSleipnir	= Category::getInstance( c_szSleipnir );
#endif // USE_LOG4CPP_STUB

}

#ifndef ORTHOLOGY_H
#define ORTHOLOGY_H

#include "orthologyi.h"

namespace libBioUtils {

class COrthology : COrthologyImpl {
public:
	bool Open( std::istream& );
};

}

#endif // ORTHOLOGY_H

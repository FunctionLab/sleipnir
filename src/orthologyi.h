/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#ifndef ORTHOLOGYI_H
#define ORTHOLOGYI_H

#include <map>

#include "file.h"

namespace Sleipnir {

    class CGene;

    class CGenome;

    class COrthologyImpl : protected CFile {
    protected:
        static const char c_cOrgSep = '|';

        typedef std::map<CGene *, size_t> TMapGeneI;

        ~COrthologyImpl();

        void Reset();

        std::vector <std::string> m_vecstrOrganisms;
        std::vector<CGenome *> m_vecpGenomes;
        TMapGeneI m_mapGenes;
        std::vector <std::vector<CGene *>> m_vecvecpGenes;
    };

}

#endif // ORTHOLOGYI_H

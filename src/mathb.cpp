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
#include "stdafx.h"
#include "mathb.h"

namespace Sleipnir {

    vector<double>    CMathImpl::s_vecdLogFact;

/*!
 * \brief
 * Returns the log of N factorial.
 * 
 * \param iN
 * Integer to be factorialed and logged.
 * 
 * \returns
 * log(N!)
 * 
 * \remarks
 * Calculated exactly and cached for small N, calculated using Stirling's approximation for large N.
 */
    double CMath::LogFact(size_t iN) {
        size_t i;

        if (s_vecdLogFact.empty()) {
            s_vecdLogFact.resize(c_iFactCutoff);
            s_vecdLogFact[0] = s_vecdLogFact[1] = 0;
            for (i = 2; i < s_vecdLogFact.size(); ++i)
                s_vecdLogFact[i] = s_vecdLogFact[i - 1] + log((double) i);
        }

        return ((iN >= c_iFactCutoff) ? LogFactStirling(iN) : s_vecdLogFact[iN]);
    }

    double CMathImpl::LogFactStirling(size_t iN) {

        if (iN < 3)
            return ((iN == 2) ? log(2.0) : 0);

        return (log(2.5066282746271) + ((iN + 0.5) * log((double) iN)) - iN);
    }

    double CMathImpl::LogFactRec(size_t iN) {

        if (iN < 3)
            return ((iN == 2) ? log(2.0) : 0);

        return (log((double) iN) + LogFactRec(iN - 1));
    }

}

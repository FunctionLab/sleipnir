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
#ifndef EXAMPLEI_H
#define EXAMPLEI_H

#include <utility>

#include "datapair.h"

namespace Sleipnir {

    class CExampleImpl {
    public:
        CExampleImpl();

        ~CExampleImpl();

        void Set(size_t, float, const CDataPair &, size_t);

        bool Equals(const CExampleImpl &, size_t) const;

        void Reset();

        bool IsEvidence(size_t) const;

        size_t GetDiscrete(size_t iDatum) const {

            return m_auData[iDatum].m_i;
        }

        float GetContinuous(size_t iDatum) const {

            return m_auData[iDatum].m_d;
        }

        bool IsSet() const {

            return !!m_auData;
        }

    protected:
        union UDatum {
            float m_d;
            size_t m_i;
        };

        UDatum *m_auData;
    };

}

#endif // EXAMPLEI_H

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
#ifndef COALESCEBASEI_H
#define COALESCEBASEI_H

#include <map>
#include <string>
#include <vector>

namespace Sleipnir {

    class CCoalesceSequencerBase {
    public:
        enum ESubsequence {
            ESubsequenceBegin = 0,
            ESubsequenceTotal = ESubsequenceBegin,
            ESubsequenceIntrons = ESubsequenceTotal + 1,
            ESubsequenceExons = ESubsequenceIntrons + 1,
            ESubsequenceEnd = ESubsequenceExons + 1
        };

        static const char *c_aszSubsequences[];

        static const char *GetSubsequence(ESubsequence eSubsequence) {

            return c_aszSubsequences[eSubsequence];
        }

        static ESubsequence GetSubsequence(const std::string &strSubsequence) {
            size_t i;

            for (i = 0; c_aszSubsequences[i]; ++i)
                if (strSubsequence == c_aszSubsequences[i])
                    return (ESubsequence) i;

            return ESubsequenceEnd;
        }

        size_t GetTypes() const {

            return m_vecstrTypes.size();
        }

        const std::string &GetType(size_t iType) const {

            return m_vecstrTypes[iType];
        }

        size_t GetType(const std::string &strType) const {
            TMapStrI::const_iterator iterType;

            return (((iterType = m_mapstriTypes.find(strType)) == m_mapstriTypes.end()) ? -1 :
                    iterType->second);
        }

    protected:
        typedef std::map <std::string, size_t> TMapStrI;

        TMapStrI m_mapstriTypes;
        std::vector <std::string> m_vecstrTypes;
    };

    template<class tType>
    class CCoalesceSequencer : public CCoalesceSequencerBase {
    public:
        tType &Get(size_t iType, ESubsequence eSubsequence) {

            return m_vecvecValues[iType][eSubsequence];
        }

        const tType &Get(size_t iType, ESubsequence eSubsequence) const {

            return m_vecvecValues[iType][eSubsequence];
        }

        const tType &Get(const std::string &strType, ESubsequence eSubsequence) const {

            return Get(GetType(strType), eSubsequence);
        }

        size_t AddType(const std::string &strType) {
            TMapStrI::const_iterator iterType;
            size_t iRet;

            if ((iterType = m_mapstriTypes.find(strType)) != m_mapstriTypes.end())
                return iterType->second;

            m_mapstriTypes[strType] = iRet = m_vecvecValues.size();
            m_vecstrTypes.push_back(strType);
            m_vecvecValues.push_back(std::vector<tType>());
            m_vecvecValues.back().resize(ESubsequenceEnd);

            return iRet;
        }

        size_t GetSubsequences(size_t iType) const {

            return m_vecvecValues[iType].size();
        }

        void Clear() {
            size_t i, j;

            for (i = 0; i < m_vecvecValues.size(); ++i)
                for (j = 0; j < m_vecvecValues[i].size(); ++j)
                    m_vecvecValues[i][j].Clear();
        }

    protected:
// Type by subsequence
        std::vector <std::vector<tType>> m_vecvecValues;
    };

}

#endif // COALESCEBASEI_H

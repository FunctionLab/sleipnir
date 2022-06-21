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
#ifndef HMM_H
#define HMM_H

#include "hmmi.h"

namespace Sleipnir {

/*!
 * \brief
 * Extremely simple Hidden Markov Model (HMM) implementation allowing learning and generation from the model.
 *
 * Implements a simple HMM of arbitrary degree over a fixed alphabet of single characters.  The HMM can be
 * learned from input strings and, after training, can generate random sequences probabilistically from the
 * model.
 */
    class CHMM : protected CHMMImpl {
    public:
        /*!
         * \brief
         * Initializes an empty HMM of the requested degree over the given alphabet.
         *
         * \param iDegree
         * Degree of the HMM.
         *
         * \param strAlphabet
         * Alphabet of characters encoded by the HMM.
         *
         * \remarks
         * iDegree must be greater than zero; strAlphabet must be nonempty.
         */
        void Open(size_t iDegree, const std::string &strAlphabet) {

            m_iDegree = iDegree;
            m_strAlphabet = strAlphabet;
            m_MatTransitions.Initialize(GetStates(), GetSymbols() - 1);
            m_MatTransitions.Clear();
        }

        /*!
         * \brief
         * Outputs a simple textual representation of the HMM.
         *
         * \param ostm
         * Stream to which the HMM is saved.
         */
        void Save(std::ostream &ostm) const {
            size_t i, j;

            ostm << m_iDegree << endl;
            ostm << m_strAlphabet << endl;
            for (i = 0; i < m_MatTransitions.GetRows(); ++i) {
                ostm << Decode(i);
                for (j = 0; j < m_MatTransitions.GetColumns(); ++j)
                    ostm << '\t' << m_MatTransitions.Get(i, j);
                ostm << endl;
            }
        }

        /*!
         * \brief
         * Updates the HMMs transition probabilities using the given string as training data.
         *
         * \param strData
         * String of characters from which transition probabilities are updated.
         *
         * \returns
         * True if the transitions were updated successfully; false otherwise.
         *
         * \remarks
         * Characters not in the HMM's alphabet are ignored.
         *
         * \see
         * Get
         */
        bool Add(const std::string &strData) {
            size_t i, iState;

            iState = 0;
            for (i = 0; i < strData.length(); ++i) {
                m_MatTransitions.Get(iState, Encode(strData[i]))++;
                iState = (i < m_iDegree) ? Encode(strData, i + 1) :
                         Encode(strData.substr(i - m_iDegree + 1, m_iDegree), m_iDegree);
                if (iState == -1)
                    return false;
            }

            return true;
        }

        /*!
         * \brief
         * Randomly generates the requested number of characters from the HMM's alphabet using its current
         * transition probabilities.
         *
         * \param iLength
         * Number of characters to randomly generate.
         *
         * \returns
         * The requested number of characters randomly generated using the HMM's current transition probabilities.
         *
         * \see
         * Add
         */
        std::string Get(size_t iLength) const {
            std::string strRet;
            size_t i, iState, iTotal, iCur;

            for (iState = 0; strRet.length() < iLength;) {
                for (iTotal = i = 0; i < m_MatTransitions.GetColumns(); ++i)
                    iTotal += m_MatTransitions.Get(iState, i);
                if ((iCur = rand()) == RAND_MAX)
                    iCur--;
                iCur = (size_t)(((float) iCur / (float) RAND_MAX) * iTotal);
                for (i = 0; (i + 1) < m_MatTransitions.GetColumns(); ++i) {
                    iTotal -= m_MatTransitions.Get(iState, i);
                    if (iCur >= iTotal)
                        break;
                }
                strRet += m_strAlphabet[i];
                iState = ((iState * GetSymbols()) + i + 1) % m_MatTransitions.GetRows();
            }

            return strRet;
        }

        /*!
         * \brief
         * Sets all transition probabilities in the HMM to uniform probabilities.
         */
        void SetUniform() {
            size_t i, j;

            for (i = 0; i < m_MatTransitions.GetRows(); ++i)
                for (j = 0; j < m_MatTransitions.GetColumns(); ++j)
                    m_MatTransitions.Set(i, j, 1);
        }
    };

}

#endif // HMM_H

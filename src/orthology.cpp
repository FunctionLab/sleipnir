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
#include "orthology.h"
#include "genome.h"
#include "meta.h"

namespace Sleipnir {

    COrthologyImpl::~COrthologyImpl() {

        Reset();
    }

    void COrthologyImpl::Reset() {
        size_t i;

        for (i = 0; i < m_vecpGenomes.size(); ++i)
            delete m_vecpGenomes[i];
        m_vecpGenomes.clear();
        m_vecstrOrganisms.clear();
        m_mapGenes.clear();
        m_vecvecpGenes.clear();
    }

/*!
 * \brief
 * Loads a new orthology file from the given stream.
 * 
 * \param istm
 * Stream from which orthology file is loaded.
 * 
 * \returns
 * True if the orthology was loaded successfully.
 * 
 * \remarks
 * Constructs CGenome and CGene objects internally as necessary.
 * 
 * \see
 * Save
 */
    bool COrthology::Open(std::istream &istm) {
        vector <string> vecstrLine;
        char *acBuf;
        size_t i, j;
        string strOrganism, strGene;
        CGenome *pGenome;
        CGene *pGene;

        Reset();
        acBuf = new char[c_iBufferSize];
        while (istm.peek() != EOF) {
            istm.getline(acBuf, c_iBufferSize - 1);
            vecstrLine.clear();
            CMeta::Tokenize(acBuf, vecstrLine);
            if (vecstrLine.empty())
                continue;

            m_vecvecpGenes.resize(m_vecvecpGenes.size() + 1);
            {
                vector < CGene * > &vecpGenes = m_vecvecpGenes[m_vecvecpGenes.size() - 1];

                for (i = 0; i < vecstrLine.size(); ++i) {
                    if (vecstrLine[i].length() == 0)
                        continue;
                    if ((j = vecstrLine[i].find(c_cOrgSep)) == string::npos) {
                        g_CatSleipnir().warn("COrthology::Open( ) illegal gene token: %s",
                                             vecstrLine[i].c_str());
                        continue;
                    }
                    strOrganism = vecstrLine[i].substr(0, j);
                    strGene = vecstrLine[i].substr(j + 1);
                    for (j = 0; j < m_vecstrOrganisms.size(); ++j)
                        if (strOrganism == m_vecstrOrganisms[j])
                            break;
                    if (j < m_vecpGenomes.size())
                        pGenome = m_vecpGenomes[j];
                    else {
                        m_vecpGenomes.push_back(pGenome = new CGenome());
                        m_vecstrOrganisms.push_back(strOrganism);
                    }
                    vecpGenes.push_back(pGene = &pGenome->AddGene(strGene));
                    m_mapGenes[pGene] = j;
                }
            }
        }
        delete[] acBuf;

        return true;
    }

/*!
 * \brief
 * Saves the current orthology to the given stream.
 * 
 * \param ostm
 * Stream into which orthology is saved.
 * 
 * \see
 * Open
 */
    void COrthology::Save(std::ostream &ostm) const {
        size_t i, j;

        for (i = 0; i < m_vecvecpGenes.size(); ++i) {
            for (j = 0; j < m_vecvecpGenes[i].size(); ++j)
                ostm << (j ? "\t" : "") << m_vecstrOrganisms[((COrthology *) this)->m_mapGenes[
                        m_vecvecpGenes[i][j]]] << c_cOrgSep << m_vecvecpGenes[i][j]->GetName();

            ostm << endl;
        }
    }

}

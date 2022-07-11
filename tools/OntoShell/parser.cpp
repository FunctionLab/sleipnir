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
#include "parser.h"

const char    CParser::SLocation::c_szRoot[] = "/";
const char    CParser::c_szDot[] = ".";
const char    CParser::c_szDotDot[] = "..";
const char *CParser::c_aszParsers[] = {"cat", "cd", "find", "help", "ls", "parentage", NULL};

CParser::SLocation::SLocation() : m_pOnto(NULL), m_iNode(-1) {}

string CParser::SLocation::ToString(bool fGloss) const {

    return (m_pOnto ? ((m_iNode == -1) ? m_pOnto->GetID() :
                       m_pOnto->GetID(m_iNode) + (fGloss ? '\t' + m_pOnto->GetGloss(m_iNode) : "")) :
            c_szRoot);
}

bool CParser::SLocation::IsValid() const {

    return (m_pOnto || (m_iNode == -1));
}

void CParser::SLocation::Invalidate() {

    m_pOnto = NULL;
    m_iNode = 0;
}

bool CParser::SLocation::operator==(const SLocation &sLoc) const {

    return ((m_pOnto == sLoc.m_pOnto) && (m_iNode == sLoc.m_iNode));
}

const char *CParser::GetCommand(size_t iCommand) {

    return c_aszParsers[iCommand];
}

bool CParser::IsRooted(const string &strLoc) {

    return (strLoc.length() && (strLoc[0] == c_cSep));
}

bool CParser::SplitLocation(const string &strLoc, vector <string> &vecstrPath) {
    size_t i;
    string strCur;

    if (!strLoc.size())
        return true;

    i = 0;
    while (true) {
        strCur.clear();
        for (; (i < strLoc.size()) && (strLoc[i] == c_cSep); ++i);
        if (i >= strLoc.size())
            break;
        for (; (i < strLoc.size()) && (strLoc[i] != c_cSep); ++i)
            strCur += strLoc[i];
        vecstrPath.push_back(strCur);
    }

    return true;
}

CParser::SLocation CParser::GetLocation(const vector<const IOntology *> &vecpOntologies,
                                        const string &strLoc, bool fLast, const SLocation *psLoc) {
    SLocation sRet;
    vector <string> vecstrPath;
    size_t i;

    if (!IsRooted(strLoc) && psLoc)
        sRet = *psLoc;
    if (!SplitLocation(strLoc, vecstrPath)) {
        sRet.Invalidate();
        return sRet;
    }

    if (!fLast) {
        if (vecstrPath.empty())
            return sRet;
        if (strLoc[strLoc.size() - 1] != c_cSep)
            vecstrPath.resize(vecstrPath.size() - 1);
    }
    for (i = 0; i < vecstrPath.size(); ++i)
        if (!MoveLocation(sRet, vecstrPath[i], vecpOntologies)) {
            sRet.Invalidate();
            break;
        }

    return sRet;
}

bool CParser::MoveLocation(SLocation &sLoc, const string &strPath,
                           const vector<const IOntology *> &vecpOntologies) {
    size_t i, iNode;

    if (strPath == c_szDot)
        return true;

    if (sLoc.m_pOnto) {
        if (strPath == c_szDotDot) {
            if (sLoc.m_iNode == -1) {
                sLoc.m_pOnto = NULL;
                return true;
            }
            switch (sLoc.m_pOnto->GetParents(sLoc.m_iNode)) {
                case 0:
                    sLoc.m_iNode = -1;
                    return true;

                case 1:
                    sLoc.m_iNode = sLoc.m_pOnto->GetParent(sLoc.m_iNode, 0);
                    return true;
            }
        } else if ((iNode = sLoc.m_pOnto->GetNode(strPath)) != -1) {
            sLoc.m_iNode = iNode;
            return true;
        }
    } else
        for (i = 0; i < vecpOntologies.size(); ++i)
            if (strPath == vecpOntologies[i]->GetID()) {
                sLoc.m_pOnto = vecpOntologies[i];
                sLoc.m_iNode = -1;
                return true;
            }

    return false;
}

void CParser::CollectGenes(const vector <SLocation> &vecLocations,
                           TSetPGenes &setpGenes) {
    size_t i, j;

    for (i = 0; i < vecLocations.size(); ++i) {
        const SLocation &sLoc = vecLocations[i];

        if (sLoc.m_iNode != -1)
            for (j = 0; j < sLoc.m_pOnto->GetGenes(sLoc.m_iNode); ++j)
                setpGenes.insert(&sLoc.m_pOnto->GetGene(sLoc.m_iNode, j));
    }
}

CParser::CParser(const IOntology **apOntologies, const CGenome &Genome) :
        m_Genome(Genome) {
    size_t i;

    if (apOntologies)
        for (i = 0; apOntologies[i]; ++i)
            m_vecpOntologies.push_back(apOntologies[i]);
}

size_t CParser::GetOntologies() const {

    return m_vecpOntologies.size();
}

const IOntology *CParser::GetOntology(size_t iOnto) const {

    return m_vecpOntologies[iOnto];
}

const CGenome &CParser::GetGenome() const {

    return m_Genome;
}

bool CParser::Recurse(SLocation sLoc, bool fRecursive, bool fZeroes,
                      vector <SLocation> &vecVisited) const {
    const IOntology *pOnto;
    size_t i, j;
    bool fOK;

    if (!sLoc.IsValid())
        return false;
    if (!fZeroes) {
        fOK = false;
        if (sLoc.m_pOnto)
            fOK = !!((sLoc.m_iNode == -1) ? m_Genome.CountGenes(sLoc.m_pOnto) :
                     sLoc.m_pOnto->GetGenes(sLoc.m_iNode, true));
        else
            for (i = 0; i < m_vecpOntologies.size(); ++i)
                if (m_Genome.CountGenes(m_vecpOntologies[i])) {
                    fOK = true;
                    break;
                }
        if (!fOK)
            return true;
    }
    for (i = 0; i < vecVisited.size(); ++i)
        if (vecVisited[i] == sLoc)
            return true;
    vecVisited.push_back(sLoc);
    if (!fRecursive)
        return true;

    if ((pOnto = sLoc.m_pOnto)) {
        if (sLoc.m_iNode == -1) {
            for (i = 0; i < pOnto->GetNodes(); ++i)
                if (!pOnto->GetParents(i)) {
                    sLoc.m_iNode = i;
                    Recurse(sLoc, fRecursive, fZeroes, vecVisited);
                }
        } else {
            j = sLoc.m_iNode;
            for (i = 0; i < pOnto->GetChildren(j); ++i) {
                sLoc.m_iNode = pOnto->GetChild(j, i);
                Recurse(sLoc, fRecursive, fZeroes, vecVisited);
            }
        }
    } else
        for (i = 0; i < m_vecpOntologies.size(); ++i) {
            sLoc.m_pOnto = m_vecpOntologies[i];
            Recurse(sLoc, fRecursive, fZeroes, vecVisited);
        }

    return true;
}

struct SSortFind {
    bool operator()(const STermFound &sOne, const STermFound &sTwo) {

        return (sOne.m_dP < sTwo.m_dP);
    }
};

void CParser::TermFinder(const CGenes &Genes, float dP, const CGenes &GenesBkg, bool fBonferroni, bool fSibs,
                         bool fBackground, vector <size_t> &veciOnto, vector <STermFound> &vecsTerms) const {
    size_t i, j;
    vector <STermFound> vecsCur;
    SSortFind sSort;

    veciOnto.resize(m_vecpOntologies.size());
    for (i = 0; i < m_vecpOntologies.size(); ++i) {
        vecsCur.clear();
        m_vecpOntologies[i]->TermFinder(Genes, vecsCur, fBonferroni,
                                        fSibs, fBackground, dP, GenesBkg.GetGenes() ? &GenesBkg : NULL);
        sort(vecsCur.begin(), vecsCur.end(), sSort);
        for (j = 0; j < vecsCur.size(); ++j)
            vecsTerms.push_back(vecsCur[j]);
        veciOnto[i] = vecsTerms.size();
    }
}

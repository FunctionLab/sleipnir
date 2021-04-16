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
* Changes made Jun 2012 by Qian Zhu
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "stdafx.h"
#include "database.h"
#include "meta.h"
#include "bayesnetint.h"
#include "datapair.h"
#include <filesystem>

namespace Sleipnir {

    const char    CDatabaseImpl::c_acDAB[] = ".dab";
    const char    CDatabaseImpl::c_acQDAB[] = ".qdab";
    const char    CDatabaseImpl::c_acExtension[] = ".db";

///////////////////////////////////////////////////////////////////////////////
// CDatabaselet
///////////////////////////////////////////////////////////////////////////////

    CDatabaselet::CDatabaselet(bool useNibble) {
        m_useNibble = useNibble;
        m_pmutx = new pthread_mutex_t();
        pthread_mutex_init(m_pmutx, NULL);
    }

    CDatabaselet::~CDatabaselet() {

        pthread_mutex_destroy(m_pmutx);
        delete m_pmutx;
        if (m_fstm.is_open()) {
            m_fstm.close();
        }
    }

/* original method for initializing databaselets, including writing header + pre-allocation */
// TODO (GW) - rename to Create() or indicate file will be written
    bool CDatabaselet::Open(const std::string &strFile, const std::vector <std::string> &vecstrGenes,
                            uint32_t iGenes, uint32_t iDatasets) {
        uint32_t iSize;
        size_t i;
        char *acFiller;

        m_fstm.clear();
        /* Open with overwriting */
        m_fstm.open(strFile.c_str(), ios_base::in | ios_base::out | ios_base::binary | ios_base::trunc);

        if (!m_fstm.is_open()) {
            g_CatSleipnir().error("CDatabaselet::Open( %s, %u, %u ) open failed", strFile.c_str(), iGenes,
                                  iDatasets);
            return false;
        }

        m_iGenes = iGenes;
        m_iDatasets = iDatasets;
        m_vecstrGenes.resize(vecstrGenes.size());
        copy(vecstrGenes.begin(), vecstrGenes.end(), m_vecstrGenes.begin());

        //allocate space
        m_fstm.write((char *) &m_iHeader, sizeof(m_iHeader));
        m_fstm.write((char *) &m_iGenes, sizeof(m_iGenes));
        m_fstm.write((char *) &m_iDatasets, sizeof(m_iDatasets));

        iSize = m_vecstrGenes.size();

        m_fstm.write((char *) &iSize, sizeof(iSize));

        //write gene-name for only the genes in the databaselets
        m_iHeader = sizeof(m_iHeader) + sizeof(m_iGenes) + sizeof(m_iDatasets) + sizeof(iSize);
        for (i = 0; i < m_vecstrGenes.size(); ++i) {
            m_fstm.write(m_vecstrGenes[i].c_str(), m_vecstrGenes[i].size() + 1);
            m_iHeader += m_vecstrGenes[i].size() + 1;
        }

        m_fstm.seekp(0, ios_base::beg);
        m_fstm.write((char *) &m_iHeader, sizeof(m_iHeader));

        //pre-allocations
        m_fstm.seekp(m_iHeader, ios_base::beg);
        acFiller = new char[GetSizeGene()];
        memset(acFiller, -1, GetSizeGene());
        for (i = 0; i < m_vecstrGenes.size(); ++i) {
            m_fstm.write(acFiller, GetSizeGene());
        }
        delete[] acFiller;
        SetFile(strFile);

        return true;

    }

    void CDatabaselet::Write(char *Qt, const size_t &iSize, const size_t offset) {
        m_fstm.seekg(m_iHeader + offset, ios_base::beg);
        m_fstm.write(Qt, iSize);
    }

/* simply opens the file without overwriting */
    bool CDatabaselet::OpenNoOverwrite() {
        m_fstm.clear();
        /* open without overwriting */
        m_fstm.open(strFileName.c_str(), ios_base::in | ios_base::out | ios_base::binary);

        if (!m_fstm.is_open()) {
            g_CatSleipnir().error("CDatabaselet::Open( %s ) open failed", strFileName.c_str());
            return false;
        }
        return true;
    }

    bool CDatabaselet::OpenWrite(unsigned char bValue, size_t iOffset, ENibbles eNibbles,
                                 unsigned char *abImage) {
        unsigned char b;

        if (m_useNibble == 0) {
            eNibbles = ENibblesBoth;
        }

        if (abImage)
            iOffset -= m_iHeader;
        if (eNibbles != ENibblesBoth) {
            if (abImage)
                b = abImage[iOffset];
            else {
                m_fstm.seekg(iOffset, ios_base::beg);
                b = m_fstm.get();
            }
        }

        switch (eNibbles) {
            case ENibblesLow:
                b = (bValue & 0xF) | (b & 0xF0);
                break;

            case ENibblesHigh:
                b = (b & 0xF) | (bValue << 4);
                break;

            case ENibblesBoth:
                b = bValue;
                break;
        }
        if (abImage)
            abImage[iOffset] = b;
        else {
            m_fstm.seekp(iOffset, ios_base::beg);
            m_fstm.put(b);
        }

        return true;
    }

    /* original file writing method */
    // vecData - is a vector of matrices, one matrix per input dab file in this group
    // iBaseGenes - is the row offset within the matrices vecData of the primary genes in this DB file
    // iBaseDatasets - is the global datdaset index of the first dataset contained in vedData
    // fBuffer - boolean whether to buffer data in memory or write directly to the file
    bool CDatabaselet::Open(const vector <CCompactFullMatrix> &vecData, size_t iBaseGenes, size_t iBaseDatasets,
                            bool fBuffer) {
        unsigned char *abImage;
        size_t iSize, iDatum, iGeneOne, iGeneTwo;
        unsigned char bOne, bTwo;

        if (fBuffer) {
            //iBaseGenes: gene id of first gene in each databaselet
            //iDataset: dataset id
            //printf("Number: %d %d %d %d\n", GetSizeGene(), GetSizePair(), iBaseGenes, iBaseDatasets);
            abImage = new unsigned char[iSize = (GetSizeGene() * m_vecstrGenes.size())];
            m_fstm.seekg(m_iHeader, ios_base::beg);
            m_fstm.read((char *) abImage, iSize);
        } else
            abImage = nullptr;

        //vecData: # of genes in databaselet x # of genes user's list
        // The gene values are ouput to the DB file in pairs (for consecutive datasets)
        //   to allow writing two nibbles at once to the same byte. This is because
        //   dataset values (for gene-pairs) are written consecutively with each other, 
        //   i.e. each row represents a gene-pair and each entry the value for a different dataset.
        
        // if the BaseDataset index is odd then we need to write the first one singleton so 
        //   the rest will align two at a time. This happens when the input block size is odd
        //   (i.e. num Dabs to handle at once is odd) with every other call to this function.
        if (iBaseDatasets % 2) {
            //iGeneOne: iterate over all genes in this databaselet (# of genes in each databaselet)
            for (iGeneOne = 0; iGeneOne < GetGenes(); ++iGeneOne) {
                //iGeneTwo: iterate overall genes in user's gene list
                for (iGeneTwo = 0; iGeneTwo < vecData[0].GetColumns(); ++iGeneTwo) {
                    //bOne, get the value of the gene located at the position (iBaseGene + iGeneOne, iGeneTwo)
                    if ((bOne = vecData[0].Get(iBaseGenes + iGeneOne, iGeneTwo))) {
                        //Offset is: m_iHeader + (GetSizeGene() * iOne) + (GetSizePair() * iTwo) + iDataset (for byte case)
                        // ENibblesHigh writes only to the high nibble of the byte
                        OpenWrite(bOne - 1, GetOffset(iGeneOne, iGeneTwo, iBaseDatasets), ENibblesHigh, abImage);
                    }
                }
            }
        }

        // Loop through the rest of the datasets two at a time to facilitate writting nibbles together
        for (iDatum = (iBaseDatasets % 2); (iDatum + 1) < vecData.size(); iDatum += 2) {
            for (iGeneOne = 0; iGeneOne < GetGenes(); ++iGeneOne) {
                for (iGeneTwo = 0; iGeneTwo < vecData[iDatum].GetColumns(); ++iGeneTwo) {
                    bOne = vecData[iDatum].Get(iBaseGenes + iGeneOne, iGeneTwo);
                    bTwo = vecData[iDatum + 1].Get(iBaseGenes + iGeneOne, iGeneTwo);
                    if (!(bOne || bTwo))
                        continue;
                    bOne -= 1;
                    bTwo -= 1;
                    if (m_useNibble) {
                        OpenWrite((bOne & 0xF) | (bTwo << 4), GetOffset(iGeneOne, iGeneTwo, iBaseDatasets +
                                                                                            iDatum), ENibblesBoth,
                                  abImage);
                    } else {
                        // ENibblesBoth writes only to both the high and low nibble of the byte together
                        OpenWrite(bOne, GetOffset(iGeneOne, iGeneTwo, iBaseDatasets + iDatum), ENibblesBoth,
                                  abImage);
                        OpenWrite(bTwo, GetOffset(iGeneOne, iGeneTwo, iBaseDatasets + iDatum + 1), ENibblesBoth,
                                  abImage);
                    }
                }
            }
        }

        // if number of datasets in a block is odd then we need to write the final
        //   dataset values. In this case iDatum will be one less than vecData.size()
        if (iDatum < vecData.size()) {
            for (iGeneOne = 0; iGeneOne < GetGenes(); ++iGeneOne) {
                for (iGeneTwo = 0; iGeneTwo < vecData[iDatum].GetColumns(); ++iGeneTwo) {
                    if ((bOne = vecData[iDatum].Get(iBaseGenes + iGeneOne, iGeneTwo))) {
                        // ENibblesLow writes only to the low nibble of the byte
                        OpenWrite(bOne - 1, GetOffset(iGeneOne, iGeneTwo, iBaseDatasets + iDatum), ENibblesLow,
                                  abImage);
                    }
                }
            }
        }
        if (fBuffer) {
            m_fstm.seekp(m_iHeader, ios_base::beg);
            m_fstm.write((char *) abImage, iSize);
            delete[] abImage;
        }

        return true;
    }

    bool CDatabase::GetGene(const string &strGene, vector<unsigned char> &vecbData) const {
        size_t iGene = this->GetGene(strGene);
        if (iGene == -1) {
            cerr << "Gene not found" << endl;
            return false;
        }
        m_vecpDBs[iGene % m_vecpDBs.size()]->Get(iGene / m_vecpDBs.size(), vecbData);
        return true;
    }

    bool CDatabase::GetGene(const size_t &iGene, vector<unsigned char> &vecbData) const {
        m_vecpDBs[iGene % m_vecpDBs.size()]->Get(iGene / m_vecpDBs.size(), vecbData);
        return true;
    }

/* mainly used by SeekMiner */
    bool CDatabaselet::Get(size_t iOne, vector<unsigned char> &vecbData) {
        size_t iSize;
        size_t i, j;
        size_t offset1 = GetOffset(iOne);
        if (this->m_fstm.is_open()) {
            cerr << "file is already opened" << endl;
        } else {
            m_fstm.clear();
            m_fstm.open(this->strFileName.c_str(), ios_base::binary | ios_base::in);
        }

        this->m_fstm.seekg(offset1, ios_base::beg);
        unsigned char *abImage = (unsigned char *) malloc(iSize = GetSizeGene());
        this->m_fstm.read((char *) abImage, iSize);

        m_fstm.close();

        if (m_useNibble == false) {
            vecbData.clear();
            vecbData.resize(iSize);
            for (i = 0; i < iSize; i++) {
                vecbData[i] = abImage[i];
            }
        } else {
            vecbData.clear();
            vecbData.resize(m_iDatasets * m_iGenes);

            for (j = 0; j < m_iGenes; j++) {
                // Loop through each dataset value for the gene pair
                for (i = 0; i < GetSizePair(); i++) {
                    size_t offset = GetOffset(iOne, j) - m_iHeader;  // This should be moved to outer loop
                    unsigned char b = abImage[offset + i];
                    unsigned char bValue = -1;
                    if ((bValue = (b & 0xF)) == 0xF) {
                        bValue = -1;
                    }
                    vecbData[j * m_iDatasets + 2 * i] = bValue;

                    if ((bValue = ((b >> 4) & 0xF)) == 0xF) {
                        bValue = -1;
                    }

                    if ((2 * i + 1) == m_iDatasets) {
                        break;
                    }
                    vecbData[j * m_iDatasets + (2 * i) + 1] = bValue;
                }
            }
        }
        free(abImage);
        return true;
    }

    bool CDatabaselet::Get(size_t iOne, size_t iTwo,
                           vector<unsigned char> &vecbData, unsigned char *charImage) {
        size_t offset = GetOffset(iOne, iTwo) - m_iHeader;
        return Get(offset, vecbData, charImage);
    }

    bool CDatabaselet::Get(size_t offset, vector<unsigned char> &vecbData, unsigned char *charImage) {
        size_t i;

        if (this->m_useNibble == false) {
            vecbData.clear();
            vecbData.resize(GetSizePair());

            for (i = 0; i < vecbData.size(); i++) {
                vecbData[i] = charImage[offset + i];
            }
        } else {
            vecbData.clear();
            vecbData.resize(m_iDatasets);


            for (i = 0; i < GetSizePair(); i++) {
                unsigned char b = charImage[offset + i];
                unsigned char bValue = -1;
                if ((bValue = (b & 0xF)) == 0xF) {
                    bValue = -1;
                }
                vecbData[2 * i] = bValue;

                if ((bValue = ((b >> 4) & 0xF)) == 0xF) {
                    bValue = -1;
                }

                if ((2 * i + 1) == m_iDatasets) {
                    break;
                }
                vecbData[(2 * i) + 1] = bValue;
            }
        }

        return true;
    }

/*	static function, combine multiple databaselets (that share the same genes, ie m_vecStrGenes),
	and output result to a single file, or output one-gene per file (if databaselet contains multiple genes)
 	bSplit: whether or not to output one-gene per file
	Works for both nibble and byte
*/
    bool CDatabaselet::Combine(std::vector<CDatabaselet *> &vecDatabaselet,
                               std::string strOutDirectory, vector <string> &vecstrGenes, bool bSplit) {

        /* for checking on consistency of databaselets */
        bool bIsConsistent = true;
        bool fUseNibble;

        size_t i, j;
        uint32_t iGenes, iDatasets;

        CDatabaselet *first = vecDatabaselet[0];
        fUseNibble = first->m_useNibble;

        iGenes = first->GetGenes();
        iDatasets = first->GetDatasets();

        vector <string> vecGenes;
        vecGenes.resize(iGenes);

        for (i = 0; i < iGenes; i++) {
            vecGenes[i] = first->GetGene(i);
        }

        for (i = 1; bIsConsistent && i < vecDatabaselet.size(); i++) {
            if (iGenes != vecDatabaselet[i]->GetGenes() || fUseNibble != vecDatabaselet[i]->m_useNibble) {
                bIsConsistent = false;
                break;
            }
            for (j = 0; j < iGenes; j++) {
                if (vecGenes[j] != vecDatabaselet[i]->GetGene(j)) {
                    bIsConsistent = false;
                    break;
                }
            }
            iDatasets += vecDatabaselet[i]->GetDatasets();
        }

        if (!bIsConsistent) {
            cerr << "Databaselets are not consistent!" << endl;
            return false;
        }

        /* adjustable parameter */
        int BUFFER = 80;

        if (BUFFER > iGenes) {
            BUFFER = iGenes;
        }


        vector <size_t> sizes;
        int numTimes = iGenes / BUFFER;
        for (i = 0; i < numTimes; i++) {
            sizes.push_back(BUFFER);
        }
        if (iGenes % BUFFER > 0) {
            numTimes++;
            sizes.push_back(iGenes % BUFFER);
        }

        size_t bb = 0;
        size_t sofar = 0;
        fprintf(stderr, "Number of times: %d\n", numTimes);
        for (bb = 0; bb < numTimes; sofar += sizes[bb], bb++) {
            fprintf(stderr, "bb sizes[bb] sofar: %zu %lu %zu\n", bb, sizes[bb], sofar);

            fprintf(stderr, "Started allocated memory %lu\n", bb);
            /* load all Databaselets into memory, for efficiency */
            auto **charImages =
                    (unsigned char **) malloc(vecDatabaselet.size() * sizeof(unsigned char *));
            size_t totSize = 0;
            size_t numPairs = first->m_iGenes * sizes[bb];

            for (i = 0; i < vecDatabaselet.size(); i++) {
                totSize += vecDatabaselet[i]->GetSizePair() * numPairs;
            }
            charImages[0] = (unsigned char *) malloc(totSize * sizeof(unsigned char));
            for (i = 1; i < vecDatabaselet.size(); i++) {
                charImages[i] = charImages[i - 1] + vecDatabaselet[i - 1]->GetSizePair() * numPairs;
            }

            fprintf(stderr, "Finished allocated memory\n");

            fprintf(stderr, "Started reading into memory\n");

            /* read databaselet into charImages */
            for (i = 0; i < vecDatabaselet.size(); i++) {
                CDatabaselet *current = vecDatabaselet[i];
                if (!current->m_fstm.is_open()) {
                    cerr << "CDatabaselet is not opened. Opening: " << current->strFileName << endl;
                    current->m_fstm.open(current->strFileName.c_str(), ios_base::binary | ios_base::in);
                }
                current->m_fstm.seekg(current->m_iHeader + current->GetSizePair() * first->m_iGenes * sofar,
                                      ios_base::beg);
                current->m_fstm.read((char *) charImages[i], current->GetSizePair() * numPairs);
            }

            fprintf(stderr, "Finished reading into memory\n");

            //Global gene name mapping
            map <string, size_t> mapstrintGenes;
            for (i = 0; i < vecstrGenes.size(); i++) {
                mapstrintGenes[vecstrGenes[i]] = i;
            }

            char acNumber[16];

            /* splitting to one gene per file after combine */
            if (bSplit) {

                for (i = 0; i < sizes[bb]; i++) {
                    /* open a new Databaselet containing only one gene */
                    string thisGene = first->GetGene(sofar + i);
                    size_t iGeneID = mapstrintGenes[thisGene];
                    sprintf(acNumber, "%08lu", iGeneID);

                    string path = strOutDirectory + "/" + acNumber + ".db";
                    vector <string> vecstrThisGene;
                    vecstrThisGene.push_back(thisGene);

                    //fprintf(stderr, "Setting up gene %d in %d\n", i, bb);

                    CDatabaselet DBS(first->m_useNibble);
                    DBS.Open(path.c_str(), vecstrThisGene, first->m_iGenes, iDatasets);

                    /* Create a new Databaselet */
                    size_t iSize;
                    unsigned char *abImage = (unsigned char *)
                            malloc(iSize = (DBS.GetSizeGene() * DBS.m_vecstrGenes.size()));
                    size_t iDatum;
                    size_t iGeneOne, iGeneTwo;
                    size_t offset2, offset3;
                    iGeneOne = i;

                    //fprintf(stderr, "Finished setting up abImage gene %d in %d\n", i, bb);

                    if (first->m_useNibble == false) {
                        /* m_iGenes is all the genes in the genome */
                        for (iGeneTwo = 0; iGeneTwo < first->m_iGenes; ++iGeneTwo) {
                            offset2 = DBS.GetSizePair() * iGeneTwo;
                            int totalSum = 0;
                            for (iDatum = 0; iDatum < vecDatabaselet.size(); iDatum++) {
                                vector<unsigned char> vc;
                                CDatabaselet *current = vecDatabaselet[iDatum];
                                size_t offset_c = current->GetSizeGene() * i + current->GetSizePair() * iGeneTwo;
                                current->Get(offset_c, vc, charImages[iDatum]);
                                offset3 = offset2 + totalSum;
                                for (j = 0; j < vc.size(); j++) {
                                    abImage[offset3 + j] = vc[j];
                                }
                                totalSum += vc.size();
                            }
                        }
                    } else {
                        size_t j;
                        unsigned char *abImage2 = (unsigned char *) malloc(iDatasets);

                        /* m_iGenes is all the genes in the genome */
                        for (iGeneTwo = 0; iGeneTwo < first->m_iGenes; ++iGeneTwo) {
                            offset2 = DBS.GetSizePair() * iGeneTwo;
                            int totalSum = 0;
                            for (iDatum = 0; iDatum < vecDatabaselet.size(); iDatum++) {
                                vector<unsigned char> vc;
                                CDatabaselet *current = vecDatabaselet[iDatum];
                                size_t offset_c = current->GetSizeGene() * i + current->GetSizePair() * iGeneTwo;
                                current->Get(offset_c, vc, charImages[iDatum]);
                                offset3 = totalSum;
                                for (j = 0; j < vc.size(); j++) {
                                    abImage2[offset3 + j] = vc[j];
                                }
                                totalSum += vc.size();
                            }
                            for (j = 0; j + 1 < iDatasets; j += 2) {
                                abImage[offset2 + j / 2] = (abImage2[j] & 0xF) | (abImage2[j + 1] << 4);
                            }
                            if (j < iDatasets) {
                                unsigned char bValue = abImage2[iDatasets - 1];
                                unsigned char b = 255;
                                abImage[offset2 + j / 2] = (bValue & 0xF) | (b & 0xF0);
                            }
                        }

                        free(abImage2);
                    }

                    //fprintf(stderr, "Writing abImage gene %d in %d\n", i, bb);

                    /* close fstream */
                    if (!DBS.m_fstm.is_open()) {
                        cerr << "CDatabaselet is not opened. Opening..." << endl;
                        DBS.m_fstm.open(DBS.strFileName.c_str(), ios_base::binary | ios_base::in);
                    }
                    DBS.m_fstm.seekp(DBS.m_iHeader, ios_base::beg);
                    DBS.m_fstm.write((char *) abImage, iSize);
                    DBS.m_fstm.close();
                    free(abImage);

                    fprintf(stderr, "Finished writing abImage gene %lu (of %lu) in %lu (of %lu)\n", i,
                            sizes[bb], bb, sizes.size());

                }

                /* do not split, just combine into one file */
            } else {

                vector <string> strTok;
                CMeta::Tokenize(first->strFileName.c_str(), strTok, "/");
                string path = strOutDirectory + "/" + strTok[strTok.size() - 1];

                CDatabaselet DBS(first->m_useNibble);
                if (bb == 0) {
                    DBS.Open(path.c_str(), first->m_vecstrGenes, first->m_iGenes, iDatasets);
                } else {
                    DBS.SetFile(path);
                    DBS.OpenNoOverwrite();
                }

                size_t iDatum;
                size_t iSize;
                auto *abImage = (unsigned char *)
                        malloc(iSize = (DBS.GetSizeGene() * sizes[bb]));
                size_t iGeneOne, iGeneTwo;
                size_t offset1, offset2, offset3;

                fprintf(stderr, "GetSizeGene() %zu\n", DBS.GetSizeGene());
                fprintf(stderr, "GetSizePair() %zu\n", DBS.GetSizePair());
                fprintf(stderr, "first->m_iGenes %d\n", first->m_iGenes);
                if (!first->m_useNibble) {
                    for (iGeneOne = 0; iGeneOne < sizes[bb]; ++iGeneOne) {
                        offset1 = DBS.GetSizeGene() * iGeneOne;
                        for (iGeneTwo = 0; iGeneTwo < first->m_iGenes; ++iGeneTwo) {
                            offset2 = DBS.GetSizePair() * iGeneTwo;
                            int totalSum = 0;
                            for (iDatum = 0; iDatum < vecDatabaselet.size(); iDatum++) {
                                vector<unsigned char> vc;
                                CDatabaselet *current = vecDatabaselet[iDatum];
                                //size_t offset_c = current->GetSizeGene() * iGeneOne + current->GetSizePair() * iGeneTwo;
                                //current->Get( offset1 + offset2, vc, charImages[iDatum]);
                                //current->Get( iGeneOne, iGeneTwo, vc, charImages[iDatum]);
                                current->Get(current->GetSizeGene() * iGeneOne + current->GetSizePair() * iGeneTwo, vc,
                                             charImages[iDatum]);
                                offset3 = offset1 + offset2 + totalSum;
                                for (j = 0; j < vc.size(); j++) {
                                    abImage[offset3 + j] = vc[j];
                                }
                                totalSum += vc.size();
                            }
                        }
                    }
                } else {
                    size_t j;
                    unsigned char *abImage2 = (unsigned char *)
                            malloc(DBS.m_iDatasets);
                    /* m_iGenes is all the genes in the genome */
                    for (iGeneOne = 0; iGeneOne < sizes[bb]; ++iGeneOne) {
                        offset1 = DBS.GetSizeGene() * iGeneOne;
                        for (iGeneTwo = 0; iGeneTwo < first->m_iGenes; ++iGeneTwo) {
                            offset2 = DBS.GetSizePair() * iGeneTwo;
                            int totalSum = 0;
                            for (iDatum = 0; iDatum < vecDatabaselet.size(); iDatum++) {
                                vector<unsigned char> vc;
                                CDatabaselet *current = vecDatabaselet[iDatum];
                                current->Get(current->GetSizeGene() * iGeneOne + current->GetSizePair() * iGeneTwo, vc,
                                             charImages[iDatum]);
                                //current->Get( offset1 + offset2, vc, charImages[iDatum]);
                                //current->Get( iGeneOne, iGeneTwo, vc, charImages[iDatum]);
                                offset3 = totalSum;
                                for (j = 0; j < vc.size(); j++) {
                                    abImage2[offset3 + j] = vc[j];
                                }
                                totalSum += vc.size();
                            }
                            for (j = 0; j + 1 < iDatasets; j += 2) {
                                abImage[offset1 + offset2 + j / 2] = (abImage2[j] & 0xF) | (abImage2[j + 1] << 4);
                            }
                            if (j < iDatasets) {
                                unsigned char bValue = abImage2[iDatasets - 1];
                                unsigned char b = 255;
                                abImage[offset1 + offset2 + j / 2] = (bValue & 0xF) | (b & 0xF0);
                            }
                        }
                    }
                    free(abImage2);
                }

                /* close the databaselet */
                if (!DBS.m_fstm.is_open()) {
                    cerr << "CDatabaselet is not opened." << endl;
                    DBS.m_fstm.open(DBS.strFileName.c_str(), ios_base::binary | ios_base::in);
                }
                //DBS.m_fstm.seekp( DBS.m_iHeader, ios_base::beg );
                DBS.m_fstm.seekg(DBS.m_iHeader + DBS.GetSizePair() * first->m_iGenes * sofar, ios_base::beg);
                DBS.m_fstm.write((char *) abImage, iSize);
                DBS.m_fstm.close();
                free(abImage);
            }

            free(charImages[0]);
            free(charImages);

        }

        return true;
    }

    bool CDatabaselet::Get(size_t iOne, size_t iTwo, vector<unsigned char> &vecbData) const {
        size_t i;

        i = vecbData.size();
        vecbData.resize(i + GetSizePair());
        pthread_mutex_lock(m_pmutx);

        m_fstm.seekg(GetOffset(iOne, iTwo), ios_base::beg);
        m_fstm.read((char *) &vecbData[i], GetSizePair());
        pthread_mutex_unlock(m_pmutx);

        return true;
    }

    bool CDatabaselet::Get(size_t iGene, vector<unsigned char> &vecbData, bool fReplace) const {
        size_t i;

        i = fReplace ? 0 : vecbData.size();
        vecbData.resize(i + GetSizeGene());
        pthread_mutex_lock(m_pmutx);

        m_fstm.seekg(GetOffset(iGene), ios_base::beg);
        m_fstm.read((char *) &vecbData[i], GetSizeGene());
        pthread_mutex_unlock(m_pmutx);

        return true;
    }

    bool CDatabaselet::Get(size_t iGene, const vector <size_t> &veciGenes, vector<unsigned char> &vecbData,
                           bool fReplace) const {
        size_t i, iOffset;

        iOffset = fReplace ? 0 : vecbData.size();
        vecbData.resize(iOffset + (veciGenes.size() * GetSizePair()));
        pthread_mutex_lock(m_pmutx);
        for (i = 0; i < veciGenes.size(); ++i, iOffset += GetSizePair()) {
            m_fstm.seekg(GetOffset(iGene, veciGenes[i]), ios_base::beg);
            m_fstm.read((char *) &vecbData[iOffset], GetSizePair());
        }
        pthread_mutex_unlock(m_pmutx);

        return true;
    }

    bool CDatabaselet::Open(const std::string &strFile) {
        uint32_t iSize;
        char *acBuffer;
        char *pc;
        size_t i;

        m_fstm.clear();
        m_fstm.open(strFile.c_str(), ios_base::binary | ios_base::in);
        if (!m_fstm.is_open()) {
            g_CatSleipnir().error("CDatabaselet::Open( %s ) open failed", strFile.c_str());
            return false;
        }

        m_fstm.read((char *) &m_iHeader, sizeof(m_iHeader));
        m_fstm.read((char *) &m_iGenes, sizeof(m_iGenes));
        m_fstm.read((char *) &m_iDatasets, sizeof(m_iDatasets));
        m_fstm.read((char *) &iSize, sizeof(iSize));

        acBuffer = new char[m_iHeader];
        m_fstm.read(acBuffer, m_iHeader - sizeof(m_iHeader) - sizeof(m_iGenes) - sizeof(m_iDatasets) - sizeof(iSize));
        m_vecstrGenes.resize(iSize);
        for (i = 0, pc = acBuffer; i < m_vecstrGenes.size(); pc += m_vecstrGenes[i++].length() + 1)
            m_vecstrGenes[i] = pc;
        delete[] acBuffer;

        SetFile(strFile);

        return true;
    }

///////////////////////////////////////////////////////////////////////////////
// CDatabase
///////////////////////////////////////////////////////////////////////////////

    bool CDatabase::Open(const std::vector <std::string> &vecstrGenes, const std::string &strInputDirectory,
                         const IBayesNet *pBayesNet, const std::string &strOutputDirectory, size_t iFiles,
                         const map <string, size_t> &mapstriZeros) {
        vector <string> vecstrNodes, vecstrSubset;
        size_t i, j;
        char acNumber[16];
        string strFile;

        if (!pBayesNet) {
            g_CatSleipnir().error("CDatabase::Open( %s, %d ) null Bayes net", strOutputDirectory.c_str(), iFiles);
            return false;
        }

        Clear();
        pBayesNet->GetNodes(vecstrNodes);
        for (i = 1; i < vecstrNodes.size(); ++i)
            vecstrNodes[i - 1] = strInputDirectory + '/' + vecstrNodes[i];

        if (vecstrNodes.size())
            vecstrNodes.resize(vecstrNodes.size() - 1);
        m_vecpDBs.resize(iFiles);

        int iNumFilesOpen = 1000;
        for (i = 0; i < m_vecpDBs.size(); ++i) {
            m_vecpDBs[i] = new CDatabaselet(m_useNibble);
        }

        size_t k;
        for (i = 0; i < m_vecpDBs.size(); ++i) {
            if (i % iNumFilesOpen == 0 && i > 0) {
                for (k = 0; k < iNumFilesOpen; k++) {
                    m_vecpDBs[i - k - 1]->CloseFile();
                }
            }
            vecstrSubset.clear();
            for (j = i; j < vecstrGenes.size(); j += m_vecpDBs.size())
                vecstrSubset.push_back(vecstrGenes[j]);
#pragma warning(disable : 4996)
            sprintf(acNumber, "%08lu", i);
#pragma warning(default : 4996)
            strFile = strOutputDirectory + '/' + acNumber + c_acExtension;
            if (!(i % 100))
                g_CatSleipnir().notice("CDatabase::Open( %s, %d ) initializing file %d/%d",
                                       strOutputDirectory.c_str(), iFiles, i, m_vecpDBs.size());
            if (!m_vecpDBs[i]->Open(strFile, vecstrSubset, vecstrGenes.size(), vecstrNodes.size())) {
                g_CatSleipnir().error("CDatabase::Open( %s, %d ) could not open file %s",
                                      strOutputDirectory.c_str(), iFiles, strFile.c_str());
                return false;
            }
        }

        for (i = 0; i < m_vecpDBs.size(); ++i) {
            m_vecpDBs[i]->CloseFile();
        }

        for (i = 0; i < vecstrGenes.size(); ++i)
            m_mapstriGenes[m_vecpDBs[i % m_vecpDBs.size()]->GetGene(i / m_vecpDBs.size())] = i;

        // Rename AddDataFromDabFiles
        return CDatabaseImpl::Open(vecstrGenes, vecstrNodes, mapstriZeros);
    }

    /* Version of Open() that takes a list of datasets as input. Key method */
    // This version is used by Data2DB to create the database
    // Rename to CreateDBFiles
    bool CDatabase::Open(const std::vector <std::string> &vecstrGenes,
                                  const std::vector <std::string> &vecstrDatasets,
                                  const std::string &strInputDirectory,
                                  const std::string &strOutputDirectory,
                                  size_t iFiles,
                                  const map <string, size_t> &mapstriZeros) {

        vector <string> vecstrNodes, vecstrSubset;
        size_t i, j;
        char acNumber[16];
        string strFile;

        Clear();

        // Create a list of the dab files that are the inputs to create the DB
        vecstrNodes.resize(vecstrDatasets.size());
        for (i = 0; i < vecstrDatasets.size(); i++) {
            // dataset names are the prefix of the dab file name
            vecstrNodes[i] = strInputDirectory + '/' + vecstrDatasets[i];
        }

        // Number of DB files to create is iFiles
        m_vecpDBs.resize(iFiles);
        int iNumFilesOpen = 1000;  // TODO - make configurable

        for (i = 0; i < m_vecpDBs.size(); ++i) {
            m_vecpDBs[i] = new CDatabaselet(m_useNibble);
        }

        size_t k;

        // vecDBs is a vector of CDatabaselets each corresponding to one DB file
        // Loop through each DB file that will be built
        for (i = 0; i < m_vecpDBs.size(); ++i) { //block size (such as 1000)
            if (i % iNumFilesOpen == 0 && i > 0) {
                // We're at some increment of the allowed open files (i.e. 1000)
                // Close the previous 1000 files so the next 1000 can be opened.
                for (k = 0; k < iNumFilesOpen; k++) {
                    m_vecpDBs[i - k - 1]->CloseFile();
                }
            }
            // Create the subset of genes in this DB file
            vecstrSubset.clear();
            for (j = i; j < vecstrGenes.size(); j += m_vecpDBs.size()) {
                // The first 1000 genes assigned round robin to the first 1000 files
                // and by every increment of 1000 there after (or numDB files)
                vecstrSubset.push_back(vecstrGenes[j]); //contains index for 1000, 2000, 3000th genes
            }
            sprintf(acNumber, "%08lu", i);
            strFile = strOutputDirectory + '/' + acNumber + c_acExtension;

            if (!(i % 100))
                g_CatSleipnir().notice("CDatabase::Open( %s, %d ) initializing file %d/%d",
                                       strOutputDirectory.c_str(), iFiles, i, m_vecpDBs.size());
            // Open the new DB file and create the header
            // vecstrSubset gives gene_ids, num genes in total, num datasets
            if (!m_vecpDBs[i]->Open(strFile, vecstrSubset, vecstrGenes.size(), vecstrNodes.size())) {
                g_CatSleipnir().error("CDatabase::Open( %s, %d ) could not open file %s",
                                      strOutputDirectory.c_str(), iFiles, strFile.c_str());
                return false;
            }
        }

        for (i = 0; i < m_vecpDBs.size(); ++i) {
            m_vecpDBs[i]->CloseFile();
        }

        for (i = 0; i < vecstrGenes.size(); ++i) {
            // DB file number
            size_t db_idx = i % m_vecpDBs.size();
            // gene index within a DB file
            size_t gene_idx = i / m_vecpDBs.size();
            std::string gene_name = m_vecpDBs[db_idx]->GetGene(gene_idx);
            m_mapstriGenes[gene_name] = i;
        }

        // Load the data from the input data files (Nodes) into the new CDatabaselet DB files
        // Rename call AddDataFromDabFiles
        return CDatabaseImpl::Open(vecstrGenes, vecstrNodes, mapstriZeros);
    }

    /* the key Open() method for Data2DB conversion */
    // Creates the output DB file content from the input Dab files
    // vecstrGenes - list of all genes
    // vecstrFiles - list of all Dab files to build DB from
    // Rename AddDataFromDabFiles
    bool CDatabaseImpl::Open(const std::vector <std::string> &vecstrGenes,
                                         const std::vector <std::string> &vecstrFiles,
                                         const map <string, size_t> &mapstriZeros) {
        size_t i, j, k, iOne, iTwo, iOutBlock, iOutBase, iOutOffset, iInBlock, iInBase, iInOffset;
        vector <size_t> veciGenes;
        float d;
        map<string, size_t>::const_iterator iterZero;

        veciGenes.resize(vecstrGenes.size());
        // m_iBlockOut is number of output DB files to process at once - default to all
        iOutBlock = (m_iBlockOut == -1) ? m_vecpDBs.size() : m_iBlockOut;
        // m_iBlockIn is number of input DAB files to process at once - default to all
        iInBlock = (m_iBlockIn == -1) ? vecstrFiles.size() : m_iBlockIn;

        int iNumFilesOpen = 1000;

        /* blocking parameter, iOutBase: number of databaselets to process at a time */
        for (iOutBase = 0; iOutBase < m_vecpDBs.size(); iOutBase += iOutBlock) {
            vector <string> vecstrMyGenes;  // names of genes in this set (block) of output DB files
            vector <size_t> veciMyGenes;  // Dab indices of genes in this block output DB file

            // Iterate through the current set (block) of output databaselets
            for (iOutOffset = 0; (iOutOffset < iOutBlock) && ((iOutBase + iOutOffset) < m_vecpDBs.size());
                 ++iOutOffset) {
                const CDatabaselet &DB = *m_vecpDBs[iOutBase + iOutOffset];
                
                // Gather all the gene names in the output set of DBs
                for (i = 0; i < DB.GetGenes(); ++i)
                    vecstrMyGenes.push_back(DB.GetGene(i));
            }

            veciMyGenes.resize(vecstrMyGenes.size());

            // Iterate through the input DAB files to get the correlation arrays for
            //  the genes in the output set of DB files.
            // Remember that each DAB file will be a dataset entry in the DB and there
            //  will be a gene correlation array for each gene in the DB file for each
            //  dataset. 
            // The values will be stored in vecData[dataset] where each vecData[dataset]
            //  will have InBlock primary genes and totalGene correlations per primary gene.
            for (iInBase = 0; iInBase < vecstrFiles.size(); iInBase += iInBlock) {
                vector <CCompactFullMatrix> vecData;
                vecData.resize(((iInBase + iInBlock) > vecstrFiles.size()) ?
                               (vecstrFiles.size() - iInBase) : iInBlock);
                for (iInOffset = 0; iInOffset < vecData.size(); ++iInOffset) {
                    CDataPair Dat;

                    if (!Dat.Open((vecstrFiles[iInBase + iInOffset] + c_acDAB).c_str(),
                                  false, m_fMemmap)) {
                        if (!Dat.Open((vecstrFiles[iInBase + iInOffset] + c_acQDAB).c_str(), false, m_fMemmap)) {
                            g_CatSleipnir().error("CDatabaseImpl::Open( ) could not open %s",
                                                  (vecstrFiles[iInBase + iInOffset] + c_acDAB).c_str());
                            g_CatSleipnir().error("CDatabaseImpl::Open( ) could not open %s",
                                                  (vecstrFiles[iInBase + iInOffset] + c_acQDAB).c_str());
                            return false;
                        }
                    }

                    // Loop through genes in output DB file and get corresponding gene indexes from the Dab file
                    for (i = 0; i < veciMyGenes.size(); ++i)
                        veciMyGenes[i] = Dat.GetGene(vecstrMyGenes[i]);
                    // Loop through all genes in the collection and get corresponding indexes in the Dab file
                    for (i = 0; i < veciGenes.size(); ++i)
                        veciGenes[i] = Dat.GetGene(vecstrGenes[i]);

                    if (m_useNibble) {
                        vecData[iInOffset].Initialize(veciMyGenes.size(), veciGenes.size(), 16, true);
                    } else {
                        vecData[iInOffset].Initialize(veciMyGenes.size(), veciGenes.size(), (unsigned char) 256, true);
                    }

                    string strName = CMeta::Filename(
                            CMeta::Deextension(CMeta::Basename(vecstrFiles[iInBase + iInOffset].c_str())));
                    size_t iZero = ((iterZero = mapstriZeros.find(strName)) == mapstriZeros.end()) ? -1
                                                                                                   : iterZero->second;
                    //#pragma omp parallel for \
				shared(Dat, veciGenes, veciMyGenes, vecData) \
				private(j, i) \
				schedule(static)

                    // Loop through the genes for this output DB file (stored as Dab indices)
                    for (i = 0; i < veciMyGenes.size(); i++) {
                        size_t iOne = veciMyGenes[i];
                        // Loop through all genes in the collection (as Dab indices)
                        for (j = 0; j < veciGenes.size(); j++) {
                            size_t iTwo = veciGenes[j];
                            size_t iVal = -1;
                            // Get the correlation value between the two genes
                            if (iOne != -1 && iTwo != -1) {
                                if (!CMeta::IsNaN(d = Dat.Get(iOne, iTwo)))
                                    iVal = Dat.Quantize(d);
                                if (iVal == -1)
                                    iVal = iZero;
                            }
                            // Write the correlation value to vecData for this inpub Dab
                            //  in preparation for writing the DB file
                            if (iVal != -1)
                                vecData[iInOffset].Set(i, j, iVal + 1);
                        }
                    }


                }

                // Output the data collected from the current set of input Dabs to all the output DB files
                for (i = iOutOffset = 0; (iOutOffset < iOutBlock) && ((iOutBase + iOutOffset) <
                                                                      m_vecpDBs.size()); ++iOutOffset) {
                    CDatabaselet &DB = *m_vecpDBs[iOutBase + iOutOffset];

                    /* close files if too many file handles opened */
                    /*if(iOutOffset>0 && (iOutOffset % iNumFilesOpen==0 ||
                        (iOutOffset + 1) >= iOutBlock || (iOutBase + iOutOffset + 1) >= m_vecpDBs.size() ) ){
                        for(k=0; k<iNumFilesOpen; k++){
                            if(iOutOffset + iOutBase - k - 1 < 0){
                                break;
                            }
                            m_vecpDBs[iOutOffset + iOutBase - 1 - k]->CloseFile();
                        }
                    }*/

                    DB.OpenNoOverwrite();

                    if (!(iOutOffset % 100))
                        cerr << "Processing offset " << iOutOffset << '/' << iOutBlock << endl;
                    //iInBase, for B=2, iInBase = 0, 2, 4 (total 5 datasets)
                    //i = 0 18 36 ...
                    // Write the vecData to the current output DB:
                    //  - vecData is a set of matrices, one per input dataset with gene correlations
                    //     for the primary genes of the output DB files
                    //  - i is the index of the first row of primary genes in the vecData matrices
                    //     corresponding to the current output DB. Since all output DBs have the
                    //     same number of primary genes in them we increment i by that geneCount
                    //     with each for loop iteration.
                    //  - iInBase is the base index number of the input datasets, 
                    //     the len of vecData will give the other index values
                    //  - m_fBuffer - whether to buffer channges in memory before writing to DB file 
                    
                    if (!DB.Open(vecData, i, iInBase, m_fBuffer)) {
                        return false;
                    }

                    i += DB.GetGenes();

                    DB.CloseFile();
                }
            }
        }

        return true;
    }

    bool CDatabase::Open(const std::string &strInputDirectory) {
        size_t i, j;
        vector <string> vecstrFiles;
        string strFile;

        FOR_EACH_DIRECTORY_FILE(strInputDirectory, strFile)
            if (!CMeta::IsExtension(strFile, c_acExtension))
                continue;

            i = atoi(CMeta::Deextension(CMeta::Basename(strFile.c_str())).c_str());
            if (vecstrFiles.size() <= i)
                vecstrFiles.resize(i + 1);
            else if (vecstrFiles[i].length() != 0) {
                g_CatSleipnir().error("CDatabase::Open( %s ) duplicate file: %s (%d)", strInputDirectory.c_str(),
                                      strFile.c_str(), i);
                return false;
            }
            vecstrFiles[i] = strInputDirectory + '/' + strFile; }

        Clear();
        m_vecpDBs.resize(vecstrFiles.size());
        for (i = 0; i < m_vecpDBs.size(); ++i) {
            m_vecpDBs[i] = new CDatabaselet(m_useNibble);
            if (!m_vecpDBs[i]->Open(vecstrFiles[i]))
                return false;
            for (j = 0; j < m_vecpDBs[i]->GetGenes(); ++j)
                m_mapstriGenes[m_vecpDBs[i]->GetGene(j)] = (j * m_vecpDBs.size()) + i;
        }

        return true;
    }

    bool CDatabaselet::Set(uint32_t &iGenes, uint32_t &iDatasets, vector <string> &vecstrSubsetGenes) {
        size_t i;
        m_vecstrGenes.clear();
        m_vecstrGenes.resize(vecstrSubsetGenes.size());
        for (i = 0; i < vecstrSubsetGenes.size(); i++) {
            m_vecstrGenes[i] = vecstrSubsetGenes[i];
        }

        m_iGenes = iGenes;
        m_iDatasets = iDatasets;
        uint32_t iSize = m_vecstrGenes.size();
        m_iHeader = sizeof(m_iHeader) + sizeof(m_iGenes) + sizeof(m_iDatasets) + sizeof(iSize);
        for (i = 0; i < m_vecstrGenes.size(); ++i) {
            m_iHeader += m_vecstrGenes[i].size() + 1;
        }

        return true;
    }

//Create a copy of current CDatabase collection that has X number of CDatabaselets
    bool CDatabase::Reorganize(const char *dest_db_dir, const size_t &num_db) {
        int dest_db = num_db;
        size_t i, j, l;
        const char c_acExtension[] = ".db";
        char acNumber[16];

        vector <string> vecstrG;
        vecstrG.resize(m_mapstriGenes.size());

        for (map<string, size_t>::iterator iter = m_mapstriGenes.begin();
             iter != m_mapstriGenes.end(); iter++) {
            string first = iter->first;
            size_t second = iter->second;
            vecstrG[second] = first;
        }

        for (i = 0; i < dest_db; i++) {
            vector <string> vecstrSubset;
            for (j = i; j < vecstrG.size(); j += dest_db)
                vecstrSubset.push_back(vecstrG[j]);
            //size of this db
            size_t iSize = GetGenes() * GetDatasets() * vecstrSubset.size();

            unsigned char *Qt = (unsigned char *) malloc(iSize);
            int tot = 0;
            for (j = 0; j < vecstrSubset.size(); j++) {
                size_t k = m_mapstriGenes.find(vecstrSubset[j])->second;
                vector<unsigned char> Qi;
                if (!GetGene(k, Qi)) {
                    cerr << "Gene error" << endl;
                    continue;
                }
                for (l = 0; l < Qi.size(); l++)
                    Qt[tot + l] = Qi[l];
                tot += Qi.size();
            }
            sprintf(acNumber, "%08lu", i);
            string dest_dir = dest_db_dir;
            string strFile = dest_dir + "/" + acNumber + c_acExtension;

            CDatabaselet DBS(false);
            DBS.Open(strFile.c_str(), vecstrSubset, vecstrG.size(),
                     GetDatasets());
            DBS.Write((char *) Qt, iSize, 0);
            free(Qt);
        }
        return true;
    }

//For SeekMiner
    bool CDatabase::Open(const string &strDBDirectory,
                         const vector <string> &vecstrGenes, const size_t &iDatasets, const size_t &iNumDBs) {
        return CDatabase::Open(strDBDirectory.c_str(), vecstrGenes, iDatasets, iNumDBs);
    }

//For SeekMiner
    bool CDatabase::Open(const char *db_dir, const vector <string> &vecstrGenes,
                         const size_t &iDatasets, const size_t &iNumDBs) {
        size_t i, j, k;
        Clear();
        m_vecpDBs.resize(iNumDBs);
        char acNumber[16];

        for (i = 0; i < m_vecpDBs.size(); ++i)
            m_vecpDBs[i] = new CDatabaselet(m_useNibble);

        string strDBDirectory = db_dir;

        for (i = 0; i < iNumDBs; i++) {
            vector <string> vecstrSubset;
            vecstrSubset.clear();
            for (j = i; j < vecstrGenes.size(); j += m_vecpDBs.size())
                vecstrSubset.push_back(vecstrGenes[j]); //contains index for 1000, 2000, 3000th genes
            sprintf(acNumber, "%08lu", i);
            string strFile = strDBDirectory + '/' + acNumber + c_acExtension;
            if (!std::filesystem::exists(strFile)) {
                return false;
            }
            uint32_t iGenes = vecstrGenes.size();
            uint32_t iDset = iDatasets;
            m_vecpDBs[i]->Set(iGenes, iDset, vecstrSubset);
            m_vecpDBs[i]->SetFile(strFile);
        }

        for (i = 0; i < vecstrGenes.size(); ++i) {
            string s = m_vecpDBs[i % m_vecpDBs.size()]->GetGene(i / m_vecpDBs.size());
            string file = m_vecpDBs[i % m_vecpDBs.size()]->GetFile();
            //fprintf(stderr, "%s\t%s\t%d\n", s.c_str(), file.c_str(), i);
            m_mapstriGenes[m_vecpDBs[i % m_vecpDBs.size()]->GetGene(i / m_vecpDBs.size())] = i;
        }
        //fprintf(stderr, "Done\n"); getchar();

        return true;
    }

}

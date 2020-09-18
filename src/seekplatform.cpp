/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a part of SEEK (Search-based exploration of expression compendium)
* which is authored and maintained by: Qian Zhu (qzhu@princeton.edu)
*
* If you use this file, please cite the following publication:
* Qian Zhu, Aaron K Wong, Arjun Krishnan, Miriam R Aure, Alicja Tadych, 
* Ran Zhang, David C Corney, Casey S Greene, Lars A Bongo, 
* Vessela N Kristensen, Moses Charikar, Kai Li & Olga G Troyanskaya
* "Targeted exploration and analysis of large cross-platform human 
* transcriptomic compendia" Nat Methods (2015)
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library for development, or use any other Sleipnir executable
* tools, please also cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "seekplatform.h"
#include "seekreader.h"
#include <filesystem>

namespace Sleipnir {

    CSeekPlatform::CSeekPlatform() {
        m_iNumGenes = 0;
        m_vecfPlatformAvg.clear();
        m_vecfPlatformStdev.clear();
        m_vecPlatformCount.clear();
        m_strPlatformName = "";
    }

    CSeekPlatform::~CSeekPlatform() {
        m_vecfPlatformAvg.clear();
        m_vecfPlatformStdev.clear();
        m_vecPlatformCount.clear();
    }

    void CSeekPlatform::Copy(const CSeekPlatform &pl) {
        m_iNumGenes = pl.m_iNumGenes;
        m_strPlatformName = pl.m_strPlatformName;
        m_vecfPlatformAvg.resize(pl.m_vecfPlatformAvg.size());
        m_vecfPlatformStdev.resize(pl.m_vecfPlatformStdev.size());
        m_vecPlatformCount.resize(pl.m_vecPlatformCount.size());
        copy(pl.m_vecfPlatformAvg.begin(), pl.m_vecfPlatformAvg.end(),
             m_vecfPlatformAvg.begin());
        copy(pl.m_vecfPlatformStdev.begin(), pl.m_vecfPlatformStdev.end(),
             m_vecfPlatformStdev.begin());
        copy(pl.m_vecPlatformCount.begin(), pl.m_vecPlatformCount.end(),
             m_vecPlatformCount.begin());
    }

    void CSeekPlatform::InitializePlatform(const utype &numGenes,
                                           const string &strPlatformName) {
        m_iNumGenes = numGenes;
        CSeekTools::InitVector(m_vecfPlatformAvg, numGenes, (float) 0);
        CSeekTools::InitVector(m_vecfPlatformStdev, numGenes, (float) 0);
        CSeekTools::InitVector(m_vecPlatformCount, numGenes, (uint32_t) 0);
        m_strPlatformName = strPlatformName;
    }

    void CSeekPlatform::SetPlatformAvg(const utype &i, const float &val) {
        m_vecfPlatformAvg[i] = val;
    }

    void CSeekPlatform::SetPlatformStdev(const utype &i, const float &val) {
        m_vecfPlatformStdev[i] = val;
    }

    void CSeekPlatform::SetPlatformCount(const utype &i, const uint32_t &val)
    {
        m_vecPlatformCount[i] = val;
    }

    float CSeekPlatform::GetPlatformAvg(const utype &i) const {
        return m_vecfPlatformAvg[i];
    }

    float CSeekPlatform::GetPlatformStdev(const utype &i) const {
        return m_vecfPlatformStdev[i];
    }

    uint32_t CSeekPlatform::GetPlatformCount(const utype &i) const
    {
        return m_vecPlatformCount[i];
    }

    void CSeekPlatform::ResetPlatform() {
        m_iNumGenes = 0;
        m_vecfPlatformAvg.clear();
        m_vecfPlatformStdev.clear();
        m_vecPlatformCount.clear();
        m_strPlatformName = "";
    }

    // ---------- SeekPlatforms Class Implementation ------------
    // SeekPlatforms (multi-platform class) implementation

    SeekPlatforms::SeekPlatforms(size_t numPlatforms, size_t numGenes)
    {
        initialize(numPlatforms, numGenes);
    }

    void SeekPlatforms::initialize(size_t numPlatforms, size_t numGenes)
    {
        m_numPlatforms = numPlatforms;
        m_numGenes = numGenes;
        m_includeCounts = true;

        platformAvgMatrix.Initialize(numPlatforms, numGenes);
        platformStdevMatrix.Initialize(numPlatforms, numGenes);
        platformCountMatrix.Initialize(numPlatforms, numGenes);

        for (int i = 0; i < numPlatforms; i++)
        {
            for (int j = 0; j < numGenes; j++)
            {
                platformAvgMatrix.Set(i, j, CMeta::GetNaN());
                platformStdevMatrix.Set(i, j, CMeta::GetNaN());
                platformCountMatrix.Set(i, j, (uint32_t) 0);
            }
        }

        m_cseek_platforms.clear();
        m_platformNames.clear();
        m_mapPlatformNameToOrderID.clear();

        m_cseek_platforms.resize(numPlatforms);
        m_platformNames.resize(numPlatforms);
    }

    void SeekPlatforms::copy(SeekPlatforms &srcPlatforms)
    {
        // copy numGenes, numPlatforms, bincludeCounts
        m_numPlatforms = srcPlatforms.getNumPlatforms();
        m_numGenes = srcPlatforms.getNumGenes();
        m_includeCounts = srcPlatforms.bIncludesCounts();

        // copy name vector
        vector<string> &srcNames = srcPlatforms.getPlatformNames();
        m_platformNames.resize(srcNames.size());
        m_platformNames = srcNames;
        // alternate method would be: copy(srcNames.begin(), srcNames.end(), m_platformNames.begin());

        // copy name->id map
        map<string, utype> srcMap = srcPlatforms.getPlatformMap();
        m_mapPlatformNameToOrderID.insert(srcMap.begin(), srcMap.end());

        // copy m_cseek_platforms
        vector<CSeekPlatform> &srcCSeekPlatforms = srcPlatforms.getCSeekPlatforms();
        m_cseek_platforms.clear();
        m_cseek_platforms.resize(m_numPlatforms);
        for (int i = 0; i < m_numPlatforms; i++) {
            m_cseek_platforms[i].Copy(srcCSeekPlatforms[i]);
        }

        // copy platform matrices
        platformAvgMatrix.Copy(srcPlatforms.platformAvgMatrix);
        platformStdevMatrix.Copy(srcPlatforms.platformStdevMatrix);
        platformCountMatrix.Copy(srcPlatforms.platformCountMatrix);
    }

    void SeekPlatforms::setCSeekPlatformData()
    {
        uint32_t numPlatforms = platformAvgMatrix.GetRows();
        uint32_t numGenes = platformAvgMatrix.GetColumns();

        assert(numPlatforms > 0);
        assert(numGenes > 0);
        assert(numGenes == m_numGenes);
        assert(numPlatforms == m_numPlatforms);
        assert(m_platformNames.size() == numPlatforms);
        assert(m_mapPlatformNameToOrderID.size() == numPlatforms);

        m_cseek_platforms.resize(numPlatforms);

        for (int platIdx = 0; platIdx < numPlatforms; platIdx++)
        {
            m_cseek_platforms[platIdx].InitializePlatform(numGenes, m_platformNames[platIdx]);
            for (int geneIdx = 0; geneIdx < numGenes; geneIdx++)
            {
                m_cseek_platforms[platIdx].SetPlatformAvg(geneIdx, platformAvgMatrix.Get(platIdx, geneIdx));
                m_cseek_platforms[platIdx].SetPlatformStdev(geneIdx, platformStdevMatrix.Get(platIdx, geneIdx));
                if (m_includeCounts)
                {
                    m_cseek_platforms[platIdx].SetPlatformCount(geneIdx, platformCountMatrix.Get(platIdx, geneIdx));
                }
            }
        }
    }

    void SeekPlatforms::loadPlatformDataFromFiles(string platformDir)
    {
        const int lineSize = 1024;
        // char *szPlatformDir = platformDir.c_str();

        string strPlatformDirectory = platformDir;
        string strAvgFile = platformDir + "/" +
                            "all_platforms.gplatavg";
        string strStdevFile = platformDir + "/" +
                              "all_platforms.gplatstdev";
        string strCountFile = platformDir + "/" +
                              "all_platforms.gplatcount";
        string strPlatformOrderFile = platformDir + "/" +
                                      "all_platforms.gplatorder";

        // Check existence of the gplatcount file, if it doesn't exist then skip loading counts
        if (!filesystem::exists(strCountFile)) {
            m_includeCounts = false;
        }

        platformAvgMatrix.Open(strAvgFile.c_str());
        platformStdevMatrix.Open(strStdevFile.c_str());
        if (m_includeCounts)
        {
            platformCountMatrix.Open(strCountFile.c_str());
        }

        m_cseek_platforms.clear();
        m_cseek_platforms.resize(platformAvgMatrix.GetRows());
        utype i, j;

        m_platformNames.clear();
        m_mapPlatformNameToOrderID.clear();
        ifstream ifsm;
        ifsm.open(strPlatformOrderFile.c_str());
        char acBuffer[lineSize];
        utype c_iBuffer = lineSize;
        i = 0;
        while (!ifsm.eof())
        {
            ifsm.getline(acBuffer, c_iBuffer - 1);
            if (acBuffer[0] == 0)
                break;
            acBuffer[c_iBuffer - 1] = 0;
            m_platformNames.push_back(acBuffer);
            m_mapPlatformNameToOrderID[acBuffer] = i;
            i++;
        }
        m_platformNames.resize(m_platformNames.size());
        ifsm.close();

        m_numPlatforms = platformAvgMatrix.GetRows();
        m_numGenes = platformAvgMatrix.GetColumns();
        assert(m_numPlatforms == m_mapPlatformNameToOrderID.size());

        setCSeekPlatformData();
    }


    void SeekPlatforms::savePlatformDataToFiles(string platformDir)
    {
        const char *szPlatformDir = platformDir.c_str();
        char outFile[1024];

        sprintf(outFile, "%s/all_platforms.gplatavg", szPlatformDir);
        platformAvgMatrix.Save(outFile);
        sprintf(outFile, "%s/all_platforms.gplatstdev", szPlatformDir);
        platformStdevMatrix.Save(outFile);
        sprintf(outFile, "%s/all_platforms.gplatcount", szPlatformDir);
        platformCountMatrix.Save(outFile);
        sprintf(outFile, "%s/all_platforms.gplatorder", szPlatformDir);
        ofstream outfile;
        outfile.open(outFile);
        for (int i = 0; i < m_platformNames.size(); i++)
        {
            outfile << m_platformNames[i] << "\n";
        }
        outfile.close();
    }


    bool SeekPlatforms::combineWithPlatform(SeekPlatforms &seekPlatforms2)
    {
        /* Combine the platform per-gene mean and stdev from two databases
         * Notes on combining standard deviations of two groups:
         * n1 = num samples in group1
         * n2 = num samples of group2
         * A1 = mean of group1
         * A2 = mean of group2
         * V1 = variance of group1
         * V2 = variance of group2
         * Find the combined mean as:
         * Ac = (n1 * A1 + n2 * A2) / (n1 + n2)
         * Let:
         * d1 = A1 - Ac
         * d2 = A2 - Ac
         * Find the combined variance as
         * Vc = n1 * (V1 + d1^2) + n2 * (V2 + d2^2) / (n1 + n2)
         * Cite 19 Recommendations 19th Aug, 2015 Sangita C.Patil(Birajdar) MAEER’s Arts, Commerce and Science College
         * Cochrane Handbook for Systematic Reviews of Interventions, 7.7.3.8, is the following:
         * Being N the number in the group, M the mean and SD the Standard Deviation, of groups 1 and 2...
         * For Mean: (N1M1 + N2M2) / (N1 + N2)
         * For SD, it is the square root of: {(N1 - 1) SD1^2 + (N2 – 1) SD2^2 + [(N1N2 / N1 + N2) (M1^2 + M2^2 – 2M1M2)]} / (N1 + N2 – 1)
        */

        if (m_numPlatforms == 0) {
            // original seekPlatforms is empty, initialize it with the new num of genes
            initialize(0, seekPlatforms2.getNumGenes());
        } else if (seekPlatforms2.getNumPlatforms() == 0) {
            // the merging seekplatform is empty, intialize with num genes
            seekPlatforms2.initialize(0, m_numGenes);
        } else {
            // If both the seekPlatforms are non-empty,
            // then the two platforms must have same number and set of genes
            assert(m_numGenes == seekPlatforms2.getNumGenes());
        }

        vector<string> &namesPlatforms2 = seekPlatforms2.getPlatformNames();
        map<string, utype> mapPlatforms2 = seekPlatforms2.getPlatformMap();
        uint32_t numPlatforms2 = seekPlatforms2.getNumPlatforms();
        assert(namesPlatforms2.size() == numPlatforms2);
        assert(mapPlatforms2.size() == numPlatforms2);

        // Divide platforms2 into two sets
        //  a. platforms in 2 that aren't in 1 -- append them to 1
        //  b. platforms in both 2 and 1 -- combine them
        vector<string> newPlatformNames;
        vector<string> mutualPlatformNames;
        for (int i=0; i < numPlatforms2; i++) {
            string &platName = namesPlatforms2[i];
            if (m_mapPlatformNameToOrderID.find(platName) == m_mapPlatformNameToOrderID.end()) {
                // platform name not found - this is a new platform
                newPlatformNames.push_back(platName);
            } else {
                mutualPlatformNames.push_back(platName);
            }
        }

        // *** Append New Platforms ***
        // loop through new platform names and add them to existing name vector and map
        // append data rows from the new platforms to the existing data matrices
        uint32_t numNewPlatforms = newPlatformNames.size();
        if (numNewPlatforms > 0) {
            bool doIncludeCounts = m_includeCounts && seekPlatforms2.bIncludesCounts();
            platformAvgMatrix.AddRows(numNewPlatforms, true);
            platformStdevMatrix.AddRows(numNewPlatforms, true);
            platformCountMatrix.AddRows(numNewPlatforms, true);
            m_numPlatforms += numNewPlatforms;
            for (int i = 0; i < numNewPlatforms; i++) {
                string &newName = newPlatformNames[i];
                uint32_t prevOrderId = mapPlatforms2[newName];
                uint32_t newOrderId = m_mapPlatformNameToOrderID.size();
                m_platformNames.push_back(newName);
                m_mapPlatformNameToOrderID[newName] = newOrderId;
                platformAvgMatrix.Set(newOrderId, seekPlatforms2.platformAvgMatrix.Get(prevOrderId));
                platformStdevMatrix.Set(newOrderId, seekPlatforms2.platformStdevMatrix.Get(prevOrderId));
                if (doIncludeCounts) {
                  platformCountMatrix.Set(newOrderId, seekPlatforms2.platformCountMatrix.Get(prevOrderId));
                }
            }
        }

        // *** Combine Mutual Platforms ***
        // loop through mutual platform names and combine the statistics
        uint32_t numMutualPlatforms = mutualPlatformNames.size();
        // We need sample counts for both platforms in order to combine them
        // assert sample counts exist for both
        assert(numMutualPlatforms == 0 || bIncludesCounts() == true);
        assert(numMutualPlatforms == 0 || seekPlatforms2.bIncludesCounts() == true);
        for (int i=0; i< numMutualPlatforms; i++) {
            string &platName = mutualPlatformNames[i];
            // Get the orderId for the current platform name
            uint32_t platId1 = m_mapPlatformNameToOrderID[platName];
            uint32_t platId2 = mapPlatforms2[platName]; 
            // For a given platform, get the platform data for the original and new
            // combine the averages
            float *avgVals1 = platformAvgMatrix.Get(platId1);
            float *avgVals2 = seekPlatforms2.platformAvgMatrix.Get(platId2);
            float *stdevVals1 = platformStdevMatrix.Get(platId1);
            float *stdevVals2 = seekPlatforms2.platformStdevMatrix.Get(platId2);
            uint32_t *countVals1 = platformCountMatrix.Get(platId1);
            uint32_t *countVals2 = seekPlatforms2.platformCountMatrix.Get(platId2);
            for (int gidx=0; gidx<m_numGenes; gidx++) {
                uint32_t count1 = countVals1[gidx];
                uint32_t count2 = countVals2[gidx];
                double avg1 = avgVals1[gidx];
                double avg2 = avgVals2[gidx];
                double var1 = pow(stdevVals1[gidx], 2);
                double var2 = pow(stdevVals2[gidx], 2);

                // Calculate combined average
                // Ac = (n1 * A1 + n2 * A2) / (n1 + n2)
                double newAvg = (avg1 * count1 + avg2 * count2) / (count1 + count2);
                platformAvgMatrix.Set(platId1, gidx, newAvg);
                // printf("Mutual Avg: A1 %0.2f, A2 %0.2f, C1 %d, C2 %d, Ac %0.2f\n", avg1, avg2, count1, count2, newAvg);

                // Calculate combined stdev
                // d1 = A1 - Ac
                double diffAvg1 = avg1 - newAvg;
                // d2 = A2 - Ac
                double diffAvg2 = avg2 - newAvg;
                // Calc the combined variance as
                // Vc = (n1 * (V1 + d1^2) + n2 * (V2 + d2^2)) / (n1 + n2)
                double newVar = (count1 * (var1 + pow(diffAvg1, 2)) + count2 * (var2 + pow(diffAvg2, 2))) /
                               (count1 + count2);
                double newStdev = sqrt(newVar);
                platformStdevMatrix.Set(platId1, gidx, newStdev);

                // Set the combined count
                platformCountMatrix.Set(platId1, gidx, count1 + count2);
            }
        }

        // Copy the matrix data into the CSeekPlatform objects
        setCSeekPlatformData();

        return true;
    }
}


/* *******  Utility Functions ******** */

// Given two directories containing platform statistics files, load and combine the
//  files stats and save into the output directory
void combinePlatformStatistics(string platDir1, string platDir2, string outDir) 
{
    Sleipnir::SeekPlatforms platform1;
    Sleipnir::SeekPlatforms platform2;

    // Load the two platforms
    platform1.loadPlatformDataFromFiles(platDir1);
    platform2.loadPlatformDataFromFiles(platDir2);

    // Combine platform2 into platform1
    platform1.combineWithPlatform(platform2);

    // Save the combined data in platform1 out to files
    platform1.savePlatformDataToFiles(outDir);
}

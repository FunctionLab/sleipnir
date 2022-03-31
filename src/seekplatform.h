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
#ifndef SEEKPLATFORM_H
#define SEEKPLATFORM_H

#include "seekbasic.h"
#include "fullmatrix.h"

namespace Sleipnir {

/*!
 * \brief
 * Representation of a microarray platform that is used by Seek
 *
 * Contains the gene \a correlation average and standard deviation for a given platform
 *
 * For each gene \a g in each dataset \a d, we calculate the average \a correlation of \a g
 * to all the genes. Then average it across
 * all the datasets in the platform.
 * The result is the platform's gene-\a correlation average, or \c PlatformAvg for short.
 *
 * The standard deviation or \c PlatformStdev measures the spread of \a correlation
 * across all datasets of the given platform.
 *
 * The purpose of the platform's average and standard deviation are to reduce potential biases
 * that might be caused by platform specific \a correlation distributions.
 * \remarks
 * The word \a correlation refers to z-score transformed, standardized Pearson correlation.
 */
    class CSeekPlatform {
    public:
        /*!
         * \brief Constructor
         */
        CSeekPlatform();

        /*!
         * \brief Destructor
         */
        ~CSeekPlatform();

        /*!
         * \brief Initialize the platform
         *
         * \param numGenes
         * The number of genes covered by the platform
         *
         * \param strPlatformName
         * Assign a name to the platform
         */
        void InitializePlatform(const utype &, const string &);

        /*!
         * \brief Set the platform \a correlation average for a particular gene
         * \param i Gene index
         * \param val The average \a correlation for the gene
         */
        void SetPlatformAvg(const utype &, const float &);

        /*!
         * \brief Set the platform standard deviation of \a correlation for a given gene
         * \param i Gene index
         * \param val The standard deviation
         */
        void SetPlatformStdev(const utype &, const float &);

        /*!
         * \brief Set the platform count of \a samples for a given gene
         * \param i Gene index
         * \param val The count
         */
        void SetPlatformCount(const utype &, const uint32_t &);

        /*!
         * \brief Get the platform-wide \a correlation average for a given gene
         * \param i Gene index
         * \return The platform-wide average
         */
        float GetPlatformAvg(const utype &) const;

        /*!
         * \brief Get the platform-wide standard deviation of \a correlation for a given gene
         * \param i Gene index
         * \return The platform-wide standard deviation
         */
        float GetPlatformStdev(const utype &) const;

        /*!
         * \brief Get the platform-wide count of \a samples for a given gene
         * \param i Gene index
         * \return The platform-wide count
         */
        uint32_t GetPlatformCount(const utype &) const;

        /*!
         * \brief Reset
         */
        void ResetPlatform();

        /*!
         * Create a copy from a given platform
         * \param pl The given platform
         */
        void Copy(const CSeekPlatform &);

    private:
        vector<float> m_vecfPlatformAvg;
        vector<float> m_vecfPlatformStdev;
        vector<uint32_t> m_vecPlatformCount;
        string m_strPlatformName;
        utype m_iNumGenes;
    };


    /*!
    * \brief
    * SeekPlatforms stores the collection of SeekPlatform objects (one per plaftorm).
    * This makes it easy to load and store all the platform data in one enclosing class.
    * 
    * Notes on why this change:
    * We want to separate the seekplatform functionality into one class/file for unit testing.
    * The code to compute, load and save seekplatform files uses fullMatrices, while the search 
    * algorithm uses vector<CSeekPlatforms>. In this SeekPlatorms class we will combine both 
    * of these together. The fullMatrices are public so the various parts of code can easily
    * set them. The function setCSeekPlatformData() copies the fullMatrices into the 
    * vector<CSeekPlatforms> format which can then be retrieved and used by the search algorithm.
    * With this change we will now pass around a SeekPlatorms object rather than pass around
    * differet individual components as was done before such as the platformNames, platformMap, 
    * and vector<CSeekPlatforms>. We also add a combineWithPlatform() function to merge platform
    * stats together. This is helpful not only when merging in new datasets, but also when 
    * using multiple databases for a search.
    * 
    */
    class SeekPlatforms {
    public:
        SeekPlatforms() { initialize(0, 0); }
        SeekPlatforms(size_t numPlatforms, size_t numGenes);
        void initialize(size_t numPlatforms, size_t numGenes);
        void clear() { initialize(0, 0); }
        void copy(const SeekPlatforms &srcPlatform);
        void setCSeekPlatformData();
        void loadPlatformDataFromFiles(string platformDir);
        void savePlatformDataToFiles(string platformDir);
        void setPlatformNames(vector<string> &platformNames) {m_platformNames = platformNames;}
        void setPlatformNameMap(map<string, utype> &mapstriPlatform) { m_mapPlatformNameToOrderID = mapstriPlatform; }
        const vector<CSeekPlatform> &getCSeekPlatforms() const { return m_cseek_platforms; }
        const map<string, utype> &getPlatformMap() const { return m_mapPlatformNameToOrderID; }
        const vector<string> &getPlatformNames() const { return m_platformNames; }
        uint32_t getNumPlatforms() const { return m_numPlatforms; }
        uint32_t getNumGenes() const { return m_numGenes; }
        bool bIncludesCounts() const { return m_includeCounts; }
        bool combineWithPlatform(SeekPlatforms &seekPlatforms2);

        // Public Member Variables
        CFullMatrix<float> platformAvgMatrix;
        CFullMatrix<float> platformStdevMatrix;
        CFullMatrix<uint32_t> platformCountMatrix;
    private: 
        // Private Member Variables
        vector<CSeekPlatform> m_cseek_platforms;  // array of individual platform objects
        vector<string> m_platformNames;
        map<string, utype> m_mapPlatformNameToOrderID;
        uint32_t m_numGenes;
        uint32_t m_numPlatforms;
        bool m_includeCounts = true; // TODO set based on file existence
    };

}

/* *******  Utility Functions ******** */

void combinePlatformStatistics(string platDir1, string platDir2, string outDir);


#endif

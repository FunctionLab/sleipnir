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
#ifndef SEEKPLATFORM_H
#define SEEKPLATFORM_H

#include "seekbasic.h"
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
class CSeekPlatform{
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
	string m_strPlatformName;
	utype m_iNumGenes;
};

}
#endif

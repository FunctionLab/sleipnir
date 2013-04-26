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
 * A microarray platform that is used by Seek
 *
 * Contains the gene correlation average and standard deviation for a given platform
 *
 * For each gene \a g in each dataset \a d, we calculate the average correlation of \a g
 * to all the genes. Then calculate the average of \a g's average correlation across
 * all the datasets in the platform.
 * The result is the platform's gene-correlation average, or \c PlatformAvg for short.
 *
 * The standard deviation or \c PlatformStdev measures the spread of the gene-correlation
 * across all datasets of the given platform.
 *
 * The purpose of the platform's average and standard deviation are to reduce potential biases
 * that might be caused by platform specific correlation distributions.
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
	void InitializePlatform(const ushort &, const string &);

	/*!
	 * \brief Set the platform gene-correlation average
	 */
	void SetPlatformAvg(const ushort &, const float &);

	/*!
	 * \brief Set the platform gene-correlation standard deviation
	 */
	void SetPlatformStdev(const ushort &, const float &);

	/*!
	 * \brief Get the platform gene-correlation average
	 */
	float GetPlatformAvg(const ushort &) const;

	/*!
	 * \brief Get the platform gene-correlation standard deviation
	 */
	float GetPlatformStdev(const ushort &) const;

	/*!
	 * \brief Reset
	 */
	void ResetPlatform();

	/*!
	 * Create a copy from a given platform
	 */
	void Copy(const CSeekPlatform &);

private:
	vector<float> m_vecfPlatformAvg;
	vector<float> m_vecfPlatformStdev;
	string m_strPlatformName;
	ushort m_iNumGenes;
};

}
#endif

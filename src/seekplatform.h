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

#include "stdafx.h"

namespace Sleipnir {

class CSeekPlatform{
public:
	CSeekPlatform();
	~CSeekPlatform();

	void InitializePlatform(const size_t &, string &);
	void SetPlatformAvg(const size_t &, float);
	void SetPlatformStdev(const size_t &, float);
	float GetPlatformAvg(const size_t &);
	float GetPlatformStdev(const size_t &);
	void ResetPlatform();

private:
	vector<float> m_vecfPlatformAvg;
	vector<float> m_vecfPlatformStdev;
	size_t m_iPlatformID;
	string m_strPlatformName;
	size_t m_iNumGenes;
};

}
#endif

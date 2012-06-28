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
#include "seekplatform.h"
#include "seekreader.h"

namespace Sleipnir {

CSeekPlatform::CSeekPlatform(){
	m_iNumGenes = 0;
	m_vecfPlatformAvg.clear();
	m_vecfPlatformStdev.clear();
	m_strPlatformName = "";
}

CSeekPlatform::~CSeekPlatform(){
}

void CSeekPlatform::InitializePlatform(const size_t &numGenes, string &strPlatformName){
	m_iNumGenes = numGenes;
	CSeekTools::InitVector(m_vecfPlatformAvg, numGenes, (float) 0);
	CSeekTools::InitVector(m_vecfPlatformStdev, numGenes, (float) 0);
	m_strPlatformName = strPlatformName;
}

void CSeekPlatform::SetPlatformAvg(const size_t &i, float val){
	m_vecfPlatformAvg[i] = val;
}
	
void CSeekPlatform::SetPlatformStdev(const size_t &i, float val){
	m_vecfPlatformStdev[i] = val;
}
	
float CSeekPlatform::GetPlatformAvg(const size_t &i){
	return m_vecfPlatformAvg[i];
}

float CSeekPlatform::GetPlatformStdev(const size_t &i){
	return m_vecfPlatformStdev[i];
}

void CSeekPlatform::ResetPlatform(){
	m_iNumGenes = 0;
	m_vecfPlatformAvg.clear();
	m_vecfPlatformStdev.clear();
	m_strPlatformName = "";
}

}

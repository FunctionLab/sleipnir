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
#ifndef PCLI_H
#define PCLI_H

#include <map>
#include <set>
#include <vector>

#include "meta.h"
#include "file.h"
#include "fullmatrix.h"

namespace Sleipnir {

class CPCLImpl: protected CFile {
protected:
	static const size_t c_iSkip = 2;
	static const char c_szEWEIGHT[];
	static const char c_szGENE[];
	static const char c_szGID[];
	static const char c_szGWEIGHT[];
	static const char c_szNAME[];
	static const char c_szOne[];
	static const char c_szExtension[];
	static const char c_szBinExtension[];
	static const char c_szDabExtension[];

	typedef std::vector<std::string> TVecStr;
	typedef std::set<size_t> TSetI;
	typedef std::map<std::string, size_t> TMapStrI;

	static size_t MedianMultiplesBin(float dValue, float dAve, float dStd,
			size_t iBins, float dBinSize) {
		size_t iRet;
		int i;

		i = (int) (0.5 + ((dValue - dAve) / dStd / dBinSize));
		iRet = iBins / 2;
		iRet = ((i < 0) && ((size_t) -i > iRet)) ? 0 : min(iBins, i + iRet);
		// cerr << dValue << '\t' << dAve << '\t' << dStd << '\t' << i << '\t' << iRet << endl;

		return iRet;
	}

	static void MedianMultiplesSmooth(float dPower,
			std::vector<float>& vecdValues) {
		static const size_t c_iRadius = 40;
		std::vector<float> vecdTmp;
		size_t i, j, k;
		float d, dSum;

		vecdTmp.resize(vecdValues.size());
		for (dSum = 0, i = 0; i < vecdTmp.size(); ++i)
			for (j = (max(i, c_iRadius) - c_iRadius); j < min(vecdTmp.size(), i
					+ c_iRadius); ++j) {
				k = max(i, j) - min(i, j);
				vecdTmp[i] += (d = (vecdValues[j]
						/ (1 + pow((float) k, dPower))));
				dSum += d;
			}
		for (i = 0; i < vecdValues.size(); ++i)
			vecdValues[i] = vecdTmp[i] / dSum;
	}

	CPCLImpl(bool fHeader) :
		m_fHeader(fHeader) {
	}
	~CPCLImpl();

	bool OpenExperiments(std::istream&, size_t, string&);
	bool OpenGene(std::istream&, std::vector<float>&, string&);
	void Reset();
	void MedianMultiplesMapped(const std::vector<std::vector<size_t> >&,
			std::vector<float>&);
	bool OpenHelper();
	bool OpenMemmap(const unsigned char*);

	void SetGene(size_t iGene, const std::string& strGene) {

		m_mapstriGenes.erase(m_vecstrGenes[iGene]);
		m_mapstriGenes[m_vecstrGenes[iGene] = strGene] = iGene;
	}

	CDataMatrix m_Data;
	TVecStr m_vecstrGenes;
	TVecStr m_vecstrExperiments;
	TVecStr m_vecstrFeatures;
	std::vector<TVecStr> m_vecvecstrFeatures;
	TSetI m_setiGenes;
	bool m_fHeader;
	TMapStrI m_mapstriGenes;

	// Memory mapped back end
	unsigned char* m_abData;
	size_t m_iData;
	HANDLE m_hndlData;
	float** m_aadData;
};

}

#endif // PCLI_H

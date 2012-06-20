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
#ifndef SEEKREADER_H
#define SEEKREADER_H

#include "seekmap.h"
#include "stdafx.h"
#include "datapair.h"
#include "seekdataset.h"
#include "database.h"

namespace Sleipnir {

class CSeekTools{
public:
	/* binary */
	template<class tType>
	static bool ReadArray(const char *fileName, vector<tType> &vData){
		FILE *f = fopen(fileName, "rb");
		if(f==NULL){
			cerr << "File not found" << endl;
			return false;
		}
		size_t iSize;
		int ret;
		ret = fread((char*) (&iSize), 1, sizeof(iSize), f);
		vData.clear();
		vData.resize(iSize);
		tType *m_Data = (tType*)malloc(iSize*sizeof(tType));
		ret = fread((char*)m_Data, 1, iSize*sizeof(tType), f);
		size_t i;
		for(i=0; i<iSize; i++){
			vData[i] = m_Data[i];
		}
		free(m_Data);
		fclose(f);
		return true;
	}

	/* binary */
	template<class tType>
	static bool WriteArray(const char *fileName, vector<tType> &vData){
		FILE *f = fopen(fileName, "wb");
		if(f==NULL){
			cerr << "File not found" << endl;
			return false;
		}
		size_t i;
		tType *m_Data = (tType*)malloc(vData.size()*sizeof(tType));
		for(i=0; i<vData.size(); i++){
			m_Data[i] = vData[i];
		}
		size_t iSize = vData.size();
		fwrite((char*) (&iSize), 1, sizeof(iSize), f);
		fwrite((char*) (m_Data), 1, iSize*sizeof(tType), f);
		free(m_Data);
		fclose(f);
		return true;
	}

	template<class tType>
	static bool InitVector(vector<tType> &vData, size_t iSize, tType tValue){
		size_t i;
		vData.clear();
		vData.resize(iSize);
		for(i=0; i<iSize; i++){
			vData[i] = tValue;
		}
		return true;
	}

	static bool CreatePresenceVector(vector<int> &srcData, vector<char> &destData, size_t iSize);
	static bool LoadDatabase(CDatabase &DB, string &strInputDirectory, 
	string &strPrepInputDirectory, vector<char> &cQuery, 
	vector<string> &vecstrQuery, vector<string> &vecstrDatasets, 
	vector<CSeekDataset*> &vc);

};


}
#endif

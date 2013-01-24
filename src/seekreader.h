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

#include "seekbasic.h"
#include "seekmap.h"
#include "datapair.h"
#include "seekdataset.h"
#include "seekplatform.h"
#include "database.h"
#include <sstream>
#include "seeknetwork.h"

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

		//do not change type
		size_t iSize;

		ushort ret;
		ret = fread((char*) (&iSize), 1, sizeof(iSize), f);
		vData.clear();
		vData.resize(iSize);
		tType *m_Data = (tType*)malloc(iSize*sizeof(tType));
		ret = fread((char*)m_Data, 1, iSize*sizeof(tType), f);
		typename vector<tType>::iterator iter;
		tType *mp;
		for(iter=vData.begin(), mp=&m_Data[0]; iter!=vData.end();
			iter++, mp++){
			*iter = *mp;
		}
		free(m_Data);
		fclose(f);
		return true;
	}

	/* binary */
	template<class tType>
	static bool WriteArray(const char *fileName, const vector<tType> &vData){
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
		//do not change type
		size_t iSize = vData.size();
		fwrite((char*) (&iSize), 1, sizeof(iSize), f);
		fwrite((char*) (m_Data), 1, iSize*sizeof(tType), f);
		free(m_Data);
		fclose(f);
		return true;
	}

	template<class tType>
	static bool WriteArrayText(const char *fileName,
		const vector<tType> &vData){
		ofstream outfile;
		outfile.open(fileName);
		size_t i;
		for(i=0; i<vData.size()-1; i++){
			outfile << vData[i] << " ";
		}
		outfile << vData[vData.size()-1] << endl;
		outfile.close();
		return true;
	}

	template<class tType>
	static bool Write2DArrayText(const char *fileName,
		const vector<vector<tType> > &vData){
		ofstream outfile;
		outfile.open(fileName);
		size_t i,j;
		for(j=0; j<vData.size(); j++){
			for(i=0; i<vData[j].size()-1; i++){
				outfile << vData[j][i] << " ";
			}
			outfile << vData[j][vData[j].size()-1] << endl;
		}
		outfile.close();
		return true;
	}

	template<class tType>
	static bool InitVector(vector<tType> &vData, const ushort &iSize,
		const tType &tValue) {
		vData.clear();
		vData.resize(iSize);
		fill(vData.begin(), vData.end(), tValue);
		return true;
	}

	template<class tType>
	static bool InitVector(vector<tType> &vData, const ushort &iSize) {
		vData.clear();
		vData.resize(iSize);
		return true;
	}

	template<class tType>
	static tType** Init2DArray(const size_t &iSize1, const size_t &iSize2,
		const tType &tValue){
		tType **f = (tType**)malloc(iSize1*sizeof(tType*));
		f[0] = (tType*)malloc(iSize1*iSize2*sizeof(tType));
		/*tType **itF = &f[1];
		tType **itLast = &f[0] + iSize1;
		for(; itF!=itLast; itF++){
			*itF = *(itF - 1) + iSize2;
		}
		tType *itVal = &f[0][0];
		tType *itValLast = &f[iSize1-1][iSize2-1] + 1;
		for(; itVal!=itValLast; itVal++){
			*itVal = tValue;
		}*/
		int i, j;
		for(i=1; i<iSize1; i++){
			f[i] = f[i-1] + iSize2;
		}
		for(i=0; i<iSize1; i++){
			for(j=0; j<iSize2; j++){
				f[i][j] = tValue;
			}
		}
		return f;
	}

	template<class tType>
	static void Free2DArray(tType** f){
		free(f[0]);
		free(f);
	}

	static bool IsNaN(const ushort &);

	static string ConvertInt(const int &);

	static bool CreatePresenceVector(const vector<ushort> &, vector<char> &,
		const ushort &);

	static bool ReadDatabaselets(const CDatabase &, 
		const vector< vector<string> > &, vector<CSeekDataset*> &, 
		//network mode options
		const int&, const bool&);

	static bool LoadDatabase(const CDatabase &, const string &,
		const string &, const string &,
		const vector<string> &, const map<string, string> &,
		const map<string, ushort> &, vector<CSeekPlatform> &,
		vector<CSeekDataset*> &);

	static bool LoadDatabase(const CDatabase &, const char *,
		const char *, const char *,
		const vector<string> &, const map<string, string> &,
		const map<string, ushort> &, vector<CSeekPlatform> &,
		vector<CSeekDataset*> &);

	static bool LoadDatabase(const CDatabase &, vector<CSeekDataset*>&,
		const vector<CSeekDataset*>&, vector<CSeekPlatform>&, 
		const vector<CSeekPlatform>&, const vector<string>&, 
		const map<string,string>&, const map<string,ushort>&);

	static bool ReadPlatforms(const string &strPlatformDirectory,
		vector<CSeekPlatform> &plat, vector<string> &vecstrPlatforms,
		map<string, ushort> &mapstriPlatforms);
	static bool ReadPlatforms(const char *plat_dir,
			vector<CSeekPlatform> &plat, vector<string> &vecstrPlatforms,
			map<string, ushort> &mapstriPlatforms);

	static bool ReadListOneColumn(const string &strFile,
		vector<string> &vecstrList, CSeekStrIntMap &mapstriList);
	static bool ReadListOneColumn(const char *file,
			vector<string> &vecstrList, CSeekStrIntMap &mapstriList);

	static bool ReadListTwoColumns(const string &strFile,
		vector<string> &list1, vector<string> &list2);
	static bool ReadListTwoColumns(const char *file,
		vector<string> &list1, vector<string> &list2);

	static bool ReadMultipleQueries(const string &strFile,
		vector< vector<string> > &qList);
	static bool ReadMultipleQueries(const char *file,
			vector< vector<string> > &qList);

	static bool ReadMultiGeneOneLine(const string &strFile,
		vector<string> &list1, const int lineSize = 1024);
	static bool ReadMultiGeneOneLine(const char *file,
			vector<string> &list1, const int lineSize = 1024);

	static bool ReadListOneColumn(const string &strFile,
		vector<string> &vecstrList);
	static bool ReadListOneColumn(const char *file,
			vector<string> &vecstrList);

	static bool ReadQuantFile(const string &strFile, vector<float> &quant);
	static bool ReadQuantFile(const char *file, vector<float> &quant);


};


}
#endif
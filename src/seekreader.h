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

/*!
 * \brief A suite of file I/O and general purpose tools that are used by Seek
 *
 * These tools are critical for initializing the search parameters, and are highly beneficial to the routine
 * manipulations of vectors and files.
 *
 * Some examples of these tools include:
 * \li Reading different search setting files, such as the dataset mapping, the gene mapping,
 * the list of queries, the correlation discretization file, etc.
 * \li Loading a set of CDatabaselet from the given directory
 * \li Reading a binary file that contains a vector of standard type elements (\c int, \c char, \c string,
 * \c float)
 * \li Reading a text file that contains a table with one or two columns
 * \li Setting all elements of a given vector to the given value
 * \li Writing a vector to a file (in binary or in text format)
 */
class CSeekTools{
public:
	/* binary */
	/*!
	 * \brief Read an array from a given binary file
	 *
	 * \param fileName The file name
	 * \param vData The destination array
	 * \return True if the reading is successful
	 * \remark This function reads an one-dimensional array. The binary file needs to be
	 * organized as follows:
	 * 1) The first field is the size of the array, \a N (\c size_t).
	 * 2) The second field is a set of \a N elements.
	 */
	template<class tType>
	static bool ReadArray(const char *fileName, vector<tType> &vData){
		FILE *f = fopen(fileName, "rb");
		if(f==NULL){
			cerr << "File not found" << endl;
			return false;
		}

		//do not change type
		size_t iSize;

		utype ret;
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
	/*!
	 * \brief Write an array in binary format
	 *
	 * \param fileName The file name
	 * \param vData The source array
	 * \return True if the reading is successful
	 * \remark This function writes an one-dimensional array. It will write the array in the following
	 * way:
	 * 1) The first field is the size of the array, \a N (\c size_t).
	 * 2) The second field is the \a N elements in the array.
	 */
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

	/*!
	 * \brief Write an array in text format
	 *
	 * \param fileName The file name
	 * \param vData The source array
	 * \return True if the reading is successful
	 * \remark This function writes an one-dimensional array in the text format. The array elements
	 * are separated by spaces.
	 */
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

	/*!
	 * \brief Write a two-dimensional array in text format
	 *
	 * \param fileName The file name
	 * \param vData The source array
	 * \return True if the reading is successful
	 * \remark This function writes an two-dimensional array in the text format.
	 * The rows are separated by new lines. The elements in a row are separated by spaces.
	 */
	template<class tType>
	static bool Write2DArrayText(const char *fileName,
		const vector<vector<tType> > &vData){
		ofstream outfile;
		outfile.open(fileName);
		size_t i,j;
		for(j=0; j<vData.size(); j++){
			if(vData[j].size()==0){
				outfile << "None" << endl;
				continue;
			}
			for(i=0; i<vData[j].size()-1; i++){
				outfile << vData[j][i] << " ";
			}
			outfile << vData[j][vData[j].size()-1] << endl;
		}
		outfile.close();
		return true;
	}

	/*!
	 * \brief Initialize a vector with a given value
	 *
	 * \param vData The source vector
	 * \param iSize The number of elements that the vector should contain
	 * \param tValue The value
	 *
	 * Resizes the source vector to the given size, then sets all elements in the vector
	 * to the given value.
	 */
	template<class tType>
	static bool InitVector(vector<tType> &vData, const utype &iSize,
		const tType &tValue) {
		vData.clear();
		vData.resize(iSize);
		fill(vData.begin(), vData.end(), tValue);
		return true;
	}

	/*!
	 * \brief Initialize a vector
	 *
	 * \param vData The source vector
	 * \param iSize The number of elements that the vector should contain
	 *
	 * Resizes the source vector to the given size.
	 */
	template<class tType>
	static bool InitVector(vector<tType> &vData, const utype &iSize) {
		vData.clear();
		vData.resize(iSize);
		return true;
	}

	/*!
	 * \brief Initialize a two-dimensional array with the given size and value
	 *
	 * \param iSize1 The first dimension size
	 * \param iSize2 The second dimension size
	 * \param tValue The value
	 *
	 * Creates a two-dimensional array of the given dimension, then populates
	 * it with the given value.
	 */
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

	/*!
	 * \brief Free a two-dimensional array
	 * \param f The two-dimensional array
	 */
	template<class tType>
	static void Free2DArray(tType** f){
		free(f[0]);
		free(f);
	}

	/*!
	 * \brief Checks if a \c utype value is invalid
	 * \param v The value to be checked
	 * A \c utype value is invalid if it is maximum (65535).
	 */
	static bool IsNaN(const utype &);

	/*!
	 * \brief Return the NaN value as a utype
	 */
	static utype GetNaN();

	/*!
	 * \brief Converts an integer to a string
	 * \param number The given integer number
	 * \return The string
	 */
	static string ConvertInt(const int &);

	/*!
	 * \brief Read a set of CDatabaselet from CDatabase instance
	 *
	 * \param DB The CDatabase instance
	 * \param vecstrAllQuery The list of queries
	 * \param vc A vector of datasets (the output)
	 * \param iClient If the network mode is enabled, the client's socket
	 * \param bNetwork If true, the network mode is enabled
	 *
	 * \remarks
	 * A CDatabaselet stores the correlation of a given gene, \a g, to all other genes in all of the
	 * datasets. In order to perform the coexpression search, the CDatabaselet's corresponding
	 * to the query genes need to be loaded from disk. This function reads the
	 * CDatabaselet files corresponding to the query genes.
	 *
	 * \remarks
	 * Once the CDatabaselet for a query gene has been read, the next step that this function
	 * performs is allotting the correlations to their corresponding datasets (CSeekDataset).
	 *
	 * \remarks
	 * The network mode is used to relay status messages between the server and the client.
	 *
	 * \remarks
	 * Assumes that the CSeekTools::LoadDatabase() has been called.
	 */
	static bool ReadDatabaselets(const vector<CDatabase*>&,
		const size_t&, const size_t&,
		const vector<vector<string> >&,
		vector<CSeekDataset*>&,
		const map<string,utype> &,
		const vector<vector<string> > &, const map<string,utype> &,
		//network mode options
		const int&, const bool&);

	/*!
	 * \brief Read the search setting files and load the CDatabase
	 *
	 * Performs the following search initializing operations:
	 * \li Reads the gene presence files \c *.gpres, the gene averages \c *.gavg,
	 * the gene variances \c *.gvar, and each dataset's correlation average and
	 * standard deviation \c *.sinfo.
	 * \li Reads the dataset-platform mapping file
	 * \li Initializes the vector of CSeekDataset with each dataset's gene-presence, gene-averages
	 *
	 * \param DB The CDatabase instance
	 * \param strPrepInputDirectory The prep directory which contains the \c *.gavg and \c *.gpres files
	 * \param strGvarInputDirectory The directory that contains the gene variance files \c *.gvar
	 * \param strSinfoInputDirectory The directory that contains the \c *.sinfo files
	 * \param vecstrDatasets The dataset definition
	 * \param mapstrstrDatasetPlatform The dataset-platform mapping
	 * \param mapstriPlatform Platform name-platform ID mapping
	 * \param vp The vector of CSeekPlatform
	 * \param vc The vector of CSeekDataset, the output
	 *
	 */
	static bool LoadDatabase(const vector<CDatabase*>&,
		const size_t&, const size_t&,
		const vector<CSeekDBSetting*>&,
		const vector<string>&,
		const map<string,string>&,
		const map<string,utype>&, vector<CSeekPlatform>&,
		vector<CSeekDataset*>&, const vector<vector<string> >&,
		const map<string,utype>&,
		const bool=false, const bool=false);

	/*!
	 * \brief Load a CDatabase by copying from an existing instance
	 *
	 * Copies the vector of initialized CSeekDataset to a new vector.
	 * Copies the vector of initialized CSeekPlatform to a new vector.
	 *
	 * \param DB The CDatabase
	 * \param vc The destination dataset vector
	 * \param vc_src The source dataset vector
	 * \param vp The destination platform vector
	 * \param vp_src The source platform vector
	 * \param vecstrDatasets The dataset definition
	 * \param mapstrstrDatasetPlatform The dataset-platform mapping
	 * \param mapstriPlatform Platform name-platform ID mapping
	 */
	static bool LoadDatabase(
		const vector<CDatabase*>&, const size_t&, const size_t&,
		vector<CSeekDataset*>&,
		const vector<CSeekDataset*>&, vector<CSeekPlatform>&, 
		const vector<CSeekPlatform>&, const vector<string>&, 
		const map<string,string>&, const map<string,utype>&);

	/*!
	 * \brief Read the platforms
	 *
	 * Reading the platforms mainly involves reading the correlation average and
	 * the correlation standard deviation for each platform in the database.
	 * The purpose is to correct the platform specific biases on the correlation values.
	 *
	 * \param strPlatformDirectory The directory that contains the platform average and standard deviation files
	 * \param plat The output
	 * \param vecstrPlatforms The platform names
	 * \param mapstriPlatform The platform name - platform ID mapping
	 * \param lineSize The maximum characters per line in the file (default 1024)
	 */
	static bool ReadPlatforms(const string &strPlatformDirectory,
		vector<CSeekPlatform> &plat, vector<string> &vecstrPlatforms,
		map<string, utype> &mapstriPlatforms, const int lineSize = 1024);

	/*!
	 * \brief Read the platforms
	 *
	 * This is the same as the previous CSeekTools::ReadPlatforms() declaration, except that the
	 * accepted string arguments are of the type \c const \c char \c *.
	 */
	static bool ReadPlatforms(const char *plat_dir,
		vector<CSeekPlatform> &plat, vector<string> &vecstrPlatforms,
		map<string, utype> &mapstriPlatforms, const int lineSize = 1024);

	/*!
	 * \brief Read a table with one column
	 *
	 * Outputs the lines in the table as a vector of strings
	 *
	 * \param strFile The file name
	 * \param vecstrList The output
	 * \param mapstriList Mapping the line to its line number
	 * \param lineSize The maximum characters per line in the file (default 1024)
	 */
	static bool ReadListOneColumn(const string &strFile,
		vector<string> &vecstrList, CSeekStrIntMap &mapstriList, const int lineSize = 1024);

	/*!
	 * \brief Read a table with one column
	 *
	 * This is the same as the previous CSeekTools::ReadListOneColumn() declaration, except that the
	 * accepted string arguments are of the type const char*.
	 */
	static bool ReadListOneColumn(const char *file,
		vector<string> &vecstrList, CSeekStrIntMap &mapstriList, const int lineSize = 1024);

	/*!
	 * \brief Read a table with one column
	 *
	 * Same as the previous CSeekTools::ReadListOneColumn() declaration, except that this does not
	 * generate the line to line number mapping.
	 */
	static bool ReadListOneColumn(const string &strFile,
		vector<string> &vecstrList, const int lineSize = 1024);

	/*!
	 * \brief Read a table with one column
	 *
	 * Same as the previous CSeekTools::ReadListOneColumn() declaration, except that this does not
	 * generate the line to line number mapping, and accepts the file name as const char *.
	 */
	static bool ReadListOneColumn(const char *file,
		vector<string> &vecstrList, const int lineSize = 1024);

	/*!
	 * \brief Read a table with two columns
	 *
	 * \param strFile The file name
	 * \param list1 The column 1 output
	 * \param list2 The column 2 output
	 * \param lineSize The maximum characters per line in the file (default 1024)
	 */
	static bool ReadListTwoColumns(const string &strFile,
		vector<string> &list1, vector<string> &list2, const int lineSize = 1024);

	/*!
	 * \brief Read a table with two columns
	 *
	 * This is the same as the previous CSeekTools::ReadListTwoColumns() declaration, except that the
	 * accepted string arguments are of the type const char *.
	 */
	static bool ReadListTwoColumns(const char *file,
		vector<string> &list1, vector<string> &list2, const int lineSize = 1024);

	/*!
	 * \brief Read a list of queries
	 *
	 * A query is specified as a set of gene names delimited by spaces.
	 * A query occupies one line in the file.
	 * \param strFile The file name
	 * \param qList The output
	 * \param lineSize The maximum characters per line in the file (default 1024)
	 */
	static bool ReadMultipleQueries(const string &strFile,
		vector< vector<string> > &qList, const int lineSize = 1024);

	/*!
	 * \brief Read a list of queries
	 *
	 * Same as the previous CSeekTools::ReadMultipleQueries() declaration, except that this function
	 * accepts the string argument as a const char *.
	 */
	static bool ReadMultipleQueries(const char *file,
		vector< vector<string> > &qList, const int lineSize = 1024);

	/*!
	 * \brief Read just one gene-set line
	 *
	 * Reads the first line in the file. The line contains a set of gene names delimited by spaces.
	 * The output is a vector of strings representing the genes in that line.
	 *
	 * \param strFile The file name
	 * \param list1 The output
	 * \param lineSize The maximum characters per line in the file (default 1024)
	 */
	static bool ReadMultiGeneOneLine(const string &strFile,
		vector<string> &list1, const int lineSize = 1024);

	/*!
	 * \brief Read just one gene-set line
	 *
	 * Same as the previous CSeekTools::ReadMultiGeneOneLine() except that the accepted string argument
	 * is of the type const char *.
	 */
	static bool ReadMultiGeneOneLine(const char *file,
		vector<string> &list1, const int lineSize = 1024);

	/*!
	 * \brief Read the correlation discretization
	 *
	 * Specifies how the correlations should be binned. The file contains the bin boundaries separated by spaces.
	 * \param strFile The file name
	 * \param quant The output
	 * \param lineSize The maximum characters per line in the file (default 5000)
	 */
	static bool ReadQuantFile(const string &strFile, vector<float> &quant, const int lineSize = 5000);

	/*!
	 * \brief Read the correlation discretization
	 *
	 * Same as the previous CSeekTools::ReadQuantFile() except that the accepted string argument is of the type
	 * const char *.
	 */
	static bool ReadQuantFile(const char *file, vector<float> &quant, const int lineSize = 5000);


};


}
#endif

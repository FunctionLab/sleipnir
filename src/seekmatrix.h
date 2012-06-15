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
#ifndef SEEKMATRIX_H
#define SEEKMATRIX_H

#include "stdafx.h"

namespace Sleipnir {

class CSeekMatrixTools {
public:
	template<class tType>
	static tType** Init(size_t iRow, size_t iColumn, tType tValue){
		tType **m_Data = (tType**)malloc(iRow*sizeof(tType*));
		m_Data[0] = (tType*)malloc(iRow*iColumn*sizeof(tType));
		size_t i, j;
		for(i=1; i<iRow; i++){
			m_Data[i] = m_Data[i-1] + iColumn;
		}
		for(i=0; i<iRow; i++){
			for(j=0; j<iColumn; j++){
				m_Data[i][j] = tValue;
			}
		}
		return m_Data;
	}

	template<class tType>
	static void Free(tType **m_Data){
		free(m_Data[0]);
		free(m_Data);
	}

};

template<class tType>
class CSeekMatrix{
public:
	CSeekMatrix(size_t iRow, size_t iColumn, tType tValue){
		m_Data = (tType**)malloc(iRow*sizeof(tType*));
		m_Data[0] = (tType*)malloc(iRow*iColumn*sizeof(tType));
		size_t i, j;
		for(i=1; i<iRow; i++){
			m_Data[i] = m_Data[i-1] + iColumn;
		}
		for(i=0; i<iRow; i++){
			for(j=0; j<iColumn; j++){
				m_Data[i][j] = tValue;
			}
		}
		m_iRow = iRow;
		m_iColumn = iColumn;
		m_cIsMatrix = true;
		m_cIsCompacted = true;
		m_veciRowSize.clear();
	}

	CSeekMatrix(size_t iRow){
		m_Data = (tType**)malloc(iRow*sizeof(tType*));
		size_t i;
		m_veciRowSize.clear();
		m_veciRowSize.resize(iRow);
		for(i=0; i<iRow; i++){
			m_Data[i] = NULL;
			m_veciRowSize[i] = 0;
		}
		m_iRow = iRow;
		m_cIsMatrix = false;
		m_cIsCompacted = false;
	}

	void InitializeRow(size_t atX, size_t iSize, tType tValue){
		m_Data[atX] = (tType*)malloc(iSize*sizeof(tType));
		size_t i;
		for(i=0; i<iSize; i++){
			m_Data[atX][i] = tValue;
		}
		m_veciRowSize[atX] = iSize;
		m_cIsMatrix = false;
		m_cIsCompacted = false;
	}

	size_t GetElements(){
		size_t i;
		size_t iTot = 0;
		for(i=0; i<m_iRow; i++){
			iTot+=m_veciRowSize[i];
		}
		return iTot;
	}

	void Clear(){
		size_t i,j;
		for(i=0; i<m_iRow; i++){
			for(j=0; j<m_iColumn; j++){
				m_Data[i][j] = 0;
			}
		}
	}

	bool Compact(){
		if(m_cIsMatrix==true || m_cIsCompacted==true){
			return true;
		}

		size_t i,j;

		size_t iSize = GetElements();
		tType **m2 = (tType**)malloc(m_iRow*sizeof(tType*));
		m2[0] = (tType*)malloc(iSize*sizeof(tType));
		for(j=0; j<m_iRow; j++){
			m2[j] = NULL;
		}

		bool isFirst = true;
		tType *prev = NULL;
		int prev_id = 0;
		for(i=0; i<m_iRow; i++){
			if(m_Data[i]==NULL) continue;
			if(isFirst==true){
				isFirst = false;
				m2[i] = (tType*)malloc(iSize*sizeof(tType));
				prev = m2[i];
				prev_id = i;
			}else{
				m2[i] = prev + m_veciRowSize[prev_id];
				prev = m2[i];
				prev_id = i;
			}
		}

		for(i=0; i<m_iRow; i++){
			if(m_Data[i]==NULL) continue;
			for(j=0; j<m_veciRowSize[i]; j++){
				m2[i][j] = m_Data[i][j];
			}
		}

		for(i=0; i<m_iRow; i++){
			if(m_Data[i]==NULL) continue;
			free(m_Data[i]);
		}
		free(m_Data);

		m_Data = m2;
		m_cIsCompacted = true;

		return true;
	}

	void Free(){
		if(m_cIsMatrix==true){
			free(m_Data[0]);
			free(m_Data);
		}else if(m_cIsCompacted==true){
			int i = -1;
			for(i=0; i<m_iRow; i++){
				if(m_Data[i]==NULL) continue;
				break;
			}
			if(i!=-1 && i<m_iRow){
				free(m_Data[i]);
				free(m_Data);
			}
		}else{
			int i = -1;
			for(i=0; i<m_iRow; i++){
				if(m_Data[i]==NULL) continue;
				free(m_Data[i]);
			}
			free(m_Data);
		}
	}

	~CSeekMatrix(){
		Free();
	}

	tType Get(size_t iRow, size_t iColumn){
		return m_Data[iRow][iColumn];
	}

	tType* GetRow(size_t iRow){
		return m_Data[iRow];
	}

	size_t GetRowSize(size_t iRow){
		return m_veciRowSize[iRow];
	}

	void Set(size_t iRow, size_t iColumn, tType tValue){
		m_Data[iRow][iColumn] = tValue;
	}

private:
	tType **m_Data;
	size_t m_iRow;
	size_t m_iColumn;
	bool m_cIsMatrix;
	bool m_cIsCompacted;
	vector<size_t> m_veciRowSize;
};


}

#endif // SEEKMATRIX_H

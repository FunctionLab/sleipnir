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
#include "seekmap.h"

namespace Sleipnir {

/*
 * Presence Data Structure
 */
CSeekPresence::CSeekPresence(size_t i){
	p = (char*)malloc(i);
	memset(p, 0, i);
	iSize = i;
}

CSeekPresence::CSeekPresence(char *cP, size_t i){
	p = (char*)malloc(i);
	memcpy(p, cP, i);
	iSize = i;
}

CSeekPresence::CSeekPresence(CSeekPresence& cP){
	p = (char*)malloc(cP.iSize);
	memcpy(p, cP.p, cP.iSize);
	iSize = cP.iSize;
}

CSeekPresence::~CSeekPresence(){
	free(p);
	iSize = 0;
}

bool CSeekPresence::Check(size_t i){
	if(p[i]==0){
		return false;
	}
	return true;
}

void CSeekPresence::Set(size_t i){
	p[i] = 1;
}

void CSeekPresence::Clear(size_t i){
	p[i] = 0;
}

void CSeekPresence::Clear(){
	memset(p, 0, iSize);
}

size_t CSeekPresence::GetSize(){
	return iSize;
}

/*
 * IntIntMap Data Structure
 */
CSeekIntIntMap::CSeekIntIntMap(size_t iSize){
	m_iF = (int*)malloc(iSize * sizeof(int));
	m_iR = (int*)malloc(iSize * sizeof(int));
	m_iSize = iSize;
	Clear();
}

CSeekIntIntMap::CSeekIntIntMap(CSeekPresence &cP, bool bReverse){
	m_iSize = cP.GetSize();
	m_iF = (int*)malloc(m_iSize * sizeof(int));
	m_iR = (int*)malloc(m_iSize * sizeof(int));
	Clear();
	Reset(cP, bReverse);
}

CSeekIntIntMap::~CSeekIntIntMap(){
	free(m_iF);
	free(m_iR);
	m_iNumSet = 0;
	m_iSize = 0;
}

int CSeekIntIntMap::GetForward(int i){
	return m_iF[i];
}

int CSeekIntIntMap::GetReverse(int i){
	return m_iR[i];
}

void CSeekIntIntMap::Add(int i){
	int j = m_iNumSet;
	m_iF[i] = j;
	m_iR[j] = i;
	m_iNumSet++;
}

void CSeekIntIntMap::Clear(){
	int i;
	for(i=0; i<m_iSize; i++){
		m_iF[i] = -1;
		m_iR[i] = -1;
	}
	m_iNumSet = 0;
}

void CSeekIntIntMap::Reset(CSeekPresence &cP, bool bReverse){
	int i;
	if(bReverse==false){
		int j = 0;
		for(i=0; i<m_iSize; i++){
			if(cP.Check(i)==true){
				Add(i);
			}
		}
	}else{
		int j = 0;
		for(i=0; i<m_iSize; i++){
			if(cP.Check(i)==false){
				Add(i);
			}
		}
	}
}

/*
 * StrIntMap Data Structure
 */
CSeekStrIntMap::CSeekStrIntMap(){
	m_mapstrint.clear();
	m_mapintstr.clear();
}

CSeekStrIntMap::~CSeekStrIntMap(){}

void CSeekStrIntMap::Set(string s, int i){
	m_mapstrint[s] = i;
	m_mapintstr[i] = s;
}

int CSeekStrIntMap::Get(string s){
	return m_mapstrint[s];
}

string CSeekStrIntMap::Get(int i){
	return m_mapintstr[i];
}

size_t CSeekStrIntMap::GetSize(){
	return m_mapintstr.size();
}

vector<string>& CSeekStrIntMap::GetAllString(){
	vector<string> vecStr;
	vecStr.clear();
	vecStr.resize(GetSize());
	map<string, int>::iterator	iter;
	size_t i = 0;
	for(iter = m_mapstrint.begin(); iter!=m_mapstrint.end(); iter++){
		vecStr[i] = iter->first;
	}
	return vecStr;
}

vector<int>& CSeekStrIntMap::GetAllInteger(){
	vector<int> vecInt;
	vecInt.clear();
	vecInt.resize(GetSize());
	map<int, string>::iterator	iter;
	size_t i = 0;
	for(iter = m_mapintstr.begin(); iter!=m_mapintstr.end(); iter++){
		vecInt[i] = iter->first;
	}
	return vecInt;
}


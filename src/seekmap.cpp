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

CSeekPresence::CSeekPresence(CSeekPresence* cP){
	p = (char*)malloc(cP->iSize);
	memcpy(p, cP->p, cP->iSize);
	iSize = cP->iSize;
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
void CSeekIntIntMap::Initialize(size_t iSize){
	m_iF = (int*)malloc(iSize * sizeof(int));
	m_iR = (int*)malloc(iSize * sizeof(int));
	m_iSize = iSize;
	Clear();
}
CSeekIntIntMap::CSeekIntIntMap(size_t iSize){
	Initialize(iSize);
}

CSeekIntIntMap::CSeekIntIntMap(CSeekPresence *cP, bool bReverse){
	Initialize(cP->GetSize());
	Reset(cP, bReverse);
}

CSeekIntIntMap::CSeekIntIntMap(vector<char> &cP, bool bReverse){
	Initialize(cP.size());
	Reset(cP, bReverse);
}


CSeekIntIntMap::CSeekIntIntMap(char *cP, size_t iSize, bool bReverse){
	Initialize(iSize);
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

int CSeekIntIntMap::GetNumSet(){
	return m_iNumSet;
}

void CSeekIntIntMap::Reset(CSeekPresence *cP, bool bReverse){
	int i;
	if(bReverse==false){
		int j = 0;
		for(i=0; i<m_iSize; i++){
			if(cP->Check(i)==true){
				Add(i);
			}
		}
	}else{
		int j = 0;
		for(i=0; i<m_iSize; i++){
			if(cP->Check(i)==false){
				Add(i);
			}
		}
	}
}

void CSeekIntIntMap::Reset(char *cP, bool bReverse){
	int i;
	if(bReverse==false){
		int j = 0;
		for(i=0; i<m_iSize; i++){
			if(cP[i]==1){
				Add(i);
			}
		}
	}else{
		int j = 0;
		for(i=0; i<m_iSize; i++){
			if(cP[i]==0){
				Add(i);
			}
		}
	}
}

void CSeekIntIntMap::Reset(vector<char> &cP, bool bReverse){
	int i;
	if(bReverse==false){
		int j = 0;
		for(i=0; i<m_iSize; i++){
			if(cP[i]==1){
				Add(i);
			}
		}
	}else{
		int j = 0;
		for(i=0; i<m_iSize; i++){
			if(cP[i]==0){
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

void CSeekStrIntMap::Clear(){
	m_mapstrint.clear();
	m_mapintstr.clear();
}

void CSeekStrIntMap::SetAll(vector<string> &s){
	Clear();
	size_t i = 0;
	for(i=0; i<s.size(); i++){
		m_mapstrint[s[i]] = i;
		m_mapintstr[i] = s[i];
	}
}

void CSeekStrIntMap::Set(string s, size_t i){
	m_mapstrint[s] = i;
	m_mapintstr[i] = s;
}

map<string, size_t>& CSeekStrIntMap::GetMapForward(){
	return m_mapstrint;
}

map<size_t, string>& CSeekStrIntMap::GetMapReverse(){
	return m_mapintstr;
}


int CSeekStrIntMap::Get(string s){
	return m_mapstrint[s];
}

string CSeekStrIntMap::Get(size_t i){
	return m_mapintstr[i];
}

size_t CSeekStrIntMap::GetSize(){
	return m_mapintstr.size();
}

vector<string> CSeekStrIntMap::GetAllString(){
	vector<string> vecStr;
	vecStr.clear();
	vecStr.resize(GetSize());
	map<string, size_t>::iterator	iter;
	size_t i = 0;
	for(iter = m_mapstrint.begin(); iter!=m_mapstrint.end(); iter++){
		vecStr[i] = iter->first;
		i++;
	}
	return vecStr;
}

vector<size_t> CSeekStrIntMap::GetAllInteger(){
	vector<size_t> vecInt;
	vecInt.clear();
	vecInt.resize(GetSize());
	map<size_t, string>::iterator	iter;
	size_t i = 0;
	for(iter = m_mapintstr.begin(); iter!=m_mapintstr.end(); iter++){
		vecInt[i] = iter->first;
		i++;
	}
	return vecInt;
}
}


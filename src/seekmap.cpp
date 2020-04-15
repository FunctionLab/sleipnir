/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a part of SEEK (Search-based exploration of expression compendium)
* which is authored and maintained by: Qian Zhu (qzhu@princeton.edu)
*
* If you use this file, please cite the following publication:
* Qian Zhu, Aaron K Wong, Arjun Krishnan, Miriam R Aure, Alicja Tadych, 
* Ran Zhang, David C Corney, Casey S Greene, Lars A Bongo, 
* Vessela N Kristensen, Moses Charikar, Kai Li & Olga G Troyanskaya
* "Targeted exploration and analysis of large cross-platform human 
* transcriptomic compendia" Nat Methods (2015)
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library for development, or use any other Sleipnir executable
* tools, please also cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "seekmap.h"

namespace Sleipnir {

/*
 * IntIntMap Data Structure
 */
void CSeekIntIntMap::Initialize(const utype &iSize){
	m_iF.resize(iSize);
	m_iR.resize(iSize);
	m_iSize = iSize;
	Clear();
}
CSeekIntIntMap::CSeekIntIntMap(const utype &iSize){
	Initialize(iSize);
}

const vector<utype>& CSeekIntIntMap::GetAllForward() const{
	return m_iF;
}

const vector<utype>& CSeekIntIntMap::GetAllReverse() const{
	return m_iR;
}

CSeekIntIntMap::CSeekIntIntMap(const vector<char> &cP, const bool bReverse){
	Initialize(cP.size());
	Reset(cP, bReverse);
}


CSeekIntIntMap::CSeekIntIntMap(const char *cP, const utype &iSize,
	const bool bReverse){
	Initialize(iSize);
	Reset(cP, bReverse);
}

CSeekIntIntMap::CSeekIntIntMap(CSeekIntIntMap *a){
	m_iNumSet = a->m_iNumSet;
	m_iSize = a->m_iSize;	
	m_iF.resize(a->m_iF.size());
	m_iR.resize(a->m_iR.size());
	copy(a->m_iF.begin(), a->m_iF.end(), m_iF.begin());
	copy(a->m_iR.begin(), a->m_iR.end(), m_iR.begin());
	m_iterR = m_iR.begin() + m_iNumSet;
}

CSeekIntIntMap::~CSeekIntIntMap(){
	m_iF.clear();
	m_iR.clear();
	m_iNumSet = 0;
	m_iSize = 0;
}

utype CSeekIntIntMap::GetForward(const utype &i) const{
	return m_iF[i];
}

utype CSeekIntIntMap::GetReverse(const utype &i) const{
	return m_iR[i];
}

void CSeekIntIntMap::Add(const utype &i){
	m_iF[i] = m_iNumSet;
	*m_iterR = i;
	m_iterR++;
	m_iNumSet++;
}

void CSeekIntIntMap::Clear(){
	vector<utype>::iterator iterF = m_iF.begin();
	vector<utype>::iterator iterR = m_iR.begin();
	for(; iterF!=m_iF.end(); iterF++, iterR++){
		*iterF = -1;
		*iterR = -1;
	}
	m_iNumSet = 0;
	m_iterR = m_iR.begin();
}

utype CSeekIntIntMap::GetNumSet() const{
	return m_iNumSet;
}

utype CSeekIntIntMap::GetSize() const{
	return m_iSize;
}

void CSeekIntIntMap::Reset(const char *cP, const bool bReverse){
	utype i;
	if(bReverse==false){
		for(i=0; i<m_iSize; i++){
			if(cP[i]==1){
				Add(i);
			}
		}
	}else{
		for(i=0; i<m_iSize; i++){
			if(cP[i]==0){
				Add(i);
			}
		}
	}
}

void CSeekIntIntMap::Reset(const vector<char> &cP, const bool bReverse){
	utype i;
	if(bReverse==false){
		for(i=0; i<m_iSize; i++){
			if(cP[i]==1){
				Add(i);
			}
		}
	}else{
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

CSeekStrIntMap::~CSeekStrIntMap(){
	m_mapstrint.clear();
	m_mapintstr.clear();
}

void CSeekStrIntMap::Clear(){
	m_mapstrint.clear();
	m_mapintstr.clear();
}

void CSeekStrIntMap::SetAll(const vector<string> &s){
	Clear();
	utype i = 0;
	for(i=0; i<s.size(); i++){
		m_mapstrint[s[i]] = i;
		m_mapintstr[i] = s[i];
	}
}

void CSeekStrIntMap::Set(const string &s, const utype &i){
	m_mapstrint[s] = i;
	m_mapintstr[i] = s;
}

map<string, utype>& CSeekStrIntMap::GetMapForward(){
	return m_mapstrint;
}

map<utype, string>& CSeekStrIntMap::GetMapReverse(){
	return m_mapintstr;
}


utype CSeekStrIntMap::Get(const string &s) const{
	map<string, utype>::const_iterator	iter = m_mapstrint.find(s);
	return iter->second;
}

string CSeekStrIntMap::Get(const utype &i) const{
	map<utype, string>::const_iterator	iter = m_mapintstr.find(i);
	return iter->second;
}

utype CSeekStrIntMap::GetSize() const{
	return m_mapintstr.size();
}

vector<string> CSeekStrIntMap::GetAllString() const{
	vector<string> vecStr;
	vecStr.clear();
	vecStr.resize(GetSize());
	map<string, utype>::const_iterator	iter = m_mapstrint.begin();
	vector<string>::iterator iterV = vecStr.begin();
	for(; iter!=m_mapstrint.end(); iter++, iterV++)
		*iterV = iter->first;
	return vecStr;
}

vector<utype> CSeekStrIntMap::GetAllInteger() const{
	vector<utype> vecInt;
	vecInt.clear();
	vecInt.resize(GetSize());
	map<utype, string>::const_iterator	iter = m_mapintstr.begin();
	vector<utype>::iterator iterV = vecInt.begin();
	for(; iter!=m_mapintstr.end(); iter++, iterV++)
		*iterV = iter->first;
	return vecInt;
}
}


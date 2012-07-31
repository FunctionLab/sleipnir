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
#include "seekmap.h"

namespace Sleipnir {

/*
 * IntIntMap Data Structure
 */
void CSeekIntIntMap::Initialize(const ushort &iSize){
	m_iF.resize(iSize);
	m_iR.resize(iSize);
	m_iSize = iSize;
	Clear();
	m_iterR = m_iR.begin();
}
CSeekIntIntMap::CSeekIntIntMap(const ushort &iSize){
	Initialize(iSize);
}

const vector<ushort>& CSeekIntIntMap::GetAllForward() const{
	return m_iF;
}

const vector<ushort>& CSeekIntIntMap::GetAllReverse() const{
	return m_iR;
}

CSeekIntIntMap::CSeekIntIntMap(const vector<char> &cP, const bool bReverse){
	Initialize(cP.size());
	Reset(cP, bReverse);
}


CSeekIntIntMap::CSeekIntIntMap(const char *cP, const ushort &iSize,
	const bool bReverse){
	Initialize(iSize);
	Reset(cP, bReverse);
}


CSeekIntIntMap::~CSeekIntIntMap(){
	m_iF.clear();
	m_iR.clear();
	m_iNumSet = 0;
	m_iSize = 0;
}

ushort CSeekIntIntMap::GetForward(const ushort &i) const{
	return m_iF[i];
}

ushort CSeekIntIntMap::GetReverse(const ushort &i) const{
	return m_iR[i];
}

void CSeekIntIntMap::Add(const ushort &i){
	m_iF[i] = m_iNumSet;
	*m_iterR = i;
	m_iterR++;
	m_iNumSet++;
}

void CSeekIntIntMap::Clear(){
	vector<ushort>::iterator iterF = m_iF.begin();
	vector<ushort>::iterator iterR = m_iR.begin();
	for(; iterF!=m_iF.end(); iterF++, iterR++){
		*iterF = -1;
		*iterR = -1;
	}
	m_iNumSet = 0;
}

ushort CSeekIntIntMap::GetNumSet() const{
	return m_iNumSet;
}

ushort CSeekIntIntMap::GetSize() const{
	return m_iSize;
}

void CSeekIntIntMap::Reset(const char *cP, const bool bReverse){
	ushort i;
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
	ushort i;
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
	ushort i = 0;
	for(i=0; i<s.size(); i++){
		m_mapstrint[s[i]] = i;
		m_mapintstr[i] = s[i];
	}
}

void CSeekStrIntMap::Set(const string &s, const ushort &i){
	m_mapstrint[s] = i;
	m_mapintstr[i] = s;
}

map<string, ushort>& CSeekStrIntMap::GetMapForward(){
	return m_mapstrint;
}

map<ushort, string>& CSeekStrIntMap::GetMapReverse(){
	return m_mapintstr;
}


ushort CSeekStrIntMap::Get(const string &s) const{
	map<string, ushort>::const_iterator	iter = m_mapstrint.find(s);
	return iter->second;
}

string CSeekStrIntMap::Get(const ushort &i) const{
	map<ushort, string>::const_iterator	iter = m_mapintstr.find(i);
	return iter->second;
}

ushort CSeekStrIntMap::GetSize() const{
	return m_mapintstr.size();
}

vector<string> CSeekStrIntMap::GetAllString() const{
	vector<string> vecStr;
	vecStr.clear();
	vecStr.resize(GetSize());
	map<string, ushort>::const_iterator	iter = m_mapstrint.begin();
	vector<string>::iterator iterV = vecStr.begin();
	for(; iter!=m_mapstrint.end(); iter++, iterV++)
		*iterV = iter->first;
	return vecStr;
}

vector<ushort> CSeekStrIntMap::GetAllInteger() const{
	vector<ushort> vecInt;
	vecInt.clear();
	vecInt.resize(GetSize());
	map<ushort, string>::const_iterator	iter = m_mapintstr.begin();
	vector<ushort>::iterator iterV = vecInt.begin();
	for(; iter!=m_mapintstr.end(); iter++, iterV++)
		*iterV = iter->first;
	return vecInt;
}
}


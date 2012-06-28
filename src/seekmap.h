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
#ifndef SEEKMAP_H
#define SEEKMAP_H

#include "stdafx.h"

namespace Sleipnir {

class CSeekPresence{
public:
	CSeekPresence(size_t);
	CSeekPresence(char*, size_t);
	CSeekPresence(CSeekPresence*);
	~CSeekPresence();
	void Clear();
	bool Check(size_t);
	void Set(size_t);
	void Clear(size_t);
	size_t GetSize();
private:
	char *p;
	int iSize;
};

class CSeekIntIntMap{
public:
	CSeekIntIntMap(size_t);
	CSeekIntIntMap(CSeekPresence*, bool=false);
	CSeekIntIntMap(vector<char>&, bool=false);
	CSeekIntIntMap(char*, size_t, bool=false);
	void Initialize(size_t);

	~CSeekIntIntMap();
	int GetForward(int);
	int GetReverse(int);
	void Add(int);
	void Clear();
	void Reset(CSeekPresence*, bool=false);
	void Reset(vector<char>&, bool=false);
	void Reset(char*, bool=false);
	int GetNumSet();

private:
	int *m_iF;
	int *m_iR;
	int m_iSize;
	int m_iNumSet;
};

class CSeekStrIntMap{
public:
	CSeekStrIntMap();
	~CSeekStrIntMap();
	void Clear();
	void Set(string, size_t);
	void SetAll(vector<string>&);
	int Get(string);
	map<string, size_t>& GetMapForward();
	map<size_t, string>& GetMapReverse();
	size_t GetSize();
	string Get(size_t);
	vector<string> GetAllString();
	vector<size_t> GetAllInteger();
private:
	map<string, size_t> m_mapstrint;
	map<size_t, string> m_mapintstr;
};

}
#endif

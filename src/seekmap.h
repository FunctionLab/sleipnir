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
	CSeekPresence(CSeekPresence&);
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
	CSeekIntIntMap(CSeekPresence&, bool=false);
	~CSeekIntIntMap();
	int GetForward(int);
	int GetReverse(int);
	void Add(int);
	void Clear();
	void Reset(CSeekPresence&, bool=false);

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
	void Set(string, int);
	int Get(string);
	size_t GetSize();
	string Get(int);
	vector<string>& GetAllString();
	vector<int>& GetAllInteger();
private:
	map<string, int> m_mapstrint;
	map<int, string> m_mapintstr;
};

}
#endif

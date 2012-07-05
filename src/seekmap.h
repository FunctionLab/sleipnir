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
/*
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
*/

class CSeekIntIntMap{
public:
	CSeekIntIntMap(const ushort&);
	//CSeekIntIntMap(const CSeekPresence*, const bool=false);
	CSeekIntIntMap(const vector<char>&, const bool=false);
	CSeekIntIntMap(const char*, const ushort &, const bool=false);
	void Initialize(const ushort&);

	~CSeekIntIntMap();
	ushort GetForward(const ushort &) const;
	ushort GetReverse(const ushort &) const;
	const vector<ushort>& GetAllForward() const;
	const vector<ushort>& GetAllReverse() const;

	void Add(const ushort&);
	void Clear();
	//void Reset(const CSeekPresence*, const bool=false) const;
	void Reset(const vector<char>&, const bool=false);
	void Reset(const char*, const bool=false);
	ushort GetNumSet() const;

private:
	vector<ushort> m_iF;
	vector<ushort> m_iR;
	vector<ushort>::iterator m_iterR;
	ushort m_iSize;
	ushort m_iNumSet;
};

class CSeekStrIntMap{
public:
	CSeekStrIntMap();
	~CSeekStrIntMap();
	void Clear();
	void Set(const string&, const ushort&);
	void SetAll(const vector<string>&);
	ushort Get(const string&) const;
	map<string, ushort>& GetMapForward();
	map<ushort, string>& GetMapReverse();
	ushort GetSize() const;
	string Get(const ushort &) const;
	vector<string> GetAllString() const;
	vector<ushort> GetAllInteger() const;
private:
	map<string, ushort> m_mapstrint;
	map<ushort, string> m_mapintstr;
};

}
#endif

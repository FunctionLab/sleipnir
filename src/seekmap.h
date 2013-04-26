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

#include "seekbasic.h"

namespace Sleipnir {

/*!
 * \brief An integer to integer mapping structure
 *
 * This map is used to conveniently get a set of available genes in a dataset, and
 * to quickly check if a given gene is present in the dataset.
 * Normally, a user needs to create two separate arrays for this purpose.
 * The first array is the presence vector (0 or 1) that is used to track the presence (1) or
 * absence (0) of each gene in a given dataset. But this array cannot be used to efficiently
 * get all available genes in the datasets, because it requires scanning through the entire
 * presence vector.
 * The second array contains only available genes in the dataset. But it cannot be used to
 * efficiently check whether a gene is present, because again it requires a scan-through.
 *
 * The solution is to use this int-int mapping structure. It encapsulates the two arrays that
 * we want: \a forward and \a reverse. As an example, let us consider a simple scenario with a
 * genome of 5 genes (0, 1, 2, 3, 4), and a dataset \a d contains 3 of 5 genes
 * (let's say 1, 3, 4).
 * We want to use this map structure to capture the presence of genes in \a d.
 * The \a forward array would contain: [-1, 0, -1, 1, 2], where -1 means the gene is absent,
 * otherwise the value means the \a n-th gene in the array.
 * The \a reverse array would contain: [1, 3, 4]. (ie., only the available genes).
 *
 * These arrays are automatically updated as genes are added to the map.
 *
 */
class CSeekIntIntMap{
public:
	/*!
	 * \brief Constructor
	 * \param iSize The number of genes in the gene-database
	 */
	CSeekIntIntMap(const ushort&);

	/*!
	 * \brief Constructor
	 * \param cP The gene presence vector (char type)
	 * \param bReverse When creating the map, whether or not to follow the reverse logic.
	 * (eg 1 means absent, 0 means present)
	 * \remark
	 * If \c bReverse is true, then this map captures only the absent genes in the dataset.
	 * By default, \c bReverse is false.
	 */
	CSeekIntIntMap(const vector<char>&, const bool=false);

	/*!
	 * \brief Constructor
	 * \param cP The gene presence array (char* type)
	 * \param iSize The size of the gene-presence array
	 * \param bReverse When creating the map, whether or not to follow the reverse logic.
	 * (eg 1 means absent, 0 means present)
	 * \remark
	 * If \c bReverse is true, then this map captures only the absent genes in the dataset.
	 * By default, \c bReverse is false.
	 */
	CSeekIntIntMap(const char*, const ushort &, const bool=false);

	/*!
	 * \brief Copy constructor
	 */
	CSeekIntIntMap(CSeekIntIntMap*);

	/*!
	 * \brief Helper function that is used by constructor
	 */
	void Initialize(const ushort&);

	/*!
	 * \brief Destructor
	 */
	~CSeekIntIntMap();

	/*!
	 * \brief Get an element from the \a forward array
	 * \param i element index
	 */
	ushort GetForward(const ushort &) const;

	/*!
	 * \brief Get an element from the \a reverse array
	 * \param i element index
	 */
	ushort GetReverse(const ushort &) const;

	/*!
	 * \brief Get the entire \a forward array
	 */
	const vector<ushort>& GetAllForward() const;

	/*!
	 * \brief Get the entire \a reverse array
	 */
	const vector<ushort>& GetAllReverse() const;

	/*!
	 * \brief Add an available gene to the map
	 */
	void Add(const ushort&);

	/*!
	 * \brief Clear the member arrays in the structure
	 */
	void Clear();

	/*!
	 * \brief Reset function
	 */
	void Reset(const vector<char>&, const bool=false);

	/*!
	 * \brief Reset function
	 */
	void Reset(const char*, const bool=false);

	/*!
	 * \brief Get the number of present genes that are currently contained in the map
	 */
	ushort GetNumSet() const;

	/*!
	 * \brief Get the genome size
	 */
	ushort GetSize() const;

private:
	vector<ushort> m_iF;
	vector<ushort> m_iR;
	vector<ushort>::iterator m_iterR;
	ushort m_iSize;
	ushort m_iNumSet;
};

/*!
 * \brief A string to integer mapping structure
 *
 * Adds <string, integer> pairs into the map.
 * Supports two major operations:
 * \c Get(string): returns the corresponding integer
 * \c Get(integer): returns the corresponding string
 */
class CSeekStrIntMap{
public:
	/*!
	 * \brief Constructor
	 */
	CSeekStrIntMap();
	/*!
	 * \brief Destructor
	 */
	~CSeekStrIntMap();
	/*!
	 * \brief Clear function
	 */
	void Clear();
	/*!
	 * \brief Add a pair to the map
	 * \param s The string
	 * \param i The integer in unsigned short (ushort)
	 */
	void Set(const string&, const ushort&);
	/*!
	 * \brief Add all the pairs at once
	 * \param s A vector of string
	 * \remarks The corresponding integers are the indices of the strings
	 * in the vector \c s.
	 */
	void SetAll(const vector<string>&);
	/*!
	 * \brief Get the corresponding integer for the given string
	 */
	ushort Get(const string&) const;
	/*!
	 * \brief Get the entire map with key=string, value=integer
	 */
	map<string, ushort>& GetMapForward();
	/*!
	 * \brief Get the entire map with key=integer, value=string
	 */
	map<ushort, string>& GetMapReverse();
	/*!
	 * \brief Get the genome size
	 */
	ushort GetSize() const;
	/*!
	 * \brief Get the corresponding string for the given integer
	 */
	string Get(const ushort &) const;
	/*!
	 * \brief Retrieve all the strings in the map as a vector
	 */
	vector<string> GetAllString() const;
	/*!
	 * \brief Retrieve all the integers in the map as a vector
	 */
	vector<ushort> GetAllInteger() const;
private:
	map<string, ushort> m_mapstrint;
	map<ushort, string> m_mapintstr;
};

}
#endif

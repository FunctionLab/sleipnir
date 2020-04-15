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
#ifndef SEEKMAP_H
#define SEEKMAP_H

#include "seekbasic.h"

namespace Sleipnir {

/*!
 * \brief An integer to integer mapping structure that is used to detect the presence of genes in a dataset
 *
 * This map is used to conveniently get a set of available genes in a dataset, and
 * to quickly check if a given gene is present in the dataset.
 * Normally, a user needs to create two separate arrays for this purpose:
 *
 * \li The first array is the presence vector (0 or 1) that is used to track the presence (1) or
 * absence (0) of each gene in a given dataset. But this array cannot be used to efficiently
 * get all available genes in the datasets, because it requires scanning through the entire
 * presence vector.
 *
 * \li The second array contains only available genes in the dataset. But it cannot be used to
 * efficiently check whether a gene is present, because again it requires a scan-through.
 *
 * The solution is to use this int-int mapping structure. 
 *
 * CSeekIntIntMap encapsulates the two arrays that we want: \a forward and \a reverse. 
 *
 * As an example, let us consider a simple scenario where we have a
 * genome of 5 genes {0, 1, 2, 3, 4}. We want to use CSeekIntIntMap to capture the gene-presence for a
 * dataset \a d which contains the following genes: {1, 3, 4}.
 * Then the two arrays would be:
 * \li The \a forward array: f = [-1, 0, -1, 1, 2], where -1 means that the gene is absent;
 * a non -1 value, \a n, means that the gene is present, and it is the \a n-th present gene in the array.
 * <br>
 * So: the absent genes {0, 2} => f[0] = -1, f[2] = -1. 
 * <br>
 * the present genes {1, 3, 4} => f[1] = 0, f[3] = 1, f[4] = 2.
 *
 * \li The \a reverse array: r = [1, 3, 4]. (ie., only the available genes).
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
	CSeekIntIntMap(const utype&);

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
	CSeekIntIntMap(const char*, const utype &, const bool=false);

	/*!
	 * \brief Copy constructor
	 * \param a A map instance
	 */
	CSeekIntIntMap(CSeekIntIntMap*);

	/*!
	 * \brief Helper function that is used by constructor
	 * \param iSize The genome size
	 */
	void Initialize(const utype&);

	/*!
	 * \brief Destructor
	 */
	~CSeekIntIntMap();

	/*!
	 * \brief Get an element from the \a forward array
	 * \param i Element index
	 * \return The item at the index
	 */
	utype GetForward(const utype &) const;

	/*!
	 * \brief Get an element from the \a reverse array
	 * \param i Element index
	 * \return The item at the index
	 */
	utype GetReverse(const utype &) const;

	/*!
	 * \brief Get the entire \a forward array
	 * \return The \a forward array
	 */
	const vector<utype>& GetAllForward() const;

	/*!
	 * \brief Get the entire \a reverse array
	 * \return The \a reverse array
	 */
	const vector<utype>& GetAllReverse() const;

	/*!
	 * \brief Add an available gene to the map
	 * \param i The gene ID to be added
	 * \remark The gene ID is a number between 0 and 21000 (the genome size). 
	 * It is specified by the gene ID mapping file \c gene_map.txt.
	 */
	void Add(const utype&);

	/*!
	 * \brief Clear the member arrays in the structure
	 * 
	 * This is used by the destructor.
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
	 * \return The number of genes that are present
	 */
	utype GetNumSet() const;

	/*!
	 * \brief Get the genome size
	 * \return The genome size
	 */
	utype GetSize() const;

private:
	vector<utype> m_iF;
	vector<utype> m_iR;
	vector<utype>::iterator m_iterR;
	utype m_iSize;
	utype m_iNumSet;
};

/*!
 * \brief A string to integer mapping structure
 *
 * Adds <string, integer> pairs into the map.
 *
 * Supports two major operations:
 *
 * \c Get(string): returns the corresponding integer
 *
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
	 * \param i The integer in utype
	 */
	void Set(const string&, const utype&);
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
	utype Get(const string&) const;
	/*!
	 * \brief Get the entire map with key=string, value=integer
	 */
	map<string, utype>& GetMapForward();
	/*!
	 * \brief Get the entire map with key=integer, value=string
	 */
	map<utype, string>& GetMapReverse();
	/*!
	 * \brief Get the genome size
	 */
	utype GetSize() const;
	/*!
	 * \brief Get the corresponding string for the given integer
	 */
	string Get(const utype &) const;
	/*!
	 * \brief Retrieve all the strings in the map as a vector
	 */
	vector<string> GetAllString() const;
	/*!
	 * \brief Retrieve all the integers in the map as a vector
	 */
	vector<utype> GetAllInteger() const;
private:
	map<string, utype> m_mapstrint;
	map<utype, string> m_mapintstr;
};

}
#endif

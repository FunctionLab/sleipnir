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
#ifndef SEEKQUERY_H
#define SEEKQUERY_H

#include "seekbasic.h"

namespace Sleipnir {
    
/*!
 * \brief
 * A query structure that is used by Seek
 *
 * Includes vectors for storing query genes, and utilities for partitioning query
 * genes into specified number of groups
 *
 * There are two ways to represent a query.
 *
 * \li A presence char vector, with number of elements = size of genome.
 * All elements are 0 except the elements indexed by the query genes, which have
 * a value of 1.
 *
 * \li A ushort vector, with number of elements = size of query.
 * A compact representation which only stores the query genes' ID. 
 */
class CSeekQuery{
public:
    
    /*!
     * \brief
     * Query partitioning mode
     *
     * Let us consider a query of size \a N.
     *
     * \c LEAVE_ONE_IN
     * : create \a N partitions, where each partition is one of the query genes
     *
     * \c LEAVE_ONE_OUT
     * : create \a N partitions, where each partition is everything but
     * excluding one of the query genes
     *
     * \c CUSTOM_PARTITION
     * : create \a X equally sized partitions, where each partition is
     * \a N / \a X number of genes. \a X is a user-given parameter.
     */
    enum PartitionMode{
        LEAVE_ONE_IN = 0,
        LEAVE_ONE_OUT = LEAVE_ONE_IN + 1,
        CUSTOM_PARTITION = LEAVE_ONE_OUT + 1
    };
    
	CSeekQuery();
	~CSeekQuery();

	bool Reset();

    /*!
     * \brief
     * Initialize with a 0-1 query presence vector
     *
     * \param query
     * A \c char-vector (0 or 1) that specifies the location of the query genes, 
     * which have a value of 1
     *
     * \remarks
     * In preparing the parameter query, we assume that genes have been mapped 
     * to integers between 0 to 21000, or whatever is pre-defined. Then a char-
     * vector of 21000 elements is constructed with the elements located at the
     * the query genes' ID having a value of 1, and the rest of elements being 0.
     */
	bool InitializeQuery(const vector<char>&);
    
    /*!
     * \brief
     * Initialize with a vector of query genes' ID
     *
     * \param query
     * A \c ushort-vector that stores the query genes' ID
     *
     * \param iGenes
     * The number of genes in the genome
     *
     * \remarks
     * This is an alternative way of defining query input. An example of a
     * valid \c query parameter is <201, 242, 42>, 
     * and \c iGenes is 21000 (if there are 21000 genes).
     */
	bool InitializeQuery(const vector<ushort>&, const ushort &);

    /*!
     * \brief
     * Return the number of query partitions
     */
	ushort GetNumFold() const;
    
    /*!
     * \brief
     * Return the query genes as a ushort-vector
     */
	const vector<ushort>& GetQuery() const;
    
    /*!
     * \brief
     * Return the query genes as a presence char-vector (0 or 1)
     */
	const vector<char>& GetQueryPresence() const;
    
    /*!
     * \brief
     * Get the \c i-th query partition
     *
     * \param i
     * The index of partition to return
     *
     * \remarks
     * No bound checking on \c i.
     */
	const vector<ushort>& GetCVQuery(ushort&) const;
    
    /*!
     * \brief
     * Create \c N random query partitions
     *
     * \param rnd
     * A random number generator
     * 
     * \param p
     * The partitioning mode
     *
     * \param iFold
     * The number of partitions to be created
     *
     * \remarks
     * The parameter \c iFold is not used if the partitioning mode is \c LEAVE_ONE_IN 
     * or \c LEAVE_ONE_OUT, because the number of partitions in this case is the 
     * size of the query.
     *
     */
	bool CreateCVPartitions(const gsl_rng*, \
		const enum PartitionMode &, const ushort=-1);

private:
	vector<ushort> queryGenes;
	vector<char> queryGenePresence;

	vector<ushort> *crossValGenes;
	ushort iNumFold;
	ushort iFoldSize;
	ushort qSize;

};

}
#endif

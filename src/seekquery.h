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
#ifndef SEEKQUERY_H
#define SEEKQUERY_H

#include "seekbasic.h"

namespace Sleipnir {

/*!
 * \brief
 * Representation of a query gene-set that is used by Seek
 *
 * Includes vectors for storing the query genes, and utilities for partitioning query
 * genes into specified number of groups
 *
 * There are two ways to represent a query.
 *
 * \li A presence \c char vector, with number of elements = size of genome.
 * All elements are 0 except the elements indexed by the query genes, which have
 * a value of 1.
 *
 * \li A \c utype vector, with number of elements = size of query.
 * A compact representation which only stores the query genes' ID. 
 */
    class CSeekQuery {
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
         * : create \a N partitions, where each partition is everything
         * excluding one of the query genes
         *
         * \c CUSTOM_PARTITION
         * : create \a X equally sized partitions, where each partition is
         * \a N / \a X number of genes. \a X is a user-given parameter.
         */
        enum PartitionMode {
            LEAVE_ONE_IN = 0,
            LEAVE_ONE_OUT = LEAVE_ONE_IN + 1,
            CUSTOM_PARTITION = LEAVE_ONE_OUT + 1
        };

        /*!
         * \brief Constructor
         */
        CSeekQuery();

        /*!
         * \brief Destructor
         */
        ~CSeekQuery();

        /*!
         * \brief The reset function
         */
        void Reset();

        /*!
         * \brief
         * Initialize with a 0-1 query presence vector
         *
         * \param query
         * A \c char vector (0 or 1) that specifies the location of the query genes,
         *
         * \remarks
         * The parameter \c query: based on the \c gene_map.txt, all genes are mapped
         * to integers between 0 to 21000 (or whatever the upper limit is). In the \c char
         * vector \c query, \c query[q] = 1 for \c q in \c Q. All other genes: \c query[g] = 0.
         */
        bool InitializeQuery(const vector<char> &);

        /*!
         * \brief
         * Initialize with a vector of query genes' ID
         *
         * \param query
         * A \c utype-vector that stores the query genes' ID
         *
         * \param iGenes
         * The number of genes in the genome
         *
         * \remarks
         * This is an alternative way of defining query input. An example of a
         * valid \c query parameter is <201, 242, 42>,
         * and \c iGenes is 21000 (if there are 21000 genes).
         */
        bool InitializeQuery(const vector <utype> &, const utype &);

        /*!
         * \brief Get the number of query partitions
         * \return The number of query partitions
         */
        utype GetNumFold() const;

        /*!
         * \brief Get the query genes as a vector
         * \return The query genes as a \c utype vector
         */
        const vector <utype> &GetQuery() const;

        /*!
         * \brief Get the query presence as a \c char vector
         * \return The query genes as a presence \c char vector
         */
        const vector<char> &GetQueryPresence() const;

        /*!
         * \brief
         * Get the \c i-th query partition
         *
         * \param i
         * The index of partition to return
         *
         * \return The vector of genes in \c i-th partition
         * \remarks
         * No bound checking on \c i.
         */
        const vector <utype> &GetCVQuery(utype &) const;

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
        bool CreateCVPartitions(const gsl_rng *, \
        const enum PartitionMode &, const utype= -1);

    private:
        vector <utype> queryGenes;
        vector<char> queryGenePresence;

        vector <utype> *crossValGenes;
        utype iNumFold;
        utype iFoldSize;
        utype qSize;

    };

}
#endif

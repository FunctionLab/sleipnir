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
#ifndef SEEKEVALUATE_H
#define SEEKEVALUATE_H

#include "seekbasic.h"
#include "seekmap.h"

namespace Sleipnir {

    template <typename K, typename V>
    struct PairedResult {
        K key;
        V val;

        PairedResult() {}

        PairedResult(K _key, V _val) : key(_key), val(_val) {}

        bool operator<(const PairedResult &item) const {
            if (val == item.val && key == item.key) {
                return false;
            }
            if (val < item.val) {
                return false;
            } else if (val > item.val) {
                return true;
            } else if (key < item.key) {
                return false;
            } else {
                return true;
            }
        }
    };

    typedef PairedResult<string, double> StrDoublePair;

    struct AResult {
        utype i;
        utype f;

        // Note: This is implemented in a confusing way. operartor< (less than)
        //  is supposed to return true when object_val < passed_in_val. This
        //  method is implemented in reverse to achieve descending sort.
        // TODO - Change this operator to the convential and achieve descending
        //   sort by passing a greater_that operator> into the sort function.
        //   Same with AResultFloat below
        bool operator<(const AResult &val) const {
            /*if(f<=val.f){
                return false;
            }else{
                return true;
            }*/
            // Sort specifies that when elem a == b must return false
            //  or an infinite loop can result. There are cases in
            //  some sort algorithms where an element is compared to itself.
            if (f == val.f && i == val.i) {
                return false;
            }
            if (f < val.f) {
                return false;
            } else if (f > val.f) {
                return true;
            } else if (i < val.i) {
                return false;
            } else {
                return true;
            }
        }
    };

    struct Ascending {
        bool operator()(const AResult &lx, const AResult &rx) const {
            if (lx.f < rx.f) {
                return true;
            } else if (lx.f > rx.f) {
                return false;
            } else if (lx.i < rx.i) {
                return true;
            } else {
                return false;
            }
            //return lx.f < rx.f;
        }
    };


    struct AResultFloat {
        utype i;
        float f;

        bool operator<(const AResultFloat &val) const {
            /*if(f<=val.f){
                return false;
            }else{
                return true;
            }*/
            if (f == val.f && i == val.i) {
                return false;
            }
            if (f < val.f) {
                return false;
            } else if (f > val.f) {
                return true;
            } else if (i < val.i) {
                return false;
            } else {
                return true;
            }
        }
    };

    struct AscendingFloat {
        bool operator()(const AResultFloat &lx, const AResultFloat &rx) const {
            if (lx.f < rx.f) {
                return true;
            } else if (lx.f > rx.f) {
                return false;
            } else if (lx.i < rx.i) {
                return true;
            } else {
                return false;
            }
            //return lx.f < rx.f;
        }
    };

/*!
 * \brief Evaluation metrics for a rank-list given some judgment gene-set
 *
 * Provide static utility functions for evaluating a ranking of genes with the user-given gold standard
 * gene-set. The typical use of such functions is in weighting datasets. Generally speaking, each dataset
 * is weighted by how well the query genes are able to retrieve each other in the dataset.
 * It is important to pick an informative measure to evaluate the retrieval of the query genes.
 * Seek provides the choice of two evaluation metrics: Rank-Biased Precision (RBP) or Average Precision.
 *
 */
    class CSeekPerformanceMeasure {
    public:
        /*!
         * \brief Sort a gene ranking by the gene score
         *
         * \param rank The vector of gene-scores to be sorted. Gene scores are inserted to this vector based
         * on their gene IDs, which are a value from 0 to the size of the vector.
         * \param mapG The gene presence map
         * \param a The output, which is a vector of (gene ID, gene score) pairs that are sorted by score
         * \param top If \c X, sort only the top \c X elements. If 0, then sort the entire vector.
         *
         * The struct \c AResult represents a (gene ID, gene score) pair. This function sorts the vector of
         * \c AResult in the descending order of the gene score.
         */
        static bool SortRankVector(const vector <utype> &rank,
                                   const CSeekIntIntMap &mapG, vector <AResult> &a,
                                   const bool bNegativeCor,
                                   const utype top = 0);

        /*!
         * \brief Calculate the rank-biased precision for a gene ranking
         *
         * \param rate The parameter \a p in the RBP formula
         * \param rbp The calculated RBP score, the output
         * \param mask The genes in the ranking to be skipped over (typically the query genes)
         * \param gold The gold-standard genes
         * \param mapG The gene presence map. Genes that are not present in the dataset are skipped over.
         * \param sing The sorted vector of (gene ID, gene score) pairs
         * \param rank The gene-score vector
         * \param top If \c X, sort only the top \c X elements. If 0, then sort the entire vector.
         *
         * First calls the CSeekPerformanceMeasure::SortRankVector() with the arguments \c rank and \c top,
         * in order to sort the gene-scores. Then with the sorted gene-ranking returned to \c sing, it calculates
         * the rank-biased precision.
         *
         * \remarks The RBP formula is given by:
         * \f[RBP=\sum_{g \in U}{(1-p)p^{rank(g)}}\f]
         * where \f$U\f$ is the gold standard gene-set, \f$p\f$ is the emphasis on ranks,
         * \f$rank(g)\f$ is the position of \f$g\f$ in the ranking
         * \f$p\f$ is typically set to 0.95 - 0.99. The recommended value is 0.99.
         * For more information, please read (Moffat et al 2008).
         */
        /* designed specifically for a CSeekDataset */
        /* mask: the query genes which are not included in RBP calcualtion */
        static bool RankBiasedPrecision(const float &rate,
                                        const vector <utype> &rank, float &rbp,
                                        const vector<char> &mask, const vector<char> &gold,
                                        const CSeekIntIntMap &mapG, vector <AResult> *sing,
                                        const bool bNegativeCor,
                /* optional */
                                        const utype top = 0);

        /*!
         * \brief Calculate the average precision for a gene ranking
         *
         * \param rank The gene-score vector
         * \param ap The calculated average precision
         * \param mask The genes in the ranking to be skipped over (typically the query genes)
         * \param gold The gold-standard genes
         * \param mapG The gene presence map. Genes that are not present in the dataset are skipped over.
         * \param ar The sorted vector of (gene ID, gene score) pairs
         */
        static bool AveragePrecision(
                const vector <utype> &rank, float &ap,
                const vector<char> &mask, const vector<char> &gold,
                const CSeekIntIntMap &mapG, vector <AResult> *ar, const bool bNegativeCor);

    };


}
#endif

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
#ifndef ORTHOLOGY_H
#define ORTHOLOGY_H

#include "orthologyi.h"

namespace Sleipnir {

/*!
 * \brief
 * An orthology is a collection of sets, each containing zero or more orthologous genes from any organism.
 * 
 * An orthology is a collection of sets referred to as orthologous clusters.  Each such cluster is a set of
 * genes; these genes can be drawn from any number of different organisms.  Genes in an orthologous cluster
 * are assumed (by some external agency) to perform similar functions in different organisms.
 * 
 * For example, for organisms O1, O2, and O3, an orthology might consist of clusters [O1:G1, O3:G2], [O1:G1,
 * O1:G2], and [O1:G3, O2:G1, O3:G1].  That is, any cluster might contain any subset of genes from any
 * subset of the available organisms, and any gene can appear in multiple clusters.
 * 
 * An orthology is stored in a tab-delimited text file in which each line represents a cluster and each
 * tab-separated token indicates a gene and the organism it is drawn from.  This file is of the form:
 * \code
 * O1:G1	O3:G2
 * O1:G1	O1:G2
 * O1:G3	O2:G1	O3:G1
 * \endcode
 * Here, there are three organisms O1, O2, O3; O1 possesses three genes G1, G2, and G3, O2 has one gene G1,
 * and O3 has two genes G1 and G2.  The organism identifiers should be short human-recognizable organism
 * identifiers (similar to those used by KEGG), and the gene identifiers should be standard unique
 * primary name.  A snippet of a real orthology file might resemble:
 * \code
 * DME|FBgn0034427	DME|FBgn0033431	MMU|MGI:104873	CEL|R04B3.2	HSA|AGA
 * MMU|MGI:87963	HSA|AGT
 * HSA|ALAD	SCE|YGL040C	MMU|MGI:96853	DME|FBgn0036271
 * DME|FBgn0036208	DME|FBgn0020764	HSA|GCAT	HSA|ALAS1	HSA|ALAS2	MMU|MGI:87989	CEL|T25B9.1	MMU|MGI:1349389	MMU|MGI:87990	SCE|YDR232W</pre>
 * \endcode
 * 
 * \remarks
 * A COrthology contains many CGenes sets, each containing CGene objects drawn from multiple underlying
 * CGenome objects, one per organism of interest.  The organism ID strings used in the orthology are
 * retained and mapped to the underlying genomes.
 */
    class COrthology : COrthologyImpl {
    public:
        bool Open(std::istream &istm);

        void Save(std::ostream &ostm) const;

        /*!
         * \brief
         * Return the number of clusters in the orthology.
         *
         * \returns
         * Number of clusters in the orthology.
         *
         * \see
         * GetCluster
         */
        size_t GetClusters() const {

            return m_vecvecpGenes.size();
        }

        /*!
         * \brief
         * Return the number of genes in the requested cluster.
         *
         * \param iCluster
         * Cluster from which genes should be returned.
         *
         * \returns
         * Number of genes in the given cluster.
         *
         * \remarks
         * For efficiency, no bounds checking is performed.  The given value must be smaller than GetClusters.
         */
        size_t GetGenes(size_t iCluster) const {

            return m_vecvecpGenes[iCluster].size();
        }

        /*!
         * \brief
         * Return the gene at the given index within the given cluster.
         *
         * \param iCluster
         * Cluster from which gene should be returned.
         *
         * \param iGene
         * Index of gene to retrieve from cluster.
         *
         * \returns
         * Gene at the given index within the given clusterl.
         *
         * \remarks
         * For efficiency, no bounds checking is performed.  The given values must be smaller than GetClusters
         * and GetGenes.
         */
        CGene &GetGene(size_t iCluster, size_t iGene) const {

            return *m_vecvecpGenes[iCluster][iGene];
        }

        /*!
         * \brief
         * Return the number of genomes (organisms) in the orthology.
         *
         * \returns
         * Number of genomes in the orthology.
         */
        size_t GetGenomes() const {

            return m_vecpGenomes.size();
        }

        /*!
         * \brief
         * Return the genome at the requested organism index.
         *
         * \param iOrganism
         * Index of genome to return.
         *
         * \returns
         * Genome at the requested index.
         *
         * \remarks
         * For efficiency, no bounds checking is performed.  The given value must be smaller than GetGenomes.
         */
        CGenome *GetGenome(size_t iOrganism) const {

            return m_vecpGenomes[iOrganism];
        }

        /*!
         * \brief
         * Return the human-readable string identifier of the requested organism index.
         *
         * \param iOrganism
         * Index of organism ID to return.
         *
         * \returns
         * Organism ID at the requested index.
         *
         * \remarks
         * For efficiency, no bounds checking is performed.  The given value must be smaller than GetGenomes.
         */
        const std::string &GetOrganism(size_t iOrganism) const {

            return m_vecstrOrganisms[iOrganism];
        }

        /*!
         * \brief
         * Return the index of the organism whose genome contains the given gene, or -1 if none exists.
         *
         * \param Gene
         * Gene whose organism index should be retrieved.
         *
         * \returns
         * Organism index associated with the given gene; -1 if none exists.
         *
         * \remarks
         * Gene comparison is done by pointer, so a different gene with the same primary identifier will
         * not match.  Gene/genome associations are maintained in a map by the orthology, so this
         * lookup is fast.
         */
        size_t GetOrganism(const CGene &Gene) const {
            TMapGeneI::const_iterator iterGenome;

            return (((iterGenome = m_mapGenes.find((CGene *) &Gene)) == m_mapGenes.end()) ? -1 :
                    iterGenome->second);
        }
    };

}

#endif // ORTHOLOGY_H

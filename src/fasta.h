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
#ifndef FASTA_H
#define FASTA_H

#include "fastai.h"

namespace Sleipnir {

/*!
 * \brief
 * Encapsulates a standard FASTA file or a modified ENCODE-style wiggle (WIG) file.
 * 
 * CFASTA performs efficient, disk-based indexing of large FASTA files; it also provides several
 * Sleipnir-specific extensions:
 * -# Sequences of different types associated with the same gene ID.  Unlike standard FASTA files, gene IDs
 *	assumed to contain no tab characters.  Instead, if an ID line contains tabs, the second tab-delimited
 *	column is interpreted as a sequence type.  This allows, for example, a gene sequence, upstream flank, and
 *	downstream flank to be associated in a separable manner with the same gene ID:
 * \code
 * > GENE
 * gene sequence here
 * > GENE	5
 * upstream flank here
 * > GENE	3
 * downstream flank here
 * \endcode
 * -# Introns and exons (or similar sequence subtypes) can be encoded within each sequence by alternating
 *	blocks of upper- and lower-case bases.  By default, upper-case bases are interpreted as exons and
 *	lower-case as introns, but there are no intrinsic semantics associated with the subtypes by CFASTA.  For
 *	example, the following "gene" would begin and end with exons separated by a single intron:
 * \code
 * > GENE
 * AACCGGTTacgtTTGGCCAA
 * \endcode
 * -# FASTA files akin to ENCODE wiggle (WIG) files are also supported, in which a single floating point value
 *	(instead of a nucleotide or amino acid) is associated with each position in a gene's sequence; each value
 *	must appear alone on separate lines.  This can be used to encode per-base conservation scores,
 *	nucleosome/TF occupancies, or other continuous data.  For example, the gene sequence above might be
 *	accompanied in a separate wiggle file by the scores:
 * \code
 * > GENE
 * 0.9
 * 0.1
 * 0.25
 * ...
 * \endcode
 * 
 * \remarks
 * Very large FASTA files are supported with relatively low memory usage by indexing the file on disk when
 * opened without loading any sequence data.  Sequence data is loaded on an as-needed basis by the Get
 * methods; this incurs a slight runtime penalty, but allows extremely large FASTA files to be accessed
 * efficiently with very low memory usage.
 */
    class CFASTA : CFASTAImpl {
    public:
        bool Open(const char *szFile, const std::set <std::string> &setstrTypes);

        void Save(std::ostream &ostm, size_t iWrap = 80) const;

        /*!
         * \brief
         * Retrieves all sequences (of any type) associated with the given gene index.
         *
         * \param iGene
         * Index of gene for which sequences are retrieved.
         *
         * \param vecsSequences
         * Zero or more output sequences associated with the given gene index.
         *
         * \returns
         * True if sequences were retrieved successfully; false otherwise.
         *
         * \remarks
         * Must be called only after a successful Open.  Performs zero or more seeks on disk to load the sequence
         * associated with the given gene index over all types in the FASTA file's index.
         *
         * \see
         * Open
         */
        bool Get(size_t iGene, std::vector <SFASTASequence> &vecsSequences) const {

            return CFASTAImpl::Get(iGene, &vecsSequences, NULL);
        }

        /*!
         * \brief
         * Retrieves all values (of any type) associated with the given gene index.
         *
         * \param iGene
         * Index of gene for which values are retrieved.
         *
         * \param vecsValues
         * Zero or more output values associated with the given gene index.
         *
         * \returns
         * True if values were retrieved successfully; false otherwise.
         *
         * \remarks
         * Must be called only after a successful Open.  Performs zero or more seeks on disk to load the values
         * associated with the given gene index over all types in the WIG file's index.
         *
         * \see
         * Open
         */
        bool Get(size_t iGene, std::vector <SFASTAWiggle> &vecsValues) const {

            return CFASTAImpl::Get(iGene, NULL, &vecsValues);
        }

        /*!
         * \brief
         * Opens a FASTA or WIG file and indexes the file without explicitly loading its contents.
         *
         * \param szFile
         * Path to FASTA/WIG file to open.
         *
         * \returns
         * True if file was loaded successfully; false otherwise.
         *
         * \remarks
         * Supports FASTA and WIG files as described in CFASTA.  No data is loaded on open, but an index is
         * created over all genes and types of interest; a file handle is held open, and data is loaded as needed
         * by the Get methods.
         *
         * \see
         * Save
         */
        bool Open(const char *szFile) {
            std::set <std::string> setstrDummy;

            return Open(szFile, setstrDummy);
        }

        /*!
         * \brief
         * Returns the total number of genes associated with this FASTA/WIG.
         *
         * \returns
         * Number of genes associated with this FASTA/WIG.
         */
        size_t GetGenes() const {

            return m_vecstrGenes.size();
        }

        /*!
         * \brief
         * Returns the gene ID associated with the given index.
         *
         * \param iGene
         * Gene index for which ID is returned.
         *
         * \returns
         * Gene ID associated with the given index.
         */
        const std::string &GetGene(size_t iGene) const {

            return CFASTAImpl::GetGene(iGene);
        }

        /*!
         * \brief
         * Returns the index of the given gene ID.
         *
         * \param strGene
         * Gene ID for which index is returned.
         *
         * \returns
         * Index of the given gene ID.
         */
        size_t GetGene(const std::string &strGene) const {
            TMapStrI::const_iterator iterGene;

            return (((iterGene = m_mapstriGenes.find(strGene)) == m_mapstriGenes.end()) ? -1 :
                    iterGene->second);
        }

        /*!
         * \brief
         * Returns the header line associated with the given gene index and sequence type.
         *
         * \param iGene
         * Gene index for which header is retrieved.
         *
         * \param strType
         * Sequence type for which header is retrieved.
         *
         * \returns
         * Header line associated with the given gene/type pair.
         *
         * \see
         * GetGene
         */
        const std::string GetHeader(size_t iGene, const std::string &strType) const {
            const TMapStrStr &mapstrstrHeaders = m_vecmapstrstrHeaders[iGene];
            TMapStrStr::const_iterator iterGene;

            return (((iterGene = mapstrstrHeaders.find(strType)) == mapstrstrHeaders.end()) ? "" :
                    iterGene->second);
        }

        /*!
         * \brief
         * Returns the set of sequence types indexed by this FASTA/WIG.
         *
         * \returns
         * Set of sequence types indexed by this FASTA/WIG.
         */
        const std::set <std::string> &GetTypes() const {

            return m_setstrTypes;
        }
    };

}

#endif // FASTA_H

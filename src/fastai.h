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
#ifndef FASTAI_H
#define FASTAI_H

#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <pthread.h>

#include "file.h"

namespace Sleipnir {

/*!
 * \brief
 * Base data associated with one entry in a FASTA/WIG file.
 */
struct SFASTABase {
	/*!
	 * \brief
	 * Type of the FASTA/WIG entry.
	 */
	std::string	m_strType;
};


/*!
 * \brief
 * Data associated with one sequence entry in a FASTA file.
 */
struct SFASTASequence : SFASTABase {
	/*!
	 * \brief
	 * True if the first subtype in the sequence is an intron; false otherwise.
	 */
	bool						m_fIntronFirst;
	/*!
	 * \brief
	 * Zero or more sequences of alternating subtypes within one entry in a FASTA file.
	 */
	std::vector<std::string>	m_vecstrSequences;
};

/*!
 * \brief
 * Data associated with one value entry in a WIG file.
 */
struct SFASTAWiggle : SFASTABase {
	/*!
	 * \brief
	 * Zero or more values associated with one entry in a WIG file.
	 */
	std::vector<float>	m_vecdValues;
};

class CFASTAImpl : public CFile {
protected:
	static const char	c_acComment[];
	static const char	c_acHeader[];

	typedef std::map<std::string, size_t>		TMapStrI;
	typedef std::map<std::string, std::string>	TMapStrStr;

	CFASTAImpl( );
	virtual ~CFASTAImpl( );

	bool Get( size_t, std::vector<SFASTASequence>*, std::vector<SFASTAWiggle>* ) const;
	bool Get( size_t, std::vector<SFASTASequence>&, size_t, const std::string&, SFASTASequence& ) const;
	bool Get( size_t, std::vector<SFASTAWiggle>&, size_t, SFASTAWiggle& ) const;

	const std::string& GetGene( size_t iGene ) const {

		return m_vecstrGenes[ iGene ]; }

	mutable std::ifstream		m_ifsm;
	TMapStrI					m_mapstriGenes;
	std::vector<std::string>	m_vecstrGenes;
	std::vector<TMapStrStr>		m_vecmapstrstrHeaders;
	std::vector<TMapStrI>		m_vecmapstriSequences;
	char*						m_szBuffer;
	std::set<std::string>		m_setstrTypes;
	mutable pthread_mutex_t		m_mutx;
};

}

#endif // FASTAI_H

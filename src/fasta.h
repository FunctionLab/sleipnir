#ifndef FASTA_H
#define FASTA_H

#include "fastai.h"

namespace Sleipnir {

struct SFASTASequence {
	std::string					m_strType;
	bool						m_fIntronFirst;
	std::vector<std::string>	m_vecstrSequences;
};

class CFASTA : CFASTAImpl {
public:
	bool Open( const char* szFile );
	void Save( std::ostream& ostm, size_t iWrap = 80 ) const;
	bool Get( size_t iGene, std::vector<SFASTASequence>& vecsSequences ) const;

	size_t GetGenes( ) const {

		return m_vecstrGenes.size( ); }

	const std::string& GetGene( size_t iGene ) const {

		return m_vecstrGenes[ iGene ]; }

	size_t GetGene( const std::string& strGene ) const {
		TMapStrI::const_iterator	iterGene;

		return ( ( ( iterGene = m_mapstriGenes.find( strGene ) ) == m_mapstriGenes.end( ) ) ? -1 :
			iterGene->second ); }

	const std::string GetHeader( size_t iGene, const std::string& strType ) const {
		const TMapStrStr&			mapstrstrHeaders	= m_vecmapstrstrHeaders[ iGene ];
		TMapStrStr::const_iterator	iterGene;

		return ( ( ( iterGene = mapstrstrHeaders.find( strType ) ) == mapstrstrHeaders.end( ) ) ? "" :
			iterGene->second ); }
};

}

#endif // FASTA_H

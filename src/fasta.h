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

class CFASTA : CFASTAImpl {
public:
	bool Open( const char* szFile, const std::set<std::string>& setstrTypes );
	void Save( std::ostream& ostm, size_t iWrap = 80 ) const;

	bool Get( size_t iGene, std::vector<SFASTASequence>& vecsSequences ) const {

		return CFASTAImpl::Get( iGene, &vecsSequences, NULL ); }

	bool Get( size_t iGene, std::vector<SFASTAWiggle>& vecsValues ) const {

		return CFASTAImpl::Get( iGene, NULL, &vecsValues ); }

	bool Open( const char* szFile ) {
		std::set<std::string>	setstrDummy;

		return Open( szFile, setstrDummy ); }

	size_t GetGenes( ) const {

		return m_vecstrGenes.size( ); }

	const std::string& GetGene( size_t iGene ) const {

		return CFASTAImpl::GetGene( iGene ); }

	size_t GetGene( const std::string& strGene ) const {
		TMapStrI::const_iterator	iterGene;

		return ( ( ( iterGene = m_mapstriGenes.find( strGene ) ) == m_mapstriGenes.end( ) ) ? -1 :
			iterGene->second ); }

	const std::string GetHeader( size_t iGene, const std::string& strType ) const {
		const TMapStrStr&			mapstrstrHeaders	= m_vecmapstrstrHeaders[ iGene ];
		TMapStrStr::const_iterator	iterGene;

		return ( ( ( iterGene = mapstrstrHeaders.find( strType ) ) == mapstrstrHeaders.end( ) ) ? "" :
			iterGene->second ); }

	const std::set<std::string>& GetTypes( ) const {

		return m_setstrTypes; }
};

}

#endif // FASTA_H

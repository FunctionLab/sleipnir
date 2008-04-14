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
#ifndef TRIEI_H
#define TRIEI_H

#include <vector>

namespace Sleipnir {

template<class tType> class CTrie;
template<class tType> class CTrieIteratorImpl;

template<class tType>
class CTrieImpl {
	friend class CTrieIteratorImpl<tType>;
protected:
	CTrieImpl( const tType& Undef ) : m_Undef(Undef), m_aabData(NULL) { }

	void Delete( unsigned char** aabData, size_t iDepth = 0 ) const {
		unsigned char	i;

		if( !aabData || ( iDepth >= m_vecbSizes.size( ) ) )
			return;

		for( i = 0; i < m_vecbSizes[ iDepth ]; ++i )
			Delete( (unsigned char**)aabData[ i ], iDepth + 1 );
		delete[] aabData; }

	std::vector<unsigned char>	m_vecbSizes;
	unsigned char**				m_aabData;
	const tType&				m_Undef;
};

template<class tType>
class CTrieIteratorImpl {
protected:
	CTrieIteratorImpl( const CTrie<tType>& Trie ) : m_Trie(Trie) {
		size_t				iDepth;
		unsigned char		b;
		unsigned char***	paabData;

		m_vecbPosition.resize( m_Trie.GetSize( ) );
		m_vecpaaPosition.resize( m_Trie.GetSize( ) );
		if( !m_Trie.m_aabData ) {
			for( iDepth = 0; iDepth < m_vecpaaPosition.size( ); ++iDepth )
				m_vecpaaPosition[ iDepth ] = NULL;
			return; }

		paabData = (unsigned char***)&m_Trie.m_aabData;
		for( iDepth = 0; ( iDepth + 1 ) < m_vecbPosition.size( ); ++iDepth )
			for( b = 0; b < m_Trie.GetSize( iDepth ); ++b )
				if( (*paabData)[ b ] ) {
					m_vecbPosition[ iDepth ] = b;
					m_vecpaaPosition[ iDepth ] = paabData;
					paabData = (unsigned char***)&(*paabData)[ b ];
					break; }
		for( b = 0; b < m_Trie.GetSize( iDepth ); ++b )
			if( *(tType*)&(*paabData)[ b ] != m_Trie.m_Undef ) {
				m_vecbPosition[ iDepth ] = b;
				m_vecpaaPosition[ iDepth ] = paabData;
				break; } }

	CTrieIteratorImpl( const CTrieIteratorImpl& Iter ) : m_Trie(Iter.m_Trie) {

		m_vecbPosition.resize( Iter.m_vecbPosition.size( ) );
		std::copy( Iter.m_vecbPosition.begin( ), Iter.m_vecbPosition.end( ), m_vecbPosition.begin( ) );
		m_vecpaaPosition.resize( Iter.m_vecpaaPosition.size( ) );
		std::copy( Iter.m_vecpaaPosition.begin( ), Iter.m_vecpaaPosition.end( ), m_vecpaaPosition.begin( ) ); }

	tType* GetValue( ) const {

		return (tType*)&(*m_vecpaaPosition[ m_vecpaaPosition.size( ) - 1 ])[
			m_vecbPosition[ m_vecbPosition.size( ) - 1 ] ]; }

	bool NextDown( size_t iDepth, bool fReset ) {
		unsigned char	b;

		for( ; ( iDepth + 1 ) < m_vecbPosition.size( ); ++iDepth ) {
			for( b = m_vecbPosition[ iDepth ]; b < m_Trie.GetSize( iDepth ); ++b )
				if( (*m_vecpaaPosition[ iDepth ])[ b ] )
					break;
			if( b >= m_Trie.GetSize( iDepth ) )
				return false;
			m_vecbPosition[ iDepth ] = b;
			m_vecbPosition[ iDepth + 1 ] = 0;
			m_vecpaaPosition[ iDepth + 1 ] = (unsigned char***)&(*m_vecpaaPosition[ iDepth ])[ b ]; }
		for( b = ( fReset ? 0 : ( m_vecbPosition[ iDepth ] + 1 ) ); b < m_Trie.GetSize( iDepth ); ++b )
			if( *(tType*)&(*m_vecpaaPosition[ iDepth ])[ b ] != m_Trie.m_Undef ) {
				m_vecbPosition[ iDepth ] = b;
				return true; }

		return false; }

	size_t NextUp( ) {
		size_t			iDepth;
		unsigned char	b;

		for( iDepth = m_vecbPosition.size( ) - 2; iDepth != -1; --iDepth )
			for( b = ( m_vecbPosition[ iDepth ] + 1 ); b < m_Trie.GetSize( iDepth ); ++b )
				if( (*m_vecpaaPosition[ iDepth ])[ b ] ) {
					m_vecbPosition[ iDepth ] = b;
					return iDepth; }

		m_vecpaaPosition[ m_vecpaaPosition.size( ) - 1 ] = NULL;
		return -1; }

	const CTrie<tType>&				m_Trie;
	std::vector<unsigned char>		m_vecbPosition;
	std::vector<unsigned char***>	m_vecpaaPosition;
};

}

#endif // TRIEI_H

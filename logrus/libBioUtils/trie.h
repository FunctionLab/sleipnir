#ifndef TRIE_H
#define TRIE_H

#include <vector>

#include "triei.h"

namespace libBioUtils {

template<class tType> class CTrieIterator;

template<class tType>
class CTrie : public CTrieImpl<tType> {
public:
	typedef CTrieIterator<tType>	iterator;

	CTrie( const std::vector<unsigned char>& vecbSizes, const tType& Undef ) : CTrieImpl( Undef ) {

		m_vecbSizes.resize( vecbSizes.size( ) );
		std::copy( vecbSizes.begin( ), vecbSizes.end( ), m_vecbSizes.begin( ) ); }

	CTrie( const IDataset* pData, const tType& Undef, bool fAnswers = true ) : CTrieImpl( Undef ) {
		size_t					i, j, k, iExp;
		vector<unsigned char>	vecbDatum;

		for( i = j = 0; i < pData->GetExperiments( ); ++i )
			if( !pData->IsHidden( i ) )
				j++;
		m_vecbSizes.resize( j );
		for( i = j = 0; i < pData->GetExperiments( ); ++i )
			if( !pData->IsHidden( i ) )
				m_vecbSizes[ j++ ] = (unsigned char)pData->GetBins( i ) + 1;

		vecbDatum.resize( GetSize( ) );
		for( i = 0; i < pData->GetGenes( ); ++i )
			for( j = ( i + 1 ); j < pData->GetGenes( ); ++j ) {
				if( !pData->IsExample( i, j ) || ( fAnswers && ( pData->GetDiscrete( i, j, 0 ) == -1 ) ) )
					continue;
				for( iExp = k = 0; k < pData->GetExperiments( ); ++k )
					if( !pData->IsHidden( k ) )
						vecbDatum[ iExp++ ] = (unsigned char)( pData->GetDiscrete( i, j, k ) + 1 );
				Set( vecbDatum )++; } }

	~CTrie( ) {

		Delete( m_aabData ); }

	tType& Set( const std::vector<unsigned char>& vecbData ) {
		size_t				iDepth;
		unsigned char		b;
		unsigned char***	paabData;

		for( iDepth = 0,paabData = &m_aabData; iDepth < m_vecbSizes.size( );
			paabData = (unsigned char***)&(*paabData)[ vecbData[ iDepth++ ] ] )
			if( !*paabData ) {
				*paabData = new unsigned char*[ m_vecbSizes[ iDepth ] ];
				for( b = 0; b < m_vecbSizes[ iDepth ]; ++b )
					if( ( iDepth + 1 ) == m_vecbSizes.size( ) )
						*(tType*)&(*paabData)[ b ] = m_Undef;
					else
						(*paabData)[ b ] = NULL; }

		return *(tType*)paabData; }

	const tType Get( const std::vector<unsigned char>& vecbData ) const {
		size_t			iDepth;
		unsigned char**	aabData;

		for( iDepth = 0,aabData = m_aabData; iDepth < m_vecbSizes.size( );
			aabData = (unsigned char**)aabData[ vecbData[ iDepth++ ] ] )
			if( !aabData )
				return m_Undef;

		return (tType)aabData; }

	size_t GetSize( ) const {

		return m_vecbSizes.size( ); }

	unsigned char GetSize( size_t iDepth ) const {

		return m_vecbSizes[ iDepth ]; }
};

template<class tType>
class CTrieIterator : protected CTrieIteratorImpl<tType> {
public:
	CTrieIterator( const CTrieIterator& Iter ) : CTrieIteratorImpl( Iter ) { }

	CTrieIterator( const CTrie<tType>& Trie ) : CTrieIteratorImpl( Trie ) { }

	bool IsDone( ) const {

		return ( m_vecpaaPosition.empty( ) || !m_vecpaaPosition[ m_vecpaaPosition.size( ) - 1 ] ); }

	const std::vector<unsigned char>& GetPosition( ) const {

		return m_vecbPosition; }

	const tType& Get( ) const {

		return Set( ); }

	tType& Set( ) const {

		return *GetValue( ); }

	bool Next( ) {
		size_t	iDepth;
		bool	fReset;

		if( IsDone( ) )
			return false;

		fReset = false;
		iDepth = m_vecbPosition.size( ) - 1;
		while( true ) {
			if( NextDown( iDepth, fReset ) )
				return true;
			fReset = true;
			if( ( iDepth = NextUp( ) ) == -1 )
				return false; } }
};

}

#endif // TRIE_H

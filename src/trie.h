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
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, and Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#ifndef TRIE_H
#define TRIE_H

#include "triei.h"

namespace Sleipnir {

template<class tType> class CTrieIterator;

/*!
 * \brief
 * A simple prefix tree implementation.
 * 
 * \param tType
 * Type of element contained by the trie.
 * 
 * A trie, or prefix tree, is an n-ary tree which acts as a key to value map.  Keys consist of a sequence
 * of elements (usually characters making up a string) stored efficiently in overlapping tree nodes.  This
 * trie implementation can store arbitrary objects as values and use arbitrary byte strings as keys.  At
 * the simplest, think of a trie like a hash or map dictionary that might save in memory usage and lookup
 * time for highly structured keys.
 * 
 * \remarks
 * Lookup time is linear in the length of the key, and memory usage is linear in the values and expected
 * logarithmic in the keys (worst case linear).  However, due to the intended usage of this trie
 * implementation, it is optimized for densely overlapping keys; the memory usage constant will be large
 * for sparse tries.
 * 
 * \see
 * CTrieIterator
 */
template<class tType>
class CTrie : public CTrieImpl<tType> {
public:
	/*!
	 * \brief
	 * Iterator type for trie traversal.
	 */
	typedef CTrieIterator<tType>	iterator;

	/*!
	 * \brief
	 * Construct a new trie with the specified branching factors and default value.
	 * 
	 * \param vecbSizes
	 * Branching factor (maximum number of different values) at each depth within the trie (i.e. at each
	 * position within the keys).
	 * 
	 * \param Default
	 * Default value to be returned when trie lookup fails for a given key.
	 * 
	 * \remarks
	 * The given sizes indicate the maximum value at each key position; for example, if a trie is initialized
	 * with sizes [2, 1, 3], then keys such as [0, 0, 0], [1, 0, 2], or [0, 0, 1] are all valid, but
	 * [2, 1, 3] or [2, 0, 0] are not.  This is equivalent to the maximum branching factor at each level
	 * of the trie, and it also dictates the maximum depth of the trie (i.e. length of a key).
	 */
	CTrie( const std::vector<unsigned char>& vecbSizes, const tType& Default ) : CTrieImpl<tType>( Default ) {

		m_vecbSizes.resize( vecbSizes.size( ) );
		std::copy( vecbSizes.begin( ), vecbSizes.end( ), m_vecbSizes.begin( ) ); }

	/*!
	 * \brief
	 * Construct a new trie encoding the data from the given dataset.
	 * 
	 * \param pData
	 * Dataset to be encoded as a trie.
	 * 
	 * \param Default
	 * Default value to be returned when trie lookup fails for a given key.
	 * 
	 * \param fAnswers
	 * If true, the first (0th) data file in the dataset represents a gold standard.
	 * 
	 * Constructs a trie encoding the given dataset.  Each key in the trie consists of gene pair values from
	 * each data file within the dataset, and the corresponding trie value is a count of how many times that
	 * combination of values occurs in the dataset.  If fAnswers is true, only gene pairs with a value
	 * present in the first (0th) data file are included in the trie.  This can provide a fairly compact
	 * way to store a highly redundant dataset, but it becomes inefficient rapidly if the dataset is
	 * sparse.
	 * 
	 * For example, suppose a dataset contains an answer file and two data files and spans three genes A,
	 * B, and C.  The answer file is binary (contains only 0, 1, and missing values), the first data file
	 * is binary, and the second data file can take three values.  The values present for each gene pair are:
	 * \code
	 * A	B	-1	0	-1
	 * A	C	0	0	0
	 * B	C	1	-1	2
	 * \endcode
	 */
	CTrie( const IDataset* pData, const tType& Default, bool fAnswers = true ) : CTrieImpl<tType>( Undef ) {
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

	/*!
	 * \brief
	 * Return a writable reference to the given key's value.
	 * 
	 * \param vecbData
	 * Key to which value should be written.
	 * 
	 * \returns
	 * Writable reference to the given key's trie value.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given key's values should all be smaller than
	 * the trie's branching factors at the appropriate levels.
	 * 
	 * \see
	 * Get
	 */
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

	/*!
	 * \brief
	 * Retrieve the value at the given key, or a default if none exists.
	 * 
	 * \param vecbData
	 * Key whose value should be retrieved.
	 * 
	 * \returns
	 * The value corresponding to the given key, or the trie's default value if none exists.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given key's values should all be smaller than
	 * the trie's branching factors at the appropriate levels.
	 * 
	 * \see
	 * Set
	 */
	const tType Get( const std::vector<unsigned char>& vecbData ) const {
		size_t			iDepth;
		unsigned char**	aabData;

		for( iDepth = 0,aabData = m_aabData; iDepth < m_vecbSizes.size( );
			aabData = (unsigned char**)aabData[ vecbData[ iDepth++ ] ] )
			if( !aabData )
				return m_Undef;

		return (tType)aabData; }

	/*!
	 * \brief
	 * Return maximum depth of the trie.
	 * 
	 * \returns
	 * Maximum depth of the trie.
	 */
	size_t GetSize( ) const {

		return m_vecbSizes.size( ); }

	/*!
	 * \brief
	 * Return branching factor of the requested trie level.
	 * 
	 * \param iDepth
	 * Trie level for which branching factor should be returned.
	 * 
	 * \returns
	 * Maximum branching factor (key value) for the requested depth.
	 * 
	 * \remarks
	 * For efficiency, no bounds checking is performed.  The given depth must be smaller than GetSize.
	 */
	unsigned char GetSize( size_t iDepth ) const {

		return m_vecbSizes[ iDepth ]; }
};

/*!
 * \brief
 * Iterator for inorder traversal of trie keys.
 * 
 * \param tType
 * Type of element contained by the trie.
 * 
 * A trie iterator provides inorder access to each key in the trie.  Given a trie, a new iterator will
 * subsequently provide read/write access to the value for each key in the dictionary.  Sample usage
 * might be:
 * \code
 * for( TTrieData::iterator SomeIterator( SomeTrie ); !SomeIterator.IsDone( ); SomeIterator.Next( ) ) {
 *   cout << SomeIterator.Get( ) << endl; }
 * \endcode
 * 
 * \see
 * CTrie
 */
template<class tType>
class CTrieIterator : protected CTrieIteratorImpl<tType> {
public:
	/*!
	 * \brief
	 * Copy constructor.
	 * 
	 * \param Iter
	 * Iterator whose values should be copied to the new iterator.
	 * 
	 * \remarks
	 * The newly constructed iterator will begin traversal in the same position as the given iterator.
	 */
	CTrieIterator( const CTrieIterator& Iter ) : CTrieIteratorImpl( Iter ) { }

	/*!
	 * \brief
	 * Construct a new iterator for the given trie.
	 * 
	 * \param Trie
	 * Trie whose key/value pairs should be traversed.
	 */
	CTrieIterator( const CTrie<tType>& Trie ) : CTrieIteratorImpl( Trie ) { }

	/*!
	 * \brief
	 * Return true if the iterator has completed its traversal.
	 * 
	 * \returns
	 * True if there are no more valid key/value pairs to be iterated.
	 * 
	 * \remarks
	 * Note that if IsDone returns true, Get, Set, and GetPosition are invalid.
	 */
	bool IsDone( ) const {

		return ( m_vecpaaPosition.empty( ) || !m_vecpaaPosition[ m_vecpaaPosition.size( ) - 1 ] ); }

	/*!
	 * \brief
	 * Return the key at the current iterator position.
	 * 
	 * \returns
	 * Key at the current iterator position.
	 */
	const std::vector<unsigned char>& GetPosition( ) const {

		return m_vecbPosition; }

	/*!
	 * \brief
	 * Return the value at the current iterator position.
	 * 
	 * \returns
	 * Value at the current iterator position.
	 * 
	 * \see
	 * Set
	 */
	const tType& Get( ) const {

		return Set( ); }

	/*!
	 * \brief
	 * Return a writable value at the current iterator position.
	 * 
	 * \returns
	 * Writable value at the current iterator position.
	 * 
	 * \see
	 * Get
	 */
	tType& Set( ) const {

		return *GetValue( ); }

	/*!
	 * \brief
	 * Advance the iterator to the next key/value position.
	 * 
	 * \returns
	 * True if the iterator was advanced; false if it was already at the end of the trie.
	 * 
	 * \remarks
	 * If IsDone returns true, Next will return false.
	 */
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

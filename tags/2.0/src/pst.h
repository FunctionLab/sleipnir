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
#ifndef PST_H
#define PST_H

#include "psti.h"
#include "meta.h"

namespace Sleipnir {

/*!
 * \brief
 * Represents a probabilistic suffix tree (PST) containing zero or more overlapping strings in a weighted
 * manner.
 * 
 * A probabilistic suffix tree (PST) represents a set of zero or more weighted strings in an efficient
 * manner.  For example, a PST might contain the strings:
 * \code
 * 2:AACCGG
 * 1:AACCT
 * 1:AACT
 * \endcode
 * indicating that the set contains <tt>AACCGG</tt> with probability 0.5, <tt>AACCT</tt> with probability
 * 0.25, and <tt>AACT</tt> with probability 0.25.  A PST encodes such strings as an overlapping tree of
 * characters, which is thus compact in memory and can be compared efficiently against an input string to
 * count the matches of any string in the tree.
 * 
 * \remarks
 * CPST is a relatively naive and inflexible (but efficient) implementation of PSTs that should only be used
 * for fixed alphabets containing a small number of characters.  It is currently intended for use mainly with
 * CCoalesce.
 * 
 * \see
 * CCoalesceMotifLibrary
 */
class CPST : protected CPSTImpl {
public:
	/*!
	 * \brief
	 * Initializes a new PST over an alphabet of the given size.
	 * 
	 * \param iArity
	 * Maximum number of unique characters in strings contained in the new PST.
	 */
	CPST( size_t iArity ) : CPSTImpl(iArity) { }

	void RemoveRCs( float dPenaltyGap, float dPenaltyMismatch, CPST& PSTOut ) const;

	/*!
	 * \brief
	 * Add the given string and strings from the given PST to the current PST, with an optional character
	 * offset between the two.
	 * 
	 * \param strSequence
	 * Fixed string to be added to the current PST.
	 * 
	 * \param PST
	 * PST containing strings to be added to the current PST.
	 * 
	 * \param iOffset
	 * Optional character offset between the given string and PST; if positive, the PST is added first, and if
	 * negative, the fixed string.
	 */
	void Add( const std::string& strSequence, const CPST& PST, int iOffset ) {

		if( iOffset < 0 ) {
			Add( strSequence );
			Add( PST, -iOffset ); }
		else {
			Add( PST );
			Add( strSequence, iOffset ); } }

	/*!
	 * \brief
	 * Add the given strings to the current PST, with an optional character offset between them.
	 * 
	 * \param strOne
	 * First string to be added to the current PST.
	 * 
	 * \param strTwo
	 * Second string to be added to the current PST.
	 * 
	 * \param iOffset
	 * Optional character offset between the given strings; if positive, strOne is added first, and if
	 * negative, strTwo.
	 */
	void Add( const std::string& strOne, const std::string& strTwo, int iOffset ) {
		size_t	i, j;

		if( iOffset < 0 )
			return Add( strTwo, strOne, -iOffset );
		i = CPSTImpl::Add( strOne, 0, m_sRoot );
		j = CPSTImpl::Add( strTwo, iOffset, m_sRoot );
		if( ( i = max( i, j ) ) > m_iDepth )
			m_iDepth = i; }

	/*!
	 * \brief
	 * Add the given string to the current PST, with an optional character offset relative to the PST's root.
	 * 
	 * \param strSequence
	 * Fixed string to be added to the current PST.
	 * 
	 * \param iOffset
	 * Optional character offset between the given string; if nonzero, iOffset characters along each branch of
	 * the current PST are skipped before strSequence is inserted.
	 */
	void Add( const std::string& strSequence, size_t iOffset = 0 ) {
		size_t	i;

		if( ( i = CPSTImpl::Add( strSequence, iOffset, m_sRoot ) ) > GetDepth( ) )
			m_iDepth = i; }

	/*!
	 * \brief
	 * Add the strings from the given PST to the current PST, with an optional character offset relative to
	 * the current PST's root.
	 * 
	 * \param PST
	 * PST containng strings to be added to the current PST.
	 * 
	 * \param iOffset
	 * Optional character offset applied to the given PST; if nonzero, iOffset characters along each branch of
	 * the current PST are skipped before strings from the given PST are added.
	 */
	void Add( const CPST& PST, size_t iOffset = 0 ) {
		size_t	i;

		if( ( i = CPSTImpl::Add( PST.m_sRoot, iOffset, m_sRoot ) ) > GetDepth( ) )
			m_iDepth = i; }

	/*!
	 * \brief
	 * Returns the match score of the PST against the given string, with an optional offset.
	 * 
	 * \param strTarget
	 * String against which PST is matched.
	 * 
	 * \param iOffset
	 * If nonzero, the number of characters to skip within the target string before scoring the match.
	 * 
	 * \returns
	 * Length-normalized probabilistic match score.
	 * 
	 * \remarks
	 * Score is equal to zero if no strings in the PST match; otherwise, it is the maximum probability over
	 * all matching strings normalized by the probability of matching by chance, (1/arity)^length.
	 */
	float GetMatch( const std::string& strTarget, size_t iOffset = 0 ) const {
		float	dPMatch;
		size_t	iMatched;

		iMatched = 0;
		return ( ( ( dPMatch = CPSTImpl::GetMatch( strTarget, m_sRoot, iOffset, iMatched ) ) && iMatched ) ?
			( dPMatch * pow( 1.0f / m_iArity, (float)( GetDepth( ) - iMatched ) ) ) : 0 ); }

	/*!
	 * \brief
	 * Returns the maximum length over all strings in the PST.
	 * 
	 * \returns
	 * Maximum length over all strings in the PST.
	 */
	size_t GetDepth( ) const {

		return m_iDepth; }

	/*!
	 * \brief
	 * Returns a string representation of the PST.
	 * 
	 * \returns
	 * String representation of the PST.
	 * 
	 * A PST is represented using pipe characters (|) to indicate choices, with parentheses indicating
	 * grouping and weights indicated with preceding counts delimited by colons (:).  For example, a PST
	 * containing the strings:
	 * \code
	 * 2:AACCGG
	 * 1:AACCT
	 * 1:AACT
	 * \endcode
	 * would be represented by the motif:
	 * \code
	 * AAC(3:C(2:GG)|(1:T))|(1:T)
	 * \endcode
	 */
	std::string GetMotif( ) const {

		return CPSTImpl::GetMotif( m_sRoot ); }

	/*!
	 * \brief
	 * Returns the best alignment score of the given string against all strings in the PST, using the
	 * requested penalties, maximum score, and offset.
	 * 
	 * \param strSequence
	 * String to be aligned with the current PST.
	 * 
	 * \param dPenaltyGap
	 * Alignment score penalty for gaps.
	 * 
	 * \param dPenaltyMismatch
	 * Alignment score penalty for mismatches.
	 * 
	 * \param dCutoff
	 * Maximum score before an alignment is considered inviable; used for optimization.
	 * 
	 * \param iOffset
	 * If nonzero, the number of characters to skip within the input string before scoring the alignment.
	 * 
	 * \returns
	 * Minimum alignment score of the given string with strings in the current PST.
	 * 
	 * \remarks
	 * Alignments are internally ungapped, so gap penalties are only incurred at the ends (i.e. by overhangs).
	 */
	float Align( const std::string& strSequence, float dPenaltyGap, float dPenaltyMismatch, float dCutoff,
		int& iOffset ) const {
		float	dRet;

		dRet = CPSTImpl::Align( m_sRoot, GetDepth( ), strSequence, strSequence.length( ), dPenaltyGap,
			dPenaltyMismatch, dCutoff, iOffset );
		iOffset *= -1;
		return dRet; }

	/*!
	 * \brief
	 * Returns the best alignment score of strings in the given PST against strings in the current PST, using
	 * the requested penalties, maximum score, and offset.
	 * 
	 * \param PST
	 * PST to be aligned with the current PST.
	 * 
	 * \param dPenaltyGap
	 * Alignment score penalty for gaps.
	 * 
	 * \param dPenaltyMismatch
	 * Alignment score penalty for mismatches.
	 * 
	 * \param dCutoff
	 * Maximum score before an alignment is considered inviable; used for optimization.
	 * 
	 * \param iOffset
	 * If nonzero, the number of characters to skip within the input PST before scoring the alignment.
	 * 
	 * \returns
	 * Minimum alignment score of the given PST with strings in the current PST.
	 * 
	 * \remarks
	 * Alignments are internally ungapped, so gap penalties are only incurred at the ends (i.e. by overhangs).
	 */
	float Align( const CPST& PST, float dPenaltyGap, float dPenaltyMismatch, float dCutoff,
		int& iOffset ) const {

		return CPSTImpl::Align( m_sRoot, GetDepth( ), PST.m_sRoot, PST.GetDepth( ), dPenaltyGap,
			dPenaltyMismatch, dCutoff, iOffset ); }

	/*!
	 * \brief
	 * Adds the given PST encoded in string form to the current PST.
	 * 
	 * \param strPST
	 * PST encoded in string form, or a string from the current PST's alphabet.
	 * 
	 * \returns
	 * True if the given PST string was added successfully.
	 * 
	 * \see
	 * GetMotif
	 */
	bool Open( const std::string& strPST ) {

		return ( ( m_iDepth = CPSTImpl::Open( strPST, m_sRoot ) ) != -1 ); }

	/*!
	 * \brief
	 * Retrieves a PWM approximately equivalent to the current PST.
	 * 
	 * \param MatPWM
	 * Output PWM approximating the current PST.
	 * 
	 * \param szSymbols
	 * String of symbols mapped to the indices of the current PST.
	 * 
	 * \remarks
	 * szSymbols must be of the same length as the current PST's arity; it is typically "ACGT".
	 * GetPWM is most useful if the PST first has its reverse complements removed.  Technically
	 * returns a PSSM or PFM of counts, not a PWM of continuous probabilities.
	 * 
	 * \see
	 * RemoveRCs
	 */
	bool GetPWM( CFullMatrix<uint16_t>& MatPWM, const char* szSymbols ) const {
		std::map<unsigned char, size_t>					mapciChars;
		std::vector<size_t>								veciOrder;
		std::map<unsigned char, size_t>::const_iterator	iterChar;
		size_t											i, j;

		if( strlen( szSymbols ) != m_iArity )
			return false;

		MatPWM.Initialize( m_iArity, GetDepth( ) );
		MatPWM.Clear( );
		if( !CPSTImpl::GetPWM( m_sRoot, 0, mapciChars, MatPWM ) )
			return false;

		veciOrder.resize( MatPWM.GetRows( ) );
		for( i = j = 0; i < m_iArity; ++i )
			veciOrder[ i ] = ( ( iterChar = mapciChars.find( szSymbols[ i ] ) ) == mapciChars.end( ) ) ?
				( mapciChars.size( ) + j++ ) : iterChar->second;
		CMeta::Permute( MatPWM.Get( ), veciOrder );
		return true; }

	/*!
	 * \brief
	 * Counts the number of possible discrete strings encoded by the PST.
	 * 
	 * \returns
	 * Number of distinct strings encoded by the PST.
	 */
	size_t Integrate( ) const {
		size_t	iRet;

		CPSTImpl::Integrate( m_sRoot, iRet = 0 );
		return iRet; }

	/*!
	 * \brief
	 * Simplifies the current PST by removing all subtrees occurring at relatively low frequency.
	 * 
	 * \returns
	 * True if simplification succeeded, false otherwise.
	 * 
	 * \remarks
	 * Removes all subtrees with frequency below 1/arity.
	 */
	bool Simplify( ) {

		return CPSTImpl::Simplify( 1.0f / m_iArity, m_sRoot ); }
};

}

#endif // PST_H

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
#ifndef COALESCEMOTIFS_H
#define COALESCEMOTIFS_H

#include "coalescemotifsi.h"

namespace Sleipnir {

struct SCoalesceModifierCache;

/*!
 * \brief
 * Manages a set of kmer, reverse complement, and probabilistic suffix tree motifs for CCoalesce.
 * 
 * A motif library for some small integer k consists of all kmers of length k, all reverse complement pairs
 * (RCs) over those kmers, and zero or more probabilistic suffix trees (PSTs) formed at runtime by merging
 * kmers, RCs, and other PSTs.  Each motif is represented by an atomic ID, which can be converted to and from
 * a string representation by the library or matched against any given string.  The library also provides
 * services for merging kmers/RCs/PSTs by non-gapped edit distance comparisons.  All such motifs are
 * candidates for (under)enrichment in clusters found by CCoalesce.
 * 
 * \remarks
 * PSTs in the library can be of any depth, regardless of k; merging overlapping kmers, for example, may form
 * a valid PST of length greater than k.
 * 
 * \see
 * CCoalesce
 */
class CCoalesceMotifLibrary : CCoalesceMotifLibraryImpl {
public:
	static std::string GetReverseComplement( const std::string& strKMer ) {

		return CCoalesceMotifLibraryImpl::GetReverseComplement( strKMer ); }

	/*!
	 * \brief
	 * Initializes a new motif library based on kmers of the given length.
	 * 
	 * \param iK
	 * Length of kmers underlying the motif library.
	 */
	CCoalesceMotifLibrary( size_t iK ) : CCoalesceMotifLibraryImpl( iK ) { }

	float GetMatch( const std::string& strSequence, uint32_t iMotif, size_t iOffset,
		SCoalesceModifierCache& sModifiers ) const;
	uint32_t Open( const std::string& strMotif );
	std::string GetPWM( uint32_t iMotif, float dCutoffPWMs, float dPenaltyGap, float dPenaltyMismatch,
		bool fNoRCs ) const;

	/*!
	 * \brief
	 * Returns the string representation of the given motif ID.
	 * 
	 * \param iMotif
	 * Motif ID to be returned as a string.
	 * 
	 * \returns
	 * String representation of the requested motif.
	 * 
	 * \remarks
	 * Kmers are represented as strings of k characters; RCs are represented as two kmers delimited by a
	 * pipe (|); PSTs are represented as described in CPST.
	 * 
	 * \see
	 * CPST
	 */
	std::string GetMotif( uint32_t iMotif ) const {

		return CCoalesceMotifLibraryImpl::GetMotif( iMotif ); }

	/*!
	 * \brief
	 * Returns a motif ID representing the merger of the two input motifs, which can be of any type.
	 * 
	 * \param iOne
	 * ID of first motif to be merged.
	 * 
	 * \param iTwo
	 * ID of second motif to be merged.
	 * 
	 * \param dCutoff
	 * Maximum edit distance threshhold for successful merging.
	 * 
	 * \returns
	 * -1 if the two motifs cannot be merged or have already been merged; the ID of the merged motif
	 * otherwise, which will always be a PST.
	 * 
	 * \remarks
	 * If the two input motifs are successfully merged, the resulting motif will always be a newly created
	 * PST to be managed by the library.  Minimum edit distances between the two input motifs are calculated
	 * using standard ungapped alignments and the current scoring penalties, with the minimum of all possible
	 * alignments used to score e.g. two input PSTs.
	 * 
	 * \see
	 * SetPenaltyGap | SetPenaltyMismatch
	 */
	uint32_t Merge( uint32_t iOne, uint32_t iTwo, float dCutoff ) {
		std::pair<uint32_t, uint32_t>	priiMerged;

		priiMerged.first = min( iOne, iTwo );
		priiMerged.second = max( iOne, iTwo );
		if( m_setpriiMerged.find( priiMerged ) != m_setpriiMerged.end( ) )
			return -1;
		m_setpriiMerged.insert( priiMerged );

		switch( GetType( iOne ) ) {
			case ETypeRC:
				switch( GetType( iTwo ) ) {
					case ETypeKMer:
						return MergeKMerRC( GetMotif( iTwo ), iOne, dCutoff );

					case ETypeRC:
						return MergeRCs( iOne, iTwo, dCutoff );

					case ETypePST:
						return MergeRCPST( iOne, *GetPST( iTwo ), dCutoff ); }

			case ETypePST:
				switch( GetType( iTwo ) ) {
					case ETypeKMer:
						return MergeKMerPST( GetMotif( iTwo ), *GetPST( iOne ), dCutoff );

					case ETypeRC:
						return MergeRCPST( iTwo, *GetPST( iOne ), dCutoff );

					case ETypePST:
						return MergePSTs( *GetPST( iOne ), *GetPST( iTwo ), dCutoff ); } }

		switch( GetType( iTwo ) ) {
			case ETypeRC:
				return MergeKMerRC( GetMotif( iOne ), iTwo, dCutoff );

			case ETypePST:
				return MergeKMerPST( GetMotif( iOne ), *GetPST( iTwo ), dCutoff ); }

		return MergeKMers( GetMotif( iOne ), GetMotif( iTwo ), dCutoff ); }

	uint32_t RemoveRCs( uint32_t iMotif, float dPenaltyGap, float dPenaltyMismatch ) {

		switch( GetType( iMotif ) ) {
			case ETypePST:
				return CCoalesceMotifLibraryImpl::RemoveRCs( *GetPST( iMotif ), dPenaltyGap,
					dPenaltyMismatch );

			case ETypeRC:
				return (uint32_t)m_veciRC2KMer[ iMotif - GetBaseRCs( ) ]; }

		return iMotif; }

	/*!
	 * \brief
	 * Returns a motif ID representing the "merger" of the input motif (which can be of any type) with an
	 * empty PST.
	 * 
	 * \param iMotif
	 * Motif to be "merged" into a PST.
	 * 
	 * \returns
	 * -1 if the motif cannot be merged; the ID of the new PST otherwise.
	 */
	uint32_t Merge( uint32_t iMotif ) {
		uint32_t	iRet;
		CPST*		pPST;

		if( !( pPST = CreatePST( iRet ) ) )
			return -1;
		switch( GetType( iMotif ) ) {
			case ETypeRC:
				return MergeRCPST( iMotif, *pPST, FLT_MAX );

			case ETypePST:
				return MergePSTs( *GetPST( iMotif ), *pPST, FLT_MAX ); }

		return MergeKMerPST( GetMotif( iMotif ), *pPST, FLT_MAX ); }

	float Align( uint32_t iOne, uint32_t iTwo, float dCutoff ) {

		switch( GetType( iOne ) ) {
			case ETypeRC:
				switch( GetType( iTwo ) ) {
					case ETypeKMer:
						return AlignKMerRC( GetMotif( iTwo ), iOne, dCutoff );

					case ETypeRC:
						return AlignRCs( iOne, iTwo, dCutoff );

					case ETypePST:
						return AlignRCPST( iOne, *GetPST( iTwo ), dCutoff ); }

			case ETypePST:
				switch( GetType( iTwo ) ) {
					case ETypeKMer:
						return AlignKMerPST( GetMotif( iTwo ), *GetPST( iOne ), dCutoff );

					case ETypeRC:
						return AlignRCPST( iTwo, *GetPST( iOne ), dCutoff );

					case ETypePST:
						return AlignPSTs( *GetPST( iOne ), *GetPST( iTwo ), dCutoff ); } }

		switch( GetType( iTwo ) ) {
			case ETypeRC:
				return AlignKMerRC( GetMotif( iOne ), iTwo, dCutoff );

			case ETypePST:
				return AlignKMerPST( GetMotif( iOne ), *GetPST( iTwo ), dCutoff ); }

		return AlignKMers( GetMotif( iOne ), GetMotif( iTwo ), dCutoff ); }

	/*!
	 * \brief
	 * Returns the number of motifs currently being managed by the library.
	 * 
	 * \returns
	 * Number of motifs (kmers, RCs, and PSTs) currently managed by the library.
	 * 
	 * \remarks
	 * For a given k, there will always be exactly 4^k kmers.  For odd k, there will be exactly (4^k)/2
	 * RCs; for even k, (4^k)/2 - (4^(k/2))/2.  The number of PSTs will vary over the lifetime of the
	 * library as they are (optionally) created by merging existing motifs.
	 */
	size_t GetMotifs( ) const {

// kmers plus reverse complements plus psts
		return ( GetBasePSTs( ) + GetPSTs( ) ); }

	/*!
	 * \brief
	 * Returns the underlying kmer length of the library.
	 * 
	 * \returns
	 * Length of kmers underlying the motif library.
	 */
	size_t GetK( ) const {

		return m_iK; }

	/*!
	 * \brief
	 * Calculates the kmer and RC motifs that match a given string of length k.
	 * 
	 * \param strKMer
	 * Kmer to be matched.
	 * 
	 * \param veciMotifs
	 * List to which matching kmer/RC motif IDs are appended.
	 * 
	 * \returns
	 * True if match was successful (even with no hits); false otherwise.
	 * 
	 * \remarks
	 * Input strings containing non-canonical bases (letters outside of the alphabet {A, C, G, T}) will be
	 * ignored.  Only kmer and RC motif IDs will be matched, regardless of the PSTs currently managed by the
	 * library.
	 * 
	 * \see
	 * GetMatch
	 */
	bool GetMatches( const std::string& strKMer, std::vector<uint32_t>& veciMotifs ) const {
		uint32_t	iMotif;
		size_t		iRC;

		if( IsIgnorableKMer( strKMer ) )
			return true;
// kmer
		if( ( iMotif = KMer2ID( strKMer ) ) == -1 )
			return false;
		veciMotifs.push_back( iMotif );
// reverse complement
		if( ( iRC = m_veciKMer2RC[ iMotif ] ) != -1 )
			veciMotifs.push_back( GetBaseRCs( ) + iRC );
		return true; }

	/*!
	 * \brief
	 * Sets the alignment score penalty for gaps (insertions or deletions).
	 * 
	 * \param dPenalty
	 * Alignment score penalty for gaps.
	 * 
	 * \remarks
	 * Alignments are internally ungapped, so gap penalties are only incurred at the ends (i.e. by overhangs).
	 * 
	 * \see
	 * GetPenaltyGap | SetPenaltyMismatch | Merge
	 */
	void SetPenaltyGap( float dPenalty ) {

		m_dPenaltyGap = dPenalty; }

	/*!
	 * \brief
	 * Gets the alignment score penalty for gaps (insertions or deletions).
	 * 
	 * \returns
	 * Alignment score penalty for gaps.
	 * 
	 * \remarks
	 * Alignments are internally ungapped, so gap penalties are only incurred at the ends (i.e. by overhangs).
	 * 
	 * \see
	 * SetPenaltyGap | GetPenaltyMismatch | Merge
	 */
	float GetPenaltyGap( ) const {

		return m_dPenaltyGap; }

	/*!
	 * \brief
	 * Sets the alignment score penalty for mismatches.
	 * 
	 * \param dPenalty
	 * Alignment score penalty for mismatches.
	 * 
	 * \see
	 * GetPenaltyMismatch | SetPenaltyGap | Merge
	 */
	void SetPenaltyMismatch( float dPenalty ) {

		m_dPenaltyMismatch = dPenalty; }

	/*!
	 * \brief
	 * Gets the alignment score penalty for mismatches.
	 * 
	 * \returns
	 * Alignment score penalty for mismatches.
	 * 
	 * \see
	 * SetPenaltyMismatch | GetPenaltyGap | Merge
	 */
	float GetPenaltyMismatch( ) const {

		return m_dPenaltyMismatch; }
};

}

#endif // COALESCEMOTIFS_H

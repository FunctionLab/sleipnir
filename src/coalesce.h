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
#ifndef COALESCE_H
#define COALESCE_H

#include "coalescei.h"

namespace Sleipnir {

class CCoalesceMotifLibrary : CCoalesceMotifLibraryImpl {
public:
	CCoalesceMotifLibrary( size_t iK ) : CCoalesceMotifLibraryImpl( iK ) { }

	float GetMatch( const std::string& strSequence, uint32_t iMotif, size_t iOffset,
		SCoalesceModifierCache& sModifiers ) const;

	std::string GetMotif( uint32_t iMotif ) const {

		return CCoalesceMotifLibraryImpl::GetMotif( iMotif ); }

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

	size_t GetMotifs( ) const {

// kmers plus reverse complements plus psts
		return ( GetBasePSTs( ) + GetPSTs( ) ); }

	size_t GetK( ) const {

		return m_iK; }

	bool GetMatches( const std::string& strKMer, std::vector<uint32_t>& veciMotifs ) const {
		uint32_t	iMotif;

		if( IsIgnorableKMer( strKMer ) )
			return true;
// kmer
		if( ( iMotif = KMer2ID( strKMer ) ) == -1 )
			return false;
		veciMotifs.push_back( iMotif );
// reverse complement
		veciMotifs.push_back( GetBaseRCs( ) + m_vecKMer2RC[ iMotif ] );
		return true; }

	void SetPenaltyGap( float dPenalty ) {

		m_dPenaltyGap = dPenalty; }

	float GetPenaltyGap( ) const {

		return m_dPenaltyGap; }

	void SetPenaltyMismatch( float dPenalty ) {

		m_dPenaltyMismatch = dPenalty; }

	float GetPenaltyMismatch( ) const {

		return m_dPenaltyMismatch; }
};

class CCoalesceCluster : public CCoalesceClusterImpl {
public:
	bool Initialize( const CPCL& PCL, CCoalesceCluster& Pot,
		std::set<std::pair<size_t, size_t> >& setpriiSeeds, size_t iPairs, float dPValue, size_t iThreads );
	void Subtract( CPCL& PCL ) const;
	void Subtract( vector<CCoalesceGeneScores>& vecGeneScores ) const;
	bool SelectConditions( const CPCL& PCL, const std::vector<CCoalesceImpl::SDataset>& vecsDatasets,
		const CCoalesceCluster& Pot, float dPValue );
	bool SelectMotifs( const vector<CCoalesceGeneScores>& vecGeneScores,
		const CCoalesceGroupHistograms& HistsCluster, const CCoalesceGroupHistograms& HistsPot, float dPValue,
		size_t iThreads, const CCoalesceMotifLibrary* pMotifs = NULL );
	bool SelectGenes( const CPCL& PCL, const std::vector<CCoalesceGeneScores>& vecGeneScores,
		const CCoalesceGroupHistograms& HistsCluster, const CCoalesceGroupHistograms& HistsPot,
		size_t iThreads, CCoalesceCluster& Pot, float dPValue, const CCoalesceMotifLibrary* pMotifs = NULL );
	void CalculateHistograms( const std::vector<CCoalesceGeneScores>& vecGeneScores,
		CCoalesceGroupHistograms& HistogramsCluster, CCoalesceGroupHistograms* pHistogramsPot ) const;
	bool Save( const std::string& strDirectory, size_t iID, const CPCL& PCL,
		const CCoalesceMotifLibrary* pMotifs = NULL ) const;
	void Save( std::ostream&, size_t iID, const CPCL& PCL, const CCoalesceMotifLibrary* pMotifs = NULL ) const;

	bool IsConverged( ) {

		return ( m_setiHistory.find( GetHash( ) ) != m_setiHistory.end( ) ); }

	bool IsEmpty( ) const {

		return ( m_setiGenes.empty( ) || m_setiConditions.empty( ) ); }

	void Add( size_t iGene ) {

		m_setiGenes.insert( iGene ); }

	void Snapshot( const std::vector<CCoalesceGeneScores>& vecGeneScores,
		CCoalesceGroupHistograms& Histograms ) {

		Histograms.SetTotal( vecGeneScores, GetGenes( ) );
		m_setiHistory.insert( GetHash( ) );
		CCoalesceClusterImpl::Snapshot( m_setiConditions, m_veciPrevConditions );
		CCoalesceClusterImpl::Snapshot( m_setsMotifs, m_vecsPrevMotifs );
		CCoalesceClusterImpl::Snapshot( m_setiGenes, m_veciPrevGenes ); }

	const std::set<size_t>& GetGenes( ) const {

		return m_setiGenes; }

	const std::set<size_t>& GetConditions( ) const {

		return m_setiConditions; }

	const std::set<SMotifMatch>& GetMotifs( ) const {

		return m_setsMotifs; }

	bool IsGene( size_t iGene ) const {

		return CCoalesceClusterImpl::IsGene( iGene ); }

	bool IsCondition( size_t iCondition ) const {

		return CCoalesceClusterImpl::IsCondition( iCondition ); }
};

class CCoalesce : CCoalesceImpl {
public:
	bool Cluster( const CPCL& PCL, const CFASTA& FASTA, std::vector<CCoalesceCluster>& vecClusters );

	void SetPValueCorrelation( float dPValue ) {

		m_dPValueCorrelation = dPValue; }

	float GetPValueCorrelation( ) const {

		return m_dPValueCorrelation; }

	void SetBins( size_t iBins ) {

		m_iBins = iBins; }

	size_t GetBins( ) const {

		return m_iBins; }

	float GetPValueCondition( ) const {

		return m_dPValueCondition; }

	void SetPValueCondition( float dPValue ) {

		m_dPValueCondition = dPValue; }

	float GetPValueMotif( ) const {

		return m_dPValueMotif; }

	void SetPValueMotif( float dPValue ) {

		m_dPValueMotif = dPValue; }

	float GetProbabilityGene( ) const {

		return m_dProbabilityGene; }

	void SetProbabilityGene( float dProbability ) {

		m_dProbabilityGene = dProbability; }

	bool IsOutputIntermediate( ) const {

		return !GetDirectoryIntermediate( ).empty( ); }

	const std::string& GetDirectoryIntermediate( ) const {

		return m_strDirectoryIntermediate; }

	void SetDirectoryIntermediate( const std::string& strDirectoryIntermediate ) {

		m_strDirectoryIntermediate = strDirectoryIntermediate; }

	void SetMotifs( CCoalesceMotifLibrary& Motifs ) {

		if( m_fMotifs && m_pMotifs && ( m_pMotifs != &Motifs ) )
			delete m_pMotifs;
		m_pMotifs = &Motifs; }

	const CCoalesceMotifLibrary* GetMotifs( ) const {

		return m_pMotifs; }

	size_t GetK( ) const {

		return m_iK; }

	void SetK( size_t iK ) {

		m_iK = iK; }

	size_t GetBasesPerMatch( ) const {

		return m_iBasesPerMatch; }

	void SetBasesPerMatch( size_t iBasesPerMatch ) {

		m_iBasesPerMatch = iBasesPerMatch; }

	const std::string& GetSequenceCache( ) const {

		return m_strSequenceCache; }

	void SetSequenceCache( const std::string& strSequenceCache ) {

		m_strSequenceCache = strSequenceCache; }

	float GetPValueMerge( ) const {

		return m_dPValueMerge; }

	void SetPValueMerge( float dPValue ) {

		m_dPValueMerge = dPValue; }

	float GetCutoffMerge( ) const {

		return m_dCutoffMerge; }

	void SetCutoffMerge( float dCutoff ) {

		m_dCutoffMerge = dCutoff; }

	float GetPenaltyGap( ) const {

		return m_dPenaltyGap; }

	void SetPenaltyGap( float dPenalty ) {

		m_dPenaltyGap = dPenalty; }

	float GetPenaltyMismatch( ) const {

		return m_dPenaltyMismatch; }

	void SetPenaltyMismatch( float dPenalty ) {

		m_dPenaltyMismatch = dPenalty; }

	size_t GetSizeMinimum( ) const {

		return m_iSizeMinimum; }

	void SetSizeMinimum( size_t iSizeGenes ) {

		m_iSizeMinimum = iSizeGenes; }

	size_t GetSizeMaximum( ) const {

		return m_iSizeMaximum; }

	void SetSizeMaximum( size_t iSizeMotifs ) {

		m_iSizeMaximum = iSizeMotifs; }

	void ClearDatasets( ) {

		m_vecsDatasets.clear( ); }

	bool AddDataset( const std::set<size_t>& setiDataset ) {
		size_t								i;
		std::set<size_t>::const_iterator	iterExperiment;

		if( setiDataset.empty( ) )
			return true;
		for( iterExperiment = setiDataset.begin( ); iterExperiment != setiDataset.end( ); ++iterExperiment )
			for( i = 0; i < m_vecsDatasets.size( ); ++i )
				if( m_vecsDatasets[ i ].IsCondition( *iterExperiment ) )
					return false;

		m_vecsDatasets.push_back( SDataset( setiDataset ) );
		return true; }

	void SetNumberCorrelation( size_t iPairs ) {

		m_iNumberCorrelation = iPairs; }

	size_t GetNumberCorrelation( ) const {

		return m_iNumberCorrelation; }

	void SetThreads( size_t iThreads ) {

		m_iThreads = iThreads; }

	size_t GetThreads( ) const {

		return m_iThreads; }

	void SetNucleosomes( const CFASTA& FASTANucleosomes ) {

		m_pPCLNucleosomes = NULL;
		m_pFASTANucleosomes = &FASTANucleosomes; }

	void SetNucleosomes( const CPCL& PCLNucleosomes ) {

		m_pFASTANucleosomes = NULL;
		m_pPCLNucleosomes = &PCLNucleosomes; }

	void ClearNucleosomes( ) {

		m_pFASTANucleosomes = NULL;
		m_pPCLNucleosomes = NULL; }
};

}

#endif // COALESCE_H

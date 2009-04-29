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
#ifndef COALESCECLUSTER_H
#define COALESCECLUSTER_H

#include "coalesceclusteri.h"

namespace Sleipnir {

class CCoalesceCluster : public CCoalesceClusterImpl {
public:
	bool Initialize( const CPCL& PCL, CCoalesceCluster& Pot, const std::vector<SCoalesceDataset>& vecsDatasets,
		std::set<std::pair<size_t, size_t> >& setpriiSeeds, size_t iPairs, float dPValue, size_t iThreads );
	void Subtract( CPCL& PCL, const CCoalesceCluster& Pot ) const;
	void Subtract( CCoalesceGeneScores& GeneScores ) const;
	bool SelectConditions( const CPCL& PCL, const CCoalesceCluster& Pot, size_t iThreads, float dPValue,
		float dZScore );
	bool SelectMotifs( const CCoalesceGroupHistograms& HistsCluster, const CCoalesceGroupHistograms& HistsPot,
		float dPValue, float dZScore, size_t iMaxMotifs, size_t iThreads,
		const CCoalesceMotifLibrary* pMotifs = NULL );
	bool SelectGenes( const CPCL& PCL, const CCoalesceGeneScores& GeneScores,
		const CCoalesceGroupHistograms& HistsCluster, const CCoalesceGroupHistograms& HistsPot,
		size_t iMinimum, size_t iThreads, CCoalesceCluster& Pot, float dPValue,
		const CCoalesceMotifLibrary* pMotifs = NULL );
	void CalculateHistograms( const CCoalesceGeneScores& GeneScores,
		CCoalesceGroupHistograms& HistogramsCluster, CCoalesceGroupHistograms* pHistogramsPot ) const;
	size_t Open( const std::string& strPCL, size_t iSkip, const CPCL& PCL,
		CCoalesceMotifLibrary* pMotifs = NULL );
	size_t Open( std::istream& istm, const CPCL& PCL, CCoalesceMotifLibrary* pMotifs = NULL );
	bool Open( const CHierarchy& Hierarchy, const std::vector<CCoalesceCluster>& vecClusters,
		const std::vector<std::string>& vecstrClusters, float dFraction, float dCutoff, size_t iCutoff,
		CCoalesceMotifLibrary* pMotifs = NULL );
	bool Save( const std::string& strDirectory, size_t iID, const CPCL& PCL,
		const CCoalesceMotifLibrary* pMotifs = NULL ) const;
	void Save( std::ostream&, size_t iID, const CPCL& PCL, const CCoalesceMotifLibrary* pMotifs = NULL,
		float dCutoffPWMs = 0, float dPenaltyGap = 0, float dPenaltyMismatch = 0, bool fNoRCs = false ) const;
	float GetSimilarity( const CCoalesceCluster& Cluster, size_t iGenes, size_t iDatasets ) const;
	void Snapshot( const CCoalesceGeneScores& GeneScores, CCoalesceGroupHistograms& Histograms );
	bool LabelMotifs( const CCoalesceMotifLibrary&, float dPenaltyGap, float dPenaltyMismatch,
		float dPValue );

	bool IsConverged( ) {

		return ( m_setiHistory.find( GetHash( ) ) != m_setiHistory.end( ) ); }

	bool IsEmpty( ) const {

		return ( m_setiGenes.empty( ) || m_setiDatasets.empty( ) ); }

	void SetGenes( size_t iGenes ) {
		size_t	i;

		m_setiGenes.clear( );
		for( i = 0; i < iGenes; ++i )
			m_setiGenes.insert( i ); }

	const std::set<size_t>& GetGenes( ) const {

		return CCoalesceClusterImpl::GetGenes( ); }

	const std::set<size_t>& GetDatasets( ) const {

		return m_setiDatasets; }

	const std::set<SMotifMatch>& GetMotifs( ) const {

		return m_setsMotifs; }

	bool IsGene( size_t iGene ) const {

		return CCoalesceClusterImpl::IsGene( iGene ); }

	bool IsDataset( size_t iDataset ) const {

		return ( m_setiDatasets.find( iDataset ) != m_setiDatasets.end( ) ); }

	void RemoveGenes( const std::vector<size_t>& veciGenes ) {
		size_t	i;

		for( i = 0; i < veciGenes.size( ); ++i )
			m_setiGenes.erase( veciGenes[ i ] ); }
};

}

#endif // COALESCECLUSTER_H

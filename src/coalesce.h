#ifndef COALESCE_H
#define COALESCE_H

#include "coalescei.h"

namespace Sleipnir {

struct SHistogram : SHistogramImpl {
public:
	void Initialize( size_t iSize ) {

		m_iTotal = 0;
		m_vecsBins.resize( iSize );
		std::fill( m_vecsBins.begin( ), m_vecsBins.end( ), 0 ); }

	bool Add( size_t iBin, unsigned short sCount ) {

		if( m_vecsBins.empty( ) )
			return false;
		if( iBin >= m_vecsBins.size( ) )
			iBin = m_vecsBins.size( ) - 1;

		m_vecsBins[ iBin ] += sCount;
		m_iTotal += sCount;
		return true; }

	bool Add( size_t iBin ) {

		return Add( iBin, 1 ); }

	unsigned short Get( size_t iBin ) const {

		if( m_vecsBins.empty( ) )
			return -1;
		if( iBin >= m_vecsBins.size( ) )
			iBin = m_vecsBins.size( ) - 1;

		return m_vecsBins[ iBin ]; }

	size_t GetBins( ) const {

		return Get( ).size( ); }

	const std::vector<unsigned short>& Get( ) const {

		return m_vecsBins; }

	bool Add( const SHistogram& Histogram ) {
		size_t	i;

		if( Histogram.GetBins( ) != GetBins( ) )
			return false;

		for( i = 0; i < GetBins( ); ++i )
			if( !Add( i, Histogram.Get( i ) ) )
				return false;

		return true; }

	size_t GetTotal( ) const {

		return m_iTotal; }
};

class CCoalesceCluster : public CCoalesceClusterImpl {
public:
	bool Initialize( const CPCL& PCL, CCoalesceCluster& Pot, float dPValue );
	void Subtract( CPCL& PCL ) const;
	bool SelectConditions( const CPCL& PCL, const CCoalesceCluster& Pot, float dPValue );
	bool SelectMotifs( const CCoalesceHistograms& HistsCluster, const CCoalesceHistograms& HistsPot,
		float dPValue );
	bool SelectGenes( const CPCL& PCL, const std::vector<CCoalesceHistograms>& vecHistograms,
		const CCoalesceHistograms& HistsCluster, const CCoalesceHistograms& HistsPot, CCoalesceCluster& Pot,
		float dPValue );
	void CalculateHistograms( const std::vector<CCoalesceHistograms>& vecHistograms,
		CCoalesceHistograms& Histograms ) const;
	bool Save( size_t iID, const CPCL& PCL ) const;

	bool IsConverged( ) {

		return ( CCoalesceClusterImpl::IsConverged( m_setiConditions, m_veciPrevConditions ) &&
			CCoalesceClusterImpl::IsConverged( m_setiGenes, m_veciPrevGenes ) &&
			CCoalesceClusterImpl::IsConverged( m_setiMotifs, m_veciPrevMotifs ) ); }

	bool IsEmpty( ) const {

		return ( m_setiGenes.empty( ) && m_setiConditions.empty( ) ); }

	void Add( size_t iGene ) {

		m_setiGenes.insert( iGene ); }

	void Snapshot( ) {

		CCoalesceClusterImpl::Snapshot( m_setiConditions, m_veciPrevConditions );
		CCoalesceClusterImpl::Snapshot( m_setiMotifs, m_veciPrevMotifs );
		CCoalesceClusterImpl::Snapshot( m_setiGenes, m_veciPrevGenes ); }

	const std::set<size_t>& GetGenes( ) const {

		return m_setiGenes; }

	const std::set<size_t>& GetConditions( ) const {

		return m_setiConditions; }

	const std::set<size_t>& GetMotifs( ) const {

		return m_setiMotifs; }

	bool IsGene( size_t iGene ) const {

		return CCoalesceClusterImpl::IsGene( iGene ); }

	bool IsCondition( size_t iCondition ) const {

		return CCoalesceClusterImpl::IsCondition( iCondition ); }
};

class CCoalesce : CCoalesceImpl {
public:
	static size_t GetKMer( const std::string& strKMer );

	static size_t CountKMers( size_t iK ) {

		return ( 1 << ( 2 * iK ) ); }

	bool Cluster( const CPCL& PCL, const CFASTA& FASTA, std::vector<CCoalesceCluster>& vecClusters );

	void SetK( size_t iK ) {

		m_iK = iK; }

	size_t GetK( ) const {

		return m_iK; }

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

		return m_fOutputIntermediate; }

	void SetOutputIntermediate( bool fOutputIntermediate ) {

		m_fOutputIntermediate = fOutputIntermediate; }
};

}

#endif // COALESCE_H

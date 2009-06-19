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
#ifndef COALESCECLUSTERI_H
#define COALESCECLUSTERI_H

#include "coalescestructsi.h"

namespace Sleipnir {

class CCoalesceCluster;
class CCoalesceGeneScores;
class CCoalesceGroupHistograms;
struct SMotifMatch;

class CCoalesceClusterImpl {
protected:
	typedef std::vector<std::map<std::string, std::set<SMotifMatch> > >	TVecMapStrSetSMotifs;

	static const char	c_cStar		= '*';
	static const char	c_szMotifs[];
	static const char	c_szConditions[];
	static const char	c_szGenes[];

	struct SDataset {
		const SCoalesceDataset*	m_psDataset;
		double					m_dP;
		float					m_dZ;
		std::vector<float>		m_vecdCentroid;

		size_t GetConditions( ) const {

			return m_psDataset->GetConditions( ); }

		size_t GetCondition( size_t iCondition ) const {

			return m_psDataset->GetCondition( iCondition ); }
	};

	struct SThreadCentroid {
		CCoalesceCluster*	m_pCluster;
		const CPCL&			m_PCL;

		SThreadCentroid( CCoalesceCluster* pCluster, const CPCL& PCL ) : m_pCluster(pCluster), m_PCL(PCL) { }
	};

	struct SThreadSignificantGene {
		size_t									m_iOffset;
		size_t									m_iStep;
		std::vector<bool>*						m_pvecfSignificant;
		const CPCL*								m_pPCL;
		const CCoalesceMotifLibrary*			m_pMotifs;
		const CCoalesceGeneScores*				m_pGeneScores;
		const CCoalesceGroupHistograms*			m_pHistsCluster;
		const CCoalesceGroupHistograms*			m_pHistsPot;
		const CCoalesceCluster*					m_pCluster;
		const CCoalesceCluster*					m_pPot;
		const std::vector<size_t>*				m_pveciDatasets;
		const std::vector<float>*				m_pvecdStdevs;
		float									m_dBeta;
		size_t									m_iMinimum;
		float									m_dProbability;
	};

	struct SThreadSelectMotif {
		uint32_t						m_iOffset;
		size_t							m_iStep;
		const CCoalesceMotifLibrary*	m_pMotifs;
		const CCoalesceGroupHistograms*	m_pHistsCluster;
		const CCoalesceGroupHistograms*	m_pHistsPot;
		float							m_dPValue;
		float							m_dZScore;
		const std::vector<uint32_t>*	m_pveciMotifs;
		std::vector<SMotifMatch>		m_vecsMotifs;
	};

	struct SThreadSeedPair {
		size_t										m_iOffset;
		size_t										m_iStep;
		const CPCL*									m_pPCL;
		float										m_dFraction;
		const std::set<std::pair<size_t, size_t> >*	m_psetpriiSeeds;
		double										m_dMaxCorr;
		double										m_dMinP;
		size_t										m_iOne;
		size_t										m_iTwo;
	};

	struct SThreadSelectCondition {
		size_t						m_iOffset;
		size_t						m_iStep;
		const std::vector<size_t>*	m_pveciCluster;
		const std::vector<size_t>*	m_pveciPot;
		std::vector<SDataset>*		m_pvecsDatasets;
		const CPCL*					m_pPCL;
	};

	static void* ThreadCentroid( void* );
	static void* ThreadSignificantGene( void* );
	static void* ThreadSelectMotif( void* );
	static void* ThreadSeedPair( void* );
	static void* ThreadSelectCondition( void* );
	static bool AddSignificant( const CCoalesceMotifLibrary&, uint32_t, const CCoalesceGroupHistograms&,
		const CCoalesceGroupHistograms&, float, float, std::vector<SMotifMatch>& );
	static size_t Open( const CHierarchy&, const std::vector<CCoalesceCluster>&,
		const std::vector<std::string>&, std::map<size_t, size_t>&, std::map<size_t, size_t>&,
		TVecMapStrSetSMotifs& );
	static bool OpenMotifs( CCoalesceMotifLibrary&, const CHierarchy&, const std::vector<SMotifMatch>&, float,
		std::set<SMotifMatch>& );

	template<class tType>
	static bool IsConverged( const std::set<tType>& setNew, std::vector<tType>& vecOld ) {
		size_t				i;
		std::vector<tType>	vecNew;

		if( setNew.size( ) != vecOld.size( ) )
			return false;
		Snapshot( setNew, vecNew );
		for( i = 0; i < vecNew.size( ); ++i )
			if( vecNew[ i ] != vecOld[ i ] )
				return false;

		return true; }

	template<class tType>
	static void Snapshot( const std::set<tType>& setNew, std::vector<tType>& vecOld ) {

		vecOld.resize( setNew.size( ) );
		std::copy( setNew.begin( ), setNew.end( ), vecOld.begin( ) );
		std::sort( vecOld.begin( ), vecOld.end( ) ); }

	template<class tType>
	static size_t GetHash( const std::set<tType>& set ) {
		size_t										iRet;
		typename std::set<tType>::const_iterator	iter;

		for( iRet = 0,iter = set.begin( ); iter != set.end( ); ++iter )
			iRet ^= GetHash( *iter );

		return iRet; }

	static size_t GetHash( size_t iValue ) {

		return ( iValue * ( (size_t)-1 / 20000 ) ); }

	static size_t GetHash( const SMotifMatch& sMotif ) {

		return sMotif.GetHash( ); }

	void Add( size_t, CCoalesceCluster& );
	bool AddCorrelatedGenes( const CPCL&, CCoalesceCluster&, float );
	bool AddSeedPair( const CPCL&, CCoalesceCluster&, std::set<std::pair<size_t, size_t> >&, float, float,
		size_t );
	void CalculateCentroid( const CPCL& );
	bool IsSignificant( size_t, const CPCL&, const std::vector<float>&, const CCoalesceMotifLibrary*,
		const CCoalesceGeneScores&, const CCoalesceGroupHistograms&, const CCoalesceGroupHistograms&,
		const CCoalesceCluster&, const std::vector<size_t>&, float, size_t, float ) const;
	bool CalculateProbabilityExpression( size_t, const CPCL&, const std::vector<float>&,
		const CCoalesceCluster&, const std::vector<size_t>&, bool, float&, float& ) const;
	bool CalculateProbabilityMotifs( const CCoalesceGeneScores&, size_t, const CCoalesceGroupHistograms&,
		const CCoalesceGroupHistograms&, bool, size_t, float&, float& ) const;
	bool SaveCopy( const CPCL&, const std::set<size_t>&, size_t, CPCL&, size_t, bool ) const;
	bool OpenMotifs( const std::set<SMotifMatch>&, CCoalesceMotifLibrary&, float );
	bool OpenMotifsHeuristic( const std::set<SMotifMatch>&, CCoalesceMotifLibrary&, float, size_t );

	size_t GetConditions( size_t iDataset ) const {

		if( iDataset < m_vecsDatasets.size( ) ) {
			const SDataset&	sDataset	= m_vecsDatasets[ iDataset ];

			if( sDataset.m_psDataset )
				return sDataset.m_psDataset->GetConditions( ); }

		return 1; }

	size_t GetCondition( size_t iDataset, size_t iCondition ) const {

		if( iDataset < m_vecsDatasets.size( ) ) {
			const SDataset&	sDataset	= m_vecsDatasets[ iDataset ];

			if( sDataset.m_psDataset )
				return sDataset.m_psDataset->GetCondition( iCondition ); }

		return iDataset; }

	bool IsGene( size_t iGene ) const {

		return ( m_setiGenes.find( iGene ) != m_setiGenes.end( ) ); }

	size_t GetHash( ) const {

		return ( GetHash( m_setiDatasets ) ^ GetHash( m_setiGenes ) ^ GetHash( m_setsMotifs ) ); }

	void GetConditions( std::set<size_t>& setiConditions ) const {
		set<size_t>::const_iterator	iterDataset;
		size_t						i;

		for( iterDataset = m_setiDatasets.begin( ); iterDataset != m_setiDatasets.end( ); ++iterDataset )
			for( i = 0; i < GetConditions( *iterDataset ); ++i )
				setiConditions.insert( GetCondition( *iterDataset, i ) ); }

	const std::set<size_t>& GetGenes( ) const {

		return m_setiGenes; }

	void Clear( ) {

		m_setiDatasets.clear( );
		m_setiGenes.clear( );
		m_setsMotifs.clear( );
		m_veciPrevDatasets.clear( );
		m_veciPrevGenes.clear( );
		m_vecsPrevMotifs.clear( );
		m_veciCounts.clear( );
		m_vecdCentroid.clear( );
		m_vecdStdevs.clear( );
		m_setiHistory.clear( );
		m_vecdPriors.clear( );
		m_vecsDatasets.clear( ); }

	std::set<size_t>			m_setiDatasets;
	std::set<size_t>			m_setiGenes;
	std::set<SMotifMatch>		m_setsMotifs;
	std::vector<size_t>			m_veciPrevDatasets;
	std::vector<size_t>			m_veciPrevGenes;
	std::vector<SMotifMatch>	m_vecsPrevMotifs;
	std::vector<size_t>			m_veciCounts;
	std::vector<float>			m_vecdCentroid;
	std::vector<float>			m_vecdStdevs;
	std::set<size_t>			m_setiHistory;
	std::vector<float>			m_vecdPriors;
	std::vector<SDataset>		m_vecsDatasets;
};

}

#endif // COALESCECLUSTERI_H

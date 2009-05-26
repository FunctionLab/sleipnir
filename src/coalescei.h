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
#ifndef COALESCEI_H
#define COALESCEI_H

#include <algorithm>
#include <iostream>
#include <set>
#include <sstream>

#include "coalescebasei.h"
#include "statistics.h"

namespace Sleipnir {

class CCoalesceCluster;
class CCoalesceGeneScores;
class CCoalesceGroupHistograms;
class CCoalesceMotifLibrary;
class CFASTA;
class CPCL;
struct SCoalesceDataset;
struct SCoalesceModifierCache;
struct SCoalesceModifiers;
struct SFASTASequence;
struct SMotifMatch;

template<class tValue = float, class tCount = unsigned short>
class CCoalesceHistogramSet {
public:
	double ZScore( size_t iMember, const CCoalesceHistogramSet& HistSet, double& dAveOne, double& dAverage,
		double& dZ, bool fCount = true ) const {

		return ZScore( iMember, HistSet, iMember, fCount, dAveOne, dAverage, dZ ); }

	double ZScore( size_t iOne, const CCoalesceHistogramSet& HistSet, size_t iTwo, bool fCount,
		double& dAveOne, double& dAverage, double& dZ ) const {
		static const double	c_dEpsilon	= 1e-6;
		tValue	AveOne, VarOne, AveTwo, VarTwo;
		double	dStd;
		size_t	iTotal;

		if( !GetEdges( ) || ( GetEdges( ) != HistSet.GetEdges( ) ) )
			return 1;

		Sums( iOne, AveOne, VarOne );
		HistSet.Sums( iTwo, AveTwo, VarTwo );
		iTotal = GetTotal( ) + HistSet.GetTotal( );
		dAverage = (double)( AveOne + AveTwo ) / iTotal;
		dAveOne = (double)AveOne / GetTotal( );
		dStd = sqrt( ( (double)( VarOne + VarTwo ) / iTotal ) - ( dAverage * dAverage ) );

		dZ = ( dAveOne - dAverage ) / dStd;
		return ( ( dStd > c_dEpsilon ) ? ( 2 * CStatistics::ZTest( dZ, fCount ?
			min( GetTotal( ), HistSet.GetTotal( ) ) : 1 ) ) :
			( ( fabs( dAveOne - dAverage ) < c_dEpsilon ) ? 0 : 1 ) ); }

	double CohensD( size_t iMember, const CCoalesceHistogramSet& HistSet, double& dAveOne, double& dAverage,
		double& dZ, bool fCount = true ) const {

		return CohensD( iMember, HistSet, iMember, fCount, dAveOne, dAverage, dZ ); }

	double CohensD( size_t iOne, const CCoalesceHistogramSet& HistSet, size_t iTwo, bool fCount,
		double& dAveOne, double& dAverage, double& dZ ) const {
		static const double	c_dEpsilon	= 1e-6;
		tValue	AveOne, VarOne, AveTwo, VarTwo;
		double	dAveTwo, dVarOne, dVarTwo, dStd;

		if( !GetEdges( ) || ( GetEdges( ) != HistSet.GetEdges( ) ) )
			return 1;

		Sums( iOne, AveOne, VarOne );
		HistSet.Sums( iTwo, AveTwo, VarTwo );
		dAverage = (double)( AveOne + AveTwo ) / ( GetTotal( ) + HistSet.GetTotal( ) );
		dAveOne = (double)AveOne / GetTotal( );
		dAveTwo = (double)AveTwo / HistSet.GetTotal( );
		dVarOne = max( 0.0, ( (double)VarOne / GetTotal( ) ) - ( dAveOne * dAveOne ) );
		dVarTwo = max( 0.0, ( (double)VarTwo / HistSet.GetTotal( ) ) - ( dAveTwo * dAveTwo ) );
		dStd = sqrt( ( dVarOne + dVarTwo ) / 2 );

		dZ = dStd ? ( ( dAveOne - dAveTwo ) / dStd ) : 0;
// dZ *= 1 - pow( (float)GetTotal( ) / ( GetTotal( ) + HistSet.GetTotal( ) ), 1 ); // doesn't work
// dZ *= exp( -(float)GetTotal( ) / ( GetTotal( ) + HistSet.GetTotal( ) ) ); // doesn't work
// dZ *= exp( -(float)GetTotal( ) / HistSet.GetTotal( ) ); // works pretty well
// dZ *= fabs( (float)( GetTotal( ) - HistSet.GetTotal( ) ) ) / max( GetTotal( ), HistSet.GetTotal( ) ); // works pretty well
// This prevents large clusters from blowing up the motif set
/*
		if( iOne == iTwo )
			dZ *= pow( fabs( (float)( GetTotal( ) - HistSet.GetTotal( ) ) ) /
				max( GetTotal( ), HistSet.GetTotal( ) ), max( 1.0f,
				log10f( (float)min( GetTotal( ), HistSet.GetTotal( ) ) ) ) );
//*/

		return ( ( dStd > c_dEpsilon ) ? ( 2 * CStatistics::ZTest( dZ, fCount ?
			min( GetTotal( ), HistSet.GetTotal( ) ) : 1 ) ) :
			( ( fabs( dAveOne - dAveTwo ) < c_dEpsilon ) ? 0 : 1 ) ); }

	void Initialize( size_t iMembers, const std::vector<tValue>& vecEdges ) {

		m_vecEdges.resize( vecEdges.size( ) );
		std::copy( vecEdges.begin( ), vecEdges.end( ), m_vecEdges.begin( ) );
		m_iZero = GetBin( 0 );
		SetMembers( iMembers ); }

	void Initialize( size_t iMembers, size_t iBins, tValue Step ) {
		std::vector<tValue>	vecEdges;
		size_t				i;

		vecEdges.resize( iBins );
		for( i = 0; i < vecEdges.size( ); ++i )
			vecEdges[ i ] = i * Step;

		Initialize( iMembers, vecEdges ); }

	bool Add( size_t iMember, tValue Value, tCount Count ) {
		size_t	iBin;

// lock
		if( m_vecEdges.empty( ) )
// unlock
			return false;

		SetMembers( max( GetMembers( ), iMember + 1 ) );
		if( ( iBin = GetBin( Value ) ) != m_iZero ) {
			m_vecCounts[ GetOffset( iMember ) + iBin ] += Count;
			m_vecTotal[ iMember ] += Count; }
// unlock

		return true; }

	tCount Integrate( size_t iMember, tValue Value, bool fUp ) const {
		tCount	Ret;
		size_t	iBin;

		for( Ret = 0,iBin = GetBin( Value ); iBin < GetEdges( ); iBin += ( fUp ? 1 : -1 ) )
			Ret += Get( iMember, iBin );

		return Ret; }

	tCount Get( size_t iMember, size_t iBin ) const {

		if( iBin >= GetEdges( ) )
			return 0;
		if( iBin == m_iZero )
			return ( GetTotal( ) - ( ( iMember < GetMembers( ) ) ? m_vecTotal[ iMember ] : 0 ) );

		return ( ( iMember < GetMembers( ) ) ? m_vecCounts[ GetOffset( iMember ) + iBin ] : 0 ); }

	tCount Get( size_t iMember, tValue Value ) const {

		return Get( iMember, GetBin( Value ) ); }

	size_t GetBin( tValue Value ) const {
		size_t	i;

		if( !GetEdges( ) )
			return -1;
		for( i = 0; i < GetEdges( ); ++i )
			if( Value <= GetEdge( i ) )
				break;

		return min( i, GetEdges( ) - 1 ); }

	size_t GetMembers( ) const {

		return m_vecTotal.size( ); }

	size_t GetEdges( ) const {

		return m_vecEdges.size( ); }

	tValue GetEdge( size_t iBin ) const {

		return m_vecEdges[ iBin ]; }

	std::string Save( size_t iMember ) const {
		std::ostringstream	ossm;
		size_t				i;

		if( GetEdges( ) ) {
			ossm << (size_t)GetTotal( ) << ':';
			for( i = 0; i < GetEdges( ); ++i )
				ossm << '\t' << (size_t)Get( iMember, i ); }

		return ossm.str( ); }

	void Clear( ) {

		fill( m_vecCounts.begin( ), m_vecCounts.end( ), 0 );
		fill( m_vecTotal.begin( ), m_vecTotal.end( ), 0 ); }

	tCount GetTotal( ) const {

		return m_Total; }

	void SetTotal( tCount Total ) {

		m_Total = Total; }

	const std::vector<tValue>& GetBins( ) const {

		return m_vecEdges; }

protected:
	void SetMembers( size_t iMembers ) {

		m_vecTotal.resize( iMembers );
		m_vecCounts.resize( GetOffset( GetMembers( ) ) ); }

	size_t GetOffset( size_t iMember ) const {

		return ( iMember * GetEdges( ) ); }

	void Sums( size_t iMember, tValue& Sum, tValue& SumSq ) const {
		size_t	i;
		tValue	Cur, BinLow, BinHigh;

		for( SumSq = 0,i = 0; i < GetEdges( ); ++i ) {
			BinLow = i ? BinHigh : ( GetEdge( i ) + GetEdge( i ) - GetEdge( i + 1 ) );
			BinHigh = GetEdge( i );
			if( ( i + 1 ) == GetEdges( ) )
				BinHigh += BinHigh - BinLow;
// Expected value of x^2 from BinLow to BinHigh
			Cur = ( pow( BinHigh, 3 ) - pow( BinLow, 3 ) ) / ( BinHigh - BinLow ) / 3;
			SumSq += Get( iMember, i ) * Cur; }
		for( Sum = 0,i = 0; i < GetEdges( ); ++i ) {
			BinLow = ( !i || ( ( i + 1 ) == GetEdges( ) ) ) ? GetEdge( i ) :
				( ( GetEdge( i ) + GetEdge( i - 1 ) ) / 2 );
			Cur = Get( iMember, i ) * BinLow;
			Sum += Cur; } }

	void AveVar( size_t iMember, double& dAve, double& dVar ) const {
		tValue	Ave, Var;

		Sums( iMember, Ave, Var );
		dAve = (double)Ave / GetTotal( );
		dVar = ( (double)Var / GetTotal( ) ) - ( dAve * dAve ); }

	size_t				m_iZero;
	tCount				m_Total;
	std::vector<tCount>	m_vecTotal;
	std::vector<tValue>	m_vecEdges;
	std::vector<tCount>	m_vecCounts;
};

// One score per type per subtype per gene per motif atom
class CCoalesceGeneScores : public CCoalesceSequencer<float*> {
public:
	CCoalesceGeneScores( ) : m_iMotifs(0), m_iGenes(0), m_iCapacity(0) {

		pthread_mutex_init( &m_mutx, NULL ); }

	virtual ~CCoalesceGeneScores( ) {
		size_t	iType, iSubtype;
		float*	ad;

		for( iType = 0; iType < GetTypes( ); ++iType )
			for( iSubtype = 0; iSubtype < GetSubsequences( iType ); ++iSubtype )
				if( ad = Get( iType, (ESubsequence)iSubtype ) )
					delete[] ad;
		pthread_mutex_destroy( &m_mutx ); }

	bool Add( size_t, const CCoalesceMotifLibrary&, const SFASTASequence&, SCoalesceModifierCache&,
		std::vector<std::vector<float> >&, std::vector<size_t>& );
	bool Add( size_t, const CCoalesceMotifLibrary&, const SFASTASequence&, SCoalesceModifierCache&, uint32_t,
		std::vector<float>&, std::vector<size_t>& );
	void Subtract( const SMotifMatch&, size_t );

	float* Get( size_t iType, ESubsequence eSubsequence, size_t iGene, bool fSet = false ) const {
		float*	ad;

		if( !( ad = Get( iType, eSubsequence ) ) )
			return NULL;
		ad += iGene * m_iCapacity;
		return ( ( fSet || !CMeta::IsNaN( ad[ 0 ] ) ) ? ad : NULL ); }

	float Get( size_t iType, ESubsequence eSubsequence, size_t iGene, uint32_t iMotif ) const {
		const float*	ad;

		return ( ( ad = Get( iType, eSubsequence, iGene ) ) ? ad[ iMotif ] : 0 ); }

	size_t GetMotifs( ) const {

		return m_iMotifs; }

	void SetGenes( size_t iGenes ) {

		m_iGenes = iGenes; }

	void Validate( ) const {
		size_t		iType, iSubsequence, iGene;
		uint32_t	iMotif;
		float		dCur, dTotal;

		for( iType = 0; iType < GetTypes( ); ++iType )
			for( iSubsequence = ( ESubsequenceBegin + 1 ); iSubsequence < GetSubsequences( iType );
				++iSubsequence )
				for( iGene = 0; iGene < m_iGenes; ++iGene ) {
					if( !Get( iType, (ESubsequence)iSubsequence, iGene ) )
						continue;
					for( iMotif = 0; iMotif < m_iMotifs; ++iMotif )
						if( ( dCur = Get( iType, (ESubsequence)iSubsequence, iGene, iMotif ) ) !=
							( dTotal = Get( iType, ESubsequenceTotal, iGene, iMotif ) ) )
							std::cerr << "INVALID" << '\t' << GetType( iType ) << '\t' << iSubsequence <<
								'\t' << iGene << '\t' << iMotif << '\t' << dCur << '\t' << dTotal <<
								std::endl; } }

protected:
	static const size_t	c_iLookahead	= 128;

	static bool Add( const CCoalesceMotifLibrary&, const std::string&, size_t, bool,
		std::vector<std::vector<float> >&, std::vector<size_t>&, size_t, SCoalesceModifierCache& );
	static bool Add( const CCoalesceMotifLibrary&, const std::string&, size_t, bool, uint32_t,
		std::vector<float>&, std::vector<size_t>&, size_t, SCoalesceModifierCache& );

	static void Add( ESubsequence eSubsequence, uint32_t iMotif, uint32_t iMotifs,
		vector<vector<float> >& vecvecdCounts, float dValue ) {
		std::vector<float>&	vecdCountsTotal	= vecvecdCounts[ ESubsequenceTotal ];
		std::vector<float>&	vecdCounts		= vecvecdCounts[ eSubsequence ];

		vecdCountsTotal.resize( iMotifs );
		vecdCountsTotal[ iMotif ] += dValue;
		vecdCounts.resize( iMotifs );
		vecdCounts[ iMotif ] += dValue; }

	float* Get( size_t iType, ESubsequence eSubsequence ) const {

		return CCoalesceSequencer<float*>::Get( iType, eSubsequence ); }

	void Set( size_t iType, ESubsequence eSubsequence, size_t iGene, uint32_t iMotif, float dValue,
		uint32_t iMotifs = 0 ) {
		float*	adTotal;
		size_t	i;

		pthread_mutex_lock( &m_mutx );
		Grow( iMotif, iMotifs );
		if( !( adTotal = Get( iType, eSubsequence ) ) ) {
			m_vecvecValues[ iType ][ eSubsequence ] = adTotal = new float[ m_iGenes * m_iCapacity ];
			memset( adTotal, 0, m_iGenes * m_iCapacity * sizeof(*adTotal) );
			for( i = 0; i < m_iGenes; ++i )
				adTotal[ i * m_iCapacity ] = CMeta::GetNaN( ); }
		{
			float*	ad	= Get( iType, eSubsequence, iGene, true );

			if( CMeta::IsNaN( ad[ 0 ] ) )
				ad[ 0 ] = 0;
			ad[ iMotif ] = dValue;
		}
		pthread_mutex_unlock( &m_mutx ); }

	void Grow( uint32_t iMotif, uint32_t iMotifs ) {
		size_t	iType, iSubtype, iGene, iTarget, iCapacity;
		float*	adTotal;
		float*	ad;
		float*	adFrom;
		float*	adTo;

		iTarget = max( iMotif + 1, iMotifs );
		if( iTarget <= m_iCapacity ) {
			if( iTarget > m_iMotifs )
				m_iMotifs = iTarget;
			return; }
		iCapacity = m_iCapacity;
		m_iCapacity = iTarget + c_iLookahead;
		for( iType = 0; iType < GetTypes( ); ++iType )
			for( iSubtype = 0; iSubtype < GetSubsequences( iType ); ++iSubtype ) {
				if( !( adTotal = Get( iType, (ESubsequence)iSubtype ) ) )
					continue;
				ad = new float[ m_iGenes * m_iCapacity ];
				for( iGene = 0,adFrom = adTotal,adTo = ad; iGene < m_iGenes;
					++iGene,adFrom += iCapacity,adTo += m_iCapacity ) {
					memcpy( adTo, adFrom, m_iMotifs * sizeof(*adTo) );
					memset( adTo + m_iMotifs, 0, ( m_iCapacity - m_iMotifs ) * sizeof(*adTo) ); }
				delete[] adTotal;
				m_vecvecValues[ iType ][ iSubtype ] = ad; }
		m_iMotifs = iTarget; }

	size_t			m_iGenes;
	size_t			m_iMotifs;
	size_t			m_iCapacity;
	pthread_mutex_t	m_mutx;
};

// One histogram per motif atom
class CCoalesceGroupHistograms : public CCoalesceSequencer<CCoalesceHistogramSet<> > {
public:
	CCoalesceGroupHistograms( size_t iBins, float dStep ) : m_iBins(iBins), m_dStep(dStep) { }

	bool Add( const CCoalesceGeneScores&, size_t, bool, uint32_t = -1 );
	void Save( std::ostream&, const CCoalesceMotifLibrary* ) const;
	bool IsSimilar( const CCoalesceMotifLibrary*, const SMotifMatch&, const SMotifMatch&, float ) const;

	void SetTotal( const CCoalesceGeneScores& GeneScores, const std::set<size_t>& setiGenes ) {
		std::set<size_t>::const_iterator	iterGene;
		size_t								i, iTypeUs, iTypeThem, iSubsequence;

		m_vecsTotals.resize( GetTypes( ) * ESubsequenceEnd );
		std::fill( m_vecsTotals.begin( ), m_vecsTotals.end( ), 0 );
		for( iterGene = setiGenes.begin( ); iterGene != setiGenes.end( ); ++iterGene )
			for( iTypeUs = 0; iTypeUs < GetTypes( ); ++iTypeUs ) {
				if( ( iTypeThem = GeneScores.GetType( GetType( iTypeUs ) ) ) == -1 )
					continue;
				for( iSubsequence = ESubsequenceBegin; iSubsequence < GetSubsequences( iTypeUs );
					++iSubsequence )
					if( GeneScores.Get( iTypeThem, (ESubsequence)iSubsequence, *iterGene ) )
						m_vecsTotals[ ( iTypeUs * ESubsequenceEnd ) + iSubsequence ]++; }

		for( i = iTypeUs = 0; iTypeUs < GetTypes( ); ++iTypeUs )
			for( iSubsequence = ESubsequenceBegin; iSubsequence < GetSubsequences( iTypeUs );
				++iSubsequence ) {
//				g_CatSleipnir.info( "CCoalesceGroupHistograms::SetTotal( ) type %s, subsequence %d contains %d genes with sequences",
//					GetType( iTypeUs ).c_str( ), iSubsequence, m_vecsTotals[ i ] );
				Get( iTypeUs, (ESubsequence)iSubsequence ).SetTotal( m_vecsTotals[ i++ ] ); } }

	void Validate( ) const {
		size_t		iType, iSubsequence, iEdge;
		uint32_t	iMotif;

		for( iType = 0; iType < GetTypes( ); ++iType ) {
			const CCoalesceHistogramSet<>&	HistsTotal	= Get( iType, ESubsequenceTotal );

			for( iSubsequence = ( ESubsequenceBegin + 1 ); iSubsequence < GetSubsequences( iType );
				++iSubsequence ) {
				const CCoalesceHistogramSet<>&	HistsCur	= Get( iType, (ESubsequence)iSubsequence );

				if( !HistsCur.GetTotal( ) )
					continue;
				for( iMotif = 0; iMotif < HistsCur.GetMembers( ); ++iMotif )
					for( iEdge = 0; iEdge < HistsCur.GetEdges( ); ++iEdge )
						if( HistsCur.Get( iMotif, iEdge ) != HistsTotal.Get( iMotif, iEdge ) ) {
							std::cerr << "INVALID" << '\t' << GetType( iType ) << '\t' << iSubsequence <<
								'\t' << iMotif << endl << HistsTotal.Save( iMotif ) << endl <<
								HistsCur.Save( iMotif ) << std::endl;
							break; } } } }

protected:
	uint32_t					m_iMotifs;
	size_t						m_iBins;
	float						m_dStep;
	std::vector<unsigned short>	m_vecsTotals;
};

class CCoalesceImpl {
protected:
	struct SThreadCombineMotif {
		size_t								m_iOffset;
		size_t								m_iStep;
		const std::vector<size_t>*			m_pveciPCL2FASTA;
		CCoalesceGeneScores*				m_pGeneScores;
		const CCoalesceMotifLibrary*		m_pMotifs;
		uint32_t							m_iMotif;
		const CFASTA*						m_pFASTA;
		const SCoalesceModifiers*			m_psModifiers;
	};

	static void* ThreadCombineMotif( void* );
	static void Normalize( CPCL& );

	CCoalesceImpl( ) : m_iK(7), m_dPValueCorrelation(0.05f), m_iBins(12), m_dZScoreCondition(0.5f),
		m_dProbabilityGene(0.95f), m_dZScoreMotif(0.5f), m_pMotifs(NULL), m_fMotifs(false),
		m_iBasesPerMatch(5000), m_dPValueMerge(0.05f), m_dCutoffMerge(2.5f), m_iSizeMinimum(5),
		m_iThreads(1), m_iSizeMerge(100), m_iSizeMaximum(1000), m_dPValueCondition(0.05f),
		m_dPValueMotif(0.05f), m_fNormalize(false) { }
	virtual ~CCoalesceImpl( );

	void Clear( );
	size_t GetMotifCount( ) const;
	bool CombineMotifs( const CFASTA&, const std::vector<size_t>&, SCoalesceModifiers&,
		const CCoalesceCluster&, size_t, CCoalesceGeneScores&, CCoalesceGroupHistograms&,
		CCoalesceGroupHistograms& ) const;
	bool InitializeDatasets( const CPCL& );
	bool InitializeGeneScores( const CPCL&, const CFASTA&, std::vector<size_t>&, SCoalesceModifiers&,
		CCoalesceGeneScores& );

	float							m_dPValueMerge;
	float							m_dProbabilityGene;
	float							m_dPValueCondition;
	float							m_dZScoreCondition;
	float							m_dPValueCorrelation;
	float							m_dPValueMotif;
	float							m_dZScoreMotif;
	float							m_dCutoffMerge;
	size_t							m_iNumberCorrelation;
	size_t							m_iBins;
	size_t							m_iK;
	size_t							m_iSizeMerge;
	size_t							m_iSizeMinimum;
	size_t							m_iSizeMaximum;
	size_t							m_iThreads;
	std::string						m_strDirectoryIntermediate;
	CCoalesceMotifLibrary*			m_pMotifs;
	bool							m_fMotifs;
	bool							m_fNormalize;
	size_t							m_iBasesPerMatch;
	std::vector<SCoalesceDataset>	m_vecsDatasets;
	std::vector<const CFASTA*>		m_vecpWiggles;
	std::vector<std::ostream*>		m_vecpostm;
};

}

#endif // COALESCEI_H

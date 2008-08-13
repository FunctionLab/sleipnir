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
#ifdef _MSC_VER
#include <hash_map>
#else // _MSC_VER
#include <ext/hash_map>

#define hash_value	hash<const char*>( )
#define stdext		__gnu_cxx
#endif // _MSC_VER
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "statistics.h"
#include "typesi.h"

namespace Sleipnir {

class CCoalesceCluster;
class CCoalesceMotifLibrary;
class CFASTA;
class CPCL;
struct SFASTASequence;

class CCoalesceMotifLibraryImpl {
protected:
	static const char	c_acBases[];
	static const size_t	c_iShift		= 2; // ceil( log2( ARRAYSIZE(c_acBases) ) )

	CCoalesceMotifLibraryImpl( size_t iK ) : m_iK(iK) { }

	size_t	m_iK;
};

template<class tValue = size_t, class tCount = unsigned short>
class CCoalesceHistogramSet {
public:
	CCoalesceHistogramSet( ) : m_iMembers(0) { }

	double PoissonTest( size_t iMember, const CCoalesceHistogramSet& HistSet ) const {
		double	dAveOne, dVarOne, dAveTwo, dVarTwo;

		if( !GetEdges( ) || ( GetEdges( ) != HistSet.GetEdges( ) ) )
			return 1;
		GetAveVar( iMember, dAveOne, dVarOne );
		HistSet.GetAveVar( iMember, dAveTwo, dVarTwo );

// Fake a Skellam with a normal
		return CStatistics::TTestWelch( dAveOne, dVarOne, GetTotal( iMember ), dAveTwo, dVarTwo,
			HistSet.GetTotal( iMember ) ); }
// The real thing isn't sensitive enough, or I'm doing it wrong
//		return CStatistics::SkellamPDF( 0, dAveOne, dAveTwo ); }

	double KSTest( size_t iMember, const CCoalesceHistogramSet& HistSet ) const {
		size_t	i;
		tCount	CumOne, CumTwo, SumOne, SumTwo;
		float	d, dMax;

		if( !GetEdges( ) || ( GetEdges( ) != HistSet.GetEdges( ) ) )
			return 1;
		SumOne = GetTotal( iMember );
		SumTwo = HistSet.GetTotal( iMember );
		CumOne = CumTwo = 0;
		for( dMax = 0,i = 0; i < GetEdges( ); ++i ) {
			CumOne += Get( iMember, i );
			CumTwo += HistSet.Get( iMember, i );
			if( ( d = fabs( ( (float)CumOne / SumOne ) - ( (float)CumTwo / SumTwo ) ) ) > dMax )
				dMax = d; }

		return CStatistics::PValueKolmogorovSmirnov( dMax, SumOne, SumTwo ); }

	double Chi2Test( size_t iMember, const CCoalesceHistogramSet& HistSet ) const {
		size_t	i, iDF;
		double	dC2, dTmp, dOT, dTO;
		tCount	One, Two;

		if( !GetEdges( ) || ( GetEdges( ) != HistSet.GetEdges( ) ) )
			return 1;
		iDF = GetEdges( ) - 1;
		dOT = sqrt( (double)GetTotal( iMember ) / HistSet.GetTotal( iMember ) );
		dTO = sqrt( (double)HistSet.GetTotal( iMember ) / GetTotal( iMember ) );
		for( dC2 = 0,i = 0; i < GetEdges( ); ++i ) {
			One = Get( iMember, i );
			Two = HistSet.Get( iMember, i );
			if( One || Two ) {
				dTmp = ( dTO * One ) - ( dOT * Two );
				dC2 += dTmp * dTmp / ( One + Two ); } }

		return ( 1 - CStatistics::Chi2CDF( dC2, iDF ) ); }

	void Initialize( size_t iMembers, const std::vector<tValue>& vecEdges ) {

		m_iMembers = iMembers;
		m_vecEdges.resize( vecEdges.size( ) );
		std::copy( vecEdges.begin( ), vecEdges.end( ), m_vecEdges.begin( ) ); }

	void Initialize( size_t iMembers, size_t iBins ) {
		std::vector<tValue>	vecEdges;
		size_t				i;

		vecEdges.resize( iBins );
		for( i = 0; i < vecEdges.size( ); ++i )
			vecEdges[ i ] = i;

		Initialize( iMembers, vecEdges ); }

	bool Add( size_t iMember, tValue Value, tCount Count ) {
		size_t	i;

		if( m_vecEdges.empty( ) || ( iMember >= m_iMembers ) )
			return false;
		for( i = 0; i < m_vecEdges.size( ); ++i )
			if( Value <= m_vecEdges[ i ] )
				break;
		i = min( i, GetEdges( ) - 1 );

		if( m_vecCounts.empty( ) )
			m_vecCounts.resize( GetOffset( m_iMembers ) + GetEdges( ) );
		m_vecCounts[ GetOffset( iMember ) + i ] += Count;
		if( m_vecTotal.empty( ) )
			m_vecTotal.resize( m_iMembers );
		m_vecTotal[ iMember ] += Count;

		return true; }

	tCount Get( size_t iMember, size_t iBin ) const {

		return ( ( ( iMember < m_iMembers ) && ( iBin < GetEdges( ) ) ) ?
			m_vecCounts[ GetOffset( iMember ) + iBin ] : 0 ); }

	const tCount* Get( size_t iMember ) const {

		return ( ( ( iMember < m_iMembers ) && GetEdges( ) ) ? &m_vecCounts[ GetOffset( iMember ) ] : NULL ); }

	size_t GetMembers( ) const {

		return m_iMembers; }

	size_t GetEdges( ) const {

		return m_vecEdges.size( ); }

	tValue GetEdge( size_t iBin ) {

		return m_vecEdges[ iBin ]; }

	bool Add( const CCoalesceHistogramSet& Histograms ) {
		size_t	i, j;

		if( ( Histograms.GetMembers( ) != GetMembers( ) ) ||
			( Histograms.GetEdges( ).size( ) != GetEdges( ).size( ) ) )
			return false;

		for( i = 0; i < GetMembers( ); ++i )
			for( j = 0; j < GetEdges( ).size( ); ++j )
				if( !Add( i, GetEdges( )[ j ], Histograms.Get( i, j ) ) )
					return false;

		return true; }

	tCount GetTotal( size_t iMember ) const {

		return ( ( iMember < m_vecTotal.size( ) ) ? m_vecTotal[ iMember ] : 0 ); }

	std::string Save( size_t iMember ) const {
		std::ostringstream	ossm;
		size_t				i;

		if( GetEdges( ) ) {
			ossm << (size_t)Get( iMember, 0 );
			for( i = 1; i < GetEdges( ); ++i )
				ossm << '\t' << (size_t)Get( iMember, i ); }

		return ossm.str( ); }

	void Clear( ) {

		fill( m_vecCounts.begin( ), m_vecCounts.end( ), 0 );
		fill( m_vecTotal.begin( ), m_vecTotal.end( ), 0 ); }

protected:
	size_t GetOffset( size_t iMember ) const {

		return ( iMember * GetEdges( ) ); }

	void GetAveVar( size_t iMember, double& dAve, double& dVar ) const {
		size_t	i;
		tCount	Cur, Ave, Var;

		for( Ave = Var = 0,i = 0; i < GetEdges( ); ++i ) {
			Cur = Get( iMember, i ) * i;
			Ave += Cur;
			Var += Cur * Cur; }

		dAve = (double)Ave / GetTotal( iMember );
		dVar = ( (double)Var / GetTotal( iMember ) ) - ( dAve * dAve ); }

	size_t				m_iMembers;
	std::vector<tCount>	m_vecTotal;
	std::vector<tValue>	m_vecEdges;
	std::vector<tCount>	m_vecCounts;
};

class CCoalesceSequencerBase {
public:
	enum ESubsequence {
		ESubsequenceBegin	= 0,
		ESubsequenceTotal	= ESubsequenceBegin,
		ESubsequenceIntrons	= ESubsequenceTotal + 1,
		ESubsequenceExons	= ESubsequenceIntrons + 1,
		ESubsequenceEnd		= ESubsequenceExons + 1
	};

	static const char* GetSubsequence( ESubsequence eSubsequence ) {
		static const char*	c_aszSubsequences[]	= {"total", "introns", "exons"};

		return c_aszSubsequences[ eSubsequence ]; }

	size_t GetTypes( ) const {

		return m_vecstrTypes.size( ); }

	const std::string& GetType( size_t iType ) const {

		return m_vecstrTypes[ iType ]; }

	size_t GetType( const std::string& strType ) const {
		TMapStrI::const_iterator	iterType;

		return ( ( ( iterType = m_mapstriTypes.find( strType ) ) == m_mapstriTypes.end( ) ) ? -1 :
			iterType->second ); }

protected:
	typedef std::map<std::string, size_t>	TMapStrI;

	TMapStrI					m_mapstriTypes;
	std::vector<std::string>	m_vecstrTypes;
};

template<class tType>
class CCoalesceSequencer : public CCoalesceSequencerBase {
public:
	tType& Get( size_t iType, ESubsequence eSubsequence ) {

		return m_vecvecValues[ iType ][ eSubsequence ]; }

	const tType& Get( size_t iType, ESubsequence eSubsequence ) const {

		return m_vecvecValues[ iType ][ eSubsequence ]; }

	const tType& Get( const std::string& strType, ESubsequence eSubsequence ) const {

		return Get( GetType( strType ), eSubsequence ); }

	size_t AddType( const std::string& strType ) {
		TMapStrI::const_iterator	iterType;
		size_t						iRet;

		if( ( iterType = m_mapstriTypes.find( strType ) ) != m_mapstriTypes.end( ) )
			return iterType->second;

		m_mapstriTypes[ strType ] = iRet = m_vecvecValues.size( );
		m_vecstrTypes.push_back( strType );
		m_vecvecValues.push_back( std::vector<tType>( ) );
		m_vecvecValues.back( ).resize( ESubsequenceEnd );

		return iRet; }

	size_t GetSubsequences( size_t iType ) const {

		return m_vecvecValues[ iType ].size( ); }

	void Clear( ) {
		size_t	i, j;

		for( i = 0; i < m_vecvecValues.size( ); ++i )
			for( j = 0; j < m_vecvecValues[ i ].size( ); ++j )
				m_vecvecValues[ i ][ j ].Clear( ); }

protected:
// Type by subsequence
	std::vector<std::vector<tType> >	m_vecvecValues;
};

// One score per motif atom
class CCoalesceGeneScores : public CCoalesceSequencer<std::vector<unsigned short> > {
protected:
	typedef std::vector<unsigned short>		TVecS;
	typedef std::map<uint32_t, uint32_t>	TMapII;

public:
	bool Add( CCoalesceMotifLibrary&, const SFASTASequence& );

	unsigned short Get( size_t iType, size_t iSubsequence, uint32_t iMotif ) const {
		const TMapII&			mapiiMotifs	= m_MotifMaps.Get( iType, (ESubsequence)iSubsequence );
		TMapII::const_iterator	iterMotif;

		return ( ( ( iterMotif = mapiiMotifs.find( iMotif ) ) == mapiiMotifs.end( ) ) ? 0 :
			CCoalesceSequencer<TVecS>::Get( iType, (ESubsequence)iSubsequence )[ iterMotif->second ] ); }

	bool IsEmpty( size_t iType, size_t iSubsequence ) const {

		return m_MotifMaps.Get( iType, (ESubsequence)iSubsequence ).empty( ); }

protected:
	bool Add( CCoalesceMotifLibrary&, const std::string&, size_t, bool );

	uint32_t AddMotif( size_t iType, ESubsequence eSubsequence, uint32_t iMotif ) {
		TMapII::const_iterator	iterMotif;
		size_t					i;
		uint32_t				iRet;

		for( i = m_MotifMaps.GetTypes( ); m_MotifMaps.GetTypes( ) <= iType; ++i )
			m_MotifMaps.AddType( m_vecstrTypes[ i ] );
		{
			TMapII&					mapiiMotifs	= m_MotifMaps.Get( iType, eSubsequence );

			if( ( iterMotif = mapiiMotifs.find( iMotif ) ) != mapiiMotifs.end( ) )
				return iterMotif->second;

			mapiiMotifs[ iMotif ] = iRet = (uint32_t)mapiiMotifs.size( );
			CCoalesceSequencer<TVecS>::Get( iType, eSubsequence ).push_back( 0 );
			return iRet;
		} }

	CCoalesceSequencer<TMapII>	m_MotifMaps;
};

// One histogram per motif atom
class CCoalesceGroupHistograms : public CCoalesceSequencer<CCoalesceHistogramSet<> > {
public:
	CCoalesceGroupHistograms( size_t iMotifs, size_t iBins ) : m_iMotifs(iMotifs), m_iBins(iBins) { }

	bool Add( const CCoalesceGeneScores&, bool );

	size_t GetMotifs( ) const {

		return m_iMotifs; }

protected:
	size_t	m_iMotifs;
	size_t	m_iBins;
};

struct SMotifMatch {
	SMotifMatch( ) { }

	SMotifMatch( uint32_t iMotif, const std::string& strType,
		CCoalesceSequencerBase::ESubsequence eSubsequence ) : m_iMotif(iMotif), m_strType(strType),
		m_eSubsequence(eSubsequence) { }

	std::string Save( const CCoalesceMotifLibrary* ) const;

	bool operator==( const SMotifMatch& sMotif ) const {

		return ( ( m_iMotif == sMotif.m_iMotif ) && ( m_strType == sMotif.m_strType ) &&
			( m_eSubsequence == sMotif.m_eSubsequence ) ); }

	bool operator!=( const SMotifMatch& sMotif ) const {

		return !( *this == sMotif ); }

	bool operator<( const SMotifMatch& sMotif ) const {

		if( m_iMotif == sMotif.m_iMotif ) {
			if( m_strType == sMotif.m_strType )
				return ( m_eSubsequence < sMotif.m_eSubsequence );
			return ( m_strType < sMotif.m_strType ); }

		return ( m_iMotif < sMotif.m_iMotif ); }

	size_t GetHash( ) const {
		size_t	iMotif, iType, iSubsequence;

		iMotif = m_iMotif * ( (size_t)-1 / 20000 );
		iType = stdext::hash_value( m_strType.c_str( ) ) * ( (size_t)-1 / 5 );
		iSubsequence = m_eSubsequence * ( (size_t)-1 / CCoalesceSequencerBase::ESubsequenceEnd );

		return ( iMotif ^ iType ^ iSubsequence ); }

	uint32_t								m_iMotif;
	std::string								m_strType;
	CCoalesceSequencerBase::ESubsequence	m_eSubsequence;
};

class CCoalesceClusterImpl {
protected:
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
	bool AddSeedPair( const CPCL&, CCoalesceCluster&, float );
	void CalculateCentroid( const CPCL& );
	void AddSignificant( uint32_t, const CCoalesceMotifLibrary*, const CCoalesceGroupHistograms&,
		const CCoalesceGroupHistograms&, float );
	bool IsSignificant( size_t, const CPCL&, const CCoalesceMotifLibrary*, const CCoalesceGeneScores&,
		const CCoalesceGroupHistograms&, const CCoalesceGroupHistograms&, const CCoalesceCluster&,
		float ) const;
	float CalculateProbabilityExpression( size_t, const CPCL&, const CCoalesceCluster&, bool ) const;
	float CalculateProbabilityMotifs( const CCoalesceGeneScores&, const CCoalesceGroupHistograms&,
		const CCoalesceGroupHistograms&, bool ) const;
	bool SaveCopy( const CPCL&, size_t, CPCL&, size_t, bool ) const;

	bool IsGene( size_t iGene ) const {

		return ( m_setiGenes.find( iGene ) != m_setiGenes.end( ) ); }

	bool IsCondition( size_t iCondition ) const {

		return ( m_setiConditions.find( iCondition ) != m_setiConditions.end( ) ); }

	size_t GetHash( ) const {

		return ( GetHash( m_setiConditions ) ^ GetHash( m_setiGenes ) ^ GetHash( m_setsMotifs ) ); }

	std::set<size_t>					m_setiConditions;
	std::set<size_t>					m_setiGenes;
	std::set<SMotifMatch>				m_setsMotifs;
	std::vector<size_t>					m_veciPrevConditions;
	std::vector<size_t>					m_veciPrevGenes;
	std::vector<SMotifMatch>			m_vecsPrevMotifs;
	std::vector<size_t>					m_veciCounts;
	std::vector<float>					m_vecdCentroid;
	std::vector<float>					m_vecdStdevs;
	mutable std::vector<long double>	m_vecdPIn;
	mutable std::vector<long double>	m_vecdPOut;
	std::set<size_t>					m_setiHistory;
};

class CCoalesceImpl {
protected:
	CCoalesceImpl( ) : m_iK(7), m_dPValueCorrelation(0.05f), m_iBins(12), m_dPValueCondition(0.05f),
		m_dProbabilityGene(0.95f), m_dPValueMotif(0.05f), m_pMotifs(NULL), m_fMotifs(false) { }
	virtual ~CCoalesceImpl( );

	void Clear( );
	size_t GetMotifCount( ) const;

	float					m_dProbabilityGene;
	float					m_dPValueCondition;
	float					m_dPValueCorrelation;
	float					m_dPValueMotif;
	size_t					m_iBins;
	size_t					m_iK;
	std::string				m_strDirectoryIntermediate;
	CCoalesceMotifLibrary*	m_pMotifs;
	bool					m_fMotifs;
};

}

#endif // COALESCEI_H

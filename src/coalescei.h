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
#include <fstream>
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
#include <utility>
#include <vector>

#include "pcl.h"
#include "statistics.h"
#include "typesi.h"

namespace Sleipnir {

class CCoalesceCluster;
class CCoalesceMotifLibrary;
class CFASTA;
class CPST;
struct SFASTASequence;
struct SFASTAWiggle;
struct SMotifMatch;

struct SCoalesceModifiers {

	void Initialize( const CPCL& );

	bool Add( const CFASTA* pWiggle ) {

		if( !pWiggle )
			return false;
		m_vecpWiggles.push_back( pWiggle );
		return true; }

	std::vector<const CFASTA*>			m_vecpWiggles;
	std::vector<std::vector<size_t> >	m_vecveciPCL2Wiggles;
};

struct SCoalesceModifierCache {
	SCoalesceModifierCache( const SCoalesceModifiers& Modifiers ) : m_Modifiers(Modifiers) { }

	void Get( size_t );
	void InitializeWeight( size_t, size_t );
	void AddWeight( size_t, size_t, size_t );
	void SetType( const std::string& );

	float GetWeight( ) const {
		size_t	i, iWiggles;

		for( iWiggles = i = 0; i < m_veciWiggleTypes.size( ); ++i )
			if( m_veciWiggleTypes[ i ] != -1 )
				iWiggles++;

		return ( iWiggles ? ( m_dWeight / iWiggles ) : 1 ); }

	const SCoalesceModifiers&				m_Modifiers;
	std::vector<std::vector<SFASTAWiggle> >	m_vecvecsWiggles;
	std::vector<size_t>						m_veciWiggleTypes;
	float									m_dWeight;
};

class CCoalesceMotifLibraryImpl {
protected:
	static const char	c_acBases[];
	static const size_t	c_iShift		= 2; // ceil( log2( ARRAYSIZE(c_acBases) ) )
	static const char	c_cSeparator	= '|';

	enum EType {
		ETypeKMer,
		ETypeRC,
		ETypePST
	};

	static size_t CountKMers( size_t iK ) {

		return ( 1 << ( iK << 1 ) ); }

	static size_t CountRCs( size_t iK ) {
		size_t	iKMers;

		iKMers = CountKMers( iK ) >> 1;
		return ( iKMers - ( ( iK % 2 ) ? 0 : ( CountKMers( iK >> 1 ) >> 1 ) ) ); }

	static std::string GetReverseComplement( const std::string& strKMer ) {
		std::string	strReverse;

		strReverse = strKMer;
		std::reverse( strReverse.begin( ), strReverse.end( ) );
		return GetComplement( strReverse ); }

	static std::string GetComplement( const std::string& strKMer ) {
		std::string	strRet;
		size_t		i;

		strRet.resize( strKMer.length( ) );
		for( i = 0; i < strRet.length( ); ++i )
			strRet[ i ] = GetComplement( strKMer[ i ] );

		return strRet; }

	static char GetComplement( char cBase ) {
		const char*	pc;
		size_t		i;

		if( !( pc = strchr( c_acBases, cBase ) ) )
			return cBase;
		i = pc - c_acBases;

		return c_acBases[ ( i & ~1 ) + ( 1 - ( i & 1 ) ) ]; }

	static uint32_t KMer2ID( const std::string& strKMer, bool fRC = false ) {
		size_t			i, iIndex;
		uint32_t		iRet;
		const char*		pc;
		unsigned char	c;

		for( i = iRet = 0; i < strKMer.length( ); ++i ) {
			iIndex = fRC ? ( strKMer.length( ) - i - 1 ) : i;
			if( !( pc = strchr( c_acBases, strKMer[ iIndex ] ) ) )
				return -1;
			c = (unsigned char)( pc - c_acBases );
			if( fRC )
				c = ( c & ~1 ) | ( 1 - ( c & 1 ) );
			iRet = ( iRet << c_iShift ) | c; }

		return iRet; }

	static std::string ID2KMer( uint32_t iID, size_t iK ) {
		static const size_t	c_iMask	= ( 1 << c_iShift ) - 1;
		std::string	strRet;
		size_t		i;

		strRet.resize( iK );
		for( i = 0; i < iK; ++i ) {
			strRet[ iK - i - 1 ] = c_acBases[ iID & c_iMask ];
			iID >>= c_iShift; }

		return strRet; }

	static bool IsIgnorableKMer( const std::string& strKMer ) {
		size_t	i;

		for( i = 0; i < strKMer.size( ); ++i )
			if( !strchr( c_acBases, strKMer[ i ] ) )
				return true;

		return false; }

	CCoalesceMotifLibraryImpl( size_t iK ) : m_iK(iK), m_dPenaltyGap(1), m_dPenaltyMismatch(2) {
		uint32_t	i, j, iRC;

// TODO: if I was smart, I could do this with a direct encoding...
		m_vecKMer2RC.resize( GetKMers( ) );
		m_vecRC2KMer.resize( GetRCs( ) );
		std::fill( m_vecKMer2RC.begin( ), m_vecKMer2RC.end( ), -1 );
		for( i = j = 0; i < m_vecKMer2RC.size( ); ++i )
			if( m_vecKMer2RC[ i ] == -1 ) {
				iRC = KMer2ID( ID2KMer( i, m_iK ), true );
				if( iRC != i ) {
					m_vecKMer2RC[ i ] = m_vecKMer2RC[ iRC ] = j;
					m_vecRC2KMer[ j++ ] = i; } } }

	virtual ~CCoalesceMotifLibraryImpl( );

	std::string GetMotif( uint32_t ) const;
	CPST* CreatePST( uint32_t& );
	uint32_t MergeKMers( const std::string&, const std::string&, float );
	uint32_t MergeKMerRC( const std::string&, uint32_t, float );
	uint32_t MergeKMerPST( const std::string&, const CPST&, float );
	uint32_t MergeRCs( uint32_t, uint32_t, float );
	uint32_t MergeRCPST( uint32_t, const CPST&, float );
	uint32_t MergePSTs( const CPST&, const CPST&, float );

	EType GetType( uint32_t iMotif ) const {

		if( iMotif < GetKMers( ) )
			return ETypeKMer;
		if( iMotif < GetBasePSTs( ) )
			return ETypeRC;

		return ETypePST; }

	uint32_t GetKMers( ) const {

		return (uint32_t)CountKMers( m_iK ); }

	uint32_t GetRCs( ) const {

			return (uint32_t)CountRCs( m_iK ); }

	uint32_t GetBaseRCs( ) const {

		return GetKMers( ); }

	uint32_t GetBasePSTs( ) const {

		return ( GetBaseRCs( ) + GetRCs( ) ); }

	const CPST* GetPST( uint32_t iMotif ) const {

		return m_vecpPSTs[ iMotif - GetBasePSTs( ) ]; }

	uint32_t GetPSTs( ) const {

		return (uint32_t)m_vecpPSTs.size( ); }

	float Align( const std::string& strFixed, const std::string& strMobile, float dCutoff, int& iRet ) const {
		const std::string&	strShort	= ( strFixed.length( ) < strMobile.length( ) ) ? strFixed : strMobile;
		const std::string&	strLong		= ( strFixed.length( ) < strMobile.length( ) ) ? strMobile : strFixed;
		size_t				iBegin, iEnd, iOffset, i, iLength, iMin;
		float				dRet, dCur;

		dRet = FLT_MAX;
		dCur = dCutoff / m_dPenaltyGap;
		iLength = strShort.length( ) + strLong.length( );
		iBegin = ( dCur < iLength ) ? (size_t)ceil( ( iLength - dCur ) / 2 ) : 0;
		iEnd = iLength + 1 - iBegin;
		for( iMin = 0,iOffset = iBegin; iOffset < iEnd; ++iOffset ) {
			i = ( iOffset <= strShort.length( ) ) ? iOffset : min( strShort.length( ), iLength - iOffset );
			dCur = m_dPenaltyGap * ( iLength - ( 2 * i ) );
			for( i = ( iOffset <= strShort.length( ) ) ? 0 : ( iOffset - strShort.length( ) );
				i < min( strShort.length( ), iOffset ); ++i )
				if( strShort[ i ] != strLong[ strLong.length( ) - iOffset + i ] )
					dCur += m_dPenaltyMismatch;
			if( dCur < dRet ) {
				dRet = dCur;
				iMin = iOffset; } }

		iRet = (int)strLong.length( ) - (int)iMin;
		if( strFixed.length( ) < strMobile.length( ) )
			iRet = -iRet;
		return dRet; }

	std::string GetRCOne( uint32_t iMotif ) const {

		return ID2KMer( (uint32_t)m_vecRC2KMer[ iMotif - GetBaseRCs( ) ], m_iK ); }

	float										m_dPenaltyGap;
	float										m_dPenaltyMismatch;
	size_t										m_iK;
	std::vector<uint32_t>						m_vecKMer2RC;
	std::vector<uint32_t>						m_vecRC2KMer;
	std::vector<CPST*>							m_vecpPSTs;
	std::set<std::pair<uint32_t, uint32_t> >	m_setpriiMerged;
};

template<class tValue = float, class tCount = unsigned short>
class CCoalesceHistogramSet {
public:
	CCoalesceHistogramSet( ) : m_iMembers(0) { }

	double PoissonTest( size_t iMember, const CCoalesceHistogramSet& HistSet ) const {
		double	dAveOne, dVarOne, dAveTwo, dVarTwo;

		if( !GetEdges( ) || ( GetEdges( ) != HistSet.GetEdges( ) ) )
			return 1;
		AveVar( iMember, dAveOne, dVarOne );
		HistSet.AveVar( iMember, dAveTwo, dVarTwo );

// Fake a Skellam with a normal; Student's T tends to be slightly too sensitive
		return CStatistics::TTestWelch( dAveOne, dVarOne, GetTotal( ), dAveTwo, dVarTwo,
			HistSet.GetTotal( ) ); }
// The real thing isn't sensitive enough, or I'm doing it wrong
//		return CStatistics::SkellamPDF( 0, dAveOne, dAveTwo ); }

	double ZTest( size_t iMember, const CCoalesceHistogramSet& HistSet ) const {
		double	dAverage, dZ;

		ZScore( iMember, HistSet, dAverage, dZ );
		return CStatistics::ZTest( dZ, GetTotal( ) ); }

	double ZScore( size_t iMember, const CCoalesceHistogramSet& HistSet, double& dAverage, double& dZ ) const {
		tValue	AveOne, StdOne, Ave, Std;
		double	dAveOne, dAve, dStd;

		if( !GetEdges( ) || ( GetEdges( ) != HistSet.GetEdges( ) ) )
			return 1;

		Sums( iMember, AveOne, StdOne );
		HistSet.Sums( iMember, Ave, Std );
		dAve = (double)( AveOne + Ave ) / ( GetTotal( ) + HistSet.GetTotal( ) );
		dAveOne = (double)AveOne / GetTotal( );
		dStd = sqrt( ( (double)( StdOne + Std ) / ( GetTotal( ) + HistSet.GetTotal( ) ) ) - ( dAve * dAve ) );

		dAverage = dAveOne;
		dZ = dStd ? ( ( dAveOne - dAve ) / dStd ) : 0;
		return ( dStd ? CStatistics::ZTest( dZ, GetTotal( ) ) : ( ( dAveOne == dAve ) ? 1 : 0 ) ); }

	double CohensD( size_t iMember, const CCoalesceHistogramSet& HistSet, double& dAverage, double& dZ ) const {

		return CohensD( iMember, HistSet, iMember, true, dAverage, dZ ); }

	double CohensD( size_t iOne, const CCoalesceHistogramSet& HistSet, size_t iTwo, bool fCount,
		double& dAverage, double& dZ ) const {
		static const double	c_dEpsilon	= 1e-6;
		tValue	AveOne, VarOne, AveTwo, VarTwo;
		double	dAveOne, dAveTwo, dVarOne, dVarTwo, dAve, dStd;

		if( !GetEdges( ) || ( GetEdges( ) != HistSet.GetEdges( ) ) )
			return 1;

		Sums( iOne, AveOne, VarOne );
		HistSet.Sums( iTwo, AveTwo, VarTwo );
		dAve = (double)( AveOne + AveTwo ) / ( GetTotal( ) + HistSet.GetTotal( ) );
		dAveOne = (double)AveOne / GetTotal( );
		dAveTwo = (double)AveTwo / HistSet.GetTotal( );
		dVarOne = max( 0.0, ( (double)VarOne / GetTotal( ) ) - ( dAveOne * dAveOne ) );
		dVarTwo = max( 0.0, ( (double)VarTwo / HistSet.GetTotal( ) ) - ( dAveTwo * dAveTwo ) );
		dStd = sqrt( ( dVarOne + dVarTwo ) / 2 );

		dAverage = dAveOne;
		dZ = dStd ? ( ( dAveOne - dAveTwo ) / dStd ) : 0;
// dZ *= 1 - pow( (float)GetTotal( ) / ( GetTotal( ) + HistSet.GetTotal( ) ), 1 ); // doesn't work
// dZ *= exp( -(float)GetTotal( ) / ( GetTotal( ) + HistSet.GetTotal( ) ) ); // doesn't work
// dZ *= exp( -(float)GetTotal( ) / HistSet.GetTotal( ) ); // works pretty well
// dZ *= fabs( (float)( GetTotal( ) - HistSet.GetTotal( ) ) ) / max( GetTotal( ), HistSet.GetTotal( ) ); // works pretty well
// This prevents large clusters from blowing up the motif set
		if( iOne == iTwo )
			dZ *= pow( fabs( (float)( GetTotal( ) - HistSet.GetTotal( ) ) ) /
				max( GetTotal( ), HistSet.GetTotal( ) ), max( 1.0f,
				log10f( (float)min( GetTotal( ), HistSet.GetTotal( ) ) ) ) );

		return ( ( dStd > c_dEpsilon ) ? ( 2 * CStatistics::ZTest( dZ, fCount ?
			min( GetTotal( ), HistSet.GetTotal( ) ) : 1 ) ) :
			( ( fabs( dAveOne - dAveTwo ) < c_dEpsilon ) ? 1 : 0 ) ); }

	double KSTest( size_t iMember, const CCoalesceHistogramSet& HistSet ) const {
		size_t	i;
		tCount	CumOne, CumTwo, SumOne, SumTwo;
		float	d, dMax;

		if( !GetEdges( ) || ( GetEdges( ) != HistSet.GetEdges( ) ) )
			return 1;
		SumOne = GetTotal( );
		SumTwo = HistSet.GetTotal( );
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
		if( !( GetTotal( ) && HistSet.GetTotal( ) ) )
			return ( ( GetTotal( ) == HistSet.GetTotal( ) ) ? 1 : 0 );
		iDF = GetEdges( ) - 1;
		dOT = sqrt( (double)GetTotal( ) / HistSet.GetTotal( ) );
		dTO = sqrt( (double)HistSet.GetTotal( ) / GetTotal( ) );
		for( dC2 = 0,i = 0; i < GetEdges( ); ++i ) {
			One = Get( iMember, i );
			Two = HistSet.Get( iMember, i );
			if( One || Two ) {
				dTmp = ( dTO * One ) - ( dOT * Two );
				dC2 += dTmp * dTmp / ( One + Two ); }
			else
				iDF--; }

		return ( 1 - CStatistics::Chi2CDF( dC2, iDF ) ); }

	void Initialize( size_t iMembers, const std::vector<tValue>& vecEdges ) {

		m_iMembers = iMembers;
		m_vecEdges.resize( vecEdges.size( ) );
		std::copy( vecEdges.begin( ), vecEdges.end( ), m_vecEdges.begin( ) );
		m_iZero = GetBin( 0 ); }

	void Initialize( size_t iMembers, size_t iBins, tValue Step ) {
		std::vector<tValue>	vecEdges;
		size_t				i;

		vecEdges.resize( iBins );
		for( i = 0; i < vecEdges.size( ); ++i )
			vecEdges[ i ] = i * Step;

		Initialize( iMembers, vecEdges ); }

	bool Add( size_t iMember, tValue Value, tCount Count ) {
		size_t	i;

// lock
		if( m_vecEdges.empty( ) || ( ( i = GetBin( Value ) ) == m_iZero ) )
// unlock
			return false;

		m_iMembers = max( m_iMembers, iMember + 1 );
		m_vecCounts.resize( GetOffset( m_iMembers ) );
		m_vecCounts[ GetOffset( iMember ) + i ] += Count;
		m_vecTotal.resize( m_iMembers );
		m_vecTotal[ iMember ] += Count;
// unlock

		return true; }

	tCount Integrate( size_t iMember, tValue Value, bool fUp ) const {
		tCount	Ret;
		size_t	iBin;

		for( Ret = 0,iBin = GetBin( Value ); iBin < GetEdges( ); iBin += ( fUp ? 1 : -1 ) )
			Ret += Get( iMember, iBin );

		return Ret; }

	tCount Get( size_t iMember, size_t iBin ) const {
		size_t	iOffset;

		if( ( iMember >= m_iMembers ) || ( iBin >= GetEdges( ) ) )
			return 0;

		return ( ( ( iOffset = ( GetOffset( iMember ) + iBin ) ) < m_vecCounts.size( ) ) ?
			( ( iBin == m_iZero ) ? ( GetTotal( ) - m_vecTotal[ iMember ] ) : m_vecCounts[ iOffset ] ) :
			( ( iBin == m_iZero ) ? GetTotal( ) : 0 ) ); }

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

	const tCount* Get( size_t iMember ) const {
		const tCount*	pRet;

		if( !GetEdges( ) || ( iMember >= m_iMembers ) )
			return NULL;

		pRet = &m_vecCounts[ GetOffset( iMember ) ];
		*(tCount*)pRet = Get( iMember, (size_t)0 );
		return pRet; }

	size_t GetMembers( ) const {

		return m_iMembers; }

	size_t GetEdges( ) const {

		return m_vecEdges.size( ); }

	tValue GetEdge( size_t iBin ) const {

		return m_vecEdges[ iBin ]; }

	bool Add( const CCoalesceHistogramSet& Histograms ) {
		size_t	i, j;

		if( ( Histograms.GetMembers( ) != GetMembers( ) ) ||
			( Histograms.GetEdges( ).size( ) != GetEdges( ).size( ) ) )
			return false;

		m_Total += Histograms.GetTotal( );
		for( i = 0; i < GetMembers( ); ++i )
			for( j = 1; j < GetEdges( ).size( ); ++j )
				if( !Add( i, GetEdges( )[ j ], Histograms.Get( i, j ) ) )
					return false;

		return true; }

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
	size_t GetOffset( size_t iMember ) const {

		return ( iMember * GetEdges( ) ); }

	void Sums( size_t iMember, tValue& Sum, tValue& SumSq ) const {
		size_t	i;
		tValue	Cur;

		for( Sum = SumSq = 0,i = 0; i < GetEdges( ); ++i ) {
			Cur = Get( iMember, i ) * GetEdge( i );
			Sum += Cur;
			SumSq += Cur * GetEdge( i ); } }

	void AveVar( size_t iMember, double& dAve, double& dVar ) const {
		tValue	Ave, Var;

		Sums( iMember, Ave, Var );
		dAve = (double)Ave / GetTotal( );
		dVar = ( (double)Var / GetTotal( ) ) - ( dAve * dAve ); }

	size_t				m_iMembers;
	size_t				m_iZero;
	tCount				m_Total;
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

// One score per type per subtype per gene per motif atom
class CCoalesceGeneScores : public CCoalesceSequencer<float**> {
public:
	CCoalesceGeneScores( ) : m_iMotifs(0), m_iGenes(0), m_iCapacity(0) { }

	virtual ~CCoalesceGeneScores( ) {
		size_t	iType, iSubtype, iGene;
		float**	aad;

		for( iType = 0; iType < GetTypes( ); ++iType )
			for( iSubtype = 0; iSubtype < GetSubsequences( iType ); ++iSubtype ) {
				if( !( aad = Get( iType, (ESubsequence)iSubtype ) ) )
					continue;
				for( iGene = 0; iGene < m_iGenes; ++iGene )
					if( aad[ iGene ] )
						delete[] aad[ iGene ];
				delete[] aad; } }

	bool Add( size_t, const CCoalesceMotifLibrary&, const SFASTASequence&, SCoalesceModifierCache&,
		std::vector<std::vector<float> >&, std::vector<size_t>& );
	bool Add( size_t, const CCoalesceMotifLibrary&, const SFASTASequence&, SCoalesceModifierCache&, uint32_t,
		std::vector<float>&, std::vector<size_t>& );
	void Subtract( const SMotifMatch&, size_t );

	const float* Get( size_t iType, ESubsequence eSubsequence, size_t iGene ) const {
		float**	aad;

		return ( ( aad = Get( iType, eSubsequence ) ) ? aad[ iGene ] : NULL ); }

	float Get( size_t iType, ESubsequence eSubsequence, size_t iGene, uint32_t iMotif ) const {
		const float*	ad;

		return ( ( ad = Get( iType, eSubsequence, iGene ) ) ? ad[ iMotif ] : 0 ); }

	bool IsPresent( size_t iType, ESubsequence eSubsequence ) const {

		return !!Get( iType, eSubsequence ); }

	size_t GetMotifs( ) const {

		return m_iMotifs; }

	void SetGenes( size_t iGenes ) {

		m_iGenes = iGenes; }

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

	float** Get( size_t iType, ESubsequence eSubsequence ) const {

		return CCoalesceSequencer<float**>::Get( iType, eSubsequence ); }

	void Set( size_t iType, ESubsequence eSubsequence, size_t iGene, uint32_t iMotif, float dValue,
		uint32_t iMotifs = 0 ) {
		float**	aad;
		float*	ad;

		Grow( iMotif, iMotifs );
		if( !( aad = Get( iType, eSubsequence ) ) ) {
			m_vecvecValues[ iType ][ eSubsequence ] = aad = new float*[ m_iGenes ];
			memset( aad, 0, m_iGenes * sizeof(*aad) ); }
		if( !( ad = aad[ iGene ] ) ) {
			aad[ iGene ] = ad = new float[ m_iCapacity ];
			memset( ad, 0, m_iCapacity * sizeof(*ad) ); }

		ad[ iMotif ] = dValue; }

	void Grow( uint32_t iMotif, uint32_t iMotifs ) {
		size_t	iType, iSubtype, iGene, iTarget;
		float**	aad;
		float*	ad;

		iTarget = max( iMotif + 1, iMotifs );
		if( iTarget < m_iCapacity ) {
			if( iTarget > m_iMotifs )
				m_iMotifs = iTarget;
			return; }
		m_iCapacity = iTarget + c_iLookahead;
		for( iType = 0; iType < GetTypes( ); ++iType )
			for( iSubtype = 0; iSubtype < GetSubsequences( iType ); ++iSubtype ) {
				if( !( aad = Get( iType, (ESubsequence)iSubtype ) ) )
					continue;
				for( iGene = 0; iGene < m_iGenes; ++iGene ) {
					if( !aad[ iGene ] )
						continue;
					ad = new float[ m_iCapacity ];
					memcpy( ad, aad[ iGene ], m_iMotifs * sizeof(*ad) );
					memset( ad + m_iMotifs, 0, ( m_iCapacity - m_iMotifs ) * sizeof(*ad) );
					delete[] aad[ iGene ];
					aad[ iGene ] = ad; } }
		if( iTarget > m_iMotifs )
			m_iMotifs = iTarget; }

	size_t	m_iGenes;
	size_t	m_iMotifs;
	size_t	m_iCapacity;
};

// One histogram per motif atom
class CCoalesceGroupHistograms : public CCoalesceSequencer<CCoalesceHistogramSet<> > {
public:
	CCoalesceGroupHistograms( uint32_t iMotifs, size_t iBins, float dStep ) : m_iMotifs(iMotifs),
		m_iBins(iBins), m_dStep(dStep) { }

	bool Add( const CCoalesceGeneScores&, size_t, bool, uint32_t = -1 );
	void Save( std::ostream&, const CCoalesceMotifLibrary* ) const;
	bool IsSimilar( const CCoalesceMotifLibrary*, const SMotifMatch&, const SMotifMatch&, float ) const;

	uint32_t GetMotifs( ) const {

		return m_iMotifs; }

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
					if( GeneScores.IsPresent( iTypeThem, (ESubsequence)iSubsequence ) )
						m_vecsTotals[ ( iTypeUs * ESubsequenceEnd ) + iSubsequence ]++; }

		for( i = iTypeUs = 0; iTypeUs < GetTypes( ); ++iTypeUs )
			for( iSubsequence = ESubsequenceBegin; iSubsequence < GetSubsequences( iTypeUs );
				++iSubsequence ) {
//				g_CatSleipnir.info( "CCoalesceGroupHistograms::SetTotal( ) type %s, subsequence %d contains %d genes with sequences",
//					GetType( iTypeUs ).c_str( ), iSubsequence, m_vecsTotals[ i ] );
				Get( iTypeUs, (ESubsequence)iSubsequence ).SetTotal( m_vecsTotals[ i++ ] ); } }

protected:
	uint32_t					m_iMotifs;
	size_t						m_iBins;
	float						m_dStep;
	std::vector<unsigned short>	m_vecsTotals;
};

struct SMotifMatch {
	SMotifMatch( ) { }

	SMotifMatch( uint32_t iMotif, const std::string& strType,
		CCoalesceSequencerBase::ESubsequence eSubsequence, float dZ, float dAverage ) : m_iMotif(iMotif),
		m_strType(strType), m_eSubsequence(eSubsequence), m_dZ(dZ), m_dAverage(dAverage) { }

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
	float									m_dZ;
	float									m_dAverage;
};

struct SCoalesceDataset {
	template<class tType>
	SCoalesceDataset( const tType& Conditions ) {

		m_veciConditions.resize( Conditions.size( ) );
		std::copy( Conditions.begin( ), Conditions.end( ), m_veciConditions.begin( ) ); }

	SCoalesceDataset( size_t iCondition ) {

		m_veciConditions.resize( 1 );
		m_veciConditions[ 0 ] = iCondition; }

	bool CalculateCovariance( const CPCL& );

	bool IsCondition( size_t iCondition ) const {

		return ( std::find( m_veciConditions.begin( ), m_veciConditions.end( ), iCondition ) !=
			m_veciConditions.end( ) ); }

	size_t GetConditions( ) const {

		return m_veciConditions.size( ); }

	size_t GetCondition( size_t iCondition ) const {

		return m_veciConditions[ iCondition ]; }

	std::vector<size_t>	m_veciConditions;
	CDataMatrix			m_MatSigmaChol;
	CDataMatrix			m_MatSigmaInv;
	double				m_dSigmaDetSqrt;
	std::vector<float>	m_vecdStdevs;
};

class CCoalesceClusterImpl {
protected:
	struct SDataset {
		const SCoalesceDataset*	m_psDataset;
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
		std::vector<SMotifMatch>		m_vecsMotifs;
	};

	struct SThreadSeedPair {
		size_t										m_iOffset;
		size_t										m_iStep;
		const CPCL*									m_pPCL;
		float										m_dFraction;
		const std::set<std::pair<size_t, size_t> >*	m_psetpriiSeeds;
		double										m_dMin;
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
		float						m_dPValue;
		std::vector<bool>*			m_pvecfSignificant;
	};

	static void* ThreadCentroid( void* );
	static void* ThreadSignificantGene( void* );
	static void* ThreadSelectMotif( void* );
	static void* ThreadSeedPair( void* );
	static void* ThreadSelectCondition( void* );
	static bool AddSignificant( const CCoalesceMotifLibrary&, uint32_t, const CCoalesceGroupHistograms&,
		const CCoalesceGroupHistograms&, float, std::vector<SMotifMatch>& );

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

	static void AdjustProbabilities( float dZ, double& dPIn, double& dPOut ) {
		double	dPAverage;

		dPAverage = ( dPIn + dPOut ) / 2;
		dZ = fabs( dZ );
		dPIn = ( ( dPIn * dZ ) + dPAverage ) / ( dZ + 1 );
		dPOut = ( ( dPOut * dZ ) + dPAverage ) / ( dZ + 1 ); }

	void Add( size_t, CCoalesceCluster& );
	bool AddCorrelatedGenes( const CPCL&, CCoalesceCluster&, float );
	bool AddSeedPair( const CPCL&, CCoalesceCluster&, std::set<std::pair<size_t, size_t> >&, float, float,
		size_t );
	void CalculateCentroid( const CPCL& );
	bool IsSignificant( size_t, const CPCL&, const CCoalesceMotifLibrary*, const CCoalesceGeneScores&,
		const CCoalesceGroupHistograms&, const CCoalesceGroupHistograms&, const CCoalesceCluster&,
		const std::vector<size_t>&, size_t, float ) const;
	bool CalculateProbabilityExpression( size_t, const CPCL&, const CCoalesceCluster&,
		const std::vector<size_t>&, bool, long double&, long double& ) const;
	bool CalculateProbabilityMotifs( const CCoalesceGeneScores&, size_t, const CCoalesceGroupHistograms&,
		const CCoalesceGroupHistograms&, bool, size_t, long double&, long double& ) const;
	bool SaveCopy( const CPCL&, size_t, CPCL&, size_t, bool ) const;

	const SCoalesceDataset& GetDataset( size_t iDataset ) const {

		return *m_vecsDatasets[ iDataset ].m_psDataset; }

	bool IsGene( size_t iGene ) const {

		return ( m_setiGenes.find( iGene ) != m_setiGenes.end( ) ); }

	size_t GetHash( ) const {

		return ( GetHash( m_setiDatasets ) ^ GetHash( m_setiGenes ) ^ GetHash( m_setsMotifs ) ); }

	void GetConditions( std::set<size_t>& setiConditions ) const {
		set<size_t>::const_iterator	iterDataset;
		size_t						i;

		for( iterDataset = m_setiDatasets.begin( ); iterDataset != m_setiDatasets.end( ); ++iterDataset )
			for( i = 0; i < GetDataset( *iterDataset ).GetConditions( ); ++i )
				setiConditions.insert( GetDataset( *iterDataset ).GetCondition( i ) ); }

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

	CCoalesceImpl( ) : m_iK(7), m_dPValueCorrelation(0.05f), m_iBins(12), m_dPValueCondition(0.05f),
		m_dProbabilityGene(0.95f), m_dPValueMotif(0.05f), m_pMotifs(NULL), m_fMotifs(false),
		m_iBasesPerMatch(5000), m_dPValueMerge(0.05f), m_dCutoffMerge(2.5f), m_dPenaltyGap(1),
		m_dPenaltyMismatch(2.1f), m_iSizeMinimum(5), m_iThreads(1) { }
	virtual ~CCoalesceImpl( );

	void Clear( );
	size_t GetMotifCount( ) const;
	bool CombineMotifs( const CFASTA&, const std::vector<size_t>&, SCoalesceModifiers&,
		const CCoalesceCluster&, size_t, CCoalesceGeneScores&, CCoalesceGroupHistograms&,
		CCoalesceGroupHistograms& ) const;
	bool InitializeDatasets( const CPCL& );
	bool InitializeGeneScores( const CPCL&, const CFASTA&, CPCL&, std::vector<size_t>&, SCoalesceModifiers&,
		CCoalesceGeneScores& );

	float							m_dPValueMerge;
	float							m_dProbabilityGene;
	float							m_dPValueCondition;
	float							m_dPValueCorrelation;
	float							m_dPValueMotif;
	float							m_dCutoffMerge;
	float							m_dPenaltyGap;
	float							m_dPenaltyMismatch;
	size_t							m_iNumberCorrelation;
	size_t							m_iBins;
	size_t							m_iK;
	size_t							m_iSizeMinimum;
	size_t							m_iSizeMaximum;
	size_t							m_iThreads;
	std::string						m_strDirectoryIntermediate;
	CCoalesceMotifLibrary*			m_pMotifs;
	bool							m_fMotifs;
	size_t							m_iBasesPerMatch;
	std::string						m_strSequenceCache;
	std::vector<SCoalesceDataset>	m_vecsDatasets;
	std::vector<const CFASTA*>		m_vecpWiggles;
};

}

#endif // COALESCEI_H

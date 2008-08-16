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

template<class tValue = float, class tCount = unsigned short>
class CCoalesceHistogramSet {
public:
	CCoalesceHistogramSet( ) : m_iMembers(0) { }

	double PoissonTest( size_t iMember, const CCoalesceHistogramSet& HistSet ) const {
		double	dAveOne, dVarOne, dAveTwo, dVarTwo;

		if( !GetEdges( ) || ( GetEdges( ) != HistSet.GetEdges( ) ) )
			return 1;
		GetAveVar( iMember, dAveOne, dVarOne );
		HistSet.GetAveVar( iMember, dAveTwo, dVarTwo );

// Fake a Skellam with a normal; Student's T tends to be slightly too sensitive
		return CStatistics::TTestWelch( dAveOne, dVarOne, GetTotal( ), dAveTwo, dVarTwo,
			HistSet.GetTotal( ) ); }
// The real thing isn't sensitive enough, or I'm doing it wrong
//		return CStatistics::SkellamPDF( 0, dAveOne, dAveTwo ); }

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
		std::copy( vecEdges.begin( ), vecEdges.end( ), m_vecEdges.begin( ) ); }

	void Initialize( size_t iMembers, size_t iBins, tValue Step ) {
		std::vector<tValue>	vecEdges;
		size_t				i;

		vecEdges.resize( iBins );
		for( i = 0; i < vecEdges.size( ); ++i )
			vecEdges[ i ] = i * Step;

		Initialize( iMembers, vecEdges ); }

	bool Add( size_t iMember, tValue Value, tCount Count ) {
		size_t	i;

		if( m_vecEdges.empty( ) || ( iMember >= m_iMembers ) || !( i = GetBin( Value ) ) )
			return false;

		if( m_vecCounts.empty( ) )
			m_vecCounts.resize( GetOffset( m_iMembers ) + GetEdges( ) );
		m_vecCounts[ GetOffset( iMember ) + i ] += Count;
		if( m_vecTotal.empty( ) )
			m_vecTotal.resize( m_iMembers );
		m_vecTotal[ iMember ] += Count;

		return true; }

	tCount Get( size_t iMember, size_t iBin ) const {
		size_t	iOffset;

		return ( ( ( iMember < m_iMembers ) && ( iBin < GetEdges( ) ) &&
			( ( iOffset = ( GetOffset( iMember ) + iBin ) ) < m_vecCounts.size( ) ) ) ?
			( iBin ? m_vecCounts[ iOffset ] : ( GetTotal( ) - m_vecTotal[ iMember ] ) ) : 0 ); }

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
		*pRet = Get( iMember, 0 );
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
			ossm << (size_t)Get( iMember, (size_t)0 );
			for( i = 1; i < GetEdges( ); ++i )
				ossm << '\t' << (size_t)Get( iMember, i ); }

		return ossm.str( ); }

	void Clear( ) {

		fill( m_vecCounts.begin( ), m_vecCounts.end( ), 0 );
		fill( m_vecTotal.begin( ), m_vecTotal.end( ), 0 ); }

	tCount GetTotal( ) const {

		return m_Total; }

	void SetTotal( tCount Total ) {

		m_Total = Total; }

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

		dAve = (double)Ave / GetTotal( );
		dVar = ( (double)Var / GetTotal( ) ) - ( dAve * dAve ); }

	size_t				m_iMembers;
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

// One score per motif atom
class CCoalesceGeneScores : public CCoalesceSequencer<std::vector<float> > {
protected:
	typedef std::vector<float>				TVecD;
	typedef std::map<uint32_t, uint32_t>	TMapII;

public:
	void Save( std::ofstream& ofsm ) const {
		size_t		iType, iSubsequence;
		uint32_t	iSize;

		iSize = (uint32_t)GetTypes( );
		ofsm.write( (const char*)&iSize, sizeof(iSize) );
		for( iType = 0; iType < GetTypes( ); ++iType ) {
			iSize = (uint32_t)GetType( iType ).length( );
			ofsm.write( (const char*)&iSize, sizeof(iSize) );
			ofsm.write( (const char*)GetType( iType ).c_str( ), iSize );

			for( iSubsequence = 0; iSubsequence < GetSubsequences( iType ); ++iSubsequence ) {
				const TVecD&					vecdMotifs	= Get( iType, (ESubsequence)iSubsequence );
				const std::vector<uint32_t>&	veciMotifs	= m_MotifVecs.Get( iType,
																(ESubsequence)iSubsequence );

				iSize = (uint32_t)vecdMotifs.size( );
				ofsm.write( (const char*)&iSize, sizeof(iSize) );
				if( iSize ) {
					ofsm.write( (const char*)&vecdMotifs.front( ), (std::streamsize)( vecdMotifs.size( ) *
						sizeof(vecdMotifs.front( )) ) );
					ofsm.write( (const char*)&veciMotifs.front( ), (std::streamsize)( veciMotifs.size( ) *
						sizeof(veciMotifs.front( )) ) ); } } } }

	void Open( std::ifstream& ifsm ) {
		uint32_t		iSize, iTypes, iType, iMotifs, iMotif;
		size_t			iSubsequence;
		vector<char>	veccType;

		ifsm.read( (char*)&iTypes, sizeof(iTypes) );
		for( iType = 0; iType < iTypes; ++iType ) {
			ifsm.read( (char*)&iSize, sizeof(iSize) );
			if( veccType.size( ) <= iSize )
				veccType.resize( iSize + 1 );
			ifsm.read( (char*)&veccType.front( ), iSize );
			veccType[ iSize ] = 0;
			AddType( &veccType.front( ) );
			m_MotifMaps.AddType( &veccType.front( ) );
			m_MotifVecs.AddType( &veccType.front( ) );

			for( iSubsequence = 0; iSubsequence < GetSubsequences( iType ); ++iSubsequence ) {
				TVecD&					vecdMotifs	= Get( iType, (ESubsequence)iSubsequence );
				std::vector<uint32_t>&	veciMotifs	= m_MotifVecs.Get( iType, (ESubsequence)iSubsequence );
				TMapII&					mapiiMotifs	= m_MotifMaps.Get( iType, (ESubsequence)iSubsequence );

				ifsm.read( (char*)&iMotifs, sizeof(iMotifs) );
				if( iMotifs ) {
					vecdMotifs.resize( iMotifs );
					ifsm.read( (char*)&vecdMotifs.front( ), (std::streamsize)( vecdMotifs.size( ) *
						sizeof(vecdMotifs.front( )) ) );
					veciMotifs.resize( iMotifs );
					ifsm.read( (char*)&veciMotifs.front( ), (std::streamsize)( veciMotifs.size( ) *
						sizeof(veciMotifs.front( )) ) );
					for( iMotif = 0; iMotif < veciMotifs.size( ); ++iMotif )
						mapiiMotifs[ veciMotifs[ iMotif ] ] = iMotif; } } } }

	bool Add( CCoalesceMotifLibrary&, const SFASTASequence&, std::vector<std::vector<unsigned short> >&,
		std::vector<size_t>& );

	float GetGlobal( size_t iType, size_t iSubsequence, uint32_t iMotif ) const {
		const TMapII&			mapiiMotifs	= m_MotifMaps.Get( iType, (ESubsequence)iSubsequence );
		TMapII::const_iterator	iterMotif;

		return ( ( ( iterMotif = mapiiMotifs.find( iMotif ) ) == mapiiMotifs.end( ) ) ? 0 :
			GetLocal( iType, iSubsequence, iterMotif->second ) ); }

	float GetLocal( size_t iType, size_t iSubsequence, uint32_t iMotif ) const {

		return CCoalesceSequencer<TVecD>::Get( iType, (ESubsequence)iSubsequence )[ iMotif ]; }

	uint32_t GetMotifs( size_t iType, size_t iSubsequence ) const {

		return (uint32_t)m_MotifVecs.Get( iType, (ESubsequence)iSubsequence ).size( ); }

	uint32_t GetMotif( size_t iType, size_t iSubsequence, uint32_t iMotif ) const {

		return (uint32_t)m_MotifVecs.Get( iType, (ESubsequence)iSubsequence )[ iMotif ]; }

protected:
	bool Add( CCoalesceMotifLibrary&, const std::string&, size_t, bool,
		vector<vector<unsigned short> >&, vector<size_t>& );

	uint32_t AddMotif( size_t iType, ESubsequence eSubsequence, uint32_t iMotif ) {
		TMapII::const_iterator	iterMotif;
		size_t					i;
		uint32_t				iRet;

		for( i = m_MotifMaps.GetTypes( ); m_MotifMaps.GetTypes( ) <= iType; ++i )
			m_MotifMaps.AddType( GetType( i ) );
		for( i = m_MotifVecs.GetTypes( ); m_MotifVecs.GetTypes( ) <= iType; ++i )
			m_MotifVecs.AddType( GetType( i ) );
		{
			TMapII&	mapiiMotifs	= m_MotifMaps.Get( iType, eSubsequence );

			if( ( iterMotif = mapiiMotifs.find( iMotif ) ) != mapiiMotifs.end( ) )
				return iterMotif->second;

			m_MotifVecs.Get( iType, eSubsequence ).push_back( iMotif );
			iRet = (uint32_t)mapiiMotifs.size( );
			mapiiMotifs[ iMotif ] = iRet;
			CCoalesceSequencer<TVecD>::Get( iType, eSubsequence ).push_back( 0 );
			return iRet;
		} }

	static void Add( ESubsequence eSubsequence, uint32_t iMotif, uint32_t iMotifs,
		vector<vector<unsigned short> >& vecvecsCounts ) {
		std::vector<unsigned short>&	vecsCountsTotal	= vecvecsCounts[ ESubsequenceTotal ];
		std::vector<unsigned short>&	vecsCounts		= vecvecsCounts[ eSubsequence ];

		vecsCountsTotal.resize( iMotifs );
		vecsCountsTotal[ iMotif ]++;
		vecsCounts.resize( iMotifs );
		vecsCounts[ iMotif ]++; }

	CCoalesceSequencer<TMapII>					m_MotifMaps;
	CCoalesceSequencer<std::vector<uint32_t> >	m_MotifVecs;
};

// One histogram per motif atom
class CCoalesceGroupHistograms : public CCoalesceSequencer<CCoalesceHistogramSet<> > {
public:
	CCoalesceGroupHistograms( uint32_t iMotifs, size_t iBins, float dStep ) : m_iMotifs(iMotifs),
		m_iBins(iBins), m_dStep(dStep) { }

	bool Add( const CCoalesceGeneScores&, bool );
	void Save( std::ostream&, const CCoalesceMotifLibrary* ) const;

	uint32_t GetMotifs( ) const {

		return m_iMotifs; }

	void SetTotal( const std::vector<CCoalesceGeneScores>& vecGeneScores, const std::set<size_t>& setiGenes ) {
		std::set<size_t>::const_iterator	iterGene;
		size_t								iTypeUs, iTypeThem, iSubsequence;
		unsigned short						sTotal;

		for( iTypeUs = 0; iTypeUs < GetTypes( ); ++iTypeUs )
			for( iSubsequence = ESubsequenceBegin; iSubsequence < GetSubsequences( iTypeUs ); ++iSubsequence ) {
				for( sTotal = 0,iterGene = setiGenes.begin( ); iterGene != setiGenes.end( ); ++iterGene )
					if( ( ( iTypeThem = vecGeneScores[ *iterGene ].GetType( GetType( iTypeUs ) ) ) != -1 ) &&
						!vecGeneScores[ *iterGene ].Get( iTypeThem, (ESubsequence)iSubsequence ).empty( ) )
						sTotal++;
				Get( iTypeUs, (ESubsequence)iSubsequence ).SetTotal( sTotal ); } }

protected:
	uint32_t	m_iMotifs;
	size_t		m_iBins;
	float		m_dStep;
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
	bool CalculateProbabilityExpression( size_t, const CPCL&, const CCoalesceCluster&, bool, long double&,
		long double& ) const;
	bool CalculateProbabilityMotifs( const CCoalesceGeneScores&, const CCoalesceGroupHistograms&,
		const CCoalesceGroupHistograms&, bool, long double&, long double& ) const;
	bool SaveCopy( const CPCL&, size_t, CPCL&, size_t, bool ) const;

	bool IsGene( size_t iGene ) const {

		return ( m_setiGenes.find( iGene ) != m_setiGenes.end( ) ); }

	bool IsCondition( size_t iCondition ) const {

		return ( m_setiConditions.find( iCondition ) != m_setiConditions.end( ) ); }

	size_t GetHash( ) const {

		return ( GetHash( m_setiConditions ) ^ GetHash( m_setiGenes ) ^ GetHash( m_setsMotifs ) ); }

	std::set<size_t>			m_setiConditions;
	std::set<size_t>			m_setiGenes;
	std::set<SMotifMatch>		m_setsMotifs;
	std::vector<size_t>			m_veciPrevConditions;
	std::vector<size_t>			m_veciPrevGenes;
	std::vector<SMotifMatch>	m_vecsPrevMotifs;
	std::vector<size_t>			m_veciCounts;
	std::vector<float>			m_vecdCentroid;
	std::vector<float>			m_vecdStdevs;
	std::set<size_t>			m_setiHistory;
};

class CCoalesceImpl {
protected:
	CCoalesceImpl( ) : m_iK(7), m_dPValueCorrelation(0.05f), m_iBins(12), m_dPValueCondition(0.05f),
		m_dProbabilityGene(0.95f), m_dPValueMotif(0.05f), m_pMotifs(NULL), m_fMotifs(false),
		m_iBasesPerMatch(1000) { }
	virtual ~CCoalesceImpl( );

	void Clear( );
	size_t GetMotifCount( ) const;
	void Save( const std::vector<CCoalesceGeneScores>& ) const;
	void Open( std::vector<CCoalesceGeneScores>& ) const;

	float					m_dProbabilityGene;
	float					m_dPValueCondition;
	float					m_dPValueCorrelation;
	float					m_dPValueMotif;
	size_t					m_iBins;
	size_t					m_iK;
	std::string				m_strDirectoryIntermediate;
	CCoalesceMotifLibrary*	m_pMotifs;
	bool					m_fMotifs;
	size_t					m_iBasesPerMatch;
	std::string				m_strSequenceCache;
};

}

#endif // COALESCEI_H

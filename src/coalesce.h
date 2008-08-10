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
	static size_t KMer2ID( const std::string& strKMer ) {
		size_t		i, iRet;
		const char*	pc;

		for( iRet = i = 0; i < strKMer.size( ); ++i ) {
			if( !( pc = strchr( c_acBases, strKMer[ i ] ) ) )
				return -1;
			iRet = ( iRet << c_iShift ) | ( pc - c_acBases ); }

		return iRet; }

	static std::string ID2KMer( size_t iID, size_t iK ) {
		std::string	strRet;
		size_t		i, iMask;

		iMask = ( 1 << c_iShift ) - 1;
		strRet.resize( iK );
		for( i = 0; i < iK; ++i ) {
			strRet[ iK - i - 1 ] = c_acBases[ iID & iMask ];
			iID >>= c_iShift; }

		return strRet; }

	static size_t CountKMers( size_t iK ) {

		return ( 1 << ( 2 * iK ) ); }

	static bool IsIgnorableKMer( const std::string& strKMer ) {

		return ( strKMer.find( 'N' ) != std::string::npos ); }

	CCoalesceMotifLibrary( size_t iK ) : CCoalesceMotifLibraryImpl( iK ) { }

	std::string GetMotif( size_t iMotif ) const {

		return ID2KMer( iMotif, GetK( ) ); }

	size_t GetID( const std::string& strKMer ) const {

		return KMer2ID( strKMer ); }

	size_t GetMotifs( ) const {

		return CountKMers( GetK( ) ); }

	size_t GetK( ) const {

		return m_iK; }
};

class CCoalesceCluster : public CCoalesceClusterImpl {
public:
	bool Initialize( const CPCL& PCL, CCoalesceCluster& Pot, float dPValue );
	void Subtract( CPCL& PCL ) const;
	bool SelectConditions( const CPCL& PCL, const CCoalesceCluster& Pot, float dPValue );
	bool SelectMotifs( const CCoalesceGroupHistograms& HistsCluster, const CCoalesceGroupHistograms& HistsPot,
		float dPValue, const CCoalesceMotifLibrary* pMotifs = NULL );
	bool SelectGenes( const CPCL& PCL, const std::vector<CCoalesceGeneScores>& vecGeneScores,
		const CCoalesceGroupHistograms& HistsCluster, const CCoalesceGroupHistograms& HistsPot,
		CCoalesceCluster& Pot, float dPValue, const CCoalesceMotifLibrary* pMotifs = NULL );
	void CalculateHistograms( const std::vector<CCoalesceGeneScores>& vecGeneScores,
		CCoalesceGroupHistograms& Histograms ) const;
	bool Save( const std::string& strDirectory, size_t iID, const CPCL& PCL,
		const CCoalesceMotifLibrary* pMotifs = NULL ) const;
	void Save( std::ostream&, size_t iID, const CPCL& PCL, const CCoalesceMotifLibrary* pMotifs = NULL ) const;

	bool IsConverged( ) {

// BUGBUG: this is redundant
		return ( ( m_setiHistory.find( GetHash( ) ) != m_setiHistory.end( ) ) ||
			CCoalesceClusterImpl::IsConverged( m_setiConditions, m_veciPrevConditions ) &&
			CCoalesceClusterImpl::IsConverged( m_setiGenes, m_veciPrevGenes ) &&
			CCoalesceClusterImpl::IsConverged( m_setsMotifs, m_vecsPrevMotifs ) ); }

	bool IsEmpty( ) const {

		return ( m_setiGenes.empty( ) || m_setiConditions.empty( ) ); }

	void Add( size_t iGene ) {

		m_setiGenes.insert( iGene ); }

	void Snapshot( ) {

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
};

}

#endif // COALESCE_H

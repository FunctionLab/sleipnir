#ifndef COALESCEI_H
#define COALESCEI_H

#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace Sleipnir {

class CCoalesceCluster;
class CFASTA;
class CPCL;
struct SFASTASequence;
struct SHistogram;

struct SHistogramImpl {
protected:
	size_t						m_iTotal;
	std::vector<unsigned short>	m_vecsBins;
};

class CCoalesceHistograms {
public:
	CCoalesceHistograms( size_t iMotifs, size_t iBins ) : m_iMotifs(iMotifs), m_iBins(iBins) { }

	void Add( const CCoalesceHistograms& );
	bool Add( const SFASTASequence&, size_t );

	bool IsEmpty( ) const {

		return m_mapstriTypes.empty( ); }

	size_t GetTypes( ) const {

		return m_vecvecvecsHistograms.size( ); }

	size_t GetSubsequences( size_t iType ) const {

		return m_vecvecvecsHistograms[ iType ].size( ); }

	size_t GetMotifs( ) const {

		return ( ( m_vecvecvecsHistograms.size( ) && m_vecvecvecsHistograms[ 0 ].size( ) ) ?
			m_vecvecvecsHistograms[ 0 ][ 0 ].size( ) : 0 ); }

	const SHistogram& Get( size_t iType, size_t iSubsequence, size_t iMotif ) const {

		return m_vecvecvecsHistograms[ iType ][ iSubsequence ][ iMotif ]; }

	size_t GetBins( ) const {

		return m_iBins; }

	size_t GetType( const std::string& strType ) const {
		TMapStrI::const_iterator	iterType;

		return ( ( ( iterType = m_mapstriTypes.find( strType ) ) == m_mapstriTypes.end( ) ) ? -1 :
			iterType->second ); }

	const std::string& GetType( size_t iType ) const {

		return m_vecstrTypes[ iType ]; }

protected:
	typedef std::map<std::string, size_t>	TMapStrI;

	enum ESubsequence {
		ESubsequenceIntrons	= 0,
		ESubsequenceExons	= ESubsequenceIntrons + 1,
		ESubsequenceTotal	= ESubsequenceExons + 1
	};

	bool Add( const std::string&, size_t, size_t, bool );
	size_t AddType( const std::string& );

	size_t												m_iMotifs;
	size_t												m_iBins;
	TMapStrI											m_mapstriTypes;
	std::vector<std::string>							m_vecstrTypes;
// Type by subsequence by motif
	std::vector<std::vector<std::vector<SHistogram> > >	m_vecvecvecsHistograms;
};

class CCoalesceClusterImpl {
protected:
	static bool IsConverged( std::set<size_t>& setiNew, std::vector<size_t>& veciOld ) {
		size_t				i;
		std::vector<size_t>	veciNew;

		if( setiNew.size( ) != veciOld.size( ) )
			return false;
		Snapshot( setiNew, veciNew );
		for( i = 0; i < veciNew.size( ); ++i )
			if( veciNew[ i ] != veciOld[ i ] )
				return false;

		return true; }

	static void Snapshot( std::set<size_t>& setiNew, std::vector<size_t>& veciOld ) {

		veciOld.resize( setiNew.size( ) );
		std::copy( setiNew.begin( ), setiNew.end( ), veciOld.begin( ) );
		std::sort( veciOld.begin( ), veciOld.end( ) ); }

	void Add( size_t, CCoalesceCluster& );
	bool AddCorrelatedGenes( const CPCL&, CCoalesceCluster&, float );
	bool AddSeedPair( const CPCL&, CCoalesceCluster&, float );
	void CalculateCentroid( const CPCL& );
	bool IsSignificant( size_t, const CCoalesceHistograms&, const CCoalesceHistograms&, float ) const;
	bool IsSignificant( size_t, const CPCL&, const CCoalesceHistograms&, const CCoalesceHistograms&,
		const CCoalesceHistograms&, const CCoalesceCluster&, float ) const;
	float CalculateProbabilityExpression( size_t, const CPCL&, const CCoalesceCluster&, bool ) const;
	float CalculateProbabilityMotifs( const CCoalesceHistograms&, const CCoalesceHistograms&,
		const CCoalesceHistograms&, bool ) const;
	bool SaveCopy( const CPCL&, size_t, CPCL&, size_t, bool ) const;

	bool IsGene( size_t iGene ) const {

		return ( m_setiGenes.find( iGene ) != m_setiGenes.end( ) ); }

	bool IsCondition( size_t iCondition ) const {

		return ( m_setiConditions.find( iCondition ) != m_setiConditions.end( ) ); }

	std::set<size_t>	m_setiConditions;
	std::set<size_t>	m_setiGenes;
	std::set<size_t>	m_setiMotifs;
	std::vector<size_t>	m_veciPrevConditions;
	std::vector<size_t>	m_veciPrevGenes;
	std::vector<size_t>	m_veciPrevMotifs;
	std::vector<size_t>	m_veciCounts;
	std::vector<float>	m_vecdCentroid;
	std::vector<float>	m_vecdStdevs;
};

class CCoalesceImpl {
protected:
	CCoalesceImpl( ) : m_iK(7), m_dPValueCorrelation(0.05f), m_iBins(12), m_dPValueCondition(0.05f),
		m_dProbabilityGene(0.95f), m_dPValueMotif(0.05f), m_fOutputIntermediate(false) { }

	float	m_dProbabilityGene;
	float	m_dPValueCondition;
	float	m_dPValueCorrelation;
	float	m_dPValueMotif;
	size_t	m_iBins;
	size_t	m_iK;
	bool	m_fOutputIntermediate;
};

}

#endif // COALESCEI_H

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
#include "stdafx.h"
#include "dat.h"
#include "genome.h"
#include "statistics.h"
#include "typesi.h"
#include "annotation.h"
#include "color.h"

namespace Sleipnir {

const char		CDatImpl::c_acComment[]		= "#";
const CColor&	CDatImpl::c_ColorMax		= CColor::c_Red;
const CColor&	CDatImpl::c_ColorMin		= CColor::c_Green;
const CColor&	CDatImpl::c_ColorMid		= CColor::c_Black;

static const struct {
	const char*			m_szExtension;
	const CDat::EFormat	m_eFormat;
} c_asFormats[]	= {
	{"dat",	CDat::EFormatText},
	{"das",	CDat::EFormatSparse},
	{"pcl",	CDat::EFormatPCL},
	{NULL,	CDat::EFormatBinary}
};

size_t CDatImpl::MapGene( TMapStrI& mapGenes, TVecStr& vecGenes, const std::string& strToken ) {
	TMapStrI::iterator	iterGenes;
	size_t				iRet;

	if( ( iterGenes = mapGenes.find( strToken ) ) == mapGenes.end( ) ) {
		iRet = mapGenes.size( );
		mapGenes[ strToken ] = iRet;
		vecGenes.push_back( strToken ); }
	else
		iRet = iterGenes->second;

	return iRet; }

void CDatImpl::ResizeNaN( TAF& vecf, size_t iLast ) {
	size_t	i;

	if( ( i = vecf.size( ) ) > iLast )
		return;

	vecf.resize( iLast + 1 );
	for( ; i < vecf.size( ); ++i )
		vecf[ i ] = CMeta::GetNaN( ); }

void CDatImpl::DabGene( std::istream& istm, char* acBuffer ) {
	size_t	i;

	i = 0;
	do
		istm.seekg( 1, ios_base::cur );
	while( acBuffer[ i++ ] = istm.get( ) ); }

CDatImpl::~CDatImpl( ) {

	Reset( ); }

void CDatImpl::Reset( ) {

	m_Data.Reset( );
	m_vecstrGenes.clear( );

	if( m_pPCL && m_fPCLMemory )
		delete m_pPCL;
	m_pPCL = NULL;
	if( m_pMeasure && m_fMeasureMemory )
		delete m_pMeasure;
	m_pMeasure = NULL;

	CMeta::Unmap( m_abData, m_hndlData, m_iData );
	m_abData = NULL;
	m_hndlData = 0;
	m_iData = 0;
	if( m_aadData )
		delete[] m_aadData;
	m_aadData = NULL; }

void CDatImpl::SlimCache( const CSlim& Slim, vector<vector<size_t> >& vecveciGenes ) const {
	size_t	iS, iG;

	vecveciGenes.resize( Slim.GetSlims( ) );
	for( iS = 0; iS < vecveciGenes.size( ); ++iS ) {
		vecveciGenes[ iS ].resize( Slim.GetGenes( iS ) );
		for( iG = 0; iG < vecveciGenes[ iS ].size( ); ++iG )
			vecveciGenes[ iS ][ iG ] =  GetGene( Slim.GetGene( iS, iG ).GetName( ) ); } }

/*!
 * \brief
 * Construct a CDat from the given ontology slim.
 * 
 * \param Slim
 * Set of ontology terms from which to generate a CDat.
 * 
 * \returns
 * True if CDat was generated successfully.
 * 
 * Generates a CDat over all genes within the given slim.  Gene pairs coannotated to at least one slim term
 * are given a value of 1; all other gene pairs are given a value of 0.  This is useful for rapidly generating
 * gold standards from functional catalogs.
 * 
 * \see
 * IOntology
 */
bool CDat::Open( const CSlim& Slim ) {
	vector<string>			vecstrGenes;
	size_t					iS1, iS2, iG1, iG2, iGene1, iGene2;
	vector<vector<size_t> >	vecveciGenes;

	Reset( );
	Slim.GetGeneNames( vecstrGenes );
	if( !Open( vecstrGenes ) )
		return false;

	SlimCache( Slim, vecveciGenes );
	for( iS1 = 0; iS1 < Slim.GetSlims( ); ++iS1 ) {
		g_CatSleipnir.info( "CDat::Open( ) processing slim: %s",
			Slim.GetSlim( iS1 ).c_str( ) );
		for( iG1 = 0; iG1 < Slim.GetGenes( iS1 ); ++iG1 ) {
			iGene1 = vecveciGenes[ iS1 ][ iG1 ];
			for( iG2 = ( iG1 + 1 ); iG2 < Slim.GetGenes( iS1 ); ++iG2 )
				Set( iGene1, vecveciGenes[ iS1 ][ iG2 ], 1 );
			for( iS2 = ( iS1 + 1 ); iS2 < Slim.GetSlims( ); ++iS2 )
				for( iG2 = 0; iG2 < Slim.GetGenes( iS2 ); ++iG2 ) {
					iGene2 = vecveciGenes[ iS2 ][ iG2 ];
					if( CMeta::IsNaN( Get( iGene1, iGene2 ) ) )
						Set( iGene1, iGene2, 0 ); } } }

	return true; }

/*!
 * \brief
 * Construct a CDat from the given ontology slims.
 * 
 * \param SlimPositives
 * Set of ontology terms from which to generate related gene pairs.
 * 
 * \param SlimNonnegatives
 * Set of ontology terms from which to generate agnostic gene pairs (neither related nor unrelated).
 * 
 * \returns
 * True if CDat was generated successfully.
 * 
 * Generates a CDat over all genes within the given slims.  Gene pairs coannotated to at least one positive
 * slim term are given a value of 1, other gene pairs coannotated to at least one nonnegative slim term are
 * given no value (NaN), and all other gene pairs are given a value of 0.  This is useful for rapidly
 * generating gold standards from functional catalogs.
 * 
 * \see
 * IOntology
 */
bool CDat::Open( const CSlim& SlimPositives, const CSlim& SlimNonnegatives ) {
	set<string>				setstrGenes;
	vector<string>			vecstrGenes;
	size_t					iS1, iS2, iG1, iG2, iGene1, iGene2;
	vector<vector<size_t> >	vecveciGenes;

	Reset( );
	SlimPositives.GetGeneNames( vecstrGenes );
	setstrGenes.insert( vecstrGenes.begin( ), vecstrGenes.end( ) );
	vecstrGenes.clear( );
	SlimNonnegatives.GetGeneNames( vecstrGenes );
	setstrGenes.insert( vecstrGenes.begin( ), vecstrGenes.end( ) );
	vecstrGenes.clear( );
	vecstrGenes.resize( setstrGenes.size( ) );
	copy( setstrGenes.begin( ), setstrGenes.end( ), vecstrGenes.begin( ) );
	if( !Open( vecstrGenes ) )
		return false;

	SlimCache( SlimPositives, vecveciGenes );
	for( iS1 = 0; iS1 < SlimPositives.GetSlims( ); ++iS1 ) {
		g_CatSleipnir.info( "CDat::Open( ) processing slim: %s",
			SlimPositives.GetSlim( iS1 ).c_str( ) );
		for( iG1 = 0; iG1 < SlimPositives.GetGenes( iS1 ); ++iG1 ) {
			iGene1 = vecveciGenes[ iS1 ][ iG1 ];
			for( iG2 = ( iG1 + 1 ); iG2 < SlimPositives.GetGenes( iS1 ); ++iG2 )
				Set( iGene1, vecveciGenes[ iS1 ][ iG2 ], 1 ); } }
	vecveciGenes.clear( );
	SlimCache( SlimNonnegatives, vecveciGenes );
	for( iS1 = 0; iS1 < SlimNonnegatives.GetSlims( ); ++iS1 ) {
		g_CatSleipnir.info( "CDat::Open( ) processing slim: %s",
			SlimNonnegatives.GetSlim( iS1 ).c_str( ) );
		for( iG1 = 0; iG1 < SlimNonnegatives.GetGenes( iS1 ); ++iG1 ) {
			iGene1 = vecveciGenes[ iS1 ][ iG1 ];
			for( iS2 = ( iS1 + 1 ); iS2 < SlimNonnegatives.GetSlims( ); ++iS2 )
				for( iG2 = 0; iG2 < SlimNonnegatives.GetGenes( iS2 ); ++iG2 ) {
					iGene2 = vecveciGenes[ iS2 ][ iG2 ];
					if( CMeta::IsNaN( Get( iGene1, iGene2 ) ) )
						Set( iGene1, iGene2, 0 ); } } }

	return true; }

/*!
 * \brief
 * Construct a CDat from the given known gene relationships and gene sets.
 * 
 * \param DatKnown
 * Known pairwise scores, either positive or negative as indicated.
 * 
 * \param vecpOther
 * Gene sets, either positive or nonnegative as indicated (possibly empty).
 * 
 * \param Genome
 * Genome containing all genes of interest.
 * 
 * \param fKnownNegatives
 * If true, DatKnown contains known negative gene pairs (0 scores); if false, it contains known related
 * gene pairs (1 scores).  In the former case, positives are generated from pairs coannotated to the
 * given gene sets; in the latter, negatives are generated from pairs not coannotated to the given gene
 * sets.
 * 
 * \returns
 * True if CDat was generated successfully.
 * 
 * Constructs a CDat by copying either the positive (1) or negative (0) values from DatKnown and calculating
 * negatives or positives from vecpOther.  In either case, values are calculated for all genes in the given
 * genome.
 * - If fKnownNegatives is true, negative (0) gene pairs are first copied from DatKnown.  Then any other
 * gene pair coannotated to at least one of the given gene sets is given a positive (1) score.  Remaining
 * gene pairs are given no value (NaN).
 * - If fKnownNegatives is false, positive (1) gene pairs are first copied from DatKnown.  Then any other
 * gene pair coannotated to at least one of the given gene sets is given a missing (NaN) value.  Remaining
 * gene pairs are given a negative (0) score.
 */
bool CDat::Open( const CDat& DatKnown, const vector<CGenes*>& vecpOther, const CGenome& Genome,
	bool fKnownNegatives ) {
	size_t			i, j, iOne, iTwo;
	vector<size_t>	veciGenes;
	float			d;

	Reset( );
	Open( Genome.GetGeneNames( ) );
	for( i = 0; i < vecpOther.size( ); ++i )
		OpenHelper( vecpOther[ i ], 1 );
	if( !fKnownNegatives )
		for( i = 0; i < GetGenes( ); ++i )
			for( j = ( i + 1 ); j < GetGenes( ); ++j )
				Set( i, j, CMeta::IsNaN( Get( i, j ) ) ? 0 : CMeta::GetNaN( ) );

	veciGenes.resize( DatKnown.GetGenes( ) );
	for( i = 0; i < veciGenes.size( ); ++i )
		veciGenes[ i ] = GetGene( DatKnown.GetGene( i ) );
	for( i = 0; i < DatKnown.GetGenes( ); ++i ) {
		iOne = veciGenes[ i ];
		for( j = ( i + 1 ); j < DatKnown.GetGenes( ); ++j ) {
			d = DatKnown.Get( i, j );
			if( CMeta::IsNaN( Get( iOne, iTwo = veciGenes[ j ] ) ) && ( fKnownNegatives == !d ) )
				Set( iOne, iTwo, d ); } }

	return true; }

/*!
 * \brief
 * Construct a copy of the given CDat.
 * 
 * \param Dat
 * Data to be copied.
 * 
 * \returns
 * True if the copy was successful.
 */
bool CDat::Open( const CDat& Dat ) {
	size_t	i;

	if( !Open( Dat.GetGeneNames( ) ) )
		return false;

	for( i = 0; i < GetGenes( ); ++i )
		memcpy( Get( i ), Dat.Get( i ), ( GetGenes( ) - i - 1 ) * sizeof(*Get( i )) );

	return true; }

/*!
 * \brief
 * Construct a CDat from the given gene sets.
 * 
 * \param vecpPositives
 * Set of gene sets from which to generate related gene pairs.
 * 
 * \param vecpNonnegatives
 * Set of ontology terms from which to generate agnostic gene pairs (neither related nor unrelated),
 * possibly empty.
 * 
 * \param dPValue
 * Hypergeometric p-value of overlap below which nonnegative gene set pairs are considered agnostic rather
 * than negative.
 * 
 * \param Genome
 * Genome containing all genes of interest.
 * 
 * \returns
 * True if CDat was generated successfully.
 * 
 * Generates a CDat over all genes in the given genome in which:
 * - Any gene pair coannotated to at least one positive gene set is assigned value 1.
 * - Any other gene pair coannotated to at least one nonnegative gene set is assigned no value (NaN).
 * - Any other gene pair annotated to nonnegative gene sets with hypergeometric p-value of overlap less than
 * the given cutoff is assigned no value (NaN).
 * - Any other gene pair is assigned value 0.
 * This is useful for rapidly generating gold standards from gene sets of interest: any gene pair known to
 * participate in some function of interest will be marked as related, any gene pair known to participate
 * in two unrelated functions will be marked as unrelated, and other gene pairs will be left unmarked.
 */
bool CDat::Open( const vector<CGenes*>& vecpPositives, const vector<CGenes*>& vecpNonnegatives,
	float dPValue, const CGenome& Genome ) {
	size_t			i, j, k, iOne, iTwo, iOverlap;
	float			d;
	const CGenes*	pBig;
	const CGenes*	pSmall;

	Reset( );
	Open( Genome.GetGeneNames( ) );
	for( i = 0; i < vecpPositives.size( ); ++i )
		OpenHelper( vecpPositives[ i ], 1 );
	if( vecpNonnegatives.size( ) ) {
		for( i = 0; i < vecpNonnegatives.size( ); ++i )
			OpenHelper( vecpNonnegatives[ i ], 0 );
		if( dPValue > 0 )
			for( i = 0; i < vecpPositives.size( ); ++i ) {
				iOne = vecpPositives[ i ]->GetGenes( );
				for( j = ( i + 1 ); j < vecpPositives.size( ); ++j ) {
					iTwo = vecpPositives[ j ]->GetGenes( );
					if( iOne < iTwo ) {
						pSmall = vecpPositives[ i ];
						pBig = vecpPositives[ j ]; }
					else {
						pSmall = vecpPositives[ j ];
						pBig = vecpPositives[ i ]; }
					for( iOverlap = k = 0; k < pSmall->GetGenes( ); ++k )
						if( pBig->IsGene( pSmall->GetGene( k ).GetName( ) ) )
							iOverlap++;
					if( CStatistics::HypergeometricCDF( iOverlap, iOne, iTwo, GetGenes( ) ) < dPValue )
						OpenHelper( pBig, pSmall, 0 ); } }
		for( i = 0; i < GetGenes( ); ++i )
			for( j = ( i + 1 ); j < GetGenes( ); ++j )
				if( CMeta::IsNaN( d = Get( i, j ) ) )
					Set( i, j, 0 );
				else if( !d )
					Set( i, j, CMeta::GetNaN( ) ); }

	return true; }

void CDatImpl::OpenHelper( const CGenes* pGenes, float dValue ) {
	vector<size_t>	veciGenes;
	size_t			i, j, iOne, iTwo;

	veciGenes.resize( pGenes->GetGenes( ) );
	for( i = 0; i < veciGenes.size( ); ++i )
		veciGenes[ i ] = GetGene( pGenes->GetGene( i ).GetName( ) );
	for( i = 0; i < veciGenes.size( ); ++i ) {
		iOne = veciGenes[ i ];
		for( j = ( i + 1 ); j < veciGenes.size( ); ++j ) {
			iTwo = veciGenes[ j ];
			if( CMeta::IsNaN( Get( iOne, iTwo ) ) )
				Set( iOne, iTwo, dValue ); } } }

void CDatImpl::OpenHelper( const CGenes* pOne, const CGenes* pTwo, float dValue ) {
	vector<size_t>	veciOne, veciTwo;
	size_t			i, j, iOne, iTwo;

	veciOne.resize( pOne->GetGenes( ) );
	for( i = 0; i < veciOne.size( ); ++i )
		veciOne[ i ] = GetGene( pOne->GetGene( i ).GetName( ) );
	veciTwo.resize( pTwo->GetGenes( ) );
	for( i = 0; i < veciTwo.size( ); ++i )
		veciTwo[ i ] = GetGene( pTwo->GetGene( i ).GetName( ) );
	for( i = 0; i < veciOne.size( ); ++i ) {
		iOne = veciOne[ i ];
		for( j = 0; j < veciTwo.size( ); ++j ) {
			iTwo = veciTwo[ j ];
			if( CMeta::IsNaN( Get( iOne, iTwo ) ) )
				Set( iOne, iTwo, dValue ); } } }

/*!
 * \brief
 * Open a CDat stored in the given stream with the given format, processing missing values and duplicates as
 * specified.
 * 
 * \param istm
 * Stream from which CDat is loaded.
 * 
 * \param eFormat
 * Format in which the stream should be parsed.
 * 
 * \param dDefault
 * Default value inserted for missing pairs (DAT format only).
 * 
 * \param fDuplicates
 * If true, allow duplicates (DAT format only), ignoring all but the last value for each gene pair.
 * 
 * \param iSkip
 * If the given stream contains a PCL, the number of columns to skip between the ID and experiments.
 * 
 * \param fZScore
 * If true and the given stream contains a PCL, z-score similarity measures after pairwise calculation.
 * 
 * \returns
 * True if CDat was successfully opened.
 * 
 * Opens a CDat from the given stream with the given format.
 * 
 * \remarks
 * dDefault and fDuplicates are ignored for non-DAT formats; iSkip and fZScore are ignored for non-PCL formats.
 * Specifying the format incorrectly will generally cause Bad Things.
 * 
 * \see
 * Save | CPCL
 */
bool CDat::Open( std::istream& istm, EFormat eFormat, float dDefault, bool fDuplicates, size_t iSkip,
	bool fZScore ) {

	switch( eFormat ) {
		case EFormatText:
			return OpenText( istm, dDefault, fDuplicates );

		case EFormatPCL:
			return OpenPCL( istm, iSkip, fZScore );

		case EFormatSparse:
			return OpenSparse( istm ); }

	return OpenBinary( istm ); }

bool CDatImpl::OpenPCL( std::istream& istm, size_t iSkip, bool fZScore ) {

	Reset( );
	m_pPCL = new CPCL( );
	if( !m_pPCL->Open( istm, iSkip ) )
		return false;

	m_pMeasure = new CMeasurePearNorm( );
	m_fMeasureMemory = true;
	if( fZScore ) {
		size_t	iN;
		double	dAve, dStd;

		AveStd( dAve, dStd, iN );
		delete m_pMeasure;
		m_pMeasure = new CMeasurePearNorm( dAve, dStd ); }

	return true; }

/*!
 * \brief
 * Construct a new CDat backed by the given PCL and similarity measure.
 * 
 * \param PCL
 * PCL from which CDat genes and pairwise values are drawn.
 * 
 * \param pMeasure
 * Similarity measure used to calculate pairwise scores between genes.
 * 
 * \param fMeasureMemory
 * If true, the CDat is responsible for freeing the given similarity measure.
 * 
 * \returns
 * True if CDat was generated successfully.
 * 
 * This opens a new CDat without precalculated similarity scores; the CDat retains a reference to the given
 * PCL and similarity measure and calculates pairwise scores as needed.  This can greatly reduce the amount
 * of memory required by the CDat, but it can increase runtime when specific pairwise values are requested
 * repeatedly.
 * 
 * \remarks
 * fMeasureMemory should be false if a static or stack-based similarity measure object is used; it can be
 * set to true for cloned/new allocated similarity measure objects that should be cleaned up with the CDat.
 * The given PCL is not copied and should thus not be destroyed before the current CDat.
 */
bool CDat::Open( const CPCL& PCL, const IMeasure* pMeasure, bool fMeasureMemory ) {

	Reset( );
	m_pPCL = (CPCL*)&PCL;
	m_fPCLMemory = false;
	m_pMeasure = pMeasure;
	m_fMeasureMemory = fMeasureMemory;

	return true; }

bool CDatImpl::OpenText( std::istream& istm, float dDefault, bool fDuplicates ) {
	const char*	pc;
	char*		pcTail;
	char*		acBuf;
	string		strToken, strCache, strValue;
	TMapStrI	mapGenes;
	size_t		iOne, iTwo, i;
	float		dScore;
	TAAF		vecvecfScores;

	Reset( );
	acBuf = new char[ c_iBufferSize ];
	while( istm.peek( ) != EOF ) {
		istm.getline( acBuf, c_iBufferSize - 1 );
		strToken = OpenToken( acBuf, &pc );
		if( !strToken.length( ) )
			break;
		if( strToken == c_acComment )
			continue;
		if( strToken != strCache ) {
			strCache = strToken;
			iOne = MapGene( mapGenes, m_vecstrGenes, strToken ); }

		strToken = OpenToken( pc, &pc );
		if( !strToken.length( ) ) {
			Reset( );
			delete[] acBuf;
			return false; }
		iTwo = MapGene( mapGenes, m_vecstrGenes, strToken );
		strValue = OpenToken( pc );
		if( !strValue.length( ) ) {
			if( CMeta::IsNaN( dScore = dDefault ) ) {
				Reset( );
				delete[] acBuf;
				return false; } }
		else if( !( dScore = (float)strtod( strValue.c_str( ), &pcTail ) ) &&
			( pcTail != ( strValue.c_str( ) + strValue.length( ) ) ) ) {
			Reset( );
			delete[] acBuf;
			return false; }

		i = ( ( iOne > iTwo ) ? iOne : iTwo );
		if( vecvecfScores.size( ) <= i )
			vecvecfScores.resize( i + 1 );
		ResizeNaN( vecvecfScores[ iOne ], i );
		ResizeNaN( vecvecfScores[ iTwo ], i );
		if( !CMeta::IsNaN( vecvecfScores[ iOne ][ iTwo ] ) && ( vecvecfScores[ iOne ][ iTwo ] != dScore ) ) {
			g_CatSleipnir.error( "CDatImpl::OpenText( ) duplicate genes %s, %s (%g:%g)",
				strCache.c_str( ), strToken.c_str( ), vecvecfScores[ iOne ][ iTwo ],
				dScore );
			if( !fDuplicates ) {
				Reset( );
				delete[] acBuf;
				return false; } }
		vecvecfScores[ iOne ][ iTwo ] = vecvecfScores[ iTwo ][ iOne ] = dScore; }
	delete[] acBuf;

	m_Data.Initialize( GetGenes( ) );
	for( iOne = 0; iOne < GetGenes( ); ++iOne )
		for( iTwo = ( iOne + 1 ); iTwo < GetGenes( ); ++iTwo )
			Set( iOne, iTwo, ( ( iTwo < vecvecfScores[ iOne ].size( ) ) ?
				vecvecfScores[ iOne ][ iTwo ] : CMeta::GetNaN( ) ) );

	return true; }

bool CDatImpl::OpenBinary( std::istream& istm ) {
	size_t	i;
	float*	adScores;

	if( !OpenGenes( istm, true, false ) )
		return false;
	m_Data.Initialize( GetGenes( ) );
	adScores = new float[ GetGenes( ) - 1 ];
	for( i = 0; ( i + 1 ) < GetGenes( ); ++i ) {
		istm.read( (char*)adScores, sizeof(*adScores) * ( GetGenes( ) - i - 1 ) );
		Set( i, adScores ); }
	delete[] adScores;

	return true; }

bool CDatImpl::OpenSparse( std::istream& istm ) {
	size_t		i;
	uint32_t	j;
	float		d;

	if( !OpenGenes( istm, true, false ) )
		return false;
	m_Data.Initialize( GetGenes( ) );
	for( i = 0; ( i + 1 ) < GetGenes( ); ++i ) {
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			Set( i, j, CMeta::GetNaN( ) );
		while( true ) {
			istm.read( (char*)&j, sizeof(j) );
			if( j == -1 )
				break;
			istm.read( (char*)&d, sizeof(d) );
			Set( i, j, d ); } }

	return true; }

/*!
 * \brief
 * Open the genes from a CDat stored in the given file, guessing the format from the file's extension.
 * 
 * \param szFile
 * Filename from which CDat genes are loaded.
 * 
 * \param iSkip
 * If the given file is a PCL, the number of columns to skip between the ID and experiments.
 * 
 * \returns
 * True if genes were successfully loaded.
 * 
 * Like Open, OpenGenes will guess the appropriate file format from the given filename's extension, but it
 * will load only a CDat's size and gene list from the file.  This is useful for rapidly obtaining gene
 * lists and counts from large CDats without incurring the memory/time penalty of loading or memory mapping
 * the whole file.
 * 
 * \remarks
 * Attempting to access a CDat's data after opening only its genes is a poor idea.  If the extension is not
 * recognized, DAB format is assumed.
 */
bool CDat::OpenGenes( const char* szFile, size_t iSkip ) {
	ifstream	ifsm;
	bool		fBinary;
	size_t		i;
	EFormat		eFormat;

	for( i = 0; c_asFormats[ i ].m_szExtension; ++i )
		if( !strcmp( szFile + strlen( szFile ) - strlen( c_asFormats[ i ].m_szExtension ),
			c_asFormats[ i ].m_szExtension ) )
			break;
	eFormat = c_asFormats[ i ].m_eFormat;
	fBinary = ( eFormat != EFormatText ) && ( eFormat != EFormatPCL );
	ifsm.open( szFile, fBinary ? ios_base::binary : ios_base::in );

	return OpenGenes( ifsm, fBinary, ( eFormat == EFormatPCL ) ); }

/*!
 * \brief
 * Open the genes from a CDat stored in the stream with the given format.
 * 
 * \param istm
 * Stream from which genes are loaded.
 * 
 * \param fBinary
 * If true, assume DAB format.
 * 
 * \param fPCL
 * If true, assume PCL format.
 * 
 * \returns
 * True if genes were successfully loaded.
 * 
 * Like Open, OpenGenes will load a CDat from the stream in the requested format, but it will load only the
 * CDat's size and gene list from the stream.  This is useful for rapidly obtaining gene lists and counts
 * from large CDats without incurring the memory/time penalty of loading or memory mapping the whole file.
 * 
 * \remarks
 * Attempting to access a CDat's data after opening only its genes is a poor idea.
 * 
 * \see
 * CPCL
 */
bool CDat::OpenGenes( std::istream& istm, bool fBinary, bool fPCL ) {

	return CDatImpl::OpenGenes( istm, fBinary, fPCL ); }

bool CDatImpl::OpenGenes( std::istream& istm, bool fBinary, bool fPCL ) {
	size_t		i, iToken;
	uint32_t	iCount;
	string		strToken, strCache;
	float		d;
	const char*	pc;
	char*		acBuf;

	Reset( );
	if( fPCL ) {
		m_pPCL = new CPCL( );
		if( m_pPCL->Open( istm ) ) {
			m_pMeasure = (IMeasure*)1;
			m_fMeasureMemory = false;
			return true; }
		return false; }
	acBuf = new char[ c_iBufferSize ];
	if( fBinary ) {
		istm.read( (char*)&iCount, sizeof(iCount) );
		if( iCount > c_iGeneLimit ) {
			delete[] acBuf;
			return false; }
		m_vecstrGenes.resize( iCount );
		for( i = 0; i < iCount; ++i ) {
			DabGene( istm, acBuf );
			m_vecstrGenes[ i ] = acBuf; } }
	else {
		set<string>					setstrGenes;
		set<string>::const_iterator	iterGenes;

		while( istm.peek( ) != EOF ) {
			istm.getline( acBuf, c_iBufferSize - 1 );
			for( iToken = 0; iToken < 3; ++iToken ) {
				strToken = OpenToken( acBuf, &pc );
				if( !strToken.length( ) )
					break;
				if( strToken != strCache ) {
					strCache = strToken;
					setstrGenes.insert( strToken ); }

				strToken = OpenToken( pc, &pc );
				setstrGenes.insert( strToken );
				d = (float)strtod( ( strToken = OpenToken( pc ) ).c_str( ), (char**)&pc );
				if( !d && ( ( pc - strToken.c_str( ) ) != strToken.length( ) ) ) {
					delete[] acBuf;
					return false; } } }
		m_vecstrGenes.reserve( setstrGenes.size( ) );
		for( iterGenes = setstrGenes.begin( ); iterGenes != setstrGenes.end( ); ++iterGenes )
			m_vecstrGenes.push_back( *iterGenes ); }
	delete[] acBuf;

	return true; }

/*!
 * \brief
 * Save a CDat to the given file, guessing the format from the file's extension.
 * 
 * \param szFile
 * Filename into which CDat is saved.
 * 
 * Save a CDat to the given file, guessing the format (DAT, DAB, or DAS) from the extension.  If null, the
 * CDat will be saved as a DAT to standard output.
 * 
 * \remarks
 * CDats cannot be saved to PCLs, only loaded from them.  If the extension is not recognized, DAB format is
 * assumed.
 * 
 * \see
 * Open
 */
void CDat::Save( const char* szFile ) const {
	size_t		i;
	EFormat		eFormat;
	ofstream	ofsm;

	if( !szFile ) {
		Save( cout, EFormatText );
		cout.flush( );
		return; }

	for( i = 0; c_asFormats[ i ].m_szExtension; ++i )
		if( !strcmp( szFile + strlen( szFile ) - strlen( c_asFormats[ i ].m_szExtension ),
			c_asFormats[ i ].m_szExtension ) )
			break;
	eFormat = c_asFormats[ i ].m_eFormat;
	ofsm.open( szFile, ( ( eFormat == EFormatText ) || ( eFormat == EFormatPCL ) ) ? ios_base::out :
		ios_base::binary );
	Save( ofsm, eFormat ); }

/*!
 * \brief
 * Save a CDat to the given stream in the requested format.
 * 
 * \param ostm
 * Stream into which CDat is saved.
 * 
 * \param eFormat
 * Format in which the CDat should be saved.
 * 
 * \remarks
 * The binary flag of the stream must match the intended format (text for DAT, binary for DAB/DAS).  CDats
 * cannot be saved to PCLs, only loaded from them.
 * 
 * \see
 * Open
 */
void CDat::Save( std::ostream& ostm, EFormat eFormat ) const {

	switch( eFormat ) {
		case EFormatText:
			SaveText( ostm );
			return;

		case EFormatSparse:
			SaveSparse( ostm );
			return; }

	SaveBinary( ostm ); }

void CDatImpl::SaveText( std::ostream& ostm ) const {
	size_t	i, j;
	float	d;

	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = Get( i, j ) ) )
				ostm << GetGene( i ) << '\t' << GetGene( j ) << '\t' << d << endl; }

void CDatImpl::SaveBinary( std::ostream& ostm ) const {
	size_t			i, j;
	const float*	pd;
	float			d;

	SaveGenes( ostm );
	if( m_pMeasure ) {
		for( i = 0; i < GetGenes( ); ++i )
			for( j = ( i + 1 ); j < GetGenes( ); ++j ) {
				d = Get( i, j );
				ostm.write( (char*)&d, sizeof(d) ); } }
	else
		for( i = 0; ( i + 1 ) < GetGenes( ); ++i ) {
			pd = m_Data.Get( i );
			ostm.write( (char*)pd, sizeof(*pd) * ( GetGenes( ) - i - 1 ) ); } }

void CDatImpl::SaveSparse( std::ostream& ostm ) const {
	uint32_t	i, j;
	float		d;

	SaveGenes( ostm );
	for( i = 0; i < GetGenes( ); ++i ) {
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = Get( i, j ) ) ) {
				ostm.write( (char*)&j, sizeof(j) );
				ostm.write( (char*)&d, sizeof(d) ); }
		j = -1;
		ostm.write( (char*)&j, sizeof(j) ); } }

void CDatImpl::SaveGenes( std::ostream& ostm ) const {
	size_t		i, j;
	uint32_t	iSize;
	string		strGene;

	iSize = GetGenes( );
	ostm.write( (char*)&iSize, sizeof(iSize) );
	for( i = 0; i < iSize; ++i ) {
		strGene = GetGene( i );
		for( j = 0; j < strGene.length( ); ++j ) {
			ostm.put( 0 );
			ostm.put( strGene[ j ] ); }
		ostm.put( 0 );
		ostm.put( 0 ); } }

/*!
 * \brief
 * Construct a new CDat with the given gene names.
 * 
 * \param vecstrGenes
 * Gene names and size to associate with the CDat.
 * 
 * \param fClear
 * If true, set each value of the new CDat to missing (NaN).
 * 
 * \param szFile
 * If non-null, use the given file as memory-mapped backing for the new CDat.
 * 
 * \returns
 * True if CDat was generated successfully.
 * 
 * Generates a CDat over the given genes (and of the same size), optionally initializing it to contain no
 * values or backing it with a memory-mapped file rather than physical memory.
 */
bool CDat::Open( const std::vector<std::string>& vecstrGenes, bool fClear, const char* szFile ) {
	size_t			i, j, iSize;
	unsigned char*	pb;

	Reset( );
	m_vecstrGenes.resize( vecstrGenes.size( ) );
	copy( vecstrGenes.begin( ), vecstrGenes.end( ), m_vecstrGenes.begin( ) );

	if( szFile ) {
		iSize = sizeof(uint32_t);
		for( i = 0; i < GetGenes( ); ++i )
			iSize += 2 * ( GetGene( i ).length( ) + 1 );
		iSize += CDistanceMatrix::GetSpace( GetGenes( ) );
		if( !CMeta::MapWrite( m_abData, m_hndlData, iSize, szFile ) )
			return false;
		*(uint32_t*)( pb = m_abData ) = GetGenes( );
		pb += sizeof(uint32_t);
		for( i = 0; i < GetGenes( ); ++i ) {
			const string&	strGene	= GetGene( i );

			for( j = 0; j < strGene.length( ); ++j ) {
				*pb++ = 0;
				*pb++ = strGene[ j ]; }
			*pb++ = 0;
			*pb++ = 0; }

		if( !OpenMemmap( pb ) )
			return false; }
	else
		m_Data.Initialize( GetGenes( ) );

	if( fClear )
		for( i = 0; i < m_vecstrGenes.size( ); ++i )
			for( j = ( i + 1 ); j < m_vecstrGenes.size( ); ++j )
				Set( i, j, CMeta::GetNaN( ) );

	return true; }

/*!
 * \brief
 * Construct a new CDat with the given gene names and values.
 * 
 * \param vecstrGenes
 * Gene names to associate with the CDat.
 * 
 * \param MatScores
 * Values to associate with the CDat.
 * 
 * \returns
 * True if CDat was generated successfully.
 * 
 * Generates a CDat over the given genes (and of the same size) and associate it with the given existing
 * pairwise values.
 * 
 * \remarks
 * vecstrGenes must be the same size as MatScores.  MatScores will not be copied; the new CDat will thus not
 * consume a substantial amount of memory, but MatScores must not be disposed of before the current CDat.
 */
bool CDat::Open( const std::vector<std::string>& vecstrGenes, const CDistanceMatrix& MatScores ) {
	size_t	i;

	if( vecstrGenes.size( ) != MatScores.GetSize( ) )
		return false;

	Reset( );
	m_vecstrGenes.reserve( vecstrGenes.size( ) );
	for( i = 0; i < vecstrGenes.size( ); ++i )
		m_vecstrGenes.push_back( vecstrGenes[ i ] );

	m_Data.Initialize( MatScores );

	return true; }

size_t CDatImpl::GetGene( const std::string& strGene ) const {
	size_t	i;

	if( m_pMeasure )
		return m_pPCL->GetGene( strGene );

	for( i = 0; i < GetGenes( ); ++i )
		if( m_vecstrGenes[ i ] == strGene )
			return i;

	return -1; }

/*!
 * \brief
 * Open a CDat stored in the given file, guessing the format from the file's extension.
 * 
 * \param szFile
 * Filename from which CDat is loaded.
 * 
 * \param fMemmap
 * If true, memory map file rather than allocating memory and copying its contents.
 * 
 * \param iSkip
 * If the given file is a PCL, the number of columns to skip between the ID and experiments.
 * 
 * \param fZScore
 * If true and the given file is a PCL, z-score similarity measures after pairwise calculation.
 * 
 * \returns
 * True if CDat was successfully opened.
 * 
 * Opens a CDat from the given file, guessing the file type from its extension: DAT, DAB, DAS, or PCL.
 * 
 * \remarks
 * A memory mapped CDat cannot be modified; doing so will generally cause a crash.  fMemmap is ignored
 * if the given file is not a DAB, and iSkip and fZScore are ignored if the given file is not a PCL.
 * PCLs are always converted to pairwise scores using CMeasurePearNorm.  If the extension is not
 * recognized, DAB format is assumed.
 * 
 * \see
 * Save | CPCL
 */
bool CDat::Open( const char* szFile, bool fMemmap, size_t iSkip, bool fZScore ) {
	ifstream	ifsm;
	EFormat		eFormat;
	size_t		i;

	if( !szFile )
		return Open( cin, EFormatText );

	for( i = 0; c_asFormats[ i ].m_szExtension; ++i )
		if( !strcmp( szFile + strlen( szFile ) - strlen( c_asFormats[ i ].m_szExtension ),
			c_asFormats[ i ].m_szExtension ) )
			break;
	eFormat = c_asFormats[ i ].m_eFormat;

	if( fMemmap && ( eFormat == EFormatBinary ) ) {
		Reset( );
		if( !CMeta::MapRead( m_abData, m_hndlData, m_iData, szFile ) ) {
			g_CatSleipnir.error( "CDat::Open( %s, %d ) failed memory mapping", szFile, fMemmap );
			return false; }
		return OpenHelper( ); }

	ifsm.open( szFile, ( ( eFormat == EFormatText ) || ( eFormat == EFormatPCL ) ) ? ios_base::in :
		ios_base::binary );
	if( !ifsm.is_open( ) )
		return false;
	return Open( ifsm, eFormat, (float)HUGE_VAL, false, iSkip, fZScore ); }

bool CDatImpl::OpenHelper( ) {
	unsigned char*	pb;
	size_t			i;

	m_vecstrGenes.resize( *(uint32_t*)( pb = m_abData ) );
	pb += sizeof(uint32_t);
	for( i = 0; i < GetGenes( ); ++i ) {
		string&	strGene	= m_vecstrGenes[ i ];

		while( *++pb )
			strGene += *pb++;
		pb++; }

	return OpenMemmap( pb ); }

bool CDatImpl::OpenMemmap( const unsigned char* pb ) {
	size_t	i;

	m_aadData = new float*[ GetGenes( ) - 1 ];
	m_aadData[ 0 ] = (float*)pb;
	for( i = 1; ( i + 1 ) < m_vecstrGenes.size( ); ++i )
		m_aadData[ i ] = m_aadData[ i - 1 ] + GetGenes( ) - i;
	m_Data.Initialize( GetGenes( ), (float**)m_aadData );

	return true; }

void CDatImpl::AveStd( double& dAve, double& dStd, size_t& iN, size_t iApproximate ) const {
	size_t	i, j;
	float	d;

	dAve = dStd = 0;
	if( iApproximate == -1 ) {
		for( iN = i = 0; i < GetGenes( ); ++i )
			for( j = ( i + 1 ); j < GetGenes( ); ++j )
				if( !CMeta::IsNaN( d = Get( i, j ) ) ) {
					iN++;
					dAve += d;
					dStd += d * d; } }
	else {
		size_t	iOne, iTwo;

		for( i = 0; i < iApproximate; ++i ) {
			iOne = rand( ) % GetGenes( );
			if( ( ( iTwo = rand( ) % GetGenes( ) ) == iOne ) ||
				CMeta::IsNaN( d = Get( iOne, iTwo ) ) ) {
				i--;
				continue; }
			dAve += d;
			dStd += d * d; }
		iN = i; }
	if( iN ) {
		dAve /= iN;
		dStd = sqrt( ( dStd / iN ) - ( dAve * dAve ) ); } }

void CDatImpl::NormalizeStdev( ) {
	double	d, dAve, dDev;
	size_t	i, j, iN;

	AveStd( dAve, dDev, iN );
	if( !( iN && dDev ) )
		return;
	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( (float)( d = Get( i, j ) ) ) )
				Set( i, j, (float)( ( d - dAve ) / dDev ) ); }

void CDatImpl::NormalizeSigmoid( ) {
	size_t	i, j;

	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			Set( i, j, 1.0f / ( 1 + exp( -Get( i, j ) ) ) ); }

void CDatImpl::NormalizeMinmax( ) {
	float	d, dMin, dMax;
	size_t	i, j;

	dMax = -( dMin = FLT_MAX );
	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = Get( i, j ) ) ) {
				if( d < dMin )
					dMin = d;
				if( d > dMax )
					dMax = d; }
	if( dMax -= dMin )
		for( i = 0; i < GetGenes( ); ++i )
			for( j = ( i + 1 ); j < GetGenes( ); ++j )
				Set( i, j, ( Get( i, j ) - dMin ) / dMax ); }

/*!
 * \brief
 * Replace each finite value in the CDat with one minus that value.
 * 
 * \remarks
 * Most useful in combination with CDat::Normalize ( true ).
 */
void CDat::Invert( ) {
	size_t	i, j;
	float	d;

	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = Get( i, j ) ) )
				Set( i, j, 1 - d ); }

/*!
 * \brief
 * Remove edges from the CDat based on the given gene file and filter type.
 * 
 * \param szGenes
 * File from which gene names are loaded, one per line.
 * 
 * \param eFilter
 * Way in which to use the given genes to remove edges.
 * 
 * \param iLimit
 * For EFilterPixie and EFilterHefalmp, the maximum number of genes to retain.
 * 
 * \returns
 * True if the filter was executed successfully.
 * 
 * Remove edges and nodes (by removing all incident edges) from the CDat based on one of several algorithms.
 * For details, see EFilter.
 * 
 * \remarks
 * EFilterTerm and to some degree EFilterEdge don't make a lot of sense for CDats that do not represent
 * gold standards.
 */
bool CDat::FilterGenes( const char* szGenes, EFilter eFilter, size_t iLimit ) {
	CGenome		Genome;
	CGenes		Genes( Genome );
	ifstream	ifsm;

	if( !szGenes )
		return false;

	ifsm.open( szGenes );
	if( !Genes.Open( ifsm ) )
		return false;
	FilterGenes( Genes, eFilter, iLimit );
	return true; }

/*!
 * \brief
 * Remove edges from the CDat based on the given gene set and filter type.
 * 
 * \param Genes
 * Gene set used to filter the CDat.
 * 
 * \param eFilter
 * Way in which to use the given genes to remove edges.
 * 
 * \param iLimit
 * For EFilterPixie and EFilterHefalmp, the maximum number of genes to retain.
 * 
 * \param dEdgeAggressiveness
 * For EFilterPixie and EFilterHefalmp, higher values result in more aggressive edge trimming.  NaN
 * completely skips edge trimming.
 * 
 * Remove edges and nodes (by removing all incident edges) from the CDat based on one of several algorithms.
 * For details, see EFilter.
 * 
 * \remarks
 * EFilterTerm and to some degree EFilterEdge don't make a lot of sense for CDats that do not represent
 * gold standards.
 */
void CDat::FilterGenes( const CGenes& Genes, EFilter eFilter, size_t iLimit, float dEdgeAggressiveness ) {
	size_t			i, j;
	vector<bool>	vecfGenes;

	vecfGenes.resize( GetGenes( ) );
	for( i = 0; i < Genes.GetGenes( ); ++i )
		if( ( j = GetGene( Genes.GetGene( i ).GetName( ) ) ) != -1 )
			vecfGenes[ j ] = true;

	switch( eFilter ) {
		case EFilterPixie:
		case EFilterHefalmp:
			FilterGenesGraph( Genes, vecfGenes, iLimit, dEdgeAggressiveness, eFilter == EFilterHefalmp );
			return; }

	for( i = 0; i < GetGenes( ); ++i ) {
		if( ( ( eFilter == EFilterExclude ) && vecfGenes[ i ] ) ||
			( ( eFilter == EFilterInclude ) && !vecfGenes[ i ] ) ) {
			for( j = ( i + 1 ); j < GetGenes( ); ++j )
				Set( i, j, CMeta::GetNaN( ) );
			continue; }
		if( ( eFilter == EFilterEdge ) && vecfGenes[ i ] )
			continue;
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			switch( eFilter ) {
				case EFilterInclude:
				case EFilterEdge:
					if( !vecfGenes[ j ] )
						Set( i, j, CMeta::GetNaN( ) );
					break;

				case EFilterTerm:
					if( !( vecfGenes[ i ] && vecfGenes[ j ] ) &&
						( !( vecfGenes[ i ] || vecfGenes[ j ] ) || Get( i, j ) ) )
						Set( i, j, CMeta::GetNaN( ) );
					break;

				case EFilterExclude:
					if( vecfGenes[ j ] )
						Set( i, j, CMeta::GetNaN( ) );
					break; } } }

struct SPixie {
	size_t	m_iNode;
	float	m_dScore;

	SPixie( size_t iNode, float dScore ) : m_iNode(iNode), m_dScore(dScore) { }

	bool operator<( const SPixie& sPixie ) const {

		return ( m_dScore < sPixie.m_dScore ); }
};

void CDatImpl::FilterGenesGraph( const CGenes& Genes, vector<bool>& vecfGenes, size_t iLimit,
	float dEdgeAggressiveness, bool fHefalmp ) {
	vector<float>				vecdNeighbors;
	size_t						i, j, iOne, iTwo, iMinOne, iMinTwo, iN;
	vector<size_t>				veciGenes, veciFinal, veciDegree;
	set<size_t>					setiN;
	set<size_t>::const_iterator	iterN;
	float						d, dMin, dCutoff;
	priority_queue<SPixie>		pqueNeighbors;
	bool						fDone;
	double						dAve, dDev;

	if( iLimit == -1 )
		iLimit = c_iNeighborhood;
	veciGenes.resize( Genes.GetGenes( ) );
	for( i = 0; i < veciGenes.size( ); ++i )
		veciGenes[ i ] = GetGene( Genes.GetGene( i ).GetName( ) );

	vecdNeighbors.resize( GetGenes( ) );
	fill( vecdNeighbors.begin( ), vecdNeighbors.end( ), 0.0f );
	if( fHefalmp )
		for( i = 0; i < GetGenes( ); ++i ) {
			size_t	iIn, iOut;
			float	dIn, dOut;

			if( vecfGenes[ i ] )
				continue;
			dIn = dOut = 0;
			for( iIn = j = 0; j < veciGenes.size( ); ++j ) {
				if( ( iOne = veciGenes[ j ] ) == -1 )
					continue;
				if( !CMeta::IsNaN( d = Get( i, iOne ) ) ) {
					iIn++;
					dIn += d; } }
			for( iOut = j = 0; j < GetGenes( ); ++j )
				if( !CMeta::IsNaN( d = Get( i, j ) ) ) {
					iOut++;
					dOut += d; }
			vecdNeighbors[ i ] = ( iIn && dOut ) ? ( dIn * iOut / iIn / dOut ) : 0; }
	else
		for( i = 0; i < veciGenes.size( ); ++i ) {
			if( ( iOne = veciGenes[ i ] ) == -1 )
				continue;
			for( j = 0; j < GetGenes( ); ++j ) {
				if( vecfGenes[ j ] )
					continue;
				if( !CMeta::IsNaN( d = Get( iOne, j ) ) )
					vecdNeighbors[ j ] += d; } }
	for( i = 0; i < vecdNeighbors.size( ); ++i )
		if( ( d = vecdNeighbors[ i ] ) > 0 )
			pqueNeighbors.push( SPixie( i, d ) );

	for( i = 0; i < veciGenes.size( ); ++i )
		if( ( iOne = veciGenes[ i ] ) != -1 )
			veciFinal.push_back( iOne );
	while( !pqueNeighbors.empty( ) && ( setiN.size( ) < iLimit ) ) {
		veciFinal.push_back( pqueNeighbors.top( ).m_iNode );
		setiN.insert( pqueNeighbors.top( ).m_iNode );
		pqueNeighbors.pop( ); }

	for( iterN = setiN.begin( ); iterN != setiN.end( ); ++iterN )
		vecfGenes[ *iterN ] = true;
	for( i = 0; i < GetGenes( ); ++i ) {
		if( !vecfGenes[ i ] ) {
			for( j = ( i + 1 ); j < GetGenes( ); ++j )
				Set( i, j, CMeta::GetNaN( ) );
			continue; }
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !vecfGenes[ j ] )
				Set( i, j, CMeta::GetNaN( ) ); }
	AveStd( dAve, dDev, iN );
	dCutoff = (float)( dAve + ( dEdgeAggressiveness * dDev ) );

	veciDegree.resize( veciFinal.size( ) );
	for( i = 0; i < veciDegree.size( ); ++i ) {
		iOne = veciFinal[ i ];
		for( j = ( i + 1 ); j < veciDegree.size( ); ++j ) {
			iTwo = veciFinal[ j ];
			if( !CMeta::IsNaN( Get( iOne, iTwo ) ) ) {
				veciDegree[ i ]++;
				veciDegree[ j ]++; } } }
	for( fDone = CMeta::IsNaN( dEdgeAggressiveness ); !fDone; ) {
		fDone = true;
		dMin = FLT_MAX;
		for( i = 0; i < veciFinal.size( ); ++i ) {
			iOne = veciFinal[ i ];
			for( j = ( i + 1 ); j < veciFinal.size( ); ++j ) {
				iTwo = veciFinal[ j ];
				if( !CMeta::IsNaN( d = Get( iOne, iTwo ) ) && ( d < dCutoff ) && ( d < dMin ) &&
					( veciDegree[ i ] > c_iDegree ) && ( veciDegree[ j ] > c_iDegree ) ) {
					fDone = false;
					dMin = d;
					iMinOne = i;
					iMinTwo = j; } } }
		if( !fDone ) {
			veciDegree[ iMinOne ]--;
			veciDegree[ iMinTwo ]--;
			Set( veciFinal[ iMinOne ], veciFinal[ iMinTwo ], CMeta::GetNaN( ) ); } } }

/*!
 * \brief
 * Generate a DOT-formatted graph from the CDat, suitable for processing with AT&T's Graphviz software.
 * 
 * \param ostm
 * Stream into which DOT is saved.
 * 
 * \param dCutoff
 * If finite, edge weights below this cutoff are not included in the DOT.
 * 
 * \param pGenome
 * If non-null, name synonyms from this genome are used to label gene nodes.
 * 
 * \param fUnlabeled
 * If true, do not label gene nodes and use only minimal unique identifiers when generating the DOT.
 * 
 * \param fHashes
 * If true, include # marks in hexadecimal color strings within the DOT.  This is moronic, but some DOT
 * parsers (*coughboostcough*) will randomly crash if # characters are included in the DOT.
 * 
 * \param pvecdColors
 * If non-null, contains weights between 0 and 1 interpolating node colors between cyan and yellow.
 * 
 * \param pvecdBorders
 * If non-null, contains border widths in pixels for each node.
 * 
 * SaveDOT is one of the most useful methods for visualizing the weighted graph implicit in a CDat,
 * particularly in combination with FilterGenes and EFilterPixie/EFilterHefalmp.  Calling SaveDOT will
 * generate a DOT graph file containing (optionally) each node and edge in the CDat, with edges colored by
 * weight (scaled from green for the minimum weight to red at the maximum).  Nodes can optionally be colored
 * or given varying border widths based on external information, or labeled with different gene names
 * (synonyms) than used internally by the CDat.  DOT files can be visualized by a variety of software, notably
 * the Graphviz package from AT&T.
 * 
 * \remarks
 * If given, pvecdColors and pvecdBorders must be of the same size as the CDat.
 * 
 * \see
 * SaveGDF | SaveNET | SaveMATISSE
 */
void CDat::SaveDOT( std::ostream& ostm, float dCutoff, const CGenome* pGenome, bool fUnlabeled, bool fHashes,
	const std::vector<float>* pvecdColors, const std::vector<float>* pvecdBorders ) const {
	size_t			i, j;
	float			d, dMin, dMax;
	bool			fAll, fLabel;
	vector<string>	vecstrNames;
	vector<bool>	vecfGenes;

	fAll = CMeta::IsNaN( dCutoff );
	ostm << "graph G {" << endl;
	ostm << "	pack = \"true\";" << endl;
	ostm << "	overlap = \"scale\";" << endl;
	ostm << "	splines = \"true\";" << endl;

	if( pvecdColors || pvecdBorders || !( fUnlabeled || fAll ) ) {
		vecfGenes.resize( GetGenes( ) );
		for( i = 0; i < vecfGenes.size( ); ++i )
			for( j = ( i + 1 ); j < vecfGenes.size( ); ++j )
				if( !CMeta::IsNaN( d = Get( i, j ) ) && ( fAll || ( d >= dCutoff ) ) )
					vecfGenes[ i ] = vecfGenes[ j ] = true; }

	vecstrNames.resize( GetGenes( ) );
	for( i = 0; i < vecstrNames.size( ); ++i ) {
		string	strName;

		fLabel = !fUnlabeled && ( fAll || vecfGenes[ i ] );
		vecstrNames[ i ] = CMeta::Filename( strName = GetGene( i ) );
		if( !isalpha( vecstrNames[ i ][ 0 ] ) )
			vecstrNames[ i ] = "_" + vecstrNames[ i ];
		if( pGenome && ( ( j = pGenome->GetGene( GetGene( i ) ) ) != -1 ) ) {
			const CGene&	Gene	= pGenome->GetGene( j );

			strName = Gene.GetSynonyms( ) ? Gene.GetSynonym( 0 ) : Gene.GetName( );
			if( fUnlabeled ) {
				if( strName != vecstrNames[ i ] ) {
					vecstrNames[ i ] += '_';
					vecstrNames[ i ] += Gene.GetSynonym( 0 ); } } }
		if( ( pvecdColors || pvecdBorders || fLabel ) && ( fAll || vecfGenes[ i ] ) ) {
			ostm << vecstrNames[ i ] << " [shape = \"ellipse\", style = \"filled";
			if( pvecdBorders )
				ostm << ", setlinewidth(" << (*pvecdBorders)[ i ] << ")";
			ostm << "\"";
			ostm << ", fillcolor = \"" << ( fHashes ? "#" : "" ) << ( pvecdColors ? CColor::Interpolate(
				(*pvecdColors)[ i ], CColor::c_Cyan, CColor::c_White, CColor::c_Yellow ).ToRGB( ) :
				"FFFFFF" ) << "\"";
			if( fLabel )
				ostm << ", label=\"" << strName << "\"";
			ostm << "];" << endl; } }

	ostm << endl;
	dMin = FLT_MAX;
	dMax = -FLT_MAX;
	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = Get( i, j ) ) && ( fAll || ( d >= dCutoff ) ) ) {
				if( d < dMin )
					dMin = d;
				if( d > dMax )
					dMax = d; }
	dMax -= dMin;
	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = Get( i, j ) ) && ( fAll || ( d >= dCutoff ) ) )
				ostm << vecstrNames[ i ] << " -- " << vecstrNames[ j ] << " [weight = " << d <<
					", color = \"" << ( fHashes ? "#" : "" ) << CColor::Interpolate( ( d - dMin ) / dMax,
					CColor::c_Green, CColor::c_Black, CColor::c_Red ).ToRGB( ) << "\"];" << endl;

	ostm << "}" << endl; }

/*!
 * \brief
 * Generate a GDF-formatted graph from the CDat, suitable for processing with the GUESS graph exploration
 * system.
 * 
 * \param ostm
 * Stream into which GDF is saved.
 * 
 * \param dCutoff
 * If finite, edge weights below this cutoff are not included in the GDF.
 * 
 * Calling SaveGDF will generate a GDF graph file containing (optionally) each node and edge in the CDat.
 * Nodes are labeled minimally, and edges are not given any inherent color information (although this can be
 * added later).  GDFs are visualizable using the GUESS graph exploration software, although it doesn't deal
 * well with large or highly connected graphs.
 * 
 * \remarks
 * Edge weights are scaled by 10x in order to fall roughly within GDF's expected range.
 * 
 * \see
 * SaveDOT | SaveNET | SaveMATISSE
 */
void CDat::SaveGDF( std::ostream& ostm, float dCutoff ) const {
	size_t			i, j;
	float			d;
	bool			fAll;
	vector<string>	vecstrNames;

	fAll = CMeta::IsNaN( dCutoff );
	ostm << "nodedef> name" << endl;

	vecstrNames.resize( GetGenes( ) );
	for( i = 0; i < vecstrNames.size( ); ++i )
		vecstrNames[ i ] = CMeta::Filename( GetGene( i ) );
	for( i = 0; i < GetGenes( ); ++i )
		ostm << vecstrNames[ i ] << endl;

	ostm << endl << "edgedef> node1,node2,weight" << endl;
	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = Get( i, j ) ) && ( fAll || ( d >= dCutoff ) ) )
				ostm << vecstrNames[ i ] << "," << vecstrNames[ j ] << "," << ( 10 * Get( i, j ) ) << endl; }

/*!
 * \brief
 * Generate a NET-formatted graph from the CDat.
 * 
 * \param ostm
 * Stream into which NET is saved.
 * 
 * \param dCutoff
 * If finite, edge weights below this cutoff are not included in the NET.
 * 
 * Calling SaveNET will generate a NET graph file containing (optionally) each node and edge in the CDat.
 * Nodes are labeled minimally, and edges are not given any inherent color information (although this can be
 * added later).  I honestly don't remember at all what software uses the NET format, and "net"'s 
 * impossible to Google these days (thanks, Microsoft...)
 * 
 * \see
 * SaveDOT | SaveGDF | SaveMATISSE
 */
void CDat::SaveNET( std::ostream& ostm, float dCutoff ) const {
	size_t	i, j;
	float	d;
	bool	fAll;

	fAll = CMeta::IsNaN( dCutoff );
	ostm << "*Vertices " << GetGenes( ) << endl;
	for( i = 0; i < GetGenes( ); ++i )
		ostm << ( i + 1 ) << " \"" << GetGene( i ) << '"' << endl;

	ostm << "*Edges" << endl;
	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = Get( i, j ) ) && ( fAll || ( d >= dCutoff ) ) )
				ostm << ( i + 1 ) << ' ' << ( j + 1 ) << ' ' << d << endl; }

/*!
 * \brief
 * Generate a MATISSE-formatted graph from the CDat, suitable for processing with the MATISSE software.
 * 
 * \param ostm
 * Stream into which MATISSE file is saved.
 * 
 * \param dCutoff
 * If finite, edge weights below this cutoff are not included in the MATISSE file.
 * 
 * \param pGenome
 * If non-null, name synonyms from this genome are used to label gene nodes.
 * 
 * Calling SaveMATISSE will generate a MATISSE graph file containing (optionally) each node and edge in the
 * CDat.  Nodes are labeled either minimally or using the given genome's names/sysnonyms, and edges are not
 * given any inherent color information (although this can be added later).  MATISSE files are visualizable
 * using the MATISSE software, although it doesn't deal well with large or highly connected graphs.
 * 
 * \see
 * SaveDOT | SaveGDF | SaveNET
 */
void CDat::SaveMATISSE( std::ostream& ostm, float dCutoff, const CGenome* pGenome ) const {
	size_t	i, j, k;
	float	d;
	bool	fAll;

	fAll = CMeta::IsNaN( dCutoff );
	for( i = 0; i < GetGenes( ); ++i ) {
		j = pGenome ? pGenome->GetGene( GetGene( i ) ) : -1;
		ostm << i << '\t' << GetGene( i ) << '\t' << ( ( ( j == -1 ) ||
			!pGenome->GetGene( j ).GetSynonyms( ) ) ? GetGene( i ) :
			pGenome->GetGene( j ).GetSynonym( 0 ) ) << '\t' << ( ( j == -1 ) ? "-" :
			pGenome->GetGene( j ).GetGloss( ) ) << endl; }

	for( k = i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = Get( i, j ) ) && ( fAll || ( d >= dCutoff ) ) )
				ostm << k++ << '\t' << i << '\t' << j << "	false	" << d << endl; }

struct SSorterRank {
	const CDat&	m_Dat;

	SSorterRank( const CDat& Dat ) : m_Dat(Dat) { }

	bool operator()( const pair<size_t, size_t>& prOne, const pair<size_t, size_t>& prTwo ) const {

		return ( m_Dat.Get( prOne.first, prOne.second ) < m_Dat.Get( prTwo.first, prTwo.second ) ); }

};

/*!
 * \brief
 * Replace each finite value in the CDat with its integer rank.
 * 
 * \remarks
 * This can be costly in memory/processing time for large CDats.
 */
void CDat::Rank( ) {
	vector<pair<size_t, size_t> >	vecprData;
	size_t							i, j, iRank;
	float							d, dPrev;

	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( !CMeta::IsNaN( Get( i, j ) ) )
				vecprData.push_back( pair<size_t, size_t>( i, j ) );
	sort( vecprData.begin( ), vecprData.end( ), SSorterRank( *this ) );
	for( iRank = i = 0; i < vecprData.size( ); ++i ) {
		d = Get( vecprData[ i ].first, vecprData[ i ].second );
		if( i && ( d != dPrev ) )
			iRank = i;
		dPrev = d;
		Set( vecprData[ i ].first, vecprData[ i ].second, (float)iRank ); } }

}

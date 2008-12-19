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
#include "dataset.h"
#include "bayesnetint.h"
#include "genome.h"
#include "compactmatrix.h"

namespace Sleipnir {

const char	CDataImpl::c_szDat[]	= ".dat";
const char	CDataImpl::c_szDab[]	= ".dab";

void CDataImpl::FilterGenes( IDataset* pData, const CGenes& Genes, CDat::EFilter eFilt ) {
	vector<bool>	vecfGenes;
	size_t			i, j;

	if( !Genes.GetGenes( ) )
		return;

	vecfGenes.resize( pData->GetGenes( ) );
	for( i = 0; i < vecfGenes.size( ); ++i )
		vecfGenes[ i ] = Genes.IsGene( pData->GetGene( i ) );

	for( i = 0; i < vecfGenes.size( ); ++i ) {
		if( ( ( eFilt == CDat::EFilterInclude ) && !vecfGenes[ i ] ) ||
			( ( eFilt == CDat::EFilterExclude ) && vecfGenes[ i ] ) ) {
			for( j = ( i + 1 ); j < vecfGenes.size( ); ++j )
				pData->Remove( i, j );
			continue; }
		if( ( eFilt == CDat::EFilterEdge ) && vecfGenes[ i ] )
			continue;
		for( j = ( i + 1 ); j < vecfGenes.size( ); ++j )
			switch( eFilt ) {
				case CDat::EFilterInclude:
				case CDat::EFilterEdge:
					if( !vecfGenes[ j ] )
						pData->Remove( i, j );
					break;

				case CDat::EFilterTerm:
					if( !( vecfGenes[ i ] && vecfGenes[ j ] ) &&
						( !( vecfGenes[ i ] || vecfGenes[ j ] ) || pData->GetDiscrete( i, j, 0 ) ) )
							pData->Remove( i, j );
					break;

				case CDat::EFilterExclude:
					if( vecfGenes[ j ] )
						pData->Remove( i, j );
					break; } } }

size_t CDataImpl::OpenMax( const char* szDataDir, const std::vector<std::string>& vecstrNodes,
	bool fAnswers, std::vector<std::string>& vecstrData, std::set<std::string>* psetstrGenes ) {
	size_t		i, iLength, iMap, iRet;
	string		strFile;
	ifstream	ifsm;
	CPCL		PCL;

	strFile = szDataDir;
	strFile += c_cSeparator;
	iLength = strFile.size( );

	iRet = 0;
	m_veciMapping.resize( vecstrNodes.size( ) );
	m_veciMapping[ 0 ] = fAnswers ? 0 : -1;
	iMap = fAnswers ? 1 : 0;
	for( i = 1; i < vecstrNodes.size( ); ++i ) {
		m_veciMapping[ i ] = -1;
		strFile.resize( iLength );
		strFile += vecstrNodes[ i ];
		strFile += c_szDab;
		ifsm.clear( );
		ifsm.open( strFile.c_str( ), ios_base::binary );
		if( ifsm.is_open( ) ) {
			iRet++;
			m_veciMapping[ i ] = iMap++;
			vecstrData.push_back( strFile );
			if( psetstrGenes )
				OpenGenes( ifsm, true, false, *psetstrGenes ); }
		else {
			strFile.resize( strFile.length( ) - strlen( c_szDab ) );
			strFile += c_szDat;
			ifsm.clear( );
			ifsm.open( strFile.c_str( ) );
			if( ifsm.is_open( ) ) {
				iRet++;
				m_veciMapping[ i ] = iMap++;
				vecstrData.push_back( strFile );
				if( psetstrGenes )
					OpenGenes( ifsm, false, false, *psetstrGenes ); }
			else {
				strFile.resize( strFile.length( ) - strlen( c_szDat ) );
				strFile += CPCL::GetExtension( );
				ifsm.clear( );
				ifsm.open( strFile.c_str( ) );
				if( ifsm.is_open( ) ) {
					iRet++;
					m_veciMapping[ i ] = iMap++;
					vecstrData.push_back( strFile );
					if( psetstrGenes )
						OpenGenes( ifsm, false, true, *psetstrGenes ); }
				else {
					g_CatSleipnir.info( "CDataImpl::OpenMax( %s ) assuming %s is hidden",
						szDataDir, vecstrNodes[ i ].c_str( ) );
					continue; } } }
		ifsm.close( ); }

	return iRet; }

bool CDataImpl::OpenGenes( std::istream& istm, bool fBinary, bool fPCL,
	std::set<std::string>& setstrGenes ) const {
	CDat	Dat;
	size_t	i;

	if( !Dat.OpenGenes( istm, fBinary, fPCL ) )
		return false;
	for( i = 0; i < Dat.GetGenes( ); ++i )
		setstrGenes.insert( Dat.GetGene( i ) );
	return true; }

bool CDataImpl::OpenGenes( const std::vector<std::string>& vecstrData ) {
	size_t						i;
	ifstream					ifsm;
	set<string>					setstrGenes;
	set<string>::const_iterator	iterGenes;

	m_veciMapping.resize( vecstrData.size( ) );
	m_veccQuants.resize( vecstrData.size( ) );
	for( i = 0; i < vecstrData.size( ); ++i ) {
		m_veciMapping[ i ] = i;
		ifsm.clear( );
		ifsm.open( vecstrData[ i ].c_str( ), ios_base::binary );
		if( !( ifsm.is_open( ) && OpenGenes( ifsm, true, false, setstrGenes ) ) ) {
			ifsm.close( );
			ifsm.clear( );
			ifsm.open( vecstrData[ i ].c_str( ) );
			if( !( ifsm.is_open( ) && OpenGenes( ifsm, false, false, setstrGenes ) ) ) {
				ifsm.close( );
				ifsm.clear( );
				ifsm.open( vecstrData[ i ].c_str( ) );
				if( !( ifsm.is_open( ) && OpenGenes( ifsm, false, true, setstrGenes ) ) ) {
					g_CatSleipnir.error( "CDataImpl::OpenGenes( ) failed to open: %s", vecstrData[ i ].c_str( ) );
					return false; } } }
		ifsm.close( ); }

	m_vecstrGenes.resize( setstrGenes.size( ) );
	i = 0;
	for( iterGenes = setstrGenes.begin( ); iterGenes != setstrGenes.end( ); ++iterGenes )
		m_vecstrGenes[ i++ ] = *iterGenes;

	return true; }

size_t CDataImpl::GetGene( const std::string& strGene ) const {
	size_t	i;

	for( i = 0; i < m_vecstrGenes.size( ); ++i )
		if( m_vecstrGenes[ i ] == strGene )
			return i;

	return -1; }

const unsigned char* CDataImpl::OpenBinary( const unsigned char* pbData ) {
	uint32_t	iVal;
	size_t		i;

	m_fContinuous = !!*(uint32_t*)pbData;
	pbData += sizeof(uint32_t);

	iVal = *(uint32_t*)pbData;
	pbData += sizeof(iVal);
	m_veciMapping.resize( iVal );
	for( i = 0; i < m_veciMapping.size( ); ++i ) {
		iVal = *(uint32_t*)pbData;
		m_veciMapping[ i ] = ( iVal == -1 ) ? (size_t)-1 : iVal;
		pbData += sizeof(iVal); }

	iVal = *(uint32_t*)pbData;
	pbData += 2 * sizeof(iVal);
	m_vecstrGenes.resize( iVal );
	for( i = 0; i < m_vecstrGenes.size( ); ++i ) {
		m_vecstrGenes[ i ] = (char*)pbData;
		pbData += m_vecstrGenes[ i ].length( ) + 1; }

	iVal = *(uint32_t*)pbData;
	pbData += sizeof(iVal);
	m_veccQuants.resize( iVal );
	for( i = 0; i < m_veccQuants.size( ); ++i )
		m_veccQuants[ i ] = *pbData++;

	return pbData; }

bool CDataImpl::OpenBinary( std::istream& istm ) {
	uint32_t	iVal;
	size_t		i, j;
	char*		ac;

	istm.read( (char*)&iVal, sizeof(iVal) );
	m_fContinuous = !!iVal;

	istm.read( (char*)&iVal, sizeof(iVal) );
	m_veciMapping.resize( iVal );
	for( i = 0; i < m_veciMapping.size( ); ++i ) {
		istm.read( (char*)&iVal, sizeof(iVal) );
		m_veciMapping[ i ] = ( iVal == -1 ) ? (size_t)-1 : iVal; }

	istm.read( (char*)&iVal, sizeof(iVal) );
	m_vecstrGenes.resize( iVal );
	istm.read( (char*)&iVal, sizeof(iVal) );
	ac = new char[ iVal ];
	istm.read( ac, iVal );
	for( i = j = 0; i < m_vecstrGenes.size( ); ++i ) {
		m_vecstrGenes[ i ] = ac + j;
		j += m_vecstrGenes[ i ].length( ) + 1; }
	delete[] ac;

	istm.read( (char*)&iVal, sizeof(iVal) );
	ac = new char[ iVal ];
	m_veccQuants.resize( iVal );
	istm.read( ac, iVal );
	copy( ac, ac + iVal, m_veccQuants.begin( ) );
	delete[] ac;

	return true; }

void CDataImpl::SaveBinary( std::ostream& ostm ) const {
	size_t		i;
	uint32_t	iVal;
	char*		ac;
	char		c;

	iVal = m_fContinuous;
	ostm.write( (char*)&iVal, sizeof(iVal) );

	iVal = (uint32_t)m_veciMapping.size( );
	ostm.write( (char*)&iVal, sizeof(iVal) );
	for( i = 0; i < m_veciMapping.size( ); ++i ) {
		iVal = ( m_veciMapping[ i ] == -1 ) ? -1 :
			(uint32_t)m_veciMapping[ i ];
		ostm.write( (char*)&iVal, sizeof(iVal) ); }

	iVal = (uint32_t)m_vecstrGenes.size( );
	ostm.write( (char*)&iVal, sizeof(iVal) );
	for( i = iVal = 0; i < m_vecstrGenes.size( ); ++i )
		iVal += (uint32_t)m_vecstrGenes[ i ].length( ) + 1;
	ostm.write( (char*)&iVal, sizeof(iVal) );
	for( i = c = 0; i < m_vecstrGenes.size( ); ++i ) {
		ostm.write( m_vecstrGenes[ i ].c_str( ), (streamsize)m_vecstrGenes[ i ].length( ) );
		ostm.write( &c, sizeof(c) ); }

	ac = new char[ iVal = (uint32_t)m_veccQuants.size( ) ];
	for( i = 0; i < iVal; ++i )
		ac[ i ] = m_veccQuants[ i ];
	ostm.write( (char*)&iVal, sizeof(iVal) );
	ostm.write( ac, iVal );
	delete[] ac; }

CDatasetImpl::CDatasetImpl( ) : m_apData(NULL) { }

CDatasetImpl::~CDatasetImpl( ) {

	Reset( ); }

void CDatasetImpl::Reset( ) {
	size_t	i;

	if( m_apData ) {
		for( i = 0; i < m_veccQuants.size( ); ++i )
			if( m_veccQuants[ i ] == (unsigned char)-1 )
				delete (CDistanceMatrix*)m_apData[ i ];
			else
				delete (CCompactMatrix*)m_apData[ i ];
		delete[] m_apData; } }

void CDatasetImpl::SaveBinary( std::ostream& ostm ) const {
	size_t		i;
	uint32_t	iData;

	CDataImpl::SaveBinary( ostm );
	for( i = iData = 0; i < m_veciMapping.size( ); ++i )
		if( m_veciMapping[ i ] != -1 )
			iData++;
	ostm.write( (char*)&iData, sizeof(iData) );
	for( i = 0; i < iData; ++i )
		if( m_veccQuants[ i ] == (unsigned char)-1 )
			((CDistanceMatrix*)m_apData[ i ])->Save( ostm, true );
		else
			((CCompactMatrix*)m_apData[ i ])->Save( ostm ); }

void CDatasetImpl::SaveText( std::ostream& ostm ) const {
	size_t			i, j, k;
	vector<float>	vecdValues;
	bool			fHit;

	vecdValues.resize( GetExperiments( ) );
	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j ) {
			fHit = false;
			for( k = 0; k < vecdValues.size( ); ++k )
				if( !CMeta::IsNaN( vecdValues[ k ] = GetContinuous( i, j, k ) ) )
					fHit = true;
			if( !fHit )
				continue;
			ostm << GetGene( i ) << '\t' << GetGene( j );
			for( k = 0; k < vecdValues.size( ); ++k ) {
				ostm << '\t';
				if( !CMeta::IsNaN( vecdValues[ k ] ) )
					ostm << vecdValues[ k ]; }
			ostm << endl; } }

/*!
 * \brief
 * Construct a dataset corresponding to the given Bayes net using the provided answer file and data
 * files from the given directory.
 * 
 * \param szAnswerFile
 * Answer file which will become the first node of the dataset.
 * 
 * \param szDataDirectory
 * Directory from which data files are loaded.
 * 
 * \param pBayesNet
 * Bayes nets whose nodes will correspond to files in the dataset.
 * 
 * \returns
 * True if dataset was constructed successfully.
 * 
 * Creates a dataset with nodes corresponding to the given Bayes net structure; the given answer file
 * is always inserted as the first (0th) data file, and thus corresponds to the first node in the Bayes
 * net (generally the class node predicting functional relationships).  Data is loaded continuously or
 * discretely as indicated by the Bayes net, and nodes for which a corresponding data file (i.e. one
 * with the same name followed by an appropriate CDat extension) cannot be located are marked as hidden.
 * 
 * \remarks
 * Each data file is loaded more or less as-is; continuous data files will be loaded directly into
 * memory, and discrete files are pre-discretized and stored in compact matrices.
 */
bool CDataset::Open( const char* szAnswerFile, const char* szDataDirectory, const IBayesNet* pBayesNet ) {
	CDataPair	Answers;

	return ( Answers.Open( szAnswerFile, pBayesNet->IsContinuous( 0 ) ) &&
		Open( Answers, szDataDirectory, pBayesNet ) ); }

bool CDatasetImpl::Open( const CDataPair* pAnswers, const char* szDataDir, const IBayesNet* pBayesNet ) {
	size_t						i;
	vector<string>				vecstrData, vecstrNodes;
	set<string>					setstrGenes;
	set<string>::const_iterator	iterGene;

	Reset( );
	m_fContinuous = pBayesNet->IsContinuous( );
	pBayesNet->GetNodes( vecstrNodes );
	m_veccQuants.resize( ( pAnswers ? 1 : 0 ) + OpenMax( szDataDir, vecstrNodes, !!pAnswers,
		vecstrData, &setstrGenes ) );
	if( pAnswers ) {
		m_vecstrGenes.resize( pAnswers->GetGenes( ) );
		for( i = 0; i < m_vecstrGenes.size( ); ++i )
			m_vecstrGenes[ i ] = pAnswers->GetGene( i ); }
	else {
		m_vecstrGenes.resize( setstrGenes.size( ) );
		for( i = 0,iterGene = setstrGenes.begin( ); iterGene != setstrGenes.end( ); ++i,++iterGene )
			m_vecstrGenes[ i ] = *iterGene; }
	m_apData = new void*[ m_veccQuants.size( ) ];
	if( pAnswers && !CDatasetImpl::Open( *pAnswers, 0 ) )
		return false;

	for( i = 0; i < vecstrData.size( ); ++i ) {
		CDataPair	Datum;

		if( !( Datum.Open( vecstrData[ i ].c_str( ), pBayesNet->IsContinuous( i + 1 ) ) &&
			CDatasetImpl::Open( Datum, i + ( pAnswers ? 1 : 0 ) ) ) )
			return false; }

	return true; }

bool CDatasetImpl::Open( const CDataPair& Datum, size_t iExp ) {
	vector<size_t>	veciGenes;
	size_t			i, j, iOne, iTwo;
	float			d;

	m_veccQuants[ iExp ] = Datum.IsContinuous( ) ? -1 : Datum.GetValues( );
	veciGenes.resize( m_vecstrGenes.size( ) );
	for( i = 0; i < veciGenes.size( ); ++i )
		veciGenes[ i ] = Datum.GetGene( m_vecstrGenes[ i ] );

	if( Datum.IsContinuous( ) ) {
		CDistanceMatrix*	pDatum;

		pDatum = new CDistanceMatrix( );
		pDatum->Initialize( m_vecstrGenes.size( ) );
		for( i = 0; i < pDatum->GetSize( ); ++i )
			for( j = ( i + 1 ); j < pDatum->GetSize( ); ++j )
				pDatum->Set( i, j, CMeta::GetNaN( ) );
		for( i = 0; i < veciGenes.size( ); ++i ) {
			if( ( iOne = veciGenes[ i ] ) == -1 )
				continue;
			for( j = ( i + 1 ); j < veciGenes.size( ); ++j )
				if( ( ( iTwo = veciGenes[ j ] ) != -1 ) &&
					!CMeta::IsNaN( d = Datum.Get( iOne, iTwo ) ) )
					pDatum->Set( i, j, d ); }
		m_apData[ iExp ] = pDatum; }
	else {
		CCompactMatrix*	pDatum;

		pDatum = new CCompactMatrix( );
		pDatum->Initialize( m_vecstrGenes.size( ), m_veccQuants[ iExp ] + 1, true );
		for( i = 0; i < veciGenes.size( ); ++i )
			if( ( iOne = veciGenes[ i ] ) != -1 )
				for( j = ( i + 1 ); j < veciGenes.size( ); ++j )
					if( ( ( iTwo = veciGenes[ j ] ) != -1 ) &&
						!CMeta::IsNaN( d = Datum.Get( iOne, iTwo ) ) )
						pDatum->Set( i, j, (unsigned char)( Datum.Quantize( d ) + 1 ) );
		m_apData[ iExp ] = pDatum; }

	return true; }

/*!
 * Construct a dataset corresponding to the given answer file and data files.
 * 
 * \param szAnswerFile
 * Answer file which will become the first node of the dataset.
 * 
 * \param vecstrDataFiles
 * Vector of file paths to load.
 * 
 * \returns
 * True if dataset was constructed successfully.
 * 
 * Creates a dataset with nodes corresponding to the given data files; the given answer file is inserted
 * as the first (0th) node.  All files are assumed to be continuous.
 */
bool CDataset::Open( const char* szAnswerFile, const std::vector<std::string>& vecstrDataFiles ) {
	size_t	i;

	Reset( );
	m_veccQuants.resize( 1 + vecstrDataFiles.size( ) );
	{
		CDataPair	Answers;

		if( !Answers.Open( szAnswerFile, true ) )
			return false;

		m_vecstrGenes.resize( Answers.GetGenes( ) );
		for( i = 0; i < m_vecstrGenes.size( ); ++i )
			m_vecstrGenes[ i ] = Answers.GetGene( i );
		m_apData = new void*[ m_veccQuants.size( ) ];
		if( !CDatasetImpl::Open( Answers, 0 ) )
			return false;
	}

	m_veciMapping.resize( m_veccQuants.size( ) );
	m_veciMapping[ 0 ] = 0;
	for( i = 1; i <= vecstrDataFiles.size( ); ++i ) {
		CDataPair	Datum;

		if( !Datum.Open( vecstrDataFiles[ i - 1 ].c_str( ), true ) ||
			!CDatasetImpl::Open( Datum, i ) )
			return false;
		m_veciMapping[ i ] = i; }

	return true; }

/*!
 * \brief
 * Construct a dataset corresponding to the given files.
 * 
 * \param vecstrDataFiles
 * Vector of file paths to load.
 * 
 * \returns
 * True if dataset was constructed successfully.
 * 
 * Creates a dataset with nodes corresponding to the given data files.  All files are assumed to be
 * continuous.
 */
bool CDataset::Open( const std::vector<std::string>& vecstrDataFiles ) {
	size_t	i;

	if( !OpenGenes( vecstrDataFiles ) )
		return false;
	Reset( );
	m_apData = new void*[ m_veccQuants.size( ) ];

	for( i = 0; i < vecstrDataFiles.size( ); ++i ) {
		CDataPair	Datum;

		if( !( Datum.Open( vecstrDataFiles[ i ].c_str( ), true ) &&
			CDatasetImpl::Open( Datum, i ) ) )
			return false; }

	return true; }

float CDatasetImpl::GetContinuous( size_t iY, size_t iX, size_t iNode ) const {
	size_t	iMap;

	if( ( iMap = m_veciMapping[ iNode ] ) == -1 )
		return CMeta::GetNaN( );

	return ( ( m_veccQuants[ iMap ] == (unsigned char)-1 ) ?
		((CDistanceMatrix*)m_apData[ iMap ])->Get( iY, iX ) :
		((CCompactMatrix*)m_apData[ iMap ])->Get( iY, iX ) ); }

size_t CDataset::GetDiscrete( size_t iY, size_t iX, size_t iNode ) const {
	size_t	iMap;

	if( ( iMap = m_veciMapping[ iNode ] ) == -1 )
		return -1;

	return ( ( m_veccQuants[ iMap ] == (unsigned char)-1 ) ? -1 :
		( ((CCompactMatrix*)m_apData[ iMap ])->Get( iY, iX ) - 1 ) ); }

bool CDataset::IsExample( size_t iY, size_t iX ) const {
	size_t	i;

	for( i = 0; i < m_veccQuants.size( ); ++i )
		if( m_veccQuants[ i ] == (unsigned char)-1 ) {
			if( !CMeta::IsNaN( ((CDistanceMatrix*)m_apData[ i ])->Get( iY, iX ) ) )
				return true; }
		else if( ((CCompactMatrix*)m_apData[ i ])->Get( iY, iX ) )
			return true;

	return false; }

void CDataset::Remove( size_t iY, size_t iX ) {
	size_t	i;

	for( i = 0; i < m_veccQuants.size( ); ++i )
		if( m_veccQuants[ i ] == (unsigned char)-1 )
			((CDistanceMatrix*)m_apData[ i ])->Set( iY, iX, CMeta::GetNaN( ) );
		else
			((CCompactMatrix*)m_apData[ i ])->Set( iY, iX, 0 ); }

// CDataOverlayImpl

const vector<string>& CDataOverlayImpl::GetGeneNames( ) const {

	return m_pDataset->GetGeneNames( ); }

size_t CDataOverlayImpl::GetExperiments( ) const {

	return m_pDataset->GetExperiments( ); }

size_t CDataOverlayImpl::GetGene( const std::string& strGene ) const {

	return m_pDataset->GetGene( strGene ); }

size_t CDataOverlayImpl::GetBins( size_t iExp ) const {

	return m_pDataset->GetBins( iExp ); }

size_t CDataOverlayImpl::GetGenes( ) const {

	return m_pDataset->GetGenes( ); }

bool CDataOverlayImpl::IsHidden( size_t iNode ) const {

	return m_pDataset->IsHidden( iNode ); }

size_t CDataOverlayImpl::GetDiscrete( size_t iX, size_t iY, size_t iNode ) const {

	return m_pDataset->GetDiscrete( iX, iY, iNode ); }

float CDataOverlayImpl::GetContinuous( size_t iX, size_t iY, size_t iNode ) const {

	return m_pDataset->GetContinuous( iX, iY, iNode ); }

const string& CDataOverlayImpl::GetGene( size_t iGene ) const {

	return m_pDataset->GetGene( iGene ); }

void CDataOverlayImpl::Save( ostream& ostm, bool fBinary ) const {

	m_pDataset->Save( ostm, fBinary ); }

// CDataFilter

/*!
 * \brief
 * Associates the data filter with the given dataset, gene set, and filter type.
 * 
 * \param pDataset
 * Dataset to be associated with the overlaying mask.
 * 
 * \param Genes
 * Gene set used to filter the dataset.
 * 
 * \param eFilter
 * Way in which to use the given genes to remove gene pairs.
 * 
 * \param pAnswers
 * If non-null, answer set to be used for filter types requiring answers (e.g. CDat::EFilterTerm).
 * 
 * \remarks
 * No calculation occurs during Attach; the gene set and filter type are stored, and calculations are
 * performed dynamically in IsExample.  The gene set is not copied, so the given object should not be
 * destroyed until after the filter.  If a filter types requiring an answer file is given without an
 * accompanying answer file, the filter won't crash, but results might be a little odd.
 */
void CDataFilter::Attach( const IDataset* pDataset, const CGenes& Genes, CDat::EFilter eFilter,
	const CDat* pAnswers ) {
	size_t	i;

	m_pDataset = pDataset;
	m_pGenes = &Genes;
	m_eFilter = eFilter;
	m_pAnswers = pAnswers;

	m_vecfGenes.resize( GetGenes( ) );
	for( i = 0; i < m_vecfGenes.size( ); ++i )
		m_vecfGenes[ i ] = m_pGenes->IsGene( GetGene( i ) );
	if( m_pAnswers ) {
		m_veciAnswers.resize( GetGenes( ) );
		for( i = 0; i < m_veciAnswers.size( ); ++i )
			m_veciAnswers[ i ] = m_pAnswers->GetGene( GetGene( i ) ); } }

bool CDataFilter::IsExample( size_t iY, size_t iX ) const {

	if( !m_pDataset )
		return false;
	if( !( m_pGenes && m_pGenes->GetGenes( ) ) )
		return m_pDataset->IsExample( iY, iX );

	switch( m_eFilter ) {
		case CDat::EFilterInclude:
			if( !( m_vecfGenes[ iY ] && m_vecfGenes[ iX ] ) )
				return false;
			break;

		case CDat::EFilterExclude:
			if( m_vecfGenes[ iY ] || m_vecfGenes[ iX ] )
				return false;
			break;

		case CDat::EFilterEdge:
			if( !( m_vecfGenes[ iY ] || m_vecfGenes[ iX ] ) )
				return false;
			break;

		case CDat::EFilterTerm:
			if( ( m_pAnswers && ( ( m_veciAnswers[ iY ] == -1 ) || ( m_veciAnswers[ iX ] == -1 ) ) ) ||
				( !( m_vecfGenes[ iY ] && m_vecfGenes[ iX ] ) &&
				( !( m_vecfGenes[ iY ] || m_vecfGenes[ iX ] ) || ( m_pAnswers &&
				( m_pAnswers->Get( m_veciAnswers[ iY ], m_veciAnswers[ iX ] ) > 0 ) ) ) ) )
				return false;
			break; }

	return m_pDataset->IsExample( iY, iX ); }

// CDataMask

/*!
 * \brief
 * Associates the data mask with the given dataset and randomly hides a fraction of its data.
 * 
 * \param pDataset
 * Dataset to be associated with the overlaying mask.
 * 
 * \param dFraction
 * Fraction of gene pairs (between 0 and 1) to be randomly masked.
 */
void CDataMask::AttachRandom( const IDataset* pDataset, float dFraction ) {
	size_t	i, j;

	Attach( pDataset );
	for( i = 0; i < m_Mask.GetSize( ); ++i )
		for( j = ( i + 1 ); j < m_Mask.GetSize( ); ++j )
			m_Mask.Set( i, j, ( m_Mask.Get( i, j ) && ( ( (float)rand( ) / RAND_MAX ) < dFraction ) ) ); }

/*!
 * \brief
 * Associates the data mask with the given mask's underlying dataset and reverses its mask.
 * 
 * \param DataMask
 * Mask to be reversed by the current mask.
 * 
 * This associates the current mask with the given mask's underlying dataset and generates an inverted
 * mask: all pairs hidden in the given mask are unhidden, and all unhidden pairs are hidden.  Data pairs
 * missing in the underlying dataset will return false from IsExample regardless.
 */
void CDataMask::AttachComplement( const CDataMask& DataMask ) {
	size_t	i, j;

	Attach( DataMask.m_pDataset );
	for( i = 0; i < m_Mask.GetSize( ); ++i )
		for( j = ( i + 1 ); j < m_Mask.GetSize( ); ++j )
			m_Mask.Set( i, j, ( m_Mask.Get( i, j ) && !DataMask.m_Mask.Get( i, j ) ) ); }

/*!
 * \brief
 * Associates the data mask with the given dataset.
 * 
 * \param pDataset
 * Dataset to be associated with the overlaying mask.
 */
void CDataMask::Attach( const IDataset* pDataset ) {
	size_t	i, j;

	m_pDataset = pDataset;
	m_Mask.Initialize( m_pDataset->GetGenes( ) );
	for( i = 0; i < m_Mask.GetSize( ); ++i )
		for( j = ( i + 1 ); j < m_Mask.GetSize( ); ++j )
			m_Mask.Set( i, j, m_pDataset->IsExample( i, j ) ); }

// CDataSubset

/*!
 * \brief
 * Construct a data subset corresponding to the given Bayes net using data files from the given directory;
 * loads only the requested number of genes at a time.
 * 
 * \param szDataDirectory
 * Directory from which data files are loaded.
 * 
 * \param pBayesNet
 * Bayes nets whose nodes will correspond to files in the data subset.
 * 
 * \param iGeneSize
 * Number of genes to be loaded at once; all data for these genes will be loaded.
 * 
 * \returns
 * True if data subset was constructed successfully.
 * 
 * Creates a data subset with nodes corresponding to the given Bayes net structure.  Data is loaded
 * continuously or discretely as indicated by the Bayes net, and nodes for which a corresponding data
 * file (i.e. one with the same name followed by an appropriate CDat extension) cannot be located are
 * marked as hidden.  No data is loaded during initialization; this instead occurs during Open.
 * 
 * \see
 * CDataset::Open
 */
bool CDataSubset::Initialize( const char* szDataDirectory, const IBayesNet* pBayesNet, size_t iGeneSize ) {
	size_t	i;

	m_iSize = iGeneSize;
	m_vecstrData.clear( );
	m_fContinuous = pBayesNet->IsContinuous( );
	{
		set<string>					setstrGenes;
		set<string>::const_iterator	iterGenes;
		vector<string>				vecstrNodes;

		pBayesNet->GetNodes( vecstrNodes );
		OpenMax( szDataDirectory, vecstrNodes, false, m_vecstrData, &setstrGenes );
		m_vecstrGenes.resize( setstrGenes.size( ) );
		i = 0;
		for( iterGenes = setstrGenes.begin( ); iterGenes != setstrGenes.end( ); ++iterGenes )
			m_vecstrGenes[ i++ ] = *iterGenes;
	}
	m_Examples.Initialize( m_iSize, m_vecstrGenes.size( ) );

	return true; }

/*!
 * \brief
 * Construct a data subset corresponding using the given data files; loads only the requested number of
 * genes at a time.
 * 
 * \param vecstrDataFiles
 * Vector of file paths to load.
 * 
 * \param iGeneSize
 * Number of genes to be loaded at once; all data for these genes will be loaded.
 * 
 * \returns
 * True if data subset was constructed successfully.
 * 
 * Creates a dataset with nodes corresponding to the given data files; all files are assumed to be
 * continuous.  No data is loaded during initialization; this instead occurs during Open.
 * 
 * \see
 * CDataset::Open
 */
bool CDataSubset::Initialize( const std::vector<std::string>& vecstrDataFiles, size_t iGeneSize ) {
	size_t	i;

	m_iSize = iGeneSize;
	m_vecstrData.resize( vecstrDataFiles.size( ) );
	m_veccQuants.resize( vecstrDataFiles.size( ) );
	for( i = 0; i < vecstrDataFiles.size( ); ++i )
		m_vecstrData[ i ] = vecstrDataFiles[ i ];
	m_fContinuous = true;

	if( !OpenGenes( vecstrDataFiles ) )
		return false;
	m_Examples.Initialize( m_iSize, m_vecstrGenes.size( ) );

	return true; }

/*!
 * \brief
 * Load data for the subset's configured number of genes starting at the given offset.
 * 
 * \param iGeneOffset
 * First gene ID to load into the data subset.
 * 
 * \returns
 * True if data was loaded successfully.
 * 
 * For a data subset configured with a particular size, loads all data for that many genes beginning at the
 * given index.
 * 
 * \see
 * Initialize
 */
bool CDataSubset::Open( size_t iGeneOffset ) {
	size_t	i, j;

	m_iOffset = iGeneOffset;
	for( i = 0; i < m_Examples.GetRows( ); ++i )
		for( j = 0; j < m_Examples.GetColumns( ); ++j )
			m_Examples.Get( i, j ).Reset( );

	m_iSize = ( ( m_iOffset + m_Examples.GetRows( ) ) > m_vecstrGenes.size( ) ) ?
		( m_vecstrGenes.size( ) - m_iOffset ) : m_Examples.GetRows( );
	for( i = 0; i < m_vecstrData.size( ); ++i ) {
		CDataPair	Datum;

		if( !( Datum.Open( m_vecstrData[ i ].c_str( ), m_fContinuous ) &&
			CDataSubsetImpl::Open( Datum, i ) ) )
			return false; }

	return true; }

bool CDataSubsetImpl::Open( const CDataPair& Datum, size_t iExp ) {
	vector<size_t>	veciGenes;
	size_t			i, j, iOne, iTwo;
	float			d;

	m_veccQuants[ iExp ] = Datum.IsContinuous( ) ? -1 : Datum.GetValues( );
	veciGenes.resize( m_vecstrGenes.size( ) );
	for( i = 0; i < veciGenes.size( ); ++i )
		veciGenes[ i ] = Datum.GetGene( m_vecstrGenes[ i ] );

	for( i = 0; i < m_iSize; ++i ) {
		if( ( iOne = veciGenes[ i + m_iOffset ] ) == -1 )
			continue;
		for( j = 0; j < veciGenes.size( ); ++j )
			if( ( ( iTwo = veciGenes[ j ] ) != -1 ) &&
				!CMeta::IsNaN( d = Datum.Get( iOne, iTwo ) ) )
				m_Examples.Get( i, j ).Set( iExp, d, Datum, m_vecstrData.size( ) ); }

	return true; }

}

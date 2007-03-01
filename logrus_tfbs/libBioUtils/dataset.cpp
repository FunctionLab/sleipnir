#include "stdafx.h"
#include "dataset.h"
#include "bayesnet.h"
#include "genome.h"

namespace libBioUtils {

const char	CDataImpl::c_szDat[]	= ".dat";
const char	CDataImpl::c_szDab[]	= ".dab";
const char	CDataImpl::c_szPcl[]	= ".pcl";

void CDataImpl::FilterGenes( IDataset* pData, const CGenes& Genes, CDat::EFilter eFilt ) {
	vector<bool>	vecfGenes;
	size_t			i, j;

	if( !Genes.GetGenes( ) )
		return;

	vecfGenes.resize( pData->GetGenes( ) );
	for( i = 0; i < vecfGenes.size( ); ++i )
		vecfGenes[ i ] = Genes.IsGene( pData->GetGene( i ) );

	for( i = 0; i < vecfGenes.size( ); ++i ) {
		if( vecfGenes[ i ] ) {
			if( eFilt == CDat::EFilterInclude )
				continue;
			if( eFilt == CDat::EFilterExclude ) {
				for( j = ( i + 1 ); j < vecfGenes.size( ); ++j )
					pData->Remove( i, j );
				continue; } }
		for( j = ( i + 1 ); j < vecfGenes.size( ); ++j )
			switch( eFilt ) {
				case CDat::EFilterInclude:
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

size_t CDataImpl::OpenMax( const char* szDataDir, const vector<string>& vecstrNodes,
	bool fAnswers, vector<string>& vecstrData, set<string>* psetstrGenes ) {
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
				strFile += c_szPcl;
				ifsm.clear( );
				ifsm.open( strFile.c_str( ) );
				if( ifsm.is_open( ) ) {
					iRet++;
					m_veciMapping[ i ] = iMap++;
					vecstrData.push_back( strFile );
					if( psetstrGenes )
						OpenGenes( ifsm, false, true, *psetstrGenes ); }
				else {
					g_CatBioUtils.info( "CDataImpl::OpenMax( %s ) assuming %s is hidden",
						szDataDir, vecstrNodes[ i ].c_str( ) );
					continue; } } }
		ifsm.close( ); }

	return iRet; }

bool CDataImpl::OpenGenes( istream& istm, bool fBinary, bool fPCL, set<string>& setstrGenes ) const {
	CDat	Dat;
	size_t	i;

	if( !Dat.OpenGenes( istm, fBinary, fPCL ) )
		return false;
	for( i = 0; i < Dat.GetGenes( ); ++i )
		setstrGenes.insert( Dat.GetGene( i ) );
	return true; }

bool CDataImpl::OpenGenes( const vector<string>& vecstrData ) {
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
				if( !( ifsm.is_open( ) && OpenGenes( ifsm, false, true, setstrGenes ) ) )
					return false; } }
		ifsm.close( ); }

	m_vecstrGenes.resize( setstrGenes.size( ) );
	i = 0;
	for( iterGenes = setstrGenes.begin( ); iterGenes != setstrGenes.end( ); ++iterGenes )
		m_vecstrGenes[ i++ ] = *iterGenes;

	return true; }

bool CDataImpl::IsHidden( size_t iNode ) const {

	return ( m_veciMapping[ iNode ] == -1 ); }

const string& CDataImpl::GetGene( size_t iGene ) const {

	return m_vecstrGenes[ iGene ]; }

size_t CDataImpl::GetGenes( ) const {

	return m_vecstrGenes.size( ); }

const vector<string>& CDataImpl::GetGeneNames( ) const {

	return m_vecstrGenes; }

size_t CDataImpl::GetExperiments( ) const {

	return m_veciMapping.size( ); }

size_t CDataImpl::GetGene( const string& strGene ) const {
	size_t	i;

	for( i = 0; i < m_vecstrGenes.size( ); ++i )
		if( m_vecstrGenes[ i ] == strGene )
			return i;

	return -1; }

unsigned char CDataImpl::GetBins( size_t iExp ) const {

	return m_veccQuants[ iExp ]; }

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

bool CDataImpl::OpenBinary( istream& istm ) {
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

void CDataImpl::SaveBinary( ostream& ostm ) const {
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

bool CDataset::Open( const char* szAnswers, const char* szDataDir,
	const IBayesNet* pBayesNet ) {
	CDataPair	Answers;

	return ( Answers.Open( szAnswers, pBayesNet->IsContinuous( 0 ) ) &&
		Open( Answers, szDataDir, pBayesNet ) ); }

bool CDataset::Open( const CDataPair& Answers, const char* szDataDir,
	const IBayesNet* pBayesNet ) {

	return CDatasetImpl::Open( &Answers, szDataDir, pBayesNet ); }

bool CDataset::Open( const char* szDataDir, const IBayesNet* pBayesNet ) {

	return CDatasetImpl::Open( NULL, szDataDir, pBayesNet ); }

bool CDatasetImpl::Open( const CDataPair* pAnswers, const char* szDataDir,
	const IBayesNet* pBayesNet ) {
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
						pDatum->Set( i, j, (unsigned char)( Datum.Quantify( d ) + 1 ) );
		m_apData[ iExp ] = pDatum; }

	return true; }

bool CDataset::Open( const char* szAnswers, const vector<string>& vecstrData ) {
	size_t	i;

	Reset( );
	m_veccQuants.resize( 1 + vecstrData.size( ) );
	{
		CDataPair	Answers;

		if( !Answers.Open( szAnswers, true ) )
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
	for( i = 1; i <= vecstrData.size( ); ++i ) {
		CDataPair	Datum;

		if( !Datum.Open( vecstrData[ i - 1 ].c_str( ), true ) ||
			!CDatasetImpl::Open( Datum, i ) )
			return false;
		m_veciMapping[ i ] = i; }

	return true; }

bool CDataset::Open( const vector<string>& vecstrData ) {
	size_t	i;

	if( !OpenGenes( vecstrData ) )
		return false;
	Reset( );
	m_apData = new void*[ m_veccQuants.size( ) ];

	for( i = 0; i < vecstrData.size( ); ++i ) {
		CDataPair	Datum;

		if( !( Datum.Open( vecstrData[ i ].c_str( ), true ) &&
			CDatasetImpl::Open( Datum, i ) ) )
			return false; }

	return true; }

bool CDataset::IsHidden( size_t iNode ) const {

	return CDataImpl::IsHidden( iNode ); }

float CDataset::GetContinuous( size_t iX, size_t iY, size_t iNode ) const {
	size_t	iMap;

	if( ( iMap = m_veciMapping[ iNode ] ) == -1 )
		return CMeta::GetNaN( );

	return ( ( m_veccQuants[ iMap ] == (unsigned char)-1 ) ?
		((CDistanceMatrix*)m_apData[ iMap ])->Get( iX, iY ) :
		((CCompactMatrix*)m_apData[ iMap ])->Get( iX, iY ) ); }

size_t CDataset::GetDiscrete( size_t iX, size_t iY, size_t iNode ) const {
	size_t	iMap;

	if( ( iMap = m_veciMapping[ iNode ] ) == -1 )
		return -1;

	return ( ( m_veccQuants[ iMap ] == (unsigned char)-1 ) ? -1 :
		( ((CCompactMatrix*)m_apData[ iMap ])->Get( iX, iY ) - 1 ) ); }

const string& CDataset::GetGene( size_t iGene ) const {

	return CDataImpl::GetGene( iGene ); }

size_t CDataset::GetGenes( ) const {

	return CDataImpl::GetGenes( ); }

bool CDataset::IsExample( size_t iX, size_t iY ) const {
	size_t	i;

	for( i = 0; i < m_veccQuants.size( ); ++i )
		if( m_veccQuants[ i ] == (unsigned char)-1 ) {
			if( !CMeta::IsNaN( ((CDistanceMatrix*)m_apData[ i ])->Get( iX, iY ) ) )
				return true; }
		else if( ((CCompactMatrix*)m_apData[ i ])->Get( iX, iY ) )
			return true;

	return false; }

const vector<string>& CDataset::GetGeneNames( ) const {

	return CDataImpl::GetGeneNames( ); }

size_t CDataset::GetExperiments( ) const {

	return CDataImpl::GetExperiments( ); }

size_t CDataset::GetGene( const string& strGene ) const {

	return CDataImpl::GetGene( strGene ); }

size_t CDataset::GetBins( size_t iExp ) const {

	return CDataImpl::GetBins( iExp ); }

void CDataset::Remove( size_t iX, size_t iY ) {
	size_t	i;

	for( i = 0; i < m_veccQuants.size( ); ++i )
		if( m_veccQuants[ i ] == (unsigned char)-1 )
			((CDistanceMatrix*)m_apData[ i ])->Set( iX, iY, CMeta::GetNaN( ) );
		else
			((CCompactMatrix*)m_apData[ i ])->Set( iX, iY, 0 ); }

void CDataset::FilterGenes( const CGenes& Genes, CDat::EFilter eFilt ) {

	CDataImpl::FilterGenes( this, Genes, eFilt ); }

// CDataOverlayImpl

const vector<string>& CDataOverlayImpl::GetGeneNames( ) const {

	return m_pDataset->GetGeneNames( ); }

size_t CDataOverlayImpl::GetExperiments( ) const {

	return m_pDataset->GetExperiments( ); }

size_t CDataOverlayImpl::GetGene( const string& strGene ) const {

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

// CDataFilter

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
	m_veciAnswers.resize( GetGenes( ) );
	for( i = 0; i < m_veciAnswers.size( ); ++i )
		m_veciAnswers[ i ] = m_pAnswers->GetGene( GetGene( i ) ); }

void CDataFilter::Remove( size_t iX, size_t iY ) {

	throw 0; }

bool CDataFilter::IsExample( size_t iX, size_t iY ) const {

	if( !( m_pGenes && m_pGenes->GetGenes( ) ) )
		return true;

	switch( m_eFilter ) {
		case CDat::EFilterInclude:
			if( !( m_vecfGenes[ iX ] && m_vecfGenes[ iY ] ) )
				return false;

		case CDat::EFilterExclude:
			if( m_vecfGenes[ iX ] || m_vecfGenes[ iY ] )
				return false;

		case CDat::EFilterTerm:
			if( ( ( m_veciAnswers[ iX ] == -1 ) || ( m_veciAnswers[ iY ] == -1 ) ) ||
				( !( m_vecfGenes[ iX ] && m_vecfGenes[ iY ] ) &&
				( !( m_vecfGenes[ iX ] || m_vecfGenes[ iY ] ) || ( m_pAnswers->Get( m_veciAnswers[ iX ],
				m_veciAnswers[ iY ] ) > 0 ) ) ) )
				return false; }

	return m_pDataset->IsExample( iX, iY ); }

// CDataMask

void CDataMask::AttachRandom( const IDataset* pDataset, float dFrac ) {
	size_t	i, j;

	Attach( pDataset );
	for( i = 0; i < m_Mask.GetSize( ); ++i )
		for( j = ( i + 1 ); j < m_Mask.GetSize( ); ++j )
			m_Mask.Set( i, j, ( m_Mask.Get( i, j ) && ( ( (float)rand( ) / RAND_MAX ) < dFrac ) ) ); }

void CDataMask::AttachComplement( const CDataMask& Data ) {
	size_t	i, j;

	Attach( Data.m_pDataset );
	for( i = 0; i < m_Mask.GetSize( ); ++i )
		for( j = ( i + 1 ); j < m_Mask.GetSize( ); ++j )
			m_Mask.Set( i, j, ( m_Mask.Get( i, j ) && !Data.m_Mask.Get( i, j ) ) ); }

void CDataMask::Attach( const IDataset* pDataset ) {
	size_t	i, j;

	m_pDataset = pDataset;
	m_Mask.Initialize( m_pDataset->GetGenes( ) );
	for( i = 0; i < m_Mask.GetSize( ); ++i )
		for( j = ( i + 1 ); j < m_Mask.GetSize( ); ++j )
			m_Mask.Set( i, j, m_pDataset->IsExample( i, j ) ); }

bool CDataMask::IsExample( size_t iX, size_t iY ) const {

	return m_Mask.Get( iX, iY ); }

void CDataMask::Remove( size_t iX, size_t iY ) {

	m_Mask.Set( iX, iY, false ); }

// CDataSubset

bool CDataSubset::Initialize( const char* szDataDir, const IBayesNet* pBayesNet,
	size_t iSize ) {
	size_t	i;

	m_iSize = iSize;
	m_vecstrData.clear( );
	m_fContinuous = pBayesNet->IsContinuous( );
	{
		set<string>					setstrGenes;
		set<string>::const_iterator	iterGenes;
		vector<string>				vecstrNodes;

		pBayesNet->GetNodes( vecstrNodes );
		OpenMax( szDataDir, vecstrNodes, false, m_vecstrData, &setstrGenes );
		m_vecstrGenes.resize( setstrGenes.size( ) );
		i = 0;
		for( iterGenes = setstrGenes.begin( ); iterGenes != setstrGenes.end( ); ++iterGenes )
			m_vecstrGenes[ i++ ] = *iterGenes;
	}
	m_Examples.Initialize( m_iSize, m_vecstrGenes.size( ) );

	return true; }

bool CDataSubset::Initialize( const vector<string>& vecstrData, size_t iSize ) {
	size_t	i;

	m_iSize = iSize;
	m_vecstrData.resize( vecstrData.size( ) );
	m_veccQuants.resize( vecstrData.size( ) );
	for( i = 0; i < vecstrData.size( ); ++i )
		m_vecstrData[ i ] = vecstrData[ i ];
	m_fContinuous = true;

	if( !OpenGenes( vecstrData ) )
		return false;
	m_Examples.Initialize( m_iSize, m_vecstrGenes.size( ) );

	return true; }

bool CDataSubset::Open( size_t iOffset ) {
	size_t	i, j;

	m_iOffset = iOffset;
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

bool CDataSubset::IsHidden( size_t iNode ) const {

	return CDataImpl::IsHidden( iNode ); }

size_t CDataSubset::GetDiscrete( size_t iX, size_t iY, size_t iNode ) const {
	size_t	iMap;

	if( ( iMap = m_veciMapping[ iNode ] ) == -1 )
		return -1;

	return m_Examples.Get( iX - m_iOffset, iY ).GetDiscrete( iMap ); }

float CDataSubset::GetContinuous( size_t iX, size_t iY, size_t iNode ) const {
	size_t	iMap;

	if( ( iMap = m_veciMapping[ iNode ] ) == -1 )
		return CMeta::GetNaN( );

	return m_Examples.Get( iX - m_iOffset, iY ).GetContinuous( iMap ); }

const string& CDataSubset::GetGene( size_t iGene ) const {

	return CDataImpl::GetGene( iGene ); }

size_t CDataSubset::GetGenes( ) const {

	return CDataImpl::GetGenes( ); }

bool CDataSubset::IsExample( size_t iX, size_t iY ) const {

	if( ( iX < m_iOffset ) || ( ( iX - m_iOffset ) >= m_iSize ) )
		return false;

	return m_Examples.Get( iX - m_iOffset, iY ).IsSet( ); }

const vector<string>& CDataSubset::GetGeneNames( ) const {

	return CDataImpl::GetGeneNames( ); }

size_t CDataSubset::GetExperiments( ) const {

	return CDataImpl::GetExperiments( ); }

size_t CDataSubset::GetGene( const string& strGene ) const {

	return CDataImpl::GetGene( strGene ); }

size_t CDataSubset::GetBins( size_t iExp ) const {

	return CDataImpl::GetBins( iExp ); }

void CDataSubset::Remove( size_t iX, size_t iY ) {

	m_Examples.Get( iX - m_iOffset, iY ).Reset( ); }

void CDataSubset::FilterGenes( const CGenes& Genes, CDat::EFilter eFilt ) {

	CDataImpl::FilterGenes( this, Genes, eFilt ); }

}

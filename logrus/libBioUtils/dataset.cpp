#include "stdafx.h"
#include "dataset.h"
#include "bayesnet.h"
#include "datapair.h"
#include "examplei.h"
#include "meta.h"
#include "genome.h"

namespace libBioUtils {

const char	CDataImpl::c_szDat[]	= ".dat";
const char	CDataImpl::c_szDab[]	= ".dab";

size_t CDataImpl::OpenMax( const char* szDataDir, const vector<string>& vecstrNodes,
	bool fAnswers, vector<string>& vecstrData, set<string>* psetstrGenes ) {
	size_t		i, iLength, iMap, iRet;
	string		strFile;
	ifstream	ifsm;

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
				OpenGenes( ifsm, true, *psetstrGenes ); }
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
					OpenGenes( ifsm, false, *psetstrGenes ); }
			else {
				g_CatBioUtils.info( "CDataImpl::OpenMax( %s ) assuming %s is hidden",
					szDataDir, vecstrNodes[ i ].c_str( ) );
				continue; } }
		ifsm.close( ); }

	return iRet; }

bool CDataImpl::OpenGenes( istream& istm, bool fBinary, set<string>& setstrGenes ) const {
	CDat	Dat;
	size_t	i;

	if( !Dat.OpenGenes( istm, fBinary ) )
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
		if( !( ifsm.is_open( ) && OpenGenes( ifsm, true, setstrGenes ) ) ) {
			ifsm.close( );
			ifsm.clear( );
			ifsm.open( vecstrData[ i ].c_str( ) );
			if( !( ifsm.is_open( ) && OpenGenes( ifsm, false, setstrGenes ) ) )
				return false; }
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
		m_veciMapping[ i ] = *(uint32_t*)pbData;
		pbData += sizeof(uint32_t); }

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
		m_veciMapping[ i ] = iVal; }

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
		iVal = (uint32_t)m_veciMapping[ i ];
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
	copy( m_veccQuants.begin( ), m_veccQuants.end( ), ac );
	ostm.write( (char*)&iVal, sizeof(iVal) );
	ostm.write( ac, iVal );
	delete[] ac; }

bool CDataset::Open( const char* szAnswers, const char* szDataDir,
	const IBayesNet* pBayesNet ) {
	CDataPair	Answers;

	return ( Answers.Open( szAnswers, pBayesNet->IsContinuous( ) ) &&
		Open( Answers, szDataDir, pBayesNet ) ); }

bool CDataset::Open( const CDataPair& Answers, const char* szDataDir,
	const IBayesNet* pBayesNet ) {
	size_t			i;
	vector<string>	vecstrData;

	m_fContinuous = pBayesNet->IsContinuous( );
	m_iMax = 1 + OpenMax( szDataDir, pBayesNet->GetNodes( ), true, vecstrData );
	m_veccQuants.resize( m_iMax );
	m_vecstrGenes.resize( Answers.GetGenes( ) );
	for( i = 0; i < m_vecstrGenes.size( ); ++i )
		m_vecstrGenes[ i ] = Answers.GetGene( i );
	m_Examples.Initialize( m_vecstrGenes.size( ) );
	if( !CDatasetImpl::Open( Answers, 0, m_iMax ) )
		return false;

	for( i = 0; i < vecstrData.size( ); ++i ) {
		CDataPair	Datum;

		if( !( Datum.Open( vecstrData[ i ].c_str( ), m_fContinuous ) &&
			CDatasetImpl::Open( Datum, i + 1, 0 ) ) )
			return false; }

	TrimExamples( );
	return true; }

bool CDatasetImpl::Open( const CDataPair& Datum, size_t iExp, size_t iMax ) {
	vector<size_t>	veciGenes;
	size_t			i, j, iOne, iTwo;
	float			d;

	m_veccQuants[ iExp ] = Datum.IsContinuous( ) ? -1 : Datum.GetValues( );
	veciGenes.resize( m_vecstrGenes.size( ) );
	for( i = 0; i < veciGenes.size( ); ++i )
		veciGenes[ i ] = Datum.GetGene( m_vecstrGenes[ i ] );

	for( i = 0; i < veciGenes.size( ); ++i ) {
		if( ( iOne = veciGenes[ i ] ) == -1 )
			continue;
		for( j = ( i + 1 ); j < veciGenes.size( ); ++j )
			if( ( ( iTwo = veciGenes[ j ] ) != -1 ) &&
				!CMeta::IsNaN( d = Datum.Get( iOne, iTwo ) ) )
				m_Examples.Get( i, j ).Set( iExp, d, Datum, iMax ); }

	return true; }

bool CDataset::Open( const char* szAnswers, const vector<string>& vecstrData ) {
	size_t	i;

	m_iMax = 1 + vecstrData.size( );
	m_veccQuants.resize( m_iMax );
	{
		CDataPair	Answers;

		if( !Answers.Open( szAnswers, true ) )
			return false;

		m_vecstrGenes.resize( Answers.GetGenes( ) );
		for( i = 0; i < m_vecstrGenes.size( ); ++i )
			m_vecstrGenes[ i ] = Answers.GetGene( i );
		m_Examples.Initialize( m_vecstrGenes.size( ) );
		if( !CDatasetImpl::Open( Answers, 0, m_iMax ) )
			return false;
	}

	m_veciMapping.resize( m_iMax );
	m_veciMapping[ 0 ] = 0;
	for( i = 1; i <= vecstrData.size( ); ++i ) {
		CDataPair	Datum;

		if( !Datum.Open( vecstrData[ i - 1 ].c_str( ), true ) ||
			!CDatasetImpl::Open( Datum, i, 0 ) )
			return false;
		m_veciMapping[ i ] = i; }

	TrimExamples( );
	return true; }

bool CDataset::Open( const vector<string>& vecstrData ) {
	size_t	i;

	if( !OpenGenes( vecstrData ) )
		return false;
	m_Examples.Initialize( m_vecstrGenes.size( ) );

	for( i = 0; i < vecstrData.size( ); ++i ) {
		CDataPair	Datum;

		if( !( Datum.Open( vecstrData[ i ].c_str( ), true ) &&
			CDatasetImpl::Open( Datum, i, vecstrData.size( ) ) ) )
			return false; }

	return true; }

void CDatasetImpl::TrimExamples( ) {
	size_t	i, j;

	for( i = 0; i < m_Examples.GetSize( ); ++i )
		for( j = ( i + 1 ); j < m_Examples.GetSize( ); ++j ) {
			CExampleImpl&	Example	= m_Examples.Get( i, j );

			if( !Example.IsEvidence( m_fContinuous, m_iMax ) )
				Example.Reset( ); } }

bool CDataset::IsHidden( size_t iNode ) const {

	return CDataImpl::IsHidden( iNode ); }

float CDataset::GetContinuous( size_t iX, size_t iY, size_t iNode ) const {
	size_t	iMap;

	if( ( iMap = m_veciMapping[ iNode ] ) == -1 )
		return CMeta::GetNaN( );

	return m_Examples.Get( iX, iY ).GetContinuous( iMap ); }

size_t CDataset::GetDiscrete( size_t iX, size_t iY, size_t iNode ) const {
	size_t	iMap;

	if( ( iMap = m_veciMapping[ iNode ] ) == -1 )
		return -1;

	return m_Examples.Get( iX, iY ).GetDiscrete( iMap ); }

const string& CDataset::GetGene( size_t iGene ) const {

	return CDataImpl::GetGene( iGene ); }

size_t CDataset::GetGenes( ) const {

	return CDataImpl::GetGenes( ); }

bool CDataset::IsExample( size_t iX, size_t iY ) const {

	return m_Examples.Get( iX, iY ).IsSet( ); }

const vector<string>& CDataset::GetGeneNames( ) const {

	return CDataImpl::GetGeneNames( ); }

size_t CDataset::GetExperiments( ) const {

	return CDataImpl::GetExperiments( ); }

size_t CDataset::GetGene( const string& strGene ) const {

	return CDataImpl::GetGene( strGene ); }

size_t CDataset::GetBins( size_t iExp ) const {

	return CDataImpl::GetBins( iExp ); }

void CDataset::Remove( size_t iX, size_t iY ) {

	m_Examples.Get( iX, iY ).Reset( ); }

void CDataMask::AttachRandom( const IDataset* pDataset, float dFrac ) {
	size_t	i, j;

	m_pDataset = pDataset;
	m_Mask.Initialize( m_pDataset->GetGenes( ) );
	for( i = 0; i < m_Mask.GetSize( ); ++i )
		for( j = ( i + 1 ); j < m_Mask.GetSize( ); ++j )
			m_Mask.Set( i, j, ( m_pDataset->IsExample( i, j ) &&
				( ( (float)rand( ) / RAND_MAX ) < dFrac ) ) ); }

void CDataMask::AttachComplement( const CDataMask& Data ) {
	size_t	i, j;

	m_pDataset = Data.m_pDataset;
	m_Mask.Initialize( m_pDataset->GetGenes( ) );
	for( i = 0; i < m_Mask.GetSize( ); ++i )
		for( j = ( i + 1 ); j < m_Mask.GetSize( ); ++j )
			m_Mask.Set( i, j, ( m_pDataset->IsExample( i, j ) &&
				!Data.m_Mask.Get( i, j ) ) ); }

void CDataMask::AttachGenes( const IDataset* pDataset, istream& istm ) {
	static const size_t	c_iBuffer	= 1024;
	char		szBuf[ c_iBuffer ];
	set<size_t>	setiGenes;
	size_t		i, j;
	bool		fOne;

	m_pDataset = pDataset;
	m_Mask.Initialize( m_pDataset->GetGenes( ) );
	while( istm.peek( ) != EOF ) {
		istm.getline( szBuf, c_iBuffer - 1 );
		if( ( i = GetGene( szBuf ) ) != -1 )
			setiGenes.insert( i ); }

	for( i = 0; i < m_Mask.GetSize( ); ++i ) {
		fOne = ( setiGenes.find( i ) != setiGenes.end( ) );
		for( j = ( i + 1 ); j < m_Mask.GetSize( ); ++j )
			m_Mask.Set( i, j, ( m_pDataset->IsExample( i, j ) && ( fOne ||
				( setiGenes.find( j ) != setiGenes.end( ) ) ) ) ); } }

size_t CDataMask::GetGenes( ) const {

	return m_pDataset->GetGenes( ); }

bool CDataMask::IsHidden( size_t iNode ) const {

	return m_pDataset->IsHidden( iNode ); }

size_t CDataMask::GetDiscrete( size_t iX, size_t iY, size_t iNode ) const {

	return m_pDataset->GetDiscrete( iX, iY, iNode ); }

float CDataMask::GetContinuous( size_t iX, size_t iY, size_t iNode ) const {

	return m_pDataset->GetContinuous( iX, iY, iNode ); }

const string& CDataMask::GetGene( size_t iGene ) const {

	return m_pDataset->GetGene( iGene ); }

bool CDataMask::IsExample( size_t iX, size_t iY ) const {

	return m_Mask.Get( iX, iY ); }

const vector<string>& CDataMask::GetGeneNames( ) const {

	return m_pDataset->GetGeneNames( ); }

size_t CDataMask::GetExperiments( ) const {

	return m_pDataset->GetExperiments( ); }

size_t CDataMask::GetGene( const string& strGene ) const {

	return m_pDataset->GetGene( strGene ); }

size_t CDataMask::GetBins( size_t iExp ) const {

	return m_pDataset->GetBins( iExp ); }

void CDataMask::Remove( size_t iX, size_t iY ) {

	m_Mask.Set( iX, iY, false ); }

bool CDataSubset::Initialize( const char* szDataDir, const IBayesNet* pBayesNet,
	size_t iSize ) {
	size_t	i;

	m_iSize = iSize;
	m_vecstrData.clear( );
	m_fContinuous = pBayesNet->IsContinuous( );
	{
		set<string>					setstrGenes;
		set<string>::const_iterator	iterGenes;

		OpenMax( szDataDir, pBayesNet->GetNodes( ), false, m_vecstrData, &setstrGenes );
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

CDatasetCompactImpl::CDatasetCompactImpl( ) : m_iData(0), m_aData(NULL) {

	m_fContinuous = false; }

CDatasetCompactImpl::~CDatasetCompactImpl( ) {

	if( m_aData )
		delete[] m_aData; }

bool CDatasetCompact::Open( const vector<string>& vecstrData ) {
	size_t	i;

	if( !OpenGenes( vecstrData ) )
		return false;
	if( m_aData )
		delete[] m_aData;
	m_aData = new CCompactMatrix[ m_iData = (uint32_t)vecstrData.size( ) ];

	for( i = 0; i < vecstrData.size( ); ++i ) {
		CDataPair	Datum;

		if( !( Datum.Open( vecstrData[ i ].c_str( ), false ) &&
			CDatasetCompactImpl::Open( Datum, i ) ) )
			return false; }

	return true; }

bool CDatasetCompact::Open( const CDataPair& Answers, const char* szDataDir,
	const IBayesNet* pBayesNet ) {
	size_t			i, j, k;
	vector<string>	vecstrData;

	if( pBayesNet->IsContinuous( ) )
		return false;

	m_iData = 1 + (uint32_t)OpenMax( szDataDir, pBayesNet->GetNodes( ), true, vecstrData );
	m_veccQuants.resize( m_iData );
	if( m_aData )
		delete[] m_aData;
	m_aData = new CCompactMatrix[ m_iData ];

	m_vecstrGenes.resize( Answers.GetGenes( ) );
	for( i = 0; i < m_vecstrGenes.size( ); ++i )
		m_vecstrGenes[ i ] = Answers.GetGene( i );

	if( !CDatasetCompactImpl::Open( Answers, 0 ) )
		return false;
	for( i = 0; i < vecstrData.size( ); ++i ) {
		CDataPair	Datum;

		if( !( Datum.Open( vecstrData[ i ].c_str( ), false ) &&
			CDatasetCompactImpl::Open( Datum, i + 1 ) ) )
			return false; }

	for( i = 0; i < m_vecstrGenes.size( ); ++i )
		for( j = ( i + 1 ); j < m_vecstrGenes.size( ); ++j ) {
			for( k = 1; k < m_iData; ++k )
				if( m_aData[ k ].Get( i, j ) )
					break;
			if( k >= m_iData )
				m_aData[ 0 ].Set( i, j, 0 ); }

	return true; }

bool CDatasetCompactImpl::Open( const CDataPair& Datum, size_t iExp ) {
	vector<size_t>	veciGenes;
	size_t			i, j, iOne, iTwo;
	float			d;
	CCompactMatrix&	Target	= m_aData[ iExp ];

	m_veccQuants[ iExp ] = Datum.IsContinuous( ) ? -1 : Datum.GetValues( );
	Target.Initialize( m_vecstrGenes.size( ), (unsigned char)( Datum.GetValues( ) + 1 ),
		true );
	veciGenes.resize( m_vecstrGenes.size( ) );
	for( i = 0; i < veciGenes.size( ); ++i )
		veciGenes[ i ] = Datum.GetGene( m_vecstrGenes[ i ] );

	for( i = 0; i < veciGenes.size( ); ++i )
		if( ( iOne = veciGenes[ i ] ) != -1 )
			for( j = ( i + 1 ); j < veciGenes.size( ); ++j )
				if( ( ( iTwo = veciGenes[ j ] ) != -1 ) &&
					!CMeta::IsNaN( d = Datum.Get( iOne, iTwo ) ) )
					Target.Set( i, j, (unsigned char)( Datum.Quantify( d ) + 1 ) );

	return true; }

bool CDatasetCompact::Open( const char* szDataDir, const IBayesNet* pBayesNet ) {

	return CDatasetCompactImpl::Open( szDataDir, pBayesNet ); }

bool CDatasetCompact::Open( const char* szDataDir, const IBayesNet* pBayesNet,
	const CGenes& GenesIn, const CGenes& GenesEx ) {

	return CDatasetCompactImpl::Open( szDataDir, pBayesNet, &GenesIn, &GenesEx ); }

bool CDatasetCompactImpl::Open( const char* szDataDir, const IBayesNet* pBayesNet,
	const CGenes* pGenesIn, const CGenes* pGenesEx ) {
	size_t						i;
	vector<string>				vecstrData;
	set<string>					setstrGenes;
	set<string>::const_iterator	iterGenes;

	if( pBayesNet->IsContinuous( ) )
		return false;

	m_iData = (uint32_t)OpenMax( szDataDir, pBayesNet->GetNodes( ), false, vecstrData, &setstrGenes );
	m_veccQuants.resize( m_iData );
	if( pGenesIn )
		for( i = 0; i < pGenesIn->GetGenes( ); ++i )
			setstrGenes.insert( pGenesIn->GetGene( i ).GetName( ) );
	if( pGenesEx )
		for( i = 0; i < pGenesEx->GetGenes( ); ++i )
			setstrGenes.erase( pGenesEx->GetGene( i ).GetName( ) );
	m_vecstrGenes.resize( setstrGenes.size( ) );
	for( i = 0,iterGenes = setstrGenes.begin( ); iterGenes != setstrGenes.end( );
		++iterGenes )
		m_vecstrGenes[ i++ ] = *iterGenes;

	if( m_aData )
		delete[] m_aData;
	m_aData = new CCompactMatrix[ m_iData ];

	for( i = 0; i < vecstrData.size( ); ++i ) {
		CDataPair	Datum;

		if( !( Datum.Open( vecstrData[ i ].c_str( ), false ) &&
			CDatasetCompactImpl::Open( Datum, i ) ) )
			return false; }

	if( pGenesIn )
		FilterGenes( *pGenesIn, CDat::EFilterInclude );
	if( pGenesEx )
		FilterGenes( *pGenesEx, CDat::EFilterExclude );

	return true; }

bool CDatasetCompact::FilterGenes( const char* szGenes, CDat::EFilter eFilt ) {
	ifstream	ifsm;
	CGenome		Genome;
	CGenes		Genes( Genome );

	ifsm.open( szGenes );
	if( !( ifsm.is_open( ) && Genes.Open( ifsm ) ) )
		return false;
	ifsm.close( );
	FilterGenes( Genes, eFilt );

	return true; }

void CDatasetCompact::FilterGenes( const CGenes& Genes, CDat::EFilter eFilt ) {

	CDatasetCompactImpl::FilterGenes( Genes, eFilt ); }

void CDatasetCompactImpl::FilterGenes( const CGenes& Genes, CDat::EFilter eFilt ) {
	vector<bool>	vecfGenes;
	size_t			i, j;

	if( !Genes.GetGenes( ) )
		return;

	vecfGenes.resize( GetGenes( ) );
	for( i = 0; i < vecfGenes.size( ); ++i )
		vecfGenes[ i ] = Genes.IsGene( GetGene( i ) );

	for( i = 0; i < vecfGenes.size( ); ++i ) {
		if( vecfGenes[ i ] ) {
			if( eFilt == CDat::EFilterInclude )
				continue;
			if( eFilt == CDat::EFilterExclude ) {
				for( j = ( i + 1 ); j < GetGenes( ); ++j )
					Remove( i, j );
				continue; } }
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			switch( eFilt ) {
				case CDat::EFilterInclude:
					if( !vecfGenes[ j ] )
						Remove( i, j );
					break;

				case CDat::EFilterTerm:
					if( !( vecfGenes[ i ] && vecfGenes[ j ] ) &&
						( !( vecfGenes[ i ] || vecfGenes[ j ] ) || GetDiscrete( i, j, 0 ) ) )
							Remove( i, j );
					break;

				case CDat::EFilterExclude:
					if( vecfGenes[ j ] )
						Remove( i, j );
					break; } } }

void CDatasetCompact::FilterAnswers( ) {
	size_t	i, j;

	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( IsExample( i, j ) && ( GetDiscrete( i, j, 0 ) == -1 ) )
				Remove( i, j ); }

bool CDatasetCompact::IsHidden( size_t iNode ) const {

	return CDataImpl::IsHidden( iNode ); }

size_t CDatasetCompact::GetDiscrete( size_t iX, size_t iY, size_t iNode ) const {

	return CDatasetCompactImpl::GetDiscrete( iX, iY, iNode ); }

size_t CDatasetCompactImpl::GetDiscrete( size_t iX, size_t iY, size_t iNode ) const {
	size_t	iMap;

	if( ( iMap = m_veciMapping[ iNode ] ) == -1 )
		return -1;

	return ( m_aData[ iMap ].Get( iX, iY ) - 1 ); }

float CDatasetCompact::GetContinuous( size_t iX, size_t iY, size_t iNode ) const {

	return CMeta::GetNaN( ); }

const string& CDatasetCompact::GetGene( size_t iGene ) const {

	return CDataImpl::GetGene( iGene ); }

size_t CDatasetCompact::GetGenes( ) const {

	return CDataImpl::GetGenes( ); }

bool CDatasetCompact::IsExample( size_t iX, size_t iY ) const {

	return CDatasetCompactImpl::IsExample( iX, iY ); }

bool CDatasetCompactImpl::IsExample( size_t iX, size_t iY ) const {
	size_t	i;

	for( i = 0; i < m_iData; ++i )
		if( m_aData[ i ].Get( iX, iY ) )
			return true;

	return false; }

const vector<string>& CDatasetCompact::GetGeneNames( ) const {

	return CDataImpl::GetGeneNames( ); }

size_t CDatasetCompact::GetExperiments( ) const {

	return CDataImpl::GetExperiments( ); }

size_t CDatasetCompact::GetGene( const string& strGene ) const {

	return CDataImpl::GetGene( strGene ); }

size_t CDatasetCompact::GetBins( size_t iExp ) const {

	return CDataImpl::GetBins( iExp ); }

void CDatasetCompact::Remove( size_t iX, size_t iY ) {

	CDatasetCompactImpl::Remove( iX, iY ); }

void CDatasetCompactImpl::Remove( size_t iX, size_t iY ) {
	size_t	i;

	for( i = 0; i < m_iData; ++i )
		m_aData[ i ].Set( iX, iY, 0 ); }

bool CDatasetCompactImpl::Open( const unsigned char* pbData ) {
	size_t	i;

	if( m_aData )
		delete[] m_aData;

	if( !( pbData = CDataImpl::OpenBinary( pbData ) ) )
		return false;
	m_iData = *(uint32_t*)pbData;
	pbData += sizeof(m_iData);
	m_aData = new CCompactMatrix[ m_iData ];
	for( i = 0; i < m_iData; ++i )
		if( !( pbData = m_aData[ i ].Open( pbData ) ) )
			return false;

	return true; }

bool CDatasetCompact::Open( istream& istm ) {
	size_t	i;

	if( m_aData )
		delete[] m_aData;

	if( !CDataImpl::OpenBinary( istm ) )
		return false;
	istm.read( (char*)&m_iData, sizeof(m_iData) );
	m_aData = new CCompactMatrix[ m_iData ];
	for( i = 0; i < m_iData; ++i )
		if( !m_aData[ i ].Open( istm ) )
			return false;

	return true; }

void CDatasetCompact::Save( ostream& ostm, bool fBinary ) const {

	fBinary ? SaveBinary( ostm ) : SaveText( ostm ); }

void CDatasetCompactImpl::SaveBinary( ostream& ostm ) const {
	size_t	i;

	CDataImpl::SaveBinary( ostm );
	ostm.write( (char*)&m_iData, sizeof(m_iData) );
	for( i = 0; i < m_iData; ++i )
		m_aData[ i ].Save( ostm ); }

void CDatasetCompactImpl::SaveText( ostream& ostm ) const {
	size_t	i, j, k;

	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( IsExample( i, j ) ) {
				ostm << GetGene( i ) << '\t' << GetGene( j );
				for( k = 0; k < GetExperiments( ); ++k )
					ostm << '\t' << GetDiscrete( i, j, k );
				ostm << endl; } }

const char	CDatasetCompactMap::c_szMap[]	= "CompactDataSetMap";

CDatasetCompactMap::CDatasetCompactMap( ) : m_pbData(NULL)
#ifdef _MSC_VER
	,m_hndlMap(0)
#endif // _MSC_VER
{ }

CDatasetCompactMap::~CDatasetCompactMap( ) {

#ifdef _MSC_VER
	if( m_pbData )
		UnmapViewOfFile( m_pbData );
	if( m_hndlMap )
		CloseHandle( m_hndlMap );
#else // _MSC_VER
	if( m_pbData )
		munmap( m_pbData, m_iData );
#endif // _MSC_VER
}

bool CDatasetCompactMap::Open( const char* szFile ) {
	size_t	i, j;
#ifdef _MSC_VER
	HANDLE	hFile;

	if( m_pbData )
		UnmapViewOfFile( m_pbData );
	if( m_hndlMap )
		CloseHandle( m_hndlMap );

	if( !( hFile = CreateFile( szFile, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING,
		FILE_ATTRIBUTE_READONLY, NULL ) ) )
		return false;
	if( !( m_hndlMap = CreateFileMapping( hFile, NULL, PAGE_READONLY, 0,
		(DWORD)( m_iData = GetFileSize( hFile, NULL ) ), c_szMap ) ) ) {
		CloseHandle( hFile );
		return false; }
	CloseHandle( hFile );

	if( !( m_pbData = (unsigned char*)MapViewOfFile( m_hndlMap, FILE_MAP_READ, 0, 0, 0 ) ) ) {
		CloseHandle( m_hndlMap );
		return false; }
	if( !CDatasetCompactImpl::Open( m_pbData ) ) {
		CloseHandle( m_hndlMap );
		return false; }
#else // _MSC_VER
	int			iFile;
	struct stat	sStat;

	if( m_pbData )
		munmap( m_pbData, m_iData );

	if( !( iFile = open( szFile, O_RDONLY ) ) )
		return false;
	fstat( iFile, &sStat );
	m_iData = sStat.st_size;

	if( ( m_pbData = (unsigned char*)mmap( NULL, m_iData, PROT_READ, MAP_SHARED, iFile, 0 ) ) == MAP_FAILED ) {
		g_CatBioUtils.error( "CDatasetCompactMap::Open( %s ) %s", szFile, strerror( errno ) );
		m_pbData = NULL;
		close( iFile );
		return false; }
	close( iFile );

	if( !CDatasetCompactImpl::Open( m_pbData ) ) {
		munmap( m_pbData, m_iData );
		return false; }
#endif // _MSC_VER

	m_Mask.Initialize( GetGenes( ) );
	for( i = 0; i < m_Mask.GetSize( ); ++i )
		for( j = ( i + 1 ); j < m_Mask.GetSize( ); ++j )
			m_Mask.Set( i, j, CDatasetCompact::IsExample( i, j ) );
	return true; }

void CDatasetCompactMap::Remove( size_t iX, size_t iY ) {

	m_Mask.Set( iX, iY, false ); }

bool CDatasetCompactMap::IsExample( size_t iX, size_t iY ) const {

	return m_Mask.Get( iX, iY ); }

}

#include "stdafx.h"
#include "dataset.h"
#include "bayesnet.h"
#include "genome.h"

namespace libBioUtils {

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

struct SIsGene {
	const CGenes&	m_Genes;
	bool			m_fIn;

	SIsGene( const CGenes& Genes, bool fIn ) : m_Genes(Genes), m_fIn(fIn) { }

	bool operator()( const string& strGene ) {

		return ( m_fIn == m_Genes.IsGene( strGene ) ); }
};

bool CDatasetCompact::Open( const CDataPair& Answers, const char* szDataDir,
	const IBayesNet* pBayesNet, bool fEverything ) {
	CGenome	Genome;
	CGenes	GenesIn( Genome ), GenesEx( Genome );

	return Open( Answers, szDataDir, pBayesNet, GenesIn, GenesEx, fEverything ); }

bool CDatasetCompact::Open( const CDataPair& Answers, const char* szDataDir,
	const IBayesNet* pBayesNet, const CGenes& GenesIn, const CGenes& GenesEx, bool fEverything ) {
	size_t			i;
	vector<string>	vecstrData, vecstrNodes;
	set<string>		setstrGenes;

	if( pBayesNet->IsContinuous( ) )
		return false;

	pBayesNet->GetNodes( vecstrNodes );
	m_iData = 1 + (uint32_t)OpenMax( szDataDir, vecstrNodes, true, vecstrData, fEverything ?
		&setstrGenes : NULL );
	m_veccQuants.resize( m_iData );
	if( m_aData )
		delete[] m_aData;
	m_aData = new CCompactMatrix[ m_iData ];

	if( fEverything ) {
		m_vecstrGenes.resize( setstrGenes.size( ) );
		copy( setstrGenes.begin( ), setstrGenes.end( ), m_vecstrGenes.begin( ) ); }
	else {
		m_vecstrGenes.resize( Answers.GetGenes( ) );
		for( i = 0; i < m_vecstrGenes.size( ); ++i )
			m_vecstrGenes[ i ] = Answers.GetGene( i ); }
	if( GenesIn.GetGenes( ) )
		remove_if( m_vecstrGenes.begin( ), m_vecstrGenes.end( ), SIsGene( GenesIn, false ) );
	if( GenesEx.GetGenes( ) )
		remove_if( m_vecstrGenes.begin( ), m_vecstrGenes.end( ), SIsGene( GenesEx, true ) );

	if( !CDatasetCompactImpl::Open( Answers, 0 ) )
		return false;
	for( i = 0; i < vecstrData.size( ); ++i ) {
		CDataPair	Datum;

		if( !( Datum.Open( vecstrData[ i ].c_str( ), false ) &&
			CDatasetCompactImpl::Open( Datum, i + 1 ) ) )
			return false; }

/*
	for( i = 0; i < m_vecstrGenes.size( ); ++i )
		for( j = ( i + 1 ); j < m_vecstrGenes.size( ); ++j ) {
			for( k = 1; k < m_iData; ++k )
				if( m_aData[ k ].Get( i, j ) )
					break;
			if( k >= m_iData )
				m_aData[ 0 ].Set( i, j, 0 ); }
*/

	return true; }

bool CDatasetCompact::Open( const CDataPair& Answers, const vector<string>& vecstrData,
	bool fEverything, bool fMemmap, size_t iSkip, bool fZScore ) {
	size_t	i, j, k;

	if( Answers.GetGenes( ) && Answers.IsContinuous( ) )
		return false;

	m_veciMapping.resize( m_iData = 1 + vecstrData.size( ) );
	for( i = 0; i < m_veciMapping.size( ); ++i )
		m_veciMapping[ i ] = i;
	m_veccQuants.resize( m_iData );
	if( m_aData )
		delete[] m_aData;
	m_aData = new CCompactMatrix[ m_iData ];

	if( fEverything ) {
		set<string>	setstrGenes;

		for( i = 0; i < Answers.GetGenes( ); ++i )
			setstrGenes.insert( Answers.GetGene( i ) );
		for( i = 0; i < vecstrData.size( ); ++i ) {
			CDat	Dat;

			if( !Dat.OpenGenes( vecstrData[ i ].c_str( ), iSkip ) )
					return false;
			for( j = 0; j < Dat.GetGenes( ); ++j )
				setstrGenes.insert( Dat.GetGene( j ) ); }
		m_vecstrGenes.resize( setstrGenes.size( ) );
		copy( setstrGenes.begin( ), setstrGenes.end( ), m_vecstrGenes.begin( ) ); }
	else {
		m_vecstrGenes.resize( Answers.GetGenes( ) );
		for( i = 0; i < m_vecstrGenes.size( ); ++i )
			m_vecstrGenes[ i ] = Answers.GetGene( i ); }

	if( !CDatasetCompactImpl::Open( Answers, 0 ) )
		return false;
	for( i = 0; i < vecstrData.size( ); ++i ) {
		CDataPair	Datum;

		if( !( Datum.Open( vecstrData[ i ].c_str( ), false, fMemmap, iSkip, fZScore ) &&
			CDatasetCompactImpl::Open( Datum, i + 1 ) ) )
			return false; }

	if( !fEverything && ( m_iData > 1 ) )
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
					Target.Set( i, j, (unsigned char)( Datum.Quantize( d ) + 1 ) );

	return true; }

bool CDatasetCompact::Open( const char* szDataDir, const IBayesNet* pBayesNet ) {

	return CDatasetCompactImpl::Open( szDataDir, pBayesNet ); }

bool CDatasetCompact::Open( const char* szDataDir, const IBayesNet* pBayesNet,
	const CGenes& GenesIn, const CGenes& GenesEx ) {

	if( !CDatasetCompactImpl::Open( szDataDir, pBayesNet, &GenesIn, &GenesEx ) )
		return false;
	CDataImpl::FilterGenes( this, GenesIn, CDat::EFilterInclude );
	CDataImpl::FilterGenes( this, GenesEx, CDat::EFilterExclude );

	return true; }

bool CDatasetCompactImpl::Open( const char* szDataDir, const IBayesNet* pBayesNet,
	const CGenes* pGenesIn, const CGenes* pGenesEx ) {
	size_t						i;
	vector<string>				vecstrData, vecstrNodes;
	set<string>					setstrGenes;
	set<string>::const_iterator	iterGenes;

	if( pBayesNet->IsContinuous( ) )
		return false;

	pBayesNet->GetNodes( vecstrNodes );
	m_iData = (uint32_t)OpenMax( szDataDir, vecstrNodes, false, vecstrData, &setstrGenes );
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

	return true; }

bool CDatasetCompact::FilterGenes( const char* szGenes, CDat::EFilter eFilt ) {
	ifstream	ifsm;
	CGenome		Genome;
	CGenes		Genes( Genome );

	ifsm.open( szGenes );
	if( !( ifsm.is_open( ) && Genes.Open( ifsm ) ) )
		return false;
	FilterGenes( Genes, eFilt );

	return true; }

void CDatasetCompact::FilterGenes( const CGenes& Genes, CDat::EFilter eFilt ) {

	CDataImpl::FilterGenes( this, Genes, eFilt ); }

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
	size_t	i, j, k, iVal;

	for( i = 0; i < GetGenes( ); ++i )
		for( j = ( i + 1 ); j < GetGenes( ); ++j )
			if( IsExample( i, j ) ) {
				ostm << GetGene( i ) << '\t' << GetGene( j );
				for( k = 0; k < GetExperiments( ); ++k ) {
					ostm << '\t';
					if( ( iVal = GetDiscrete( i, j, k ) ) == -1 )
						ostm << "-1";
					else
						ostm << iVal; }
				ostm << endl; } }

CDatasetCompactMap::CDatasetCompactMap( ) : m_pbData(NULL), m_hndlMap(0) { }

CDatasetCompactMap::~CDatasetCompactMap( ) {

	CMeta::Unmap( m_pbData, m_hndlMap, m_iData ); }

bool CDatasetCompactMap::Open( const char* szFile ) {
	size_t	i, j;

	CMeta::MapRead( m_pbData, m_hndlMap, m_iData, szFile );
	if( !CDatasetCompactImpl::Open( m_pbData ) ) {
		CMeta::Unmap( m_pbData, m_hndlMap, m_iData );
		return false; }

	m_Mask.Initialize( GetGenes( ) );
	for( i = 0; i < m_Mask.GetSize( ); ++i )
		for( j = ( i + 1 ); j < m_Mask.GetSize( ); ++j )
			m_Mask.Set( i, j, CDatasetCompact::IsExample( i, j ) );
	return true; }

void CDatasetCompactMap::Remove( size_t iX, size_t iY ) {

	m_Mask.Set( iX, iY, false ); }

bool CDatasetCompactMap::IsExample( size_t iX, size_t iY ) const {

	return m_Mask.Get( iX, iY ); }

bool CDatasetCompact::Open( const CGenes& GenesIn, const CGenes& GenesEx, const CDataPair& Answers,
	const vector<string>& vecstrPCLs, size_t iSkip, const IMeasure* pMeasure, const vector<float>& vecdQuants,
	const IBayesNet* pBayesNet ) {
	size_t					i, j, iPCL;
	set<string>				setstrGenes;
	set<string>::iterator	iterGene;

	g_CatBioUtils.notice( "CDatasetCompact::Open( %d ) opening PCL files",
		iSkip );

	m_veciMapping.resize( m_iData = 1 + (uint32_t)vecstrPCLs.size( ) );
	for( i = 0; i < m_veciMapping.size( ); ++i )
		m_veciMapping[ i ] = i;
	m_veccQuants.resize( m_iData );
	m_veccQuants[ 0 ] = Answers.GetValues( );
	for( i = 1; i < m_veccQuants.size( ); ++i )
		m_veccQuants[ i ] = (unsigned char)vecdQuants.size( );

	for( i = 0; i < Answers.GetGenes( ); ++i )
		setstrGenes.insert( Answers.GetGene( i ) );
	for( iPCL = 0; iPCL < vecstrPCLs.size( ); ++iPCL ) {
		ifstream	ifsm;

		ifsm.open( vecstrPCLs[ iPCL ].c_str( ) );
		if( !CDataImpl::OpenGenes( ifsm, false, true, setstrGenes ) ) {
			g_CatBioUtils.error( "CDatasetCompact::Open( %d ) could not open: %s", iSkip,
				vecstrPCLs[ iPCL ].c_str( ) );
			return false; } }
	if( GenesIn.GetGenes( ) ) {
		for( iterGene = setstrGenes.begin( ); iterGene != setstrGenes.end( ); ++iterGene )
			if( !GenesIn.IsGene( *iterGene ) )
				setstrGenes.erase( iterGene );
		for( i = 0; i < GenesIn.GetGenes( ); ++i )
			setstrGenes.insert( GenesIn.GetGene( i ).GetName( ) ); }
	if( GenesEx.GetGenes( ) )
		for( i = 0; i < GenesEx.GetGenes( ); ++i )
			setstrGenes.erase( GenesEx.GetGene( i ).GetName( ) );
	m_vecstrGenes.resize( setstrGenes.size( ) );
	copy( setstrGenes.begin( ), setstrGenes.end( ), m_vecstrGenes.begin( ) );

	if( m_aData )
		delete[] m_aData;
	m_aData = new CCompactMatrix[ m_iData ];
	if( !CDatasetCompactImpl::Open( Answers, 0 ) )
		return false;

	for( iPCL = 0; iPCL < vecstrPCLs.size( ); ++iPCL ) {
		CPCL			PCL;
		ifstream		ifsm;
		CDistanceMatrix	Dist;
		CDataPair		Datum;
		vector<size_t>	veciGenes;
		vector<string>	vecstrGenes;
		size_t			iGenes, iOne, iTwo;
		const float*	adOne;

		g_CatBioUtils.notice( "CDatasetCompact::Open( %d ) opening: %s", iSkip, vecstrPCLs[ iPCL ].c_str( ) );
		ifsm.open( vecstrPCLs[ iPCL ].c_str( ) );
		if( !PCL.Open( ifsm, iSkip ) ) {
			g_CatBioUtils.error( "CDatasetCompact::Open( %d ) could not open: %s", iSkip, vecstrPCLs[ iPCL ].c_str( ) );
			return 1; }
		if( pMeasure->IsRank( ) )
			PCL.RankTransform( );

		veciGenes.resize( PCL.GetGenes( ) );
		if( GenesIn.GetGenes( ) || GenesEx.GetGenes( ) )
			for( i = 0; i < PCL.GetGenes( ); ++i ) {
				const string&	strGene	= PCL.GetGene( i );

				if( GenesEx.GetGenes( ) && GenesEx.IsGene( strGene ) )
					veciGenes[ i ] = -1;
				else if( GenesIn.GetGenes( ) )
					veciGenes[ i ] = (unsigned int)( GenesIn.IsGene( strGene ) ? iGenes++ : -1 );
				else
					veciGenes[ i ] = (unsigned int)iGenes++;
				if( veciGenes[ i ] != -1 )
					vecstrGenes.push_back( strGene ); }
		else {
			vecstrGenes.resize( PCL.GetGenes( ) );
			copy( PCL.GetGeneNames( ).begin( ), PCL.GetGeneNames( ).end( ), vecstrGenes.begin( ) );
			for( i = 0; i < veciGenes.size( ); ++i )
				veciGenes[ i ] = i; }
		Dist.Initialize( vecstrGenes.size( ) );
		for( i = 0; i < Dist.GetSize( ); ++i )
			for( j = ( i + 1 ); j < Dist.GetSize( ); ++j )
				Dist.Set( i, j, CMeta::GetNaN( ) );
		for( i = 0; i < PCL.GetGenes( ); ++i ) {
			if( ( iOne = veciGenes[ i ] ) == -1 )
				continue;
			adOne = PCL.Get( i );
			for( j = ( i + 1 ); j < PCL.GetGenes( ); ++j )
				if( ( iTwo = veciGenes[ j ] ) != -1 )
					Dist.Set( iOne, iTwo, (float)pMeasure->Measure( adOne, PCL.GetExperiments( ), PCL.Get( j ),
						PCL.GetExperiments( ) ) ); }

		Datum.Open( vecstrGenes, Dist );
		Datum.Normalize( false );
		Datum.SetQuants( vecdQuants );
		if( !CDatasetCompactImpl::Open( Datum, iPCL + 1 ) )
			return false; }

	return true; }

}

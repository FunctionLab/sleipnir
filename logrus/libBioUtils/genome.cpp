#include "stdafx.h"
#include "genome.h"
#include "annotation.h"
#include "meta.h"

namespace libBioUtils {

CGene::CGene( const string& strName ) : CGeneImpl(strName) { }

CGeneImpl::CGeneImpl( const string& strName ) : m_strName(strName), m_iOntologies(0),
	m_apOntologies(NULL), m_apveciAnnotations(NULL), m_iSynonyms(0), m_astrSynonyms(NULL),
	m_fRNA(false), m_fDubious(false) { }

CGeneImpl::~CGeneImpl( ) {
	size_t	i;

	if( m_iOntologies ) {
		delete[] m_apOntologies;
		for( i = 0; i < m_iOntologies; ++i )
			delete m_apveciAnnotations[ i ];
		delete[] m_apveciAnnotations; }
	if( m_iSynonyms )
		delete[] m_astrSynonyms; }

CGeneImpl& CGeneImpl::operator=( const CGeneImpl& Gene ) {
	size_t	i, j;

	m_strName = Gene.m_strName;
	if( m_iSynonyms = Gene.m_iSynonyms ) {
		m_astrSynonyms = new string[ m_iSynonyms ];
		for( i = 0; i < m_iSynonyms; ++i )
			m_astrSynonyms[ i ] = Gene.m_astrSynonyms[ i ]; }
	if( m_iOntologies = Gene.m_iOntologies ) {
		m_apOntologies = new IOntology const*[ m_iOntologies ];
		m_apveciAnnotations = new vector<size_t>*[ m_iOntologies ];
		for( i = 0; i < m_iOntologies; ++i ) {
			m_apOntologies[ i ] = Gene.m_apOntologies[ i ];
			m_apveciAnnotations[ i ] = new vector<size_t>( );
			m_apveciAnnotations[ i ]->resize( Gene.m_apveciAnnotations[ i ]->size( ) );
			for( j = 0; j < m_apveciAnnotations[ i ]->size( ); ++j )
				(*m_apveciAnnotations[ i ])[ j ] = (*Gene.m_apveciAnnotations[ i ])[ j ]; } }

	return *this; }

bool CGene::AddAnnotation( const IOntology* pOntology, size_t iNode ) {
	size_t	i, iOnto;

	for( iOnto = 0; iOnto < m_iOntologies; ++iOnto )
		if( m_apOntologies[ iOnto ] == pOntology )
			break;
	if( iOnto >= m_iOntologies )
		IncrementOntologies( pOntology );

	for( i = 0; i < m_apveciAnnotations[ iOnto ]->size( ); ++i )
		if( (*m_apveciAnnotations[ iOnto ])[ i ] == iNode )
			return false;
	m_apveciAnnotations[ iOnto ]->push_back( iNode );
	return true; }

bool CGene::AddSynonym( const string& strSyn ) {
	size_t	i;
	string*	astrSynonyms;

	if( !strSyn.length( ) )
		g_CatBioUtils.warn( "CGene::AddSynonym( %s ) adding null synonym to %s",
			strSyn.c_str( ), m_strName.c_str( ) );
	if( strSyn == m_strName )
		return false;
	for( i = 0; i < m_iSynonyms; ++i )
		if( strSyn == m_astrSynonyms[ i ] )
			return false;

	astrSynonyms = new string[ ++m_iSynonyms ];
	for( i = 0; ( i + 1 ) < m_iSynonyms; ++i )
		astrSynonyms[ i ] = m_astrSynonyms[ i ];
	astrSynonyms[ i ] = strSyn;

	if( m_astrSynonyms )
		delete[] m_astrSynonyms;
	m_astrSynonyms = astrSynonyms;
	return true; }

void CGeneImpl::IncrementOntologies( const IOntology* pOntology ) {
	const IOntology**	apOntologies;
	vector<size_t>**	apveciAnnotations;

	apOntologies = new IOntology const*[ m_iOntologies + 1 ];
	if( m_apOntologies ) {
		memcpy( apOntologies, m_apOntologies, m_iOntologies * sizeof(*m_apOntologies) );
		delete[] m_apOntologies; }
	apOntologies[ m_iOntologies ] = pOntology;
	m_apOntologies = apOntologies;

	apveciAnnotations = new vector<size_t>*[ m_iOntologies + 1 ];
	if( m_apveciAnnotations ) {
		memcpy( apveciAnnotations, m_apveciAnnotations, m_iOntologies *
			sizeof(*m_apveciAnnotations) );
		delete[] m_apveciAnnotations; }
	apveciAnnotations[ m_iOntologies++ ] = new vector<size_t>;
	m_apveciAnnotations = apveciAnnotations; }

bool CGene::IsAnnotated( const IOntology* pOnto ) const {
	size_t	i;

	for( i = 0; i < m_iOntologies; ++i )
		if( pOnto == m_apOntologies[ i ] )
			return true;

	return false; }

bool CGene::IsAnnotated( const IOntology* pOnto, size_t iNode ) const {
	size_t	i, j;

	for( i = 0; i < m_iOntologies; ++i )
		if( pOnto == m_apOntologies[ i ] ) {
			for( j = 0; j < m_apveciAnnotations[ i ]->size( ); ++j )
				if( iNode == (*m_apveciAnnotations[ i ])[ j ] )
					return true;
			break; }

	return false; }


const string& CGene::GetName( ) const {

	return m_strName; }

size_t CGene::GetSynonyms( ) const {

	return m_iSynonyms; }

const string& CGene::GetSynonym( size_t iSyn ) const {

	return m_astrSynonyms[ iSyn ]; }

size_t CGene::GetOntologies( ) const {

	return m_iOntologies; }

const IOntology* CGene::GetOntology( size_t iOnto ) const {

	return m_apOntologies[ iOnto ]; }

size_t CGene::GetAnnotations( size_t iOnto ) const {

	return m_apveciAnnotations[ iOnto ]->size( ); }

size_t CGene::GetAnnotation( size_t iOnto, size_t iAnno ) const {

	return (*m_apveciAnnotations[ iOnto ])[ iAnno ]; }

void CGene::SetGloss( const string& strGloss ) {

	m_strGloss = strGloss; }

const string& CGene::GetGloss( ) const {

	return m_strGloss; }

void CGene::SetDubious( bool fDubious ) {

	m_fDubious = fDubious; }

bool CGene::GetDubious( ) const {

	return m_fDubious; }

void CGene::SetRNA( bool fRNA ) {

	m_fRNA = fRNA; }

bool CGene::GetRNA( ) const {

	return m_fRNA; }

const char	CGenomeImpl::c_szDubious[]	= "Dubious";
const char	CGenomeImpl::c_szORF[]		= "ORF";
const char*	CGenomeImpl::c_aszRNA[]		= { "ncRNA", "rRNA", "snRNA", "snoRNA", "tRNA",
	NULL };

CGenomeImpl::~CGenomeImpl( ) {
	size_t	i;

	for( i = 0; i < m_vecpGenes.size( ); ++i )
		delete m_vecpGenes[ i ]; }

bool CGenome::Open( istream& istm ) {
	static const size_t	c_iBuf	= 1024;
	char			szBuf[ c_iBuf ];
	vector<string>	vecstrLine, vecstrNames;
	size_t			i;

	while( istm.peek( ) != EOF ) {
		istm.getline( szBuf, c_iBuf - 1 );
		vecstrLine.clear( );
		CMeta::Tokenize( szBuf, vecstrLine );
		if( vecstrLine.size( ) < 16 )
			return false;

		if( vecstrLine[ 1 ] != c_szORF ) {
			for( i = 0; c_aszRNA[ i ]; ++i )
				if( vecstrLine[ 1 ] == c_aszRNA[ i ] )
					break;
			if( !c_aszRNA[ i ] )
				continue; }

		{
			CGene&	Gene	= AddGene( vecstrLine[ 3 ] );

//	1	Type		ORF, CDS, RNA, etc.
//	2	Qualifier	Dubious, Verified, etc.
//	3	ORF
//	4	Name
//	5	Aliases
//	15	DESC
			Gene.SetRNA( vecstrLine[ 1 ] != c_szORF );
			Gene.SetDubious( vecstrLine[ 2 ] == c_szDubious );
			if( vecstrLine[ 4 ].length( ) )
				AddSynonym( Gene, vecstrLine[ 4 ] );
			if( vecstrLine[ 5 ].length( ) ) {
				vecstrNames.clear( );
				CMeta::Tokenize( vecstrLine[ 5 ].c_str( ), vecstrNames, "|" );
				for( i = 0; i < vecstrNames.size( ); ++i )
					AddSynonym( Gene, vecstrNames[ i ] ); }
			Gene.SetGloss( vecstrLine[ 15 ] );
		} }

	return true; }

CGene& CGenome::AddGene( const string& strName ) {
	TMapStrI::const_iterator	iterGene;
	CGene*						pGene;

	if( ( iterGene = m_mapGenes.find( strName ) ) != m_mapGenes.end( ) )
		return GetGene( iterGene->second );

	pGene = new CGene( strName );
	m_vecpGenes.push_back( pGene );
	m_mapGenes[ strName ] = m_vecpGenes.size( ) - 1;
	return *pGene; }

CGene& CGenome::GetGene( size_t iGene ) const {

	return *m_vecpGenes[ iGene ]; }

size_t CGenome::GetGene( const string& strGene ) const {
	TMapStrI::const_iterator	iterGene;

	return ( ( ( iterGene = m_mapGenes.find( strGene ) ) == m_mapGenes.end( ) ) ? -1 :
		iterGene->second ); }

size_t CGenome::GetGenes( ) const {

	return m_vecpGenes.size( ); }

size_t CGenome::CountGenes( const IOntology* pOnto ) const {
	size_t	i, j, iRet;

	for( iRet = i = 0; i < m_vecpGenes.size( ); ++i )
		for( j = 0; j < m_vecpGenes[ i ]->GetOntologies( ); ++j )
			if( pOnto == m_vecpGenes[ i ]->GetOntology( j ) ) {
				iRet++;
				break; }

	return iRet; }

bool CGenome::AddSynonym( CGene& Gene, const string& strName ) {

	if( strName == Gene.GetName( ) )
		return false;

	m_mapGenes[ strName ] = m_mapGenes[ Gene.GetName( ) ];
	return Gene.AddSynonym( strName ); }

CGenes::CGenes( CGenome& Genome ) : CGenesImpl( Genome ) { }

CGenesImpl::CGenesImpl( CGenome& Genome ) : m_Genome(Genome) { }

bool CGenes::Open( istream& istm ) {
	static const size_t	c_iBuffer	= 1024;
	char	szBuf[ c_iBuffer ];

	if( istm.rdstate( ) != ios_base::goodbit )
		return false;

	m_mapGenes.clear( );
	m_vecpGenes.clear( );
	while( istm.peek( ) != EOF ) {
		istm.getline( szBuf, c_iBuffer - 1 );
		{
			CGene*	pGene	= &m_Genome.AddGene( szBuf );

			m_mapGenes[ pGene->GetName( ) ] = m_vecpGenes.size( );
			m_vecpGenes.push_back( pGene );
		} }

	return true; }

size_t CGenes::CountAnnotations( const IOntology* pOnto, size_t iNode, bool fKids ) const {
	size_t	i, iRet;

	for( iRet = i = 0; i < m_vecpGenes.size( ); ++i )
		if( pOnto->IsAnnotated( iNode, *m_vecpGenes[ i ], fKids ) )
			iRet++;

	return iRet; }

bool CGenes::Open( const vector<string>& vecstrGenes ) {
	size_t	i;

	m_mapGenes.clear( );
	m_vecpGenes.resize( vecstrGenes.size( ) );
	for( i = 0; i < m_vecpGenes.size( ); ++i ) {
		m_mapGenes[ vecstrGenes[ i ] ] = i;
		m_vecpGenes[ i ] = &m_Genome.AddGene( vecstrGenes[ i ] ); }

	return true; }

void CGenes::Filter( const CGenes& Genes ) {
	size_t	i, j, iSize;

	iSize = m_vecpGenes.size( );
	for( i = 0; i < Genes.GetGenes( ); ++i )
		for( j = 0; j < iSize; ++j )
			if( m_vecpGenes[ j ] == &Genes.GetGene( i ) ) {
				m_mapGenes.erase( m_vecpGenes[ j ]->GetName( ) );
				m_vecpGenes[ j ] = m_vecpGenes[ m_vecpGenes.size( ) - 1 ];
				iSize--;
				break; }
	m_vecpGenes.resize( iSize ); }

}

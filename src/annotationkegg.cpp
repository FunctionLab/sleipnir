#include "stdafx.h"
#include "annotation.h"
#include "genome.h"
#include "meta.h"

namespace Sleipnir {

COntologyKEGG::COntologyKEGG( ) {

	m_pOntology = this; }

/*!
 * \brief
 * Initializes the ontology using a KEGG "ko" file.
 * 
 * \param istm
 * Stream from which ko file is read.
 * 
 * \param Genome
 * Genome into which genes are inserted or read during annotation parsing.
 * 
 * \param strOrganism
 * Organism identifier for which annotations should be read (e.g. "SCE" or "HSA").
 * 
 * \param fSynonyms
 * If true, use the first ko synonym (when present) as the primary gene name; otherwise, use the ko file's
 * gene identifier.
 * 
 * \returns
 * True if the ontology was successfully initialized.
 * 
 * Given an input stream containing a ko file from the KEGG orthology, COntologyKEGG::Open parses both
 * the structure of KEGG (usually a flat set of non-hierarchical pathways) and the genes annotated to each
 * KEGG term.  Terms are identified by KEGG IDs (e.g. "ko00624"), and genes are identified by the ko ID or
 * first synonym as specified by the fSynonyms parameter.  Genes are retrieved from Genome if already present
 * or inserted if not; it is thus important to ensure that the proper primary gene names are used so as to
 * agree with any identifiers already present in Genome.
 * 
 * \remarks
 * For yeast, leaving fSynonyms false will use ORF IDs as the primary gene name (usually the desired
 * behavior); setting it true will use common names.  For mouse or human, leaving fSynonyms false will use
 * Entrez gene IDs as the primary gene name; setting it true will use HGNC symbols (usually the desired
 * behavior).  For worm, fSynonyms should always be false to use systematic transcript IDs as the primary
 * gene name.  For fly, leaving fSynonyms false will use Flybase gene IDs as the primary gene name (usually
 * the desired behavior); setting it true will use common names.  For other organisms, please inspect the
 * ko file.
 */
bool COntologyKEGG::Open( istream& istm, CGenome& Genome, const string& strOrganism, bool fSynonyms ) {
	SParserKEGG					sParser( istm, Genome, strOrganism, fSynonyms );
	size_t						i, j, iNode;
	TMapStrI::const_iterator	iterNode;
	vector<set<CGene*> >		vecsetpGenes;
	set<CGene*>::iterator		iterGene;

	g_CatSleipnir.info( "COntologyKEGG::Open( %s )", strOrganism.c_str( ) );
	Reset( );
	sParser.GetLine( );
	while( istm.peek( ) != EOF ) {
		if( !COntologyKEGGImpl::Open( sParser ) )
			return false;
		for( i = 0; i < sParser.m_vecstrIDs.size( ); ++i ) {
			if( ( iterNode = m_mapNodes.find( sParser.m_vecstrIDs[ i ] ) ) ==
				m_mapNodes.end( ) )
				m_mapNodes[ sParser.m_vecstrIDs[ i ] ] = iNode = m_mapNodes.size( );
			else
				iNode = iterNode->second;
			if( vecsetpGenes.size( ) <= iNode )
				vecsetpGenes.resize( iNode + 1 );
			for( j = 0; j < sParser.m_vecpGenes.size( ); ++j )
				vecsetpGenes[ iNode ].insert( sParser.m_vecpGenes[ j ] ); } }

	m_aNodes = new SNode[ m_iNodes = vecsetpGenes.size( ) ];
	for( iterNode = m_mapNodes.begin( ); iterNode != m_mapNodes.end( ); ++iterNode ) {
		i = iterNode->second;
		m_aNodes[ i ].m_strID = iterNode->first;
		m_aNodes[ i ].m_strGloss = sParser.m_mapGlosses[ iterNode->first ];
		m_aNodes[ i ].m_iGenes = vecsetpGenes[ i ].size( );
		m_aNodes[ i ].m_apGenes = new const CGene*[ m_aNodes[ i ].m_iGenes ];
		for( j = 0,iterGene = vecsetpGenes[ i ].begin( );
			iterGene != vecsetpGenes[ i ].end( ); ++j,++iterGene ) {
			(*iterGene)->AddAnnotation( this, i );
			m_aNodes[ i ].m_apGenes[ j ] = *iterGene; } }

	return true; }

const char	COntologyKEGGImpl::c_szKEGG[]		= "KEGG";
const char	COntologyKEGGImpl::c_szEntry[]		= "ENTRY";
const char	COntologyKEGGImpl::c_szName[]		= "NAME";
const char	COntologyKEGGImpl::c_szDefinition[]	= "DEFINITION";
const char	COntologyKEGGImpl::c_szClass[]		= "CLASS";
const char	COntologyKEGGImpl::c_szPath[]		= "PATH:";
const char	COntologyKEGGImpl::c_szBR[]			= "BR:";
const char	COntologyKEGGImpl::c_szDBLinks[]	= "DBLINKS";
const char	COntologyKEGGImpl::c_szGenes[]		= "GENES";
const char	COntologyKEGGImpl::c_szEnd[]		= "///";

COntologyKEGGImpl::SParserKEGG::SParserKEGG( istream& istm, CGenome& Genome, const string& strOrganism,
	bool fSynonyms ) : m_fSynonyms(fSynonyms), m_fOrganism(false), m_strOrganism(strOrganism),
	SParser( istm, Genome ) { }

void COntologyKEGGImpl::SParserKEGG::Reset( ) {

	m_vecpGenes.clear( );
	m_vecstrIDs.clear( ); }

COntologyKEGGImpl::COntologyKEGGImpl( ) : COntologyImpl( c_szKEGG ) { }

bool COntologyKEGGImpl::Open( SParserKEGG& sParser ) {

	sParser.Reset( );
	return ( OpenEntry( sParser ) && OpenName( sParser ) &&
		OpenDefinition( sParser ) && OpenClass( sParser ) &&
		OpenDBLinks( sParser ) && OpenGenes( sParser ) &&
		OpenEnd( sParser ) ); }

bool COntologyKEGGImpl::OpenEntry( SParserKEGG& sParser ) {

	return ( sParser.IsStart( c_szEntry ) && sParser.GetLine( ) ); }

bool COntologyKEGGImpl::OpenName( SParserKEGG& sParser ) {

	g_CatSleipnir.debug( "COntologyKEGGImpl::OpenName( ) %s", sParser.m_szLine );
	return ( sParser.IsStart( c_szName ) ? sParser.GetLine( ) : true ); }

bool COntologyKEGGImpl::OpenDefinition( SParserKEGG& sParser ) {

	if( !sParser.IsStart( c_szDefinition ) )
		return true;

	do
		if( !sParser.GetLine( ) )
			return false;
	while( isspace( sParser.m_szLine[ 0 ] ) );

	return true; }

bool COntologyKEGGImpl::OpenClass( SParserKEGG& sParser ) {
	size_t	i;

	if( !sParser.IsStart( c_szClass ) )
		return false;

	sParser.m_strGloss.clear( );
	i = strlen( c_szClass );
	memmove( sParser.m_szLine, sParser.m_szLine + i, strlen( sParser.m_szLine ) - i + 1 );
	sParser.m_fPathing = false;
	do
		if( !OpenGloss( sParser ) )
			return false;
	while( isspace( sParser.m_szLine[ 0 ] ) );

	return true; }

bool COntologyKEGGImpl::OpenGloss( SParserKEGG& sParser ) {
	char*			pchStartGloss;
	char*			pchEndGloss;
	char*			pchStartPath;
	char*			pchEndPath;
	vector<string>	vecstrIDs;
	size_t			i;

	for( pchStartGloss = sParser.m_szLine; isspace( *pchStartGloss ); ++pchStartGloss );
	if( ( pchEndGloss = strstr( pchStartGloss, c_szPath ) ) ||
		( pchEndGloss = strstr( pchStartGloss, c_szBR ) ) ) {
		pchStartPath = pchEndGloss + ( strncmp( pchEndGloss, c_szBR, strlen( c_szBR ) ) ?
			strlen( c_szPath ) : strlen( c_szBR ) );
		if( !( pchEndPath = strchr( pchStartPath, ']' ) ) )
			return false;
		*pchEndPath = 0;
		CMeta::Tokenize( pchStartPath, vecstrIDs, " ", true );
		for( i = 0; i < vecstrIDs.size( ); ++i )
			sParser.m_vecstrIDs.push_back( vecstrIDs[ i ] );
		if( pchEndGloss > ( pchStartGloss + 1 ) ) {
			*( pchEndGloss - 1 ) = 0;
			if( sParser.m_fPathing )
				sParser.m_strGloss.clear( );
			else if( sParser.m_strGloss.length( ) )
				sParser.m_strGloss += ' ';
			sParser.m_strGloss += pchStartGloss; }
		sParser.m_fPathing = true;
		for( i = 0; i < vecstrIDs.size( ); ++i )
			sParser.m_mapGlosses[ vecstrIDs[ i ] ] = sParser.m_strGloss; }
	else {
		if( sParser.m_fPathing ) {
			sParser.m_fPathing = false;
			sParser.m_strGloss.clear( ); }
		else if( sParser.m_strGloss.length( ) )
			sParser.m_strGloss += ' ';
		sParser.m_strGloss += pchStartGloss; }

	return sParser.GetLine( ); }

bool COntologyKEGGImpl::OpenDBLinks( SParserKEGG& sParser ) {

	if( !sParser.IsStart( c_szDBLinks ) )
		return true;

	do
		if( !sParser.GetLine( ) )
			return false;
	while( isspace( sParser.m_szLine[ 0 ] ) );

	return true; }

bool COntologyKEGGImpl::OpenGenes( SParserKEGG& sParser ) {
	size_t	i;

	if( !sParser.IsStart( c_szGenes ) )
		return true;

	i = strlen( c_szGenes );
	memmove( sParser.m_szLine, sParser.m_szLine + i, strlen( sParser.m_szLine ) - i + 1 );
	do
		if( !OpenOrganism( sParser ) )
			return false;
	while( isspace( sParser.m_szLine[ 0 ] ) );

	return true; }

bool COntologyKEGGImpl::OpenOrganism( SParserKEGG& sParser ) {
	char*	pch;
	size_t	i;

	for( pch = sParser.m_szLine; *pch && isspace( *pch ); ++pch );
	if( !*pch )
		return false;

	if( sParser.m_fOrganism ) {
		if( ( strlen( pch ) > 3 ) && ( pch[ 3 ] == ':' ) )
			sParser.m_fOrganism = false; }
	else if( !strncmp( pch, ( sParser.m_strOrganism + ':' ).c_str( ),
		i = ( sParser.m_strOrganism.length( ) + 1 ) ) ) {
		sParser.m_fOrganism = true;
		pch += i + 1; }
	if( sParser.m_fOrganism )
		while( *pch )
			pch = OpenGene( sParser, pch );

	return sParser.GetLine( ); }

char* COntologyKEGGImpl::OpenGene( SParserKEGG& sParser, char* pch ) {
	char*			pchEnd;
	char*			pchSyn;
	bool			fInc, fSyn;
	CGene*			pGene;
	string			strName;
	vector<string>	vecstrSynonyms;
	size_t			i;

	for( pchEnd = pch; *pchEnd && !isspace( *pchEnd ) && ( *pchEnd != '(' ); ++pchEnd );
	if( fInc = !!*pchEnd )
		fSyn = ( *pchEnd == '(' );
	*pchEnd = 0;
	strName = pch;

	if( fInc ) {
		++pchEnd;
		if( fSyn ) {
			pchSyn = pchEnd;
			for( ; *pchEnd != ')'; ++pchEnd );
			*(pchEnd++) = 0;
			vecstrSynonyms.push_back( pchSyn ); } }
	for( ; *pchEnd && !isspace( *pchEnd ); ++pchEnd );
	if( isspace( *pchEnd ) )
		++pchEnd;

	pGene = &sParser.m_Genome.AddGene( ( sParser.m_fSynonyms && !vecstrSynonyms.empty( ) ) ?
		vecstrSynonyms[ 0 ] : strName );
	if( sParser.m_fSynonyms )
		sParser.m_Genome.AddSynonym( *pGene, strName );
	for( i = 0; i < vecstrSynonyms.size( ); ++i )
		sParser.m_Genome.AddSynonym( *pGene, vecstrSynonyms[ i ] );
	sParser.m_vecpGenes.push_back( pGene );

	return pchEnd; }

bool COntologyKEGGImpl::OpenEnd( SParserKEGG& sParser ) {

	return ( sParser.IsStart( c_szEnd ) && sParser.GetLine( ) ); }

}

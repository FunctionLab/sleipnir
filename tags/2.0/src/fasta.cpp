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
#include "fasta.h"
#include "meta.h"

namespace Sleipnir {

const char CFASTAImpl::c_acComment[]	= "#";
const char CFASTAImpl::c_acHeader[]		= ">";

CFASTAImpl::CFASTAImpl( ) {

	m_szBuffer = new char[ c_iBufferSize ];
	pthread_mutex_init( &m_mutx, NULL ); }

CFASTAImpl::~CFASTAImpl( ) {

	pthread_mutex_destroy( &m_mutx );
	delete[] m_szBuffer; }

/*!
 * \brief
 * Opens a FASTA or WIG file and indexes the file without explicitly loading its contents.
 * 
 * \param szFile
 * Path to FASTA/WIG file to open.
 * 
 * \param setstrTypes
 * If nonempty, set of sequence types to be loaded; types not in the set are ignored.
 * 
 * \returns
 * True if file was loaded successfully; false otherwise.
 * 
 * \remarks
 * Supports FASTA and WIG files as described in CFASTA.  No data is loaded on open, but an index is created
 * over all genes and types of interest; a file handle is held open, and data is loaded as needed by the Get
 * methods.
 * 
 * \see
 * Save
 */
bool CFASTA::Open( const char* szFile, const std::set<std::string>& setstrTypes ) {
	char*				pc;
	vector<string>		vecstrLine;
	TMapStrI::iterator	iterGene;
	size_t				i, iGene;

	m_ifsm.close( );
	m_ifsm.clear( );
	m_mapstriGenes.clear( );
	m_vecstrGenes.clear( );
	m_vecmapstriSequences.clear( );
	m_vecmapstrstrHeaders.clear( );

	m_ifsm.open( szFile, ios_base::binary );
	if( !m_ifsm.is_open( ) )
		return false;

	while( !m_ifsm.eof( ) ) {
		m_ifsm.getline( m_szBuffer, c_iBufferSize );
		m_szBuffer[ c_iBufferSize - 1 ] = 0;
		if( !m_szBuffer[ 0 ] || !strchr( c_acHeader, m_szBuffer[ 0 ] ) )
			continue;
		for( pc = m_szBuffer + strlen( m_szBuffer ) - 1; ( pc >= m_szBuffer ) && strchr( "\n\r", *pc ); --pc )
			*pc = 0;
		for( pc = m_szBuffer + 1; *pc && isspace( *pc ); ++pc );
		vecstrLine.clear( );
		CMeta::Tokenize( pc, vecstrLine );
		if( vecstrLine.empty( ) ) {
			g_CatSleipnir( ).warn( "CFASTA::Open( %s ) invalid header line: %s", szFile, m_szBuffer );
			continue; }
		if( vecstrLine.size( ) < 2 )
			vecstrLine.push_back( "" );
		if( !setstrTypes.empty( ) && ( setstrTypes.find( vecstrLine[ 1 ] ) == setstrTypes.end( ) ) )
			continue;
		m_setstrTypes.insert( vecstrLine[ 1 ] );
		if( ( iterGene = m_mapstriGenes.find( vecstrLine[ 0 ] ) ) == m_mapstriGenes.end( ) ) {
			m_mapstriGenes[ vecstrLine[ 0 ] ] = iGene = m_vecstrGenes.size( );
			m_vecstrGenes.push_back( vecstrLine[ 0 ] );
			while( m_vecmapstriSequences.size( ) <= iGene )
				m_vecmapstriSequences.push_back( TMapStrI( ) ); }
		else
			iGene = iterGene->second;
		m_vecmapstriSequences[ iGene ][ vecstrLine[ 1 ] ] = m_ifsm.tellg( );
		if( vecstrLine.size( ) > 2 ) {
			string	strHeader;

			while( m_vecmapstrstrHeaders.size( ) <= iGene )
				m_vecmapstrstrHeaders.push_back( map<string, string>( ) );
			strHeader = vecstrLine[ 2 ];
			for( i = 3; i < vecstrLine.size( ); ++i )
				strHeader += '\t' + vecstrLine[ i ];
			m_vecmapstrstrHeaders[ iGene ][ vecstrLine[ 1 ] ] = strHeader; } }
	while( m_vecmapstrstrHeaders.size( ) <= iGene )
		m_vecmapstrstrHeaders.push_back( map<string, string>( ) );

	return true; }

/*!
 * \brief
 * Saves a copy of the FASTA file to the given output stream.
 * 
 * \param ostm
 * Output stream to which FASTA file is saved.
 * 
 * \param iWrap
 * If given, column at which output FASTA is linewrapped.
 * 
 * \remarks
 * Currently only supports FASTA files, not WIGs.
 * 
 * \see
 * Open
 */
void CFASTA::Save( std::ostream& ostm, size_t iWrap ) const {
	size_t	i, j, iGene;

	for( iGene = 0; iGene < GetGenes( ); ++iGene ) {
		vector<SFASTASequence>	vecsSequences;

		Get( iGene, vecsSequences );
		for( i = 0; i < vecsSequences.size( ); ++i ) {
			const SFASTASequence&	sSequence	= vecsSequences[ i ];
			string					strSequence;
			bool					fIntron;

			ostm << c_acHeader[ 0 ] << ' ' << GetGene( iGene );
			if( !sSequence.m_strType.empty( ) )
				ostm << '\t' << sSequence.m_strType;
			if( ( strSequence = GetHeader( iGene, sSequence.m_strType ) ).length( ) ) {
				if( sSequence.m_strType.empty( ) )
					ostm << '\t';
				ostm << '\t' << strSequence; }
			ostm << endl;

			for( strSequence = "",fIntron = sSequence.m_fIntronFirst,j = 0;
				j < sSequence.m_vecstrSequences.size( ); ++j,fIntron = !fIntron ) {
				string	strCur;

				strCur = sSequence.m_vecstrSequences[ j ];
				transform( strCur.begin( ), strCur.end( ), strCur.begin( ), fIntron ? ::tolower : ::toupper );
				strSequence += strCur; }
			for( j = 0; j < strSequence.length( ); j += iWrap )
				ostm << strSequence.substr( j, iWrap ) << endl; } } }

bool CFASTAImpl::Get( size_t iGene, vector<SFASTASequence>* pvecsSequences,
	vector<SFASTAWiggle>* pvecsValues ) const {

	if( iGene == -1 )
		return false;

	const TMapStrI&				mapstriSequences	= m_vecmapstriSequences[ iGene ];
	TMapStrI::const_iterator	iterGene;
	char*						pc;
	bool						fOpen;

	pthread_mutex_lock( &m_mutx );
	fOpen = m_ifsm.is_open( );
	pthread_mutex_unlock( &m_mutx );
	if( !fOpen )
		return false;
	for( iterGene = mapstriSequences.begin( ); iterGene != mapstriSequences.end( ); ++iterGene ) {
		SFASTASequence	sSequence;
		SFASTAWiggle	sValues;
		string			strSequence;

		sSequence.m_strType = sValues.m_strType = iterGene->first;
		pthread_mutex_lock( &m_mutx );
		m_ifsm.clear( );
		m_ifsm.seekg( iterGene->second );
		if( (size_t)m_ifsm.tellg( ) != iterGene->second ) {
			g_CatSleipnir( ).error( "CFASTA::Get( %d ) error parsing: %s %s at %d (%d)", iGene,
				GetGene( iGene ).c_str( ), iterGene->first.c_str( ), iterGene->second,
				(size_t)m_ifsm.tellg( ) );
			pthread_mutex_unlock( &m_mutx );
			return false; }
		while( !m_ifsm.eof( ) ) {
			m_ifsm.getline( m_szBuffer, c_iBufferSize );
			m_szBuffer[ c_iBufferSize - 1 ] = 0;
			if( !m_szBuffer[ 0 ] || strchr( c_acComment, m_szBuffer[ 0 ] ) )
				continue;
			if( strchr( c_acHeader, m_szBuffer[ 0 ] ) )
				break;
			if( pvecsValues )
				sValues.m_vecdValues.push_back( (float)atof( m_szBuffer ) );
			if( pvecsSequences ) {
				for( pc = m_szBuffer + strlen( m_szBuffer ) - 1;
					( pc >= m_szBuffer ) && strchr( "\n\r", *pc ); --pc )
					*pc = 0;
				strSequence += m_szBuffer; } }
		pthread_mutex_unlock( &m_mutx );

		if( pvecsValues && !Get( iGene, *pvecsValues, iterGene->second, sValues ) )
			return false;
		if( pvecsSequences && !Get( iGene, *pvecsSequences, iterGene->second, strSequence, sSequence ) )
			return false; }

	return true; }

bool CFASTAImpl::Get( size_t iGene, std::vector<SFASTASequence>& vecsSequences, size_t iOffset,
	const std::string& strSequence, SFASTASequence& sSequence ) const {
	size_t	iBegin, iEnd;

	if( strSequence.empty( ) ) {
		g_CatSleipnir( ).debug( "CFASTA::Get( %d ) no sequence found: %s %s at %d", iGene,
			GetGene( iGene ).c_str( ), sSequence.m_strType.c_str( ), iOffset );
		return true; }
	for( iBegin = 0; iBegin < strSequence.size( ); iBegin = iEnd ) {
		bool	fBegin;
		string	strCur;

		fBegin = !!isupper( strSequence[ iBegin ] );
		if( !iBegin )
			sSequence.m_fIntronFirst = !fBegin;
		for( iEnd = iBegin + 1; ( iEnd < strSequence.size( ) ) &&
			( fBegin == !!isupper( strSequence[ iEnd ] ) ); ++iEnd );
		strCur = strSequence.substr( iBegin, iEnd - iBegin );
		transform( strCur.begin( ), strCur.end( ), strCur.begin( ), ::toupper );
		sSequence.m_vecstrSequences.push_back( strCur ); }
	vecsSequences.push_back( sSequence );

	return true; }

bool CFASTAImpl::Get( size_t iGene, vector<SFASTAWiggle>& vecsValues, size_t iOffset,
	SFASTAWiggle& sValues ) const {

	if( sValues.m_vecdValues.empty( ) ) {
		g_CatSleipnir( ).debug( "CFASTA::Get( %d ) no values found: %s %s at %d", iGene,
			GetGene( iGene ).c_str( ), sValues.m_strType.c_str( ), iOffset );
		return true; }
	vecsValues.push_back( sValues );

	return true; }

}

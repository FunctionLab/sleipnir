#include "stdafx.h"
#include "fasta.h"
#include "meta.h"

namespace Sleipnir {

const char CFASTAImpl::c_acComment[]	= "#";
const char CFASTAImpl::c_acHeader[]		= ">";

bool CFASTA::Open( const char* szFile ) {
	const char*			pc;
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
		for( pc = m_szBuffer + 1; *pc && isspace( *pc ); ++pc );
		vecstrLine.clear( );
		CMeta::Tokenize( pc, vecstrLine );
		if( vecstrLine.empty( ) ) {
			g_CatSleipnir.warn( "CFASTA::Open( %s ) invalid header line: %s", szFile, m_szBuffer );
			continue; }
		if( vecstrLine.size( ) < 2 )
			vecstrLine.push_back( "" );
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

void CFASTA::Save( ostream& ostm, size_t iWrap ) const {
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

bool CFASTA::Get( size_t iGene, vector<SFASTASequence>& vecsSequences ) const {
	const TMapStrI&				mapstriSequences	= m_vecmapstriSequences[ iGene ];
	TMapStrI::const_iterator	iterGene;
	size_t						iBegin, iEnd;

	if( !m_ifsm.is_open( ) )
		return false;
	for( iterGene = mapstriSequences.begin( ); iterGene != mapstriSequences.end( ); ++iterGene ) {
		SFASTASequence	sSequence;
		string			strSequence;

		sSequence.m_strType = iterGene->first;
		m_ifsm.clear( );
		m_ifsm.seekg( iterGene->second );
		if( (size_t)m_ifsm.tellg( ) != iterGene->second ) {
			g_CatSleipnir.error( "CFASTA::Get( %d ) error parsing: %s %s at %d (%d)", iGene,
				GetGene( iGene ).c_str( ), iterGene->first.c_str( ), iterGene->second, (size_t)m_ifsm.tellg( ) );
			return false; }
		while( !m_ifsm.eof( ) ) {
			m_ifsm.getline( m_szBuffer, c_iBufferSize );
			m_szBuffer[ c_iBufferSize - 1 ] = 0;
			if( !m_szBuffer[ 0 ] || strchr( c_acComment, m_szBuffer[ 0 ] ) )
				continue;
			if( strchr( c_acHeader, m_szBuffer[ 0 ] ) )
				break;
			strSequence += m_szBuffer; }

		if( strSequence.empty( ) ) {
			g_CatSleipnir.error( "CFASTA::Get( %d ) no sequence found: %s %s at %d", iGene,
				GetGene( iGene ).c_str( ), iterGene->first.c_str( ), iterGene->second );
			return false; }
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
		vecsSequences.push_back( sSequence ); }

	return true; }

}

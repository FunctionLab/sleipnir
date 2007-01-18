#include "stdafx.h"
#include "orthology.h"
#include "genome.h"
#include "meta.h"

namespace libBioUtils {

COrthologyImpl::~COrthologyImpl( ) {

	Reset( ); }

void COrthologyImpl::Reset( ) {
	TMapStrGenome::iterator	iterGenome;

	for( iterGenome = m_mapGenomes.begin( ); iterGenome != m_mapGenomes.end( ); ++iterGenome )
		delete iterGenome->second;
	m_mapGenomes.clear( );
	m_mapOrthology.clear( );
	m_vecvecpGenes.clear( ); }

bool COrthology::Open( istream& istm ) {
	TMapStrGenome::iterator	iterGenome;
	vector<string>			vecstrLine;
	char*					acBuf;
	size_t					i, j;
	string					strOrganism, strGene;
	CGenome*				pGenome;

	Reset( );
	acBuf = new char[ c_iBufferSize ];
	while( istm.peek( ) != EOF ) {
		istm.getline( acBuf, c_iBufferSize - 1 );
		vecstrLine.clear( );
		CMeta::Tokenize( acBuf, vecstrLine );
		if( vecstrLine.empty( ) )
			continue;

		m_vecvecpGenes.resize( m_vecvecpGenes.size( ) + 1 );
		{
			vector<CGene*>&	vecpGenes	= m_vecvecpGenes[ m_vecvecpGenes.size( ) - 1 ];

			for( i = 0; i < vecstrLine.size( ); ++i ) {
				if( vecstrLine[ i ].length( ) == 0 )
					continue;
				if( ( j = vecstrLine[ i ].find( c_cOrgSep ) ) == string::npos ) {
					g_CatBioUtils.warn( "COrthology::Open( ) illegal gene token: %s",
						vecstrLine[ i ].c_str( ) );
					continue; }
				strOrganism = vecstrLine[ i ].substr( 0, j );
				strGene = vecstrLine[ i ].substr( j + 1 );
				if( ( iterGenome = m_mapGenomes.find( strOrganism ) ) == m_mapGenomes.end( ) )
					m_mapGenomes[ strOrganism ] = pGenome = new CGenome( );
				else
					pGenome = iterGenome->second;
				vecpGenes.push_back( &pGenome->AddGene( strGene ) ); }
		} }
	delete[] acBuf;

	return true; }

}

#include "stdafx.h"
#include "orthology.h"
#include "genome.h"
#include "meta.h"

namespace libBioUtils {

COrthologyImpl::~COrthologyImpl( ) {

	Reset( ); }

void COrthologyImpl::Reset( ) {
	size_t	i;

	for( i = 0; i < m_vecpGenomes.size( ); ++i )
		delete m_vecpGenomes[ i ];
	m_vecpGenomes.clear( );
	m_vecstrOrganisms.clear( );
	m_mapGenes.clear( );
	m_mapOrthology.clear( );
	m_vecvecpGenes.clear( ); }

bool COrthology::Open( istream& istm ) {
	vector<string>	vecstrLine;
	char*			acBuf;
	size_t			i, j;
	string			strOrganism, strGene;
	CGenome*		pGenome;
	CGene*			pGene;

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
				for( j = 0; j < m_vecstrOrganisms.size( ); ++j )
					if( strOrganism == m_vecstrOrganisms[ j ] )
						break;
				if( j < m_vecpGenomes.size( ) )
					pGenome = m_vecpGenomes[ j ];
				else
					m_vecpGenomes.push_back( pGenome = new CGenome( ) );
				vecpGenes.push_back( pGene = &pGenome->AddGene( strGene ) );
				m_mapGenes[ pGene ] = pGenome; }
		} }
	delete[] acBuf;

	return true; }

}

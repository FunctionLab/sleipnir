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
* Changes made Jun 2012 by Qian Zhu
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "stdafx.h"
#include "database.h"
#include "meta.h"
#include "bayesnetint.h"
#include "datapair.h"

namespace Sleipnir {

const char	CDatabaseImpl::c_acDAB[]		= ".dab";
const char	CDatabaseImpl::c_acQDAB[]		= ".qdab";
const char	CDatabaseImpl::c_acExtension[]	= ".db";

///////////////////////////////////////////////////////////////////////////////
// CDatabaselet
///////////////////////////////////////////////////////////////////////////////

CDatabaselet::CDatabaselet( bool useNibble) {
	m_useNibble = useNibble;
	m_pmutx = new pthread_mutex_t( );
	pthread_mutex_init( m_pmutx, NULL );
}

CDatabaselet::~CDatabaselet( ) {

	pthread_mutex_destroy( m_pmutx );
	delete m_pmutx;
	if(m_fstm.is_open()){
		m_fstm.close();
	}
}

/* original method for initializing databaselets, including writing header + pre-allocation */
bool CDatabaselet::Open( const std::string& strFile, const std::vector<std::string>& vecstrGenes,
	uint32_t iGenes, uint32_t iDatasets ) {
	uint32_t	iSize;
	size_t		i;
	char*		acFiller;

	m_fstm.clear( );
	/* Open with overwriting */
	m_fstm.open( strFile.c_str( ), ios_base::in | ios_base::out | ios_base::binary | ios_base::trunc );

	if( !m_fstm.is_open( ) ) {
		g_CatSleipnir( ).error( "CDatabaselet::Open( %s, %u, %u ) open failed", strFile.c_str( ), iGenes,
			iDatasets );
		return false; }


	m_iGenes = iGenes;
	m_iDatasets = iDatasets;
	m_vecstrGenes.resize( vecstrGenes.size( ) );
	copy( vecstrGenes.begin( ), vecstrGenes.end( ), m_vecstrGenes.begin( ) );

	//allocate space
	m_fstm.write( (char*)&m_iHeader, sizeof(m_iHeader) );
	m_fstm.write( (char*)&m_iGenes, sizeof(m_iGenes) );
	m_fstm.write( (char*)&m_iDatasets, sizeof(m_iDatasets) );

	iSize = m_vecstrGenes.size( );

	m_fstm.write((char*)&iSize, sizeof(iSize));

	//write gene-name for only the genes in the databaselets
	m_iHeader = sizeof(m_iHeader) + sizeof(m_iGenes) + sizeof(m_iDatasets) + sizeof(iSize);
	for( i = 0; i < m_vecstrGenes.size( ); ++i ) {
		m_fstm.write( m_vecstrGenes[ i ].c_str( ), m_vecstrGenes[ i ].size( ) + 1);
		m_iHeader += m_vecstrGenes[ i ].size( ) + 1;
	}

	m_fstm.seekp( 0, ios_base::beg);
	m_fstm.write( (char*)&m_iHeader, sizeof(m_iHeader) );

	//pre-allocations
	m_fstm.seekp( m_iHeader, ios_base::beg );
	acFiller = new char[ GetSizeGene( ) ];
	memset( acFiller, -1, GetSizeGene( ) );
	for( i = 0; i < m_vecstrGenes.size( ); ++i ){
		m_fstm.write( acFiller, GetSizeGene( ) );
	}
	delete[] acFiller;
	SetFile(strFile);

	return true;

}

/* simply opens the file without overwriting */
bool CDatabaselet::OpenNoOverwrite() {
	m_fstm.clear( );
	/* open without overwriting */
	m_fstm.open( strFileName.c_str( ), ios_base::in | ios_base::out | ios_base::binary);

	if( !m_fstm.is_open( ) ) {
		g_CatSleipnir( ).error( "CDatabaselet::Open( %s ) open failed", strFileName.c_str( ));
		return false;
	}
	return true;
}

bool CDatabaselet::OpenWrite( unsigned char bValue, size_t iOffset, ENibbles eNibbles,
	unsigned char* abImage ) {
	unsigned char	b;

	if(m_useNibble==0){
		eNibbles = ENibblesBoth;
	}

	if( abImage )
		iOffset -= m_iHeader;
	if( eNibbles != ENibblesBoth ) {
		if( abImage )
			b = abImage[ iOffset ];
		else {
			m_fstm.seekg( iOffset, ios_base::beg);
			b = m_fstm.get( );
		}
	}

	switch( eNibbles ) {
		case ENibblesLow:
			b = ( bValue & 0xF ) | ( b & 0xF0 );
			break;

		case ENibblesHigh:
			b = ( b & 0xF ) | ( bValue << 4 );
			break;

		case ENibblesBoth:
			b = bValue;
			break;
	}
	if( abImage )
		abImage[ iOffset ] = b;
	else {
		m_fstm.seekp( iOffset, ios_base::beg );
		m_fstm.put( b );
	}

	return true; }

/* original file writing method */
bool CDatabaselet::Open( const vector<CCompactFullMatrix>& vecData, size_t iBaseGenes, size_t iBaseDatasets, bool fBuffer ) {
	unsigned char*	abImage;
	size_t			iSize, iDatum, iGeneOne, iGeneTwo;
	unsigned char	bOne, bTwo;

	if( fBuffer ) {
		//iBaseGenes: gene id of first gene in each databaselet
		//iDataset: dataset id
		//printf("Number: %d %d %d %d\n", GetSizeGene(), GetSizePair(), iBaseGenes, iBaseDatasets);
		abImage = new unsigned char[ iSize = ( GetSizeGene( ) * m_vecstrGenes.size( ) ) ];
		m_fstm.seekg( m_iHeader, ios_base::beg );
		m_fstm.read( (char*)abImage, iSize );
	}
	else
		abImage = NULL;

	//vecData: # of genes in databaselet x # of genes user's list

	//if this is not the first dataset in the dataset block
	if( iBaseDatasets % 2 ){
		//iGeneOne: iterate over all genes in this databaselet (# of genes in each databaselet)
		for( iGeneOne = 0; iGeneOne < GetGenes( ); ++iGeneOne ){
			//iGeneTwo: iterate overall genes in user's gene list
			for( iGeneTwo = 0; iGeneTwo < vecData[ 0 ].GetColumns( ); ++iGeneTwo ){
				//bOne, get the value of the gene located at the position (iBaseGene + iGeneOne, iGeneTwo)
				if( bOne = vecData[ 0 ].Get( iBaseGenes + iGeneOne, iGeneTwo ) ){
					//Offset is: m_iHeader + (GetSizeGene() * iOne) + (GetSizePair() * iTwo) + iDataset (for byte case)
					OpenWrite( bOne - 1, GetOffset( iGeneOne, iGeneTwo, iBaseDatasets ), ENibblesHigh, abImage );
				}
			}
		}
	}

	for( iDatum = ( iBaseDatasets % 2 ); ( iDatum + 1 ) < vecData.size( ); iDatum += 2 ){
		for( iGeneOne = 0; iGeneOne < GetGenes( ); ++iGeneOne ){
			for( iGeneTwo = 0; iGeneTwo < vecData[ iDatum ].GetColumns( ); ++iGeneTwo ) {
				bOne = vecData[ iDatum ].Get( iBaseGenes + iGeneOne, iGeneTwo );
				bTwo = vecData[ iDatum + 1 ].Get( iBaseGenes + iGeneOne, iGeneTwo );
				if( !( bOne || bTwo ) )
					continue;
				bOne -= 1;
				bTwo -= 1;
				if(m_useNibble){
					OpenWrite( ( bOne & 0xF ) | ( bTwo << 4 ), GetOffset( iGeneOne, iGeneTwo, iBaseDatasets +
							iDatum ), ENibblesBoth, abImage );
				}else{
					OpenWrite( bOne, GetOffset( iGeneOne, iGeneTwo, iBaseDatasets + iDatum ), ENibblesBoth,
							abImage );
					OpenWrite( bTwo, GetOffset( iGeneOne, iGeneTwo, iBaseDatasets + iDatum + 1 ), ENibblesBoth,
							abImage );
				}
			}
		}
	}

	if( iDatum < vecData.size( ) ){
		for( iGeneOne = 0; iGeneOne < GetGenes( ); ++iGeneOne ){
			for( iGeneTwo = 0; iGeneTwo < vecData[ iDatum ].GetColumns( ); ++iGeneTwo ){
				if( bOne = vecData[ iDatum ].Get( iBaseGenes + iGeneOne, iGeneTwo ) ){
					OpenWrite( bOne - 1, GetOffset( iGeneOne, iGeneTwo, iBaseDatasets + iDatum ), ENibblesLow,
						abImage );
				}
			}
		}
	}
	if( fBuffer ) {
		m_fstm.seekp( m_iHeader, ios_base::beg );
		m_fstm.write( (char*)abImage, iSize );
		delete[] abImage;
	}

	return true; }

/* 	A faster and simpler writing method for the matrix.
	takes UcharFullMatrix
	and requires buffering to be enabled, and works only with byte output
*/
bool CDatabaselet::OpenFast( const vector<CUcharFullMatrix>& vecData, size_t iBaseGenes, size_t iBaseDatasets, bool fBuffer) {
	if(fBuffer){
		cerr << "Requires buferring to be enabled." << endl;
		return false;
	}
	if(m_useNibble){
		cerr << "Requires byte." << endl;
		return false;
	}

	unsigned char*	abImage;
	size_t			iSize, iDatum, iGeneOne, iGeneTwo;
	unsigned char	bOne, bTwo;

	abImage = new unsigned char[ iSize = ( GetSizeGene( ) * m_vecstrGenes.size( ) ) ];
	m_fstm.seekg( m_iHeader );
	m_fstm.read( (char*)abImage, iSize );

	for( iDatum = 0; iDatum  < vecData.size( ); iDatum ++ ){
		for( iGeneOne = 0; iGeneOne < GetGenes( ); ++iGeneOne ){
			size_t index = vecData[iDatum].GetGeneIndex(GetGene(iGeneOne));
			size_t iOffset = (GetSizeGene() * iGeneOne) + iBaseDatasets + iDatum;

			for( iGeneTwo = 0; iGeneTwo < vecData[iDatum].GetColumns(); ++iGeneTwo ){
				if( bOne = vecData[ iDatum].Get(index, iGeneTwo ) ){
					abImage[ iOffset + GetSizePair() * iGeneTwo ] = bOne - 1;
				}
			}
		}
	}

	m_fstm.seekp( m_iHeader );
	m_fstm.write( (char*)abImage, iSize );
	delete[] abImage;

	return true;
}

bool CDatabaselet::Get( size_t iOne, size_t iTwo,
		vector<unsigned char>& vecbData, unsigned char *charImage){
	size_t	i;
	size_t offset = GetOffset(iOne, iTwo) - m_iHeader;

	if(this->m_useNibble==false){
		vecbData.clear();
		vecbData.resize(GetSizePair());

		for(i=0; i<vecbData.size(); i++){
			vecbData[i] = charImage[offset + i];
		}
	}else{
		vecbData.clear();
		vecbData.resize(m_iDatasets);


		for(i=0; i<GetSizePair(); i++){
			unsigned char b = charImage[offset + i];
			unsigned char bValue = -1;
			if( ( bValue = ( b & 0xF ) ) == 0xF ){
				bValue = -1;
			}
			vecbData[ 2 * i ] = bValue;

			if( ( bValue = ( ( b >> 4 ) & 0xF ) ) == 0xF ){
				bValue = -1;
			}

			if((2 * i + 1)==m_iDatasets){
				break;
			}
			vecbData[ (2 * i) + 1 ] = bValue;
		}
	}

	return true;
}

/*	static function, combine multiple databaselets (that share the same genes, ie m_vecStrGenes),
	and output result to a single file, or output one-gene per file (if databaselet contains multiple genes)
 	bSplit: whether or not to output one-gene per file
	Works for both nibble and byte
*/
bool CDatabaselet::Combine(std::vector<CDatabaselet*>& vecDatabaselet,
		std::string strOutDirectory, bool bSplit){

	/* for checking on consistency of databaselets */
	bool bIsConsistent = true;
	bool fUseNibble;

	size_t i, j;
	uint32_t iGenes, iDatasets;

	CDatabaselet *first = vecDatabaselet[0];
	fUseNibble = first->m_useNibble;

	iGenes = first->GetGenes();
	iDatasets = first->GetDatasets();

	vector<string> vecGenes;
	vecGenes.resize(iGenes);

	for(i=0; i<iGenes; i++){
		vecGenes[i] = first->GetGene(i);
	}

	for(i=1; bIsConsistent && i<vecDatabaselet.size(); i++){
		if(iGenes!=vecDatabaselet[i]->GetGenes() || fUseNibble!=vecDatabaselet[i]->m_useNibble){
			bIsConsistent = false;
			break;
		}
		for(j=0; j<iGenes; j++){
			if(vecGenes[j]!=vecDatabaselet[i]->GetGene(j)){
				bIsConsistent = false;
				break;
			}
		}
		iDatasets+=vecDatabaselet[i]->GetDatasets();
	}

	if(!bIsConsistent){
		cerr << "Databaselets are not consistent!" << endl;
		return false;
	}

	/* load all Databaselets into memory, for efficiency */
	unsigned char **charImages =
			(unsigned char**)malloc(vecDatabaselet.size()*sizeof(unsigned char*));
	size_t iImageSize = iDatasets * iGenes * first->m_iGenes;
	charImages[0] = (unsigned char*)malloc(iImageSize*sizeof(unsigned char));
	for(i=1; i<vecDatabaselet.size(); i++){
		charImages[i] = charImages[i-1] + vecDatabaselet[i-1]->m_iDatasets * first->m_iGenes * iGenes;
	}

	/* read databaselet into charImages */
	for(i=0; i<vecDatabaselet.size(); i++){
		CDatabaselet *current = vecDatabaselet[i];
		if(current->m_fstm.is_open()){
			current->m_fstm.seekg(current->m_iHeader, ios_base::beg);
			current->m_fstm.read((char*) charImages[i], iImageSize);
		}else{
			cerr << "CDatabaselet is not open." << endl;
			free(charImages[0]);
			free(charImages);
			return false;
		}
	}

	/* splitting to one gene per file after combine */
	if(bSplit){

		for(i=0; i<iGenes; i++){

			/* open a new Databaselet containing only one gene */
			string thisGene = first->GetGene(i);
			string path = strOutDirectory + "/" + thisGene + ".db";
			vector<string> vecstrThisGene;
			vecstrThisGene.push_back(thisGene);

			/* Create a new Databaselet */
			size_t iSize;
			CDatabaselet DBS(first->m_useNibble);
			DBS.Open(path.c_str(), vecstrThisGene, first->m_iGenes, iDatasets);
			unsigned char *abImage = (unsigned char*)
				malloc( iSize = (DBS.GetSizeGene( ) * DBS.m_vecstrGenes.size( ) ));
			size_t iDatum;
			size_t iGeneOne, iGeneTwo;
			size_t offset2, offset3;
			iGeneOne = i;

			if(first->m_useNibble==false){
				/* m_iGenes is all the genes in the genome */
				for( iGeneTwo = 0; iGeneTwo < first->m_iGenes; ++iGeneTwo ){
					offset2 = DBS.GetSizePair()*iGeneTwo;
					int totalSum = 0;
					for( iDatum = 0; iDatum  < vecDatabaselet.size(); iDatum ++ ){
						vector<unsigned char> vc;
						CDatabaselet *current = vecDatabaselet[iDatum];
						current->Get( iGeneOne, iGeneTwo, vc, charImages[iDatum]);
						offset3 = offset2 + totalSum;
						for(j=0; j<vc.size(); j++){
							abImage[offset3 + j] = vc[j];
						}
						totalSum+=vc.size();
					}
				}
			}else{
				size_t j;
				unsigned char *abImage2 = (unsigned char*)
					malloc(iDatasets);

				/* m_iGenes is all the genes in the genome */
				for( iGeneTwo = 0; iGeneTwo < first->m_iGenes; ++iGeneTwo ){
					offset2 = DBS.GetSizePair() * iGeneTwo;
					int totalSum = 0;
					for( iDatum = 0; iDatum  < vecDatabaselet.size(); iDatum ++ ){
						vector<unsigned char> vc;
						CDatabaselet *current = vecDatabaselet[iDatum];
						current->Get( iGeneOne, iGeneTwo, vc, charImages[iDatum]);
						offset3 = totalSum;
						for(j=0; j<vc.size(); j++){
							abImage2[offset3+j] = vc[j];
						}
						totalSum+=vc.size();
					}
					for(j=0; j+1 < iDatasets; j+=2){
						abImage[offset2 + j / 2] = (abImage2[j] & 0xF) | (abImage2[j+1] << 4);
					}
					if(j<iDatasets){
						unsigned char bValue = abImage2[iDatasets - 1];
						unsigned char b = 255;
						abImage[offset2 + j / 2] = ( bValue & 0xF ) | ( b & 0xF0 );
					}
				}

				free(abImage2);
			}

			/* close fstream */
			if(DBS.m_fstm.is_open()){
				DBS.m_fstm.seekp( DBS.m_iHeader, ios_base::beg );
				DBS.m_fstm.write( (char*)abImage, iSize );
				DBS.m_fstm.close();
			}else{
				cerr << "CDatabaselet is not opened." << endl;
				free(abImage);
				free(charImages[0]);
				free(charImages);
				return false;
			}

			free(abImage);

		}

	/* do not split, just combine into one file */
	}else{

		vector<string> strTok;
		CMeta::Tokenize(first->strFileName.c_str(), strTok, "/");
		string path = strOutDirectory + "/" + strTok[strTok.size()-1];

		CDatabaselet DBS(first->m_useNibble);

		DBS.Open(path.c_str(), first->m_vecstrGenes, first->m_iGenes, iDatasets);

		size_t iDatum;
		size_t iSize;
		unsigned char *abImage = (unsigned char*)
				malloc( iSize = (DBS.GetSizeGene( ) * DBS.m_vecstrGenes.size( ) ) );
		size_t iGeneOne, iGeneTwo;
		size_t offset1, offset2, offset3;

		if(first->m_useNibble==false){
			for(iGeneOne = 0; iGeneOne < first->GetGenes(); ++iGeneOne){
				offset1 = DBS.GetSizeGene() * iGeneOne;
				for( iGeneTwo = 0; iGeneTwo < first->m_iGenes; ++iGeneTwo ){
					offset2 = DBS.GetSizePair()*iGeneTwo;
					int totalSum = 0;
					for( iDatum = 0; iDatum  < vecDatabaselet.size(); iDatum ++ ){
						vector<unsigned char> vc;
						CDatabaselet *current = vecDatabaselet[iDatum];
						current->Get( iGeneOne, iGeneTwo, vc, charImages[iDatum]);
						offset3 = offset1 + offset2 + totalSum;
						for(j=0; j<vc.size(); j++){
							abImage[offset3 + j] = vc[j];
						}
						totalSum+=vc.size();
					}
				}
			}
		}else{
			size_t j;
			unsigned char *abImage2 = (unsigned char*)
				malloc(DBS.m_iDatasets);
			/* m_iGenes is all the genes in the genome */
			for(iGeneOne = 0; iGeneOne < first->GetGenes(); ++iGeneOne){
				offset1 = DBS.GetSizeGene() * iGeneOne;
				for( iGeneTwo = 0; iGeneTwo < first->m_iGenes; ++iGeneTwo ){
					offset2 = DBS.GetSizePair()*iGeneTwo;
					int totalSum = 0;
					for( iDatum = 0; iDatum  < vecDatabaselet.size(); iDatum ++ ){
						vector<unsigned char> vc;
						CDatabaselet *current = vecDatabaselet[iDatum];
						current->Get( iGeneOne, iGeneTwo, vc, charImages[iDatum]);
						offset3 = totalSum;
						for(j=0; j<vc.size(); j++){
							abImage2[offset3 + j] = vc[j];
						}
						totalSum+=vc.size();
					}
					for(j=0; j+1 < iDatasets; j+=2){
						abImage[offset1 + offset2 + j / 2] = (abImage2[j] & 0xF) | (abImage2[j+1] << 4);
					}
					if(j<iDatasets){
						unsigned char bValue = abImage2[iDatasets - 1];
						unsigned char b = 255;
						abImage[offset1 + offset2 + j / 2] = ( bValue & 0xF ) | ( b & 0xF0 );
					}
				}
			}
			free(abImage2);
		}

		/* close the databaselet */
		if(DBS.m_fstm.is_open()){
			DBS.m_fstm.seekp( DBS.m_iHeader, ios_base::beg );
			DBS.m_fstm.write( (char*)abImage, iSize );
			DBS.m_fstm.close();
		}else{
			cerr << "CDatabaselet is not opened." << endl;
			free(abImage);
			free(charImages[0]);
			free(charImages);
			return false;
		}

		free(abImage);
	}

	free(charImages[0]);
	free(charImages);

	return true;
}

bool CDatabaselet::Get( size_t iOne, size_t iTwo, vector<unsigned char>& vecbData ) const {
	size_t	i;

	i = vecbData.size( );
	vecbData.resize( i + GetSizePair( ) );
	pthread_mutex_lock( m_pmutx );

	m_fstm.seekg( GetOffset( iOne, iTwo ) );
	m_fstm.read( (char*)&vecbData[ i ], GetSizePair( ) );
	pthread_mutex_unlock( m_pmutx );

	return true; }

bool CDatabaselet::Get( size_t iGene, vector<unsigned char>& vecbData, bool fReplace ) const {
	size_t	i;

	i = fReplace ? 0 : vecbData.size( );
	vecbData.resize( i + GetSizeGene( ) );
	pthread_mutex_lock( m_pmutx );

	m_fstm.seekg( GetOffset( iGene ) );
	m_fstm.read( (char*)&vecbData[ i ], GetSizeGene( ) );
	pthread_mutex_unlock( m_pmutx );

	return true; }

bool CDatabaselet::Get( size_t iGene, const vector<size_t>& veciGenes, vector<unsigned char>& vecbData,
	bool fReplace ) const {
	size_t	i, iOffset;

	iOffset = fReplace ? 0 : vecbData.size( );
	vecbData.resize( iOffset + ( veciGenes.size( ) * GetSizePair( ) ) );
	pthread_mutex_lock( m_pmutx );
	for( i = 0; i < veciGenes.size( ); ++i,iOffset += GetSizePair( ) ) {
		m_fstm.seekg( GetOffset( iGene, veciGenes[ i ] ) );
		m_fstm.read( (char*)&vecbData[ iOffset ], GetSizePair( ) );
	}
	pthread_mutex_unlock( m_pmutx );

	return true; }

bool CDatabaselet::Open( const std::string& strFile ) {
	uint32_t	iSize;
	char*		acBuffer;
	char*		pc;
	size_t		i;

	m_fstm.clear( );
	m_fstm.open( strFile.c_str( ), ios_base::binary | ios_base::in );
	if( !m_fstm.is_open( ) ) {
		g_CatSleipnir( ).error( "CDatabaselet::Open( %s ) open failed", strFile.c_str( ) );
		return false;
	}

	m_fstm.read( (char*)&m_iHeader, sizeof(m_iHeader) );
	m_fstm.read( (char*)&m_iGenes, sizeof(m_iGenes) );
	m_fstm.read( (char*)&m_iDatasets, sizeof(m_iDatasets) );
	m_fstm.read( (char*)&iSize, sizeof(iSize) );

	acBuffer = new char[ m_iHeader ];
	m_fstm.read( acBuffer, m_iHeader - sizeof(m_iHeader) - sizeof(m_iGenes) - sizeof(m_iDatasets) - sizeof(iSize) );
	m_vecstrGenes.resize( iSize );
	for( i = 0,pc = acBuffer; i < m_vecstrGenes.size( ); pc += m_vecstrGenes[ i++ ].length( ) + 1 )
		m_vecstrGenes[ i ] = pc;
	delete[] acBuffer;

	SetFile(strFile);

	return true; }

///////////////////////////////////////////////////////////////////////////////
// CDatabase
///////////////////////////////////////////////////////////////////////////////

bool CDatabase::Open( const std::vector<std::string>& vecstrGenes, const std::string& strInputDirectory,
	const IBayesNet* pBayesNet, const std::string& strOutputDirectory, size_t iFiles) {
	vector<string>	vecstrNodes, vecstrSubset;
	size_t			i, j;
	char			acNumber[ 16 ];
	string			strFile;

	if( !pBayesNet ) {
		g_CatSleipnir( ).error( "CDatabase::Open( %s, %d ) null Bayes net", strOutputDirectory.c_str( ), iFiles );
		return false; }

	Clear( );
	pBayesNet->GetNodes( vecstrNodes );
	for( i = 1; i < vecstrNodes.size( ); ++i )
		vecstrNodes[ i - 1 ] = strInputDirectory + '/' + vecstrNodes[ i ];

	if( vecstrNodes.size( ) )
		vecstrNodes.resize( vecstrNodes.size( ) - 1 );
	m_vecpDBs.resize( iFiles );
	for( i = 0; i < m_vecpDBs.size( ); ++i ) {
		m_vecpDBs[ i ] = new CDatabaselet( m_useNibble);
		vecstrSubset.clear( );
		for( j = i; j < vecstrGenes.size( ); j += m_vecpDBs.size( ) )
			vecstrSubset.push_back( vecstrGenes[ j ] );
#pragma warning(disable : 4996)
		sprintf( acNumber, "%08u", i );
#pragma warning(default : 4996)
		strFile = strOutputDirectory + '/' + acNumber + c_acExtension;
		if( !( i % 100 ) )
			g_CatSleipnir( ).notice( "CDatabase::Open( %s, %d ) initializing file %d/%d",
				strOutputDirectory.c_str( ), iFiles, i, m_vecpDBs.size( ) );
		if( !m_vecpDBs[ i ]->Open( strFile, vecstrSubset, vecstrGenes.size( ), vecstrNodes.size( ) ) ) {
			g_CatSleipnir( ).error( "CDatabase::Open( %s, %d ) could not open file %s",
				strOutputDirectory.c_str( ), iFiles, strFile.c_str( ) );
			return false; } }
	for( i = 0; i < vecstrGenes.size( ); ++i )
		m_mapstriGenes[ m_vecpDBs[ i % m_vecpDBs.size( ) ]->GetGene( i / m_vecpDBs.size( ) ) ] = i;

	return CDatabaseImpl::Open( vecstrGenes, vecstrNodes );
}

/* Version of Open() that takes a list of datasets as input. Key method */
bool CDatabase::Open( const std::vector<std::string>& vecstrGenes, const std::vector<std::string>& vecstrDatasets,
	const std::string& strInputDirectory, const std::string& strOutputDirectory, size_t iFiles){

	vector<string>	vecstrNodes, vecstrSubset;
	size_t			i, j;
	char			acNumber[ 16 ];
	string			strFile;

	Clear();

	vecstrNodes.resize(vecstrDatasets.size());
	for(i=0; i<vecstrDatasets.size(); i++){
		vecstrNodes[i] = strInputDirectory + '/' + vecstrDatasets[i];
	}

	m_vecpDBs.resize( iFiles );
	int iNumFilesOpen = 1000;
	for( i = 0; i < m_vecpDBs.size( ); ++i ) {
		m_vecpDBs[ i ] = new CDatabaselet( m_useNibble );
	}

	size_t k;

	for( i = 0; i < m_vecpDBs.size( ); ++i ) { //block size (such as 1000)
		if(i%iNumFilesOpen==0 && i>0){
			for(k=0; k<iNumFilesOpen; k++){
				m_vecpDBs[i-k-1]->CloseFile();
			}
		}
		vecstrSubset.clear( );
		for( j = i; j < vecstrGenes.size( ); j += m_vecpDBs.size( ) )
			vecstrSubset.push_back( vecstrGenes[ j ] ); //contains index for 1000, 2000, 3000th genes
		sprintf( acNumber, "%08u", i );
		if(iFiles>=vecstrGenes.size()){
			//if one gene per file, let databaselet filename be gene-name
			strFile = strOutputDirectory + '/' + vecstrSubset[0] + c_acExtension;
		}else{
			strFile = strOutputDirectory + '/' + acNumber + c_acExtension;
		}

		if( !( i % 100 ) )
			g_CatSleipnir( ).notice( "CDatabase::Open( %s, %d ) initializing file %d/%d",
				strOutputDirectory.c_str( ), iFiles, i, m_vecpDBs.size( ) );
		if( !m_vecpDBs[ i ]->Open( strFile, vecstrSubset, vecstrGenes.size( ), vecstrNodes.size( ) ) ) {
			g_CatSleipnir( ).error( "CDatabase::Open( %s, %d ) could not open file %s",
				strOutputDirectory.c_str( ), iFiles, strFile.c_str( ) );
			return false; } }

	for( i = 0; i < m_vecpDBs.size( ); ++i ){
		m_vecpDBs[i]->CloseFile();
	}

	for( i = 0; i < vecstrGenes.size( ); ++i )
		m_mapstriGenes[ m_vecpDBs[ i % m_vecpDBs.size( ) ]->GetGene( i / m_vecpDBs.size( ) ) ] = i;

	return CDatabaseImpl::Open( vecstrGenes, vecstrNodes ); 
}

/* the key Open() method for Data2DB conversion */
bool CDatabaseImpl::Open( const std::vector<std::string>& vecstrGenes,
	const std::vector<std::string>& vecstrFiles ) {
	size_t			i, j, k, iOne, iTwo, iOutBlock, iOutBase, iOutOffset, iInBlock, iInBase, iInOffset;
	vector<size_t>	veciGenes;
	float			d;

	/* define number of threads to concurrently process datasets */
	omp_set_num_threads(4);

	veciGenes.resize( vecstrGenes.size( ) );
	iOutBlock = ( m_iBlockOut == -1 ) ? m_vecpDBs.size( ) : m_iBlockOut;
	iInBlock = ( m_iBlockIn == -1 ) ? vecstrFiles.size( ) : m_iBlockIn;

	int iNumFilesOpen = 1000;

	/* blocking parameter, iOutBase: number of databaselets to process at a time */
	for( iOutBase = 0; iOutBase < m_vecpDBs.size( ); iOutBase += iOutBlock ) {
		vector<string>	vecstrMyGenes;
		vector<size_t>	veciMyGenes;

		for( iOutOffset = 0; ( iOutOffset < iOutBlock ) && ( ( iOutBase + iOutOffset ) < m_vecpDBs.size( ) );
			++iOutOffset ) {
			const CDatabaselet&	DB	= *m_vecpDBs[ iOutBase + iOutOffset ];

			for( i = 0; i < DB.GetGenes( ); ++i )
				vecstrMyGenes.push_back( DB.GetGene( i ) );
		}

		veciMyGenes.resize( vecstrMyGenes.size( ) );

		for( iInBase = 0; iInBase < vecstrFiles.size( ); iInBase += iInBlock ) {
			vector<CCompactFullMatrix>	vecData;
			vecData.resize( ( ( iInBase + iInBlock ) > vecstrFiles.size( ) ) ?
				( vecstrFiles.size( ) - iInBase ) : iInBlock );
			for( iInOffset = 0; iInOffset < vecData.size( ); ++iInOffset ) {
				CDataPair	Dat;

				if( !Dat.Open( (vecstrFiles[ iInBase + iInOffset ] + c_acDAB).c_str( ),
						false, m_fMemmap, 0, false, false ) ) {
				    if( !Dat.Open( (vecstrFiles[ iInBase + iInOffset ] + c_acQDAB).c_str( ), false, m_fMemmap ) ) {
						g_CatSleipnir( ).error( "CDatabaseImpl::Open( ) could not open %s",
							(vecstrFiles[ iInBase + iInOffset ] + c_acDAB).c_str( ) );
				    	g_CatSleipnir( ).error( "CDatabaseImpl::Open( ) could not open %s",
							(vecstrFiles[ iInBase + iInOffset ] + c_acQDAB).c_str( ) );
				    	return false;
				    }
				}

				for( i = 0; i < veciMyGenes.size( ); ++i )
					veciMyGenes[ i ] = Dat.GetGene( vecstrMyGenes[ i ] );
				for( i = 0; i < veciGenes.size( ); ++i )
					veciGenes[ i ] = Dat.GetGene( vecstrGenes[ i ] );

				if(m_useNibble){
					vecData[ iInOffset ].Initialize( veciMyGenes.size( ), veciGenes.size( ), 16, true );
				}else{
					vecData[ iInOffset ].Initialize( veciMyGenes.size( ), veciGenes.size( ), 256, true );
				}

				//#pragma omp parallel for \
				shared(Dat, veciGenes, veciMyGenes, vecData) \
				private(j, i) \
				schedule(static)
				for(i=0; i<veciMyGenes.size(); i++){
					size_t t = veciMyGenes[i];
					if(t==-1) continue;
					float *d_array = Dat.GetFullRow(t);
					for(j=0; j<veciGenes.size(); j++){
						size_t s = veciGenes[j];
						if(s == -1) continue;
						if(s == t) continue;
						vecData[iInOffset].Set(i,j,Dat.Quantize(d_array[s])+1);
					}
					free(d_array);
				}


			}

			for( i = iOutOffset = 0; ( iOutOffset < iOutBlock ) && ( ( iOutBase + iOutOffset ) <
				m_vecpDBs.size( ) ); ++iOutOffset ) {
				CDatabaselet&	DB	= *m_vecpDBs[ iOutBase + iOutOffset ];

				/* close files if too many file handles opened */
				if(iOutOffset>0 && (iOutOffset%iNumFilesOpen==0 || 
					(iOutBase + iOutOffset)==m_vecpDBs.size()-1)){
					for(k=0; k<iNumFilesOpen; k++){
						if(iOutOffset + iOutBase - 1 - k == 0){
							break;
						}
						m_vecpDBs[iOutOffset + iOutBase - 1 - k]->CloseFile();
					}
				}

				DB.OpenNoOverwrite();

				if( !( iOutOffset % 100 ) )
					cerr << "Processing offset " << iOutOffset << '/' << iOutBlock << endl;
				//iInBase, for B=2, iInBase = 0, 2, 4 (total 5 datasets)
				//i = 0 18 36 ...
				if( !DB.Open( vecData, i, iInBase, m_fBuffer ) ){
					return false;
				}

				i += DB.GetGenes( );
			}
		}
	}

	return true;
}

bool CDatabase::Open( const std::string& strInputDirectory ) {
	size_t			i, j;
	vector<string>	vecstrFiles;
	string			strFile;

	FOR_EACH_DIRECTORY_FILE(strInputDirectory, strFile)
		if( !CMeta::IsExtension( strFile, c_acExtension ) )
			continue;

		i = atoi( CMeta::Deextension( CMeta::Basename( strFile.c_str( ) ) ).c_str( ) );
		if( vecstrFiles.size( ) <= i )
			vecstrFiles.resize( i + 1 );
		else if( vecstrFiles[ i ].length( ) != 0 ) {
			g_CatSleipnir( ).error( "CDatabase::Open( %s ) duplicate file: %s (%d)", strInputDirectory.c_str( ),
				strFile.c_str( ), i );
			return false; }
		vecstrFiles[ i ] = strInputDirectory + '/' + strFile; }

	Clear( );
	m_vecpDBs.resize( vecstrFiles.size( ) );
	for( i = 0; i < m_vecpDBs.size( ); ++i ) {
		m_vecpDBs[ i ] = new CDatabaselet( m_useNibble);
		if( !m_vecpDBs[ i ]->Open( vecstrFiles[ i ] ) )
			return false;
		for( j = 0; j < m_vecpDBs[ i ]->GetGenes( ); ++j )
			m_mapstriGenes[ m_vecpDBs[ i ]->GetGene( j ) ] = ( j * m_vecpDBs.size( ) ) + i; }

	return true;
}




}

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
#include "cmdline.h"

static const char	c_szBick[]	= "[Bick]";
static const char	c_szBicd[]	= "[Bicd]";
static const char	c_szSamba[]	= "samba";
static const char	c_szFuzzy[]	= "fuzzy";
static const char	c_szList[]	= "list";
static const char	c_szParam[]	= "param";

int open_list( istream&, CGenome&, vector<float>&, vector<CGenes*>& );
int open_samba( istream&, CGenome&, vector<float>&, vector<CGenes*>& );
int open_param( istream&, CGenome&, vector<float>&, vector<CGenes*>& );
int open_fuzzy( istream&, size_t, CDat& );

int main( int iArgs, char** aszArgs ) {
	CDat								Dat;
	ifstream							ifsm;
	istream*							pistm;
	gengetopt_args_info					sArgs;
	CGenome								Genome;
	vector<float>						vecdWeights;
	vector<CGenes*>						vecpClusters;
	vector<string>						vecstrGenes;
	set<string>							setstrGenes;
	set<string>::const_iterator			iterGene;
	vector<vector<size_t> >				vecveciGenes;
	int									iRet;
	size_t								i, j, iCluster, iGeneOne, iGeneTwo, iOne, iTwo;
	CGenes*								pClusterOne;
	float								d, dOne;
	map<string, size_t>					mapstriGenes;
	map<string, size_t>::const_iterator	iterIndex;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg );

	if( sArgs.input_arg ) {
		ifsm.open( sArgs.input_arg );
		pistm = &ifsm; }
	else
		pistm = &cin;
	iRet = 1;
	if( !strcmp( c_szSamba, sArgs.type_arg ) )
		iRet = open_samba( *pistm, Genome, vecdWeights, vecpClusters );
	else if( !strcmp( c_szList, sArgs.type_arg ) )
		iRet = open_list( *pistm, Genome, vecdWeights, vecpClusters );
	else if( !strcmp( c_szParam, sArgs.type_arg ) )
		iRet = open_param( *pistm, Genome, vecdWeights, vecpClusters );
	else if( !strcmp( c_szFuzzy, sArgs.type_arg ) )
		iRet = open_fuzzy( *pistm, sArgs.skip_arg, Dat );
	else
		cmdline_parser_print_help( );
	if( iRet )
		return iRet;

	if( !vecpClusters.size( ) ) {
		cerr << "No clusters found!" << endl;
		return 1; }

	for( i = 0; i < vecpClusters.size( ); ++i )
		for( j = 0; j < vecpClusters[ i ]->GetGenes( ); ++j )
			setstrGenes.insert( vecpClusters[ i ]->GetGene( j ).GetName( ) );
	vecstrGenes.resize( setstrGenes.size( ) );
	for( i = 0,iterGene = setstrGenes.begin( ); iterGene != setstrGenes.end( ); ++i,++iterGene )
		vecstrGenes[ i ] = *iterGene;

	Dat.Open( vecstrGenes );
	for( i = 0; i < Dat.GetGenes( ); ++i )
		mapstriGenes[ Dat.GetGene( i ) ] = i;
	vecveciGenes.resize( vecpClusters.size( ) );
	for( i = 0; i < vecveciGenes.size( ); ++i ) {
		vecveciGenes[ i ].resize( vecpClusters[ i ]->GetGenes( ) );
		for( j = 0; j < vecveciGenes[ i ].size( ); ++j )
			vecveciGenes[ i ][ j ] = ( ( iterIndex = mapstriGenes.find(
				vecpClusters[ i ]->GetGene( j ).GetName( ) ) ) == mapstriGenes.end( ) ) ? -1 :
				iterIndex->second; }

	if( sArgs.size_flag )
		for( iCluster = 0; iCluster < vecpClusters.size( ); ++iCluster )
			vecdWeights[ iCluster ] = 1.0f / vecpClusters[ iCluster ]->GetGenes( );
	if( sArgs.counts_flag )
		for( iCluster = 0; iCluster < vecpClusters.size( ); ++iCluster ) {
			const vector<size_t>&	veciOne	= vecveciGenes[ iCluster ];

			cerr << "Processing cluster " << iCluster << "/" << vecpClusters.size( ) << endl;
			pClusterOne = vecpClusters[ iCluster ];
			dOne = vecdWeights[ iCluster ];
			for( iGeneOne = 0; iGeneOne < veciOne.size( ); ++iGeneOne ) {
				iOne = veciOne[ iGeneOne ];
				for( iGeneTwo = ( iGeneOne + 1 ); iGeneTwo < veciOne.size( ); ++iGeneTwo ) {
					float	d;

					iTwo = veciOne[ iGeneTwo ];
					d = Dat.Get( iOne, iTwo );
					Dat.Set( iOne, iTwo, ( CMeta::IsNaN( d ) ? 0 : d ) + dOne ); } } }
	else {
		vector<bool>	vecfClustered;

		vecfClustered.resize( Dat.GetGenes( ) );
		for( iGeneOne = 0; iGeneOne < vecfClustered.size( ); ++iGeneOne ) {
			const string&	strGene	= Dat.GetGene( iGeneOne );

			for( iCluster = 0; iCluster < vecpClusters.size( ); ++iCluster )
				if( vecpClusters[ iCluster ]->IsGene( strGene ) ) {
					vecfClustered[ iGeneOne ] = true;
					break; } }
		for( iGeneOne = 0; iGeneOne < Dat.GetGenes( ); ++iGeneOne ) {
			if( !vecfClustered[ iGeneOne ] )
				continue;
			for( iGeneTwo = ( iGeneOne + 1 ); iGeneTwo < Dat.GetGenes( ); ++iGeneTwo )
				if( vecfClustered[ iGeneTwo ] )
					Dat.Set( iGeneOne, iGeneTwo, 0 ); }
		for( iCluster = 0; iCluster < vecpClusters.size( ); ++iCluster ) {
			const vector<size_t>&	veciOne	= vecveciGenes[ iCluster ];

			cerr << "Processing cluster " << iCluster << "/" << vecpClusters.size( ) << endl;
			pClusterOne = vecpClusters[ iCluster ];
			dOne = vecdWeights[ iCluster ];
			for( iGeneOne = 0; iGeneOne < veciOne.size( ); ++iGeneOne ) {
				iOne = veciOne[ iGeneOne ];
				for( iGeneTwo = ( iGeneOne + 1 ); iGeneTwo < veciOne.size( ); ++iGeneTwo ) {
					iTwo = veciOne[ iGeneTwo ];
					if( CMeta::IsNaN( d = Dat.Get( iOne, iTwo ) ) || ( dOne > d ) )
						Dat.Set( iOne, iTwo, dOne ); } } } }

	Dat.Save( sArgs.output_arg );
	for( i = 0; i < vecpClusters.size( ); ++i )
		delete vecpClusters[ i ];

	return 0; }

struct SSorter {
	const vector<float>&	m_vecdWeights;

	SSorter( const vector<float>& vecdWeights ) : m_vecdWeights(vecdWeights) { }

	bool operator()( size_t iOne, size_t iTwo ) {

		return ( m_vecdWeights[ iOne ] < m_vecdWeights[ iTwo ] ); }
};

int open_samba( istream& istm, CGenome& Genome, vector<float>& vecdWeights, vector<CGenes*>& vecpClusters ) {
	static const size_t	c_iLine	= 1024;
	char				szLine[ c_iLine ];
	vector<string>		vecstrLine, vecstrGenes;
	CGenes*				pCluster;
	size_t				iCluster, iCur, i;
	vector<size_t>		veciIndices;
	vector<float>		vecdCopy;
	vector<CGenes*>		vecpCopy;

	while( !istm.eof( ) ) {
		istm.getline( szLine, c_iLine - 1 );
		if( !strcmp( szLine, c_szBick ) )
			continue;
		if( !strcmp( szLine, c_szBicd ) )
			break;
		vecstrLine.clear( );
		CMeta::Tokenize( szLine, vecstrLine );
		if( vecstrLine.size( ) != 2 ) {
			cerr << "Illegal line: " << szLine;
			return 1; }
		vecdCopy.push_back( (float)atof( vecstrLine[ 1 ].c_str( ) ) ); }

	iCluster = -1;
	while( !istm.eof( ) ) {
		istm.getline( szLine, c_iLine - 1 );
		if( !szLine[ 0 ] )
			break;
		vecstrLine.clear( );
		CMeta::Tokenize( szLine, vecstrLine );
		if( vecstrLine.size( ) != 3 ) {
			cerr << "Illegal line: " << szLine;
			return 1; }
		if( atoi( vecstrLine[ 1 ].c_str( ) ) != 1 )
			continue;
		if( ( iCur = atoi( vecstrLine[ 0 ].c_str( ) ) ) != iCluster ) {
			if( vecstrGenes.size( ) ) {
				vecpCopy.push_back( pCluster = new CGenes( Genome ) );
				pCluster->Open( vecstrGenes ); }
			vecstrGenes.clear( );
			iCluster = iCur; }
		vecstrGenes.push_back( vecstrLine[ 2 ] ); }
	if( vecstrGenes.size( ) ) {
		vecpCopy.push_back( pCluster = new CGenes( Genome ) );
		pCluster->Open( vecstrGenes ); }

	veciIndices.resize( vecdCopy.size( ) );
	for( i = 0; i < veciIndices.size( ); ++i )
		veciIndices[ i ] = i;
	sort( veciIndices.begin( ), veciIndices.end( ), SSorter( vecdCopy ) );
	for( i = 0; i < veciIndices.size( ); ++i ) {
		vecdWeights.push_back( vecdCopy[ veciIndices[ i ] ] );
		vecpClusters.push_back( vecpCopy[ veciIndices[ i ] ] ); }

	return 0; }

int open_list( istream& istm, CGenome& Genome, vector<float>& vecdWeights, vector<CGenes*>& vecpClusters ) {
	static const size_t	c_iLine	= 1024;
	char			szLine[ c_iLine ];
	vector<string>	vecstrLine, vecstrGenes;
	CGenes*			pCluster;
	size_t			iCluster, iCur;

	iCluster = -1;
	while( !istm.eof( ) ) {
		istm.getline( szLine, c_iLine - 1 );
		if( !szLine[ 0 ] )
			break;
		vecstrLine.clear( );
		CMeta::Tokenize( szLine, vecstrLine );
		if( vecstrLine.size( ) != 2 ) {
			cerr << "Illegal line: " << szLine;
			return 1; }
		if( ( iCur = atoi( vecstrLine[ 1 ].c_str( ) ) ) != iCluster ) {
			if( vecstrGenes.size( ) ) {
				vecdWeights.push_back( 1 );
				vecpClusters.push_back( pCluster = new CGenes( Genome ) );
				pCluster->Open( vecstrGenes ); }
			vecstrGenes.clear( );
			iCluster = iCur; }
		vecstrGenes.push_back( vecstrLine[ 0 ] ); }
	if( vecstrGenes.size( ) ) {
		vecdWeights.push_back( 1 );
		vecpClusters.push_back( pCluster = new CGenes( Genome ) );
		pCluster->Open( vecstrGenes ); }

	return 0; }

int open_param( istream& istm, CGenome& Genome, vector<float>& vecdWeights, vector<CGenes*>& vecpClusters ) {
	static const size_t	c_iLine	= 1024;
	char			szLine[ c_iLine ];
	vector<string>	vecstrLine, vecstrGenes;
	CGenes*			pCluster;
	size_t			iCluster, iCur;
	float			dParam;

	while( !istm.eof( ) ) {
		istm.getline( szLine, c_iLine - 1 );
		if( !szLine[ 0 ] )
			break;
		vecstrLine.clear( );
		CMeta::Tokenize( szLine, vecstrLine );
		if( vecstrLine.size( ) != 2 ) {
			cerr << "Illegal line: " << szLine;
			return 1; }
		if( vecstrLine[ 0 ] == c_szParam ) {
			if( vecstrGenes.size( ) ) {
				vecdWeights.push_back( dParam );
				vecpClusters.push_back( pCluster = new CGenes( Genome ) );
				pCluster->Open( vecstrGenes ); }
			iCluster = -1;
			vecstrGenes.clear( );
			dParam = (float)atof( vecstrLine[ 1 ].c_str( ) );
			continue; }
		if( ( iCur = atoi( vecstrLine[ 1 ].c_str( ) ) ) != iCluster ) {
			if( vecstrGenes.size( ) ) {
				vecdWeights.push_back( dParam );
				vecpClusters.push_back( pCluster = new CGenes( Genome ) );
				pCluster->Open( vecstrGenes ); }
			iCluster = iCur;
			vecstrGenes.clear( ); }
		vecstrGenes.push_back( vecstrLine[ 0 ] ); }
	if( vecstrGenes.size( ) ) {
		vecdWeights.push_back( dParam );
		vecpClusters.push_back( pCluster = new CGenes( Genome ) );
		pCluster->Open( vecstrGenes ); }

	return 0; }

int open_fuzzy( istream& istm, size_t iSkip, CDat& Dat ) {
	CPCL	PCL;
	size_t	i, j, k;
	float	dCur, dMax;

	if( !PCL.Open( istm, iSkip ) ) {
		cerr << "Could not open fuzzy PCL" << endl;
		return 1; }

	Dat.Open( PCL.GetGeneNames( ) );
	for( i = 0; i < Dat.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j ) {
			dMax = 0;
			for( k = 0; k < PCL.GetExperiments( ); ++k )
				if( ( dCur = PCL.Get( i, k ) * PCL.Get( j, k ) ) > dMax )
					dMax = dCur;
			Dat.Set( i, j, dMax ); }

	return 0; }

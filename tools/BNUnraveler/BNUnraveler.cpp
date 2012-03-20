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

static const char	c_acDab[]	= ".dab";

struct SEvaluate {
	CBayesNetSmile*			m_pBN;
	const CDataPair*		m_pDat;
	const CGenes*			m_pGenes;
	const CGenes*			m_pGenesIn;
	CDat*					m_pYes;
	CDat*					m_pNo;
	const CDataPair*		m_pAnswers;
	size_t					m_iZero;
	size_t					m_iNode;
	const vector<size_t>*	m_pveciGenes;
	bool					m_fEverything;
	const char*				m_szName;
};

void* initialize( void* );
void* evaluate( void* );
void* finalize( void* );

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info		sArgs;
	map<string,size_t>		mapZeros;
	CDataPair				Answers;
	vector<CBayesNetSmile*>	vecpBNs;
	vector<size_t>			veciZeros, veciGenes;
	CGenome					Genome;
	CGenes					GenesIn( Genome );
	size_t					i, iNode, iThread, iTerm;
	vector<CGenes*>			vecpGenes;
	CDatasetCompact			Dataset;
	vector<pthread_t>		vecpthdThreads;
	vector<string>			vecstrNodes, vecstrTmps;
	string					strFile;
	vector<CDat*>			vecpYes, vecpNo;
	vector<SEvaluate>		vecsData;

#ifdef WIN32
	pthread_win32_process_attach_np( );
#endif // WIN32

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) || !sArgs.inputs_num ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg );

	vecpBNs.resize( sArgs.inputs_num );
	for( i = 0; i < vecpBNs.size( ); ++i ) {
		vecpBNs[ i ] = new CBayesNetSmile( !!sArgs.group_flag );
		if( !vecpBNs[ i ]->Open( ( strFile = (string)sArgs.input_arg + '/' + CMeta::Deextension(
			CMeta::Basename( sArgs.inputs[ i ] ) ) + ( sArgs.xdsl_flag ? ".xdsl" : ".dsl" ) ).c_str( ) ) ) {
			cerr << "Couldn't open: " << strFile << endl;
			return 1; } }

	if( sArgs.zeros_arg ) {
		ifstream		ifsm;
		vector<string>	vecstrZeros;
		char			acLine[ 1024 ];

		ifsm.open( sArgs.zeros_arg );
		if( !ifsm.is_open( ) ) {
			cerr << "Couldn't open: " << sArgs.zeros_arg << endl;
			return 1; }
		while( !ifsm.eof( ) ) {
			ifsm.getline( acLine, ARRAYSIZE(acLine) - 1 );
			acLine[ ARRAYSIZE(acLine) - 1 ] = 0;
			vecstrZeros.clear( );
			CMeta::Tokenize( acLine, vecstrZeros );
			if( vecstrZeros.empty( ) )
				continue;
			mapZeros[ vecstrZeros[ 0 ] ] = atoi( vecstrZeros[ 1 ].c_str( ) ); } }
	vecpBNs[ 0 ]->GetNodes( vecstrNodes );
	veciZeros.resize( vecstrNodes.size( ) );
	for( i = 0; i < veciZeros.size( ); ++i ) {
		map<string,size_t>::const_iterator	iterZero;

		veciZeros[ i ] = ( ( iterZero = mapZeros.find( vecstrNodes[ i ] ) ) == mapZeros.end( ) ) ? -1 :
			iterZero->second; }

	if( sArgs.answers_arg && !Answers.Open( sArgs.answers_arg, false, !!sArgs.memmap_flag ) ) {
		cerr << "Couldn't open: " << sArgs.answers_arg << endl;
		return 1; }
	if( sArgs.genome_arg ) {
		ifstream	ifsm;
		CGenes		Genes( Genome );

		ifsm.open( sArgs.genome_arg );
		if( !Genes.Open( ifsm ) ) {
			cerr << "Couldn't open: " << sArgs.genome_arg << endl;
			return 1; } }
	else if( sArgs.answers_arg )
		for( i = 0; i < Answers.GetGenes( ); ++i )
			Genome.AddGene( Answers.GetGene( i ) );
	else {
		vector<string>	vecstrFiles;

		vecstrFiles.resize( vecstrNodes.size( ) - 1 );
		for( i = 0; i < vecstrFiles.size( ); ++i )
			vecstrFiles[ i ] = (string)sArgs.directory_arg + '/' + vecstrNodes[ i + 1 ] + c_acDab;
		if( !Dataset.OpenGenes( vecstrFiles ) ) {
			cerr << "Couldn't open: " << sArgs.directory_arg << endl;
			return 1; }
		for( i = 0; i < Dataset.GetGenes( ); ++i )
			Genome.AddGene( Dataset.GetGene( i ) ); }
	if( sArgs.genes_arg ) {
		ifstream	ifsm;

		ifsm.open( sArgs.genes_arg );
		if( !GenesIn.Open( ifsm, !sArgs.genome_arg ) ) {
			cerr << "Couldn't open: " << sArgs.genes_arg << endl;
			return 1; } }

	vecpGenes.resize( sArgs.inputs_num );
	vecpYes.resize( vecpGenes.size( ) );
	vecpNo.resize( vecpGenes.size( ) );
	vecstrTmps.resize( vecpNo.size( ) );
	for( i = 0; i < vecpGenes.size( ); ++i ) {
		ifstream	ifsm;
		char		acTemp[ L_tmpnam + 1 ];

		vecpGenes[ i ]  = new CGenes( Genome );
		ifsm.open( sArgs.inputs[ i ] );
		if( !vecpGenes[ i ]->Open( ifsm, false ) ) {
			cerr << "Couldn't open: " << sArgs.inputs[ i ] << endl;
			return 1; }
		vecpYes[ i ] = new CDat( );
		vecpYes[ i ]->Open( Genome.GetGeneNames( ), false, ( (string)sArgs.output_arg + '/' +
			CMeta::Basename( sArgs.inputs[ i ] ) + c_acDab ).c_str( ) );
		vecpNo[ i ] = new CDat( );
#pragma warning( disable : 4996 )
		vecstrTmps[ i ] = tmpnam( acTemp );
#pragma warning( default : 4996 )
		vecpNo[ i ]->Open( Genome.GetGeneNames( ), false, vecstrTmps[ i ].c_str( ) ); }

	vecsData.resize( vecpGenes.size( ) );
	vecpthdThreads.resize( vecsData.size( ) );
	for( iTerm = 0; iTerm < vecpGenes.size( ); iTerm += iThread ) {
		for( iThread = 0; ( ( sArgs.threads_arg == -1 ) || ( iThread < (size_t)sArgs.threads_arg ) ) &&
			( ( iTerm + iThread ) < vecpGenes.size( ) ); ++iThread ) {
			i = iTerm + iThread;
			vecsData[ i ].m_pBN = vecpBNs[ i ];
			vecsData[ i ].m_pGenes = vecpGenes[ i ];
			vecsData[ i ].m_pYes = vecpYes[ i ];
			vecsData[ i ].m_pNo = vecpNo[ i ];
			vecsData[ i ].m_szName = sArgs.inputs[ i ];
			if( pthread_create( &vecpthdThreads[ i ], NULL, initialize, &vecsData[ i ] ) ) {
				cerr << "Couldn't create initialization thread: " << sArgs.inputs[ i ] << endl;
				return 1; } }
		for( i = 0; i < iThread; ++i )
			pthread_join( vecpthdThreads[ iTerm + i ], NULL ); }

	veciGenes.resize( vecpYes[ 0 ]->GetGenes( ) );
	for( iNode = 1; iNode < vecstrNodes.size( ); ++iNode ) {
		CDataPair	Dat;

		if( !Dat.Open( ( strFile = ( (string)sArgs.directory_arg + '/' + vecstrNodes[ iNode ] +
			c_acDab ) ).c_str( ), false, !!sArgs.memmap_flag ) ) {
			cerr << "Couldn't open: " << strFile << endl;
			return 1; }
		for( i = 0; i < veciGenes.size( ); ++i )
			veciGenes[ i ] = Dat.GetGene( vecpYes[ 0 ]->GetGene( i ) );
		for( iTerm = 0; iTerm < vecpGenes.size( ); iTerm += iThread ) {
			for( iThread = 0; ( ( sArgs.threads_arg == -1 ) || ( iThread < (size_t)sArgs.threads_arg ) ) &&
				( ( iTerm + iThread ) < vecpGenes.size( ) ); ++iThread ) {
				i = iTerm + iThread;
				vecsData[ i ].m_pBN = vecpBNs[ i ];
				vecsData[ i ].m_pDat = &Dat;
				vecsData[ i ].m_pGenes = vecpGenes[ i ];
				vecsData[ i ].m_pGenesIn = GenesIn.GetGenes( ) ? &GenesIn : NULL;
				vecsData[ i ].m_pYes = vecpYes[ i ];
				vecsData[ i ].m_pNo = vecpNo[ i ];
				vecsData[ i ].m_pAnswers = Answers.GetGenes( ) ? &Answers : NULL;
				vecsData[ i ].m_iZero = sArgs.zero_flag ? 0 : veciZeros[ iNode ];
				vecsData[ i ].m_iNode = iNode;
				vecsData[ i ].m_pveciGenes = &veciGenes;
				vecsData[ i ].m_fEverything = !!sArgs.everything_flag;
				vecsData[ i ].m_szName = sArgs.inputs[ i ];
				if( pthread_create( &vecpthdThreads[ i ], NULL, evaluate, &vecsData[ i ] ) ) {
					cerr << "Couldn't create evaluation thread: " << sArgs.inputs[ i ] << endl;
					return 1; } }
			for( i = 0; i < iThread; ++i )
				pthread_join( vecpthdThreads[ iTerm + i ], NULL ); } }

	for( iTerm = 0; iTerm < vecpGenes.size( ); iTerm += iThread ) {
		for( iThread = 0; ( ( sArgs.threads_arg == -1 ) || ( iThread < (size_t)sArgs.threads_arg ) ) &&
			( ( iTerm + iThread ) < vecpGenes.size( ) ); ++iThread ) {
			i = iTerm + iThread;
			vecsData[ i ].m_pBN = vecpBNs[ i ];
			vecsData[ i ].m_pYes = vecpYes[ i ];
			vecsData[ i ].m_pNo = vecpNo[ i ];
			vecsData[ i ].m_szName = sArgs.inputs[ i ];
			if( pthread_create( &vecpthdThreads[ i ], NULL, finalize, &vecsData[ i ] ) ) {
				cerr << "Couldn't create finalization thread: " << sArgs.inputs[ i ] << endl;
				return 1; } }
		for( i = 0; i < iThread; ++i )
			pthread_join( vecpthdThreads[ iTerm + i ], NULL ); }

	for( i = 0; i < vecstrTmps.size( ); ++i )
		_unlink( vecstrTmps[ i ].c_str( ) );
	for( i = 0; i < vecpGenes.size( ); ++i ) {
		delete vecpBNs[ i ];
		delete vecpYes[ i ];
		delete vecpNo[ i ];
		delete vecpGenes[ i ]; }

#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; }

void* initialize( void* pData ) {
	SEvaluate*	psData;
	size_t		i;
	float*		adBuffer;

	psData = (SEvaluate*)pData;
	adBuffer = new float[ psData->m_pYes->GetGenes( ) ];
	for( i = 0; i < psData->m_pNo->GetGenes( ); ++i )
		adBuffer[ i ] = CMeta::GetNaN( );
	for( i = 0; i < psData->m_pNo->GetGenes( ); ++i ) {
		if( !( i % 1000 ) )
			cerr << "IN: " << psData->m_szName << ", " << i << endl;
		psData->m_pNo->Set( i, adBuffer ); }
	for( i = 0; i < psData->m_pYes->GetGenes( ); ++i ) {
		if( !( i % 1000 ) )
			cerr << "IY: " << psData->m_szName << ", " << i << endl;
		psData->m_pYes->Set( i, adBuffer ); }
	delete[] adBuffer;

	return NULL; }

void* evaluate( void* pData ) {
	SEvaluate*		psData;
	CDataMatrix		MatCPT;
	size_t			i, j, iOne, iTwo, iBin, iIndex;
	vector<bool>	vecfGenes, vecfGenesIn;
	bool			fTermOne, fListOne;
	float*			adYes;
	float*			adNo;
	float			dNo, dYes;

	psData = (SEvaluate*)pData;
	psData->m_pBN->GetCPT( 0, MatCPT );
	dNo = log( MatCPT.Get( 0, 0 ) );
	dYes = log( MatCPT.Get( 1, 0 ) );
	adYes = new float[ psData->m_pYes->GetGenes( ) ];
	adNo = new float[ psData->m_pNo->GetGenes( ) ];
	if( !psData->m_fEverything ) {
		vecfGenes.resize( psData->m_pYes->GetGenes( ) );
		for( i = 0; i < vecfGenes.size( ); ++i )
			vecfGenes[ i ] = psData->m_pGenes->IsGene( psData->m_pYes->GetGene( i ) ); }
	if( psData->m_pGenesIn ) {
		vecfGenesIn.resize( psData->m_pYes->GetGenes( ) );
		for( i = 0; i < vecfGenes.size( ); ++i )
			vecfGenesIn[ i ] = psData->m_pGenesIn->IsGene( psData->m_pYes->GetGene( i ) ); }
	psData->m_pBN->GetCPT( psData->m_iNode, MatCPT );
	for( i = 0; i < psData->m_pYes->GetGenes( ); ++i ) {
		if( !( i % 1000 ) )
			cerr << "C: " << psData->m_szName << ", " << i << endl;
		if( ( ( iOne = (*psData->m_pveciGenes)[ i ] ) == -1 ) && ( psData->m_iZero == -1 ) )
			continue;
		fTermOne = psData->m_fEverything || vecfGenes[ i ];
		fListOne = !psData->m_pGenesIn || vecfGenesIn[ i ];
		memcpy( adYes, psData->m_pYes->Get( i ), ( psData->m_pYes->GetGenes( ) - i - 1 ) * sizeof(*adYes) );
		memcpy( adNo, psData->m_pNo->Get( i ), ( psData->m_pNo->GetGenes( ) - i - 1 ) * sizeof(*adNo) );
		for( j = ( i + 1 ); j < psData->m_pYes->GetGenes( ); ++j ) {
			if( ( ( iTwo = (*psData->m_pveciGenes)[ j ] ) == -1 ) && ( psData->m_iZero == -1 ) )
				continue;
			if( !( fTermOne || vecfGenes[ j ] ) || !( fListOne || vecfGenesIn[ j ] ) )
				continue;
			if( psData->m_pAnswers && CMeta::IsNaN( psData->m_pAnswers->Get( i, j ) ) )
				continue;

			iBin = psData->m_pDat->Quantize( iOne, iTwo, psData->m_iZero );
			if( iBin == -1 )
				continue;

			if( CMeta::IsNaN( adYes[ iIndex = ( j - i - 1 ) ] ) ) {
				adYes[ iIndex ] = dYes;
				adNo[ iIndex ] = dNo; }
			adNo[ iIndex ] += log( MatCPT.Get( iBin, 0 ) );
			adYes[ iIndex ] += log( MatCPT.Get( iBin, 1 ) ); }
		psData->m_pNo->Set( i, adNo );
		psData->m_pYes->Set( i, adYes ); }
	delete[] adYes;
	delete[] adNo;

	return NULL; }

void* finalize( void* pData ) {
	SEvaluate*	psData;
	size_t		i, j;
	float*		adYes;
	float*		adNo;
	float		dPrior;
	CDataMatrix	MatCPT;

	psData = (SEvaluate*)pData;
	psData->m_pBN->GetCPT( 0, MatCPT );
	dPrior = MatCPT.Get( 1, 0 );
	adYes = new float[ psData->m_pYes->GetGenes( ) ];
	adNo = new float[ psData->m_pNo->GetGenes( ) ];
	for( i = 0; i < psData->m_pYes->GetGenes( ); ++i ) {
		if( !( i % 1000 ) )
			cerr << "F: " << psData->m_szName << ", " << i << endl;
		memcpy( adYes, psData->m_pYes->Get( i ), ( psData->m_pYes->GetGenes( ) - i - 1 ) * sizeof(*adYes) );
		memcpy( adNo, psData->m_pNo->Get( i ), ( psData->m_pNo->GetGenes( ) - i - 1 ) * sizeof(*adNo) );
		for( j = 0; j < ( psData->m_pYes->GetGenes( ) - i - 1 ); ++j )
			adYes[ j ] = CMeta::IsNaN( adYes[ j ] ) ? dPrior :
				(float)( 1 / ( 1 + exp( (double)adNo[ j ] - (double)adYes[ j ] ) ) );
		psData->m_pYes->Set( i, adYes ); }
	delete[] adNo;
	delete[] adYes;

	return NULL; }

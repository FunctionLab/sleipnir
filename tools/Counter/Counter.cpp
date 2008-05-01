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
static const char	c_acQuant[]	= ".quant";

typedef CFullMatrix<size_t>	CCountMatrix;

struct SLearn {
	CCountMatrix*		m_pMatCounts;
	const CGenes*		m_pGenes;
	const CDataPair*	m_pAnswers;
	const CDataPair*	m_pDat;
	size_t				m_iZero;
};

struct SEvaluate {
	const CBayesNetMinimal*	m_pBN;
	const CDataPair*		m_pDat;
	const CGenes*			m_pGenes;
	CDat*					m_pYes;
	CDat*					m_pNo;
	size_t					m_iZero;
	size_t					m_iNode;
	const vector<size_t>*	m_pveciGenes;
	bool					m_fFirst;
	string					m_strName;
};

void* learn( void* );
void* evaluate( void* );
void* finalize( void* );
int main_count( const gengetopt_args_info&, const map<string, size_t>& );
int main_xdsls( const gengetopt_args_info&, const map<string, size_t>& );
int main_inference( const gengetopt_args_info&, const map<string, size_t>&, const map<string, size_t>& );

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	map<string, size_t>	mapstriZeros, mapstriDatasets;
	int					iRet;

#ifdef WIN32
	pthread_win32_process_attach_np( );
#endif // WIN32
	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );
#ifdef SMILEXML_LIB
	EnableXdslFormat( );
#endif

	if( sArgs.zeros_arg ) {
		ifstream		ifsm;
		vector<string>	vecstrLine;
		char			acLine[ 1024 ];

		ifsm.open( sArgs.zeros_arg );
		if( !ifsm.is_open( ) ) {
			cerr << "Couldn't open: " << sArgs.zeros_arg << endl;
			return 1; }
		while( !ifsm.eof( ) ) {
			ifsm.getline( acLine, ARRAYSIZE(acLine) - 1 );
			acLine[ ARRAYSIZE(acLine) - 1 ] = 0;
			vecstrLine.clear( );
			CMeta::Tokenize( acLine, vecstrLine );
			if( vecstrLine.empty( ) )
				continue;
			mapstriZeros[ vecstrLine[ 0 ] ] = atoi( vecstrLine[ 1 ].c_str( ) ); } }

	if( sArgs.datasets_arg ) {
		ifstream		ifsm;
		char			acLine[ 1024 ];
		vector<string>	vecstrLine;

		ifsm.open( sArgs.datasets_arg );
		if( !ifsm.is_open( ) ) {
			cerr << "Could not open: " << sArgs.datasets_arg << endl;
			return 1; }
		while( !ifsm.eof( ) ) {
			ifsm.getline( acLine, ARRAYSIZE(acLine) - 1 );
			acLine[ ARRAYSIZE(acLine) - 1 ] = 0;
			if( !acLine[ 0 ] )
				continue;
			vecstrLine.clear( );
			CMeta::Tokenize( acLine, vecstrLine );
			if( vecstrLine.size( ) != 2 ) {
				cerr << "Illegal datasets line: " << acLine << endl;
				return 1; }
			mapstriDatasets[ vecstrLine[ 1 ] ] = atol( vecstrLine[ 0 ].c_str( ) ) - 1; }
		ifsm.close( ); }

	if( sArgs.answers_arg )
		iRet = main_count( sArgs, mapstriZeros );
	else if( sArgs.counts_arg )
		iRet = main_xdsls( sArgs, mapstriDatasets );
	else if( sArgs.networks_arg )
		iRet = main_inference( sArgs, mapstriZeros, mapstriDatasets );

	pthread_exit( NULL );
#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return iRet; }

int main_count( const gengetopt_args_info& sArgs, const map<string, size_t>& mapstriZeros ) {
	size_t								i, j, k, m, iTerm, iThread;
	vector<vector<CCountMatrix*>* >		vecpvecpMats;
	vector<CCountMatrix*>				vecpMatRoots;
	vector<CGenes*>						vecpGenes;
	CDataPair							Answers, Dat;
	string								strFile;
	vector<pthread_t>					vecpthdThreads;
	vector<SLearn>						vecsData;
	map<string, size_t>::const_iterator	iterZero;
	CGenome								Genome;
	vector<string>						vecstrNames;

	if( !Answers.Open( sArgs.answers_arg, false, !!sArgs.memmap_flag ) ) {
		cerr << "Couldn't open: " << sArgs.answers_arg << endl;
		return 1; }

	vecpGenes.resize( sArgs.inputs_num );
	for( i = 0; i < vecpGenes.size( ); ++i ) {
		ifstream	ifsm;

		vecpGenes[ i ]  = new CGenes( Genome );
		ifsm.open( sArgs.inputs[ i ] );
		if( !vecpGenes[ i ]->Open( ifsm ) ) {
			cerr << "Couldn't open: " << sArgs.inputs[ i ] << endl;
			return 1; } }
	if( !vecpGenes.size( ) ) {
		vecpGenes.insert( vecpGenes.begin( ), new CGenes( Genome ) );
		vecpGenes[ 0 ]->Open( Answers.GetGeneNames( ) ); }
	vecpMatRoots.resize( vecpGenes.size( ) );
	vecpthdThreads.resize( vecpMatRoots.size( ) );
	vecsData.resize( vecpthdThreads.size( ) );
	for( iTerm = 0; iTerm < vecpMatRoots.size( ); iTerm += iThread ) {
		cerr << "Learning root " << iTerm << '/' << vecpMatRoots.size( ) << endl;
		for( iThread = 0; ( ( sArgs.threads_arg == -1 ) || ( iThread < (size_t)sArgs.threads_arg ) ) &&
			( ( iTerm + iThread ) < vecpMatRoots.size( ) ); ++iThread ) {
			i = iTerm + iThread;
			vecsData[ i ].m_pMatCounts = vecpMatRoots[ i ] = new CCountMatrix( );
			vecsData[ i ].m_pDat = NULL;
			vecsData[ i ].m_pGenes = vecpGenes[ i ];
			vecsData[ i ].m_pAnswers = &Answers;
			vecsData[ i ].m_iZero = -1;
			if( pthread_create( &vecpthdThreads[ i ], NULL, learn, &vecsData[ i ] ) ) {
				cerr << "Couldn't create root thread: " << sArgs.inputs[ i ] << endl;
				return 1; } }
		for( i = 0; i < iThread; ++i )
			pthread_join( vecpthdThreads[ iTerm + i ], NULL ); }

#ifdef _MSC_VER
	HANDLE			hSearch;
	WIN32_FIND_DATA	sEntry;
	bool			fContinue;

	for( fContinue = true,hSearch = FindFirstFile( ( (string)sArgs.directory_arg + "/*.quant" ).c_str( ),
		&sEntry ); fContinue && ( hSearch != INVALID_HANDLE_VALUE );
		fContinue = !!FindNextFile( hSearch, &sEntry ) ) {
		strFile = sEntry.cFileName;
#else // _MSC_VER
	DIR*			pDir;
	struct dirent*	psEntry;

	pDir = opendir( sArgs.directory_arg );
	for( psEntry = readdir( pDir ); psEntry; psEntry = readdir( pDir ) ) {
		strFile = psEntry->d_name;
#endif // _MSC_VER
		string					strName;
		vector<CCountMatrix*>*	pvecpMatCounts;

		if( ( strFile.length( ) < ARRAYSIZE(c_acQuant) ) || ( ( i = strFile.rfind( c_acQuant ) ) !=
			( strFile.length( ) - ARRAYSIZE(c_acQuant) + 1 ) ) )
			continue;

		strName = (string)sArgs.directory_arg + "/" + strFile.substr( 0, i ) + c_acDab;
		if( !Dat.Open( strName.c_str( ), false, !!sArgs.memmap_flag ) ) {
			cerr << "Couldn't open: " << strName << endl;
			return 1; }
		cerr << "Processing: " << strName << endl;
		strName = CMeta::Filename( CMeta::Deextension( CMeta::Basename( strName.c_str( ) ) ) );
		vecstrNames.push_back( strName );
		pvecpMatCounts = new vector<CCountMatrix*>( );
		pvecpMatCounts->resize( vecpGenes.size( ) );
		for( iTerm = 0; iTerm < vecpMatRoots.size( ); iTerm += iThread ) {
			cerr << "Learning term " << iTerm << '/' << vecpMatRoots.size( ) << endl;
			for( iThread = 0; ( ( sArgs.threads_arg == -1 ) || ( iThread < (size_t)sArgs.threads_arg ) ) &&
				( ( iTerm + iThread ) < vecpMatRoots.size( ) ); ++iThread ) {
				i = iTerm + iThread;
				vecsData[ i ].m_pMatCounts = (*pvecpMatCounts)[ i ] = new CCountMatrix( );
				vecsData[ i ].m_pDat = &Dat;
				vecsData[ i ].m_pGenes = vecpGenes[ i ];
				vecsData[ i ].m_pAnswers = &Answers;
				vecsData[ i ].m_iZero = ( ( iterZero = mapstriZeros.find( strName ) ) ==
					mapstriZeros.end( ) ) ? -1 : iterZero->second;
				if( pthread_create( &vecpthdThreads[ i ], NULL, learn, &vecsData[ i ] ) ) {
					cerr << "Couldn't create root thread: " << sArgs.inputs[ i ] << endl;
					return 1; } }
			for( i = 0; i < iThread; ++i )
				pthread_join( vecpthdThreads[ iTerm + i ], NULL ); }
		vecpvecpMats.push_back( pvecpMatCounts ); }
#ifdef _MSC_VER
	FindClose( hSearch );
#else // _MSC_VER
	closedir( pDir );
#endif // _MSC_VER

	for( i = 0; i < vecpMatRoots.size( ); ++i ) {
		cout << ( sArgs.inputs ? CMeta::Deextension( CMeta::Basename( sArgs.inputs[ i ] ) ) : "global" ) <<
			'\t' << vecpvecpMats.size( ) << endl;
		for( j = 0; j < vecpMatRoots[ i ]->GetRows( ); ++j )
			cout << ( j ? "\t" : "" ) << vecpMatRoots[ i ]->Get( j, 0 );
		cout << endl;
		for( j = 0; j < vecpvecpMats.size( ); ++j ) {
			cout << vecstrNames[ j ] << endl;
			for( k = 0; k < (*vecpvecpMats[ j ])[ i ]->GetColumns( ); ++k ) {
				for( m = 0; m < (*vecpvecpMats[ j ])[ i ]->GetRows( ); ++m )
					cout << ( m ? "\t" : "" ) << (*vecpvecpMats[ j ])[ i ]->Get( m, k );
				cout << endl; } } }

	for( i = 0; i < vecpvecpMats.size( ); ++i ) {
		for( j = 0; j < vecpvecpMats[ i ]->size( ); ++j )
			delete (*vecpvecpMats[ i ])[ j ];
		delete vecpvecpMats[ i ]; }
	for( i = 0; i < vecpMatRoots.size( ); ++i ) {
		delete vecpMatRoots[ i ];
		delete vecpGenes[ i ]; }

	return 0; }

void* learn( void* pData ) {
	SLearn*			psData;
	size_t			i, j, iAnswer, iVal, iOne, iTwo;
	vector<bool>	vecfGenes;
	vector<size_t>	veciGenes;

	psData = (SLearn*)pData;
	vecfGenes.resize( psData->m_pAnswers->GetGenes( ) );
	for( i = 0; i < vecfGenes.size( ); ++i )
		vecfGenes[ i ] = psData->m_pGenes->IsGene( psData->m_pAnswers->GetGene( i ) );
	if( psData->m_pDat ) {
		psData->m_pMatCounts->Initialize( psData->m_pDat->GetValues( ), psData->m_pAnswers->GetValues( ) );
		veciGenes.resize( psData->m_pAnswers->GetGenes( ) );
		for( i = 0; i < veciGenes.size( ); ++i )
			veciGenes[ i ] = psData->m_pDat->GetGene( psData->m_pAnswers->GetGene( i ) ); }
	else
		psData->m_pMatCounts->Initialize( psData->m_pAnswers->GetValues( ), 1 );
	psData->m_pMatCounts->Clear( );
	for( i = 0; i < psData->m_pAnswers->GetGenes( ); ++i ) {
		if( psData->m_pDat )
			iOne = veciGenes[ i ];
		for( j = ( i + 1 ); j < psData->m_pAnswers->GetGenes( ); ++j ) {
			if( ( ( iAnswer = psData->m_pAnswers->Quantize( psData->m_pAnswers->Get( i, j ) ) ) == -1 ) ||
				!( vecfGenes[ i ] || vecfGenes[ j ] ) || ( ( vecfGenes[ i ] != vecfGenes[ j ] ) && iAnswer ) )
				continue;
			if( psData->m_pDat ) {
				iTwo = veciGenes[ j ];
				iVal = -1;
				if( ( iOne != -1 ) && ( iTwo != -1 ) )
					iVal = psData->m_pDat->Quantize( psData->m_pDat->Get( iOne, iTwo ) );
				if( iVal == -1 )
					iVal = psData->m_iZero;
				if( iVal == -1 )
					continue;
				psData->m_pMatCounts->Get( iVal, iAnswer )++; }
			else
				psData->m_pMatCounts->Get( iAnswer, 0 )++; } }

	pthread_exit( NULL );
	return NULL; }

int main_xdsls( const gengetopt_args_info& sArgs, const map<string, size_t>& mapstriDatasets ) {
	static const size_t	c_iBuffer	= 1024;
	char						szBuffer[ c_iBuffer ];
	string						strFile;
	CBayesNetMinimal			BNDefault;
	CBayesNetMinimal*			pBN;
	vector<CBayesNetMinimal*>	vecpBNs;
	ifstream					ifsm;
	vector<string>				vecstrLine;
	ofstream					ofsm;
	uint32_t					iSize;
	size_t						i;
	vector<float>				vecdAlphas;

	if( mapstriDatasets.empty( ) ) {
		cerr << "No datasets given" << endl;
		return 1; }
	if( sArgs.alphas_arg ) {
		vecdAlphas.resize( mapstriDatasets.size( ) );
		ifsm.clear( );
		ifsm.open( sArgs.alphas_arg );
		if( !ifsm.is_open( ) ) {
			cerr << "Could not open: " << sArgs.alphas_arg << endl;
			return 1; }
		while( !ifsm.eof( ) ) {
			ifsm.getline( szBuffer, c_iBuffer - 1 );
			szBuffer[ c_iBuffer - 1 ] = 0;
			if( !szBuffer[ 0 ] )
				continue;
			vecstrLine.clear( );
			CMeta::Tokenize( szBuffer, vecstrLine );
			if( vecstrLine.size( ) != 2 ) {
				cerr << "Illegal alphas line: " << szBuffer << endl;
				return 1; }
			vecdAlphas[ mapstriDatasets.find( vecstrLine[ 0 ] )->second ] =
				(float)atof( vecstrLine[ 1 ].c_str( ) ); }
		ifsm.close( ); }

	if( !BNDefault.OpenCounts( sArgs.default_arg, mapstriDatasets, vecdAlphas, sArgs.pseudocounts_arg ) ) {
		cerr << "Could not open default counts: " << ( sArgs.default_arg ? sArgs.default_arg : "not given" ) <<
			endl;
		return 1; }

#ifdef _MSC_VER
	HANDLE			hSearch;
	WIN32_FIND_DATA	sEntry;
	bool			fContinue;

	for( fContinue = true,hSearch = FindFirstFile( ( (string)sArgs.counts_arg + "/*.txt" ).c_str( ),
		&sEntry ); fContinue && ( hSearch != INVALID_HANDLE_VALUE );
		fContinue = !!FindNextFile( hSearch, &sEntry ) ) {
		strFile = sEntry.cFileName;
#else // _MSC_VER
	DIR*			pDir;
	struct dirent*	psEntry;

	pDir = opendir( sArgs.counts_arg );
	for( psEntry = readdir( pDir ); psEntry; psEntry = readdir( pDir ) ) {
		strFile = psEntry->d_name;
#endif // _MSC_VER

		strFile = (string)sArgs.counts_arg + '/' + strFile;
		cerr << "Processing: " << strFile << endl;
		pBN = new CBayesNetMinimal( );
		if( !pBN->OpenCounts( strFile.c_str( ), mapstriDatasets, vecdAlphas, sArgs.pseudocounts_arg,
			&BNDefault ) )
			return 1;
		vecpBNs.push_back( pBN ); }
#ifdef _MSC_VER
	FindClose( hSearch );
#else // _MSC_VER
	closedir( pDir );
#endif // _MSC_VER

	cerr << "Created " << vecpBNs.size( ) << " Bayesian classifiers" << endl;

	if( sArgs.output_arg ) {
		if( sArgs.smile_flag ) {
			CBayesNetSmile						BNSmile;
			vector<string>						vecstrNames;
			map<string, size_t>::const_iterator	iterName;

			vecstrNames.resize( mapstriDatasets.size( ) );
			for( iterName = mapstriDatasets.begin( ); iterName != mapstriDatasets.end( ); ++iterName )
				vecstrNames[ iterName->second ] = iterName->first;

			BNSmile.Open( BNDefault, vecstrNames );
			BNSmile.Save( ( (string)sArgs.output_arg + '/' + BNDefault.GetID( ) +
				( sArgs.xdsl_flag ? ".x" : "." ) + "dsl" ).c_str( ) );
			for( i = 0; i < vecpBNs.size( ); ++i ) {
				BNSmile.Open( *vecpBNs[ i ], vecstrNames );
				BNSmile.Save( ( (string)sArgs.output_arg + '/' + vecpBNs[ i ]->GetID( ) +
					( sArgs.xdsl_flag ? ".x" : "." ) + "dsl" ).c_str( ) ); } }
		else {
			ofsm.open( sArgs.output_arg, ios_base::binary );
			BNDefault.Save( ofsm );
			iSize = (uint32_t)vecpBNs.size( );
			ofsm.write( (const char*)&iSize, sizeof(iSize) );
			for( i = 0; i < vecpBNs.size( ); ++i )
				vecpBNs[ i ]->Save( ofsm );
			ofsm.close( ); } }
	for( i = 0; i < vecpBNs.size( ); ++i )
		delete vecpBNs[ i ];

	return 0; }

int main_inference( const gengetopt_args_info& sArgs, const map<string, size_t>& mapstriZeros,
	const map<string, size_t>& mapstriDatasets ) {
	map<string, size_t>::const_iterator	iterDataset;
	vector<size_t>						veciGenes;
	size_t								i, iTerm, iThread;
	vector<CGenes*>						vecpGenes;
	vector<CDat*>						vecpYes, vecpNo;
	vector<string>						vecstrTmps;
	CGenome								Genome;
	vector<SEvaluate>					vecsData;
	CBayesNetMinimal					BNDefault;
	vector<CBayesNetMinimal*>			vecpBNs;
	ifstream							ifsm;
	uint32_t							iSize;
	map<string, size_t>::const_iterator	iterZero;
	vector<pthread_t>					vecpthdThreads;
	bool								fFirst;

	ifsm.open( sArgs.networks_arg, ios_base::binary );
	if( !BNDefault.Open( ifsm ) ) {
		cerr << "Could not open: " << sArgs.networks_arg << endl;
		return 1; }
	ifsm.read( (char*)&iSize, sizeof(iSize) );
	vecpBNs.resize( iSize );
	for( i = 0; i < vecpBNs.size( ); ++i ) {
		vecpBNs[ i ] = new CBayesNetMinimal( );
		if( !vecpBNs[ i ]->Open( ifsm ) ) {
			cerr << "Could not open: " << sArgs.networks_arg << endl;
			return 1; } }
	ifsm.close( );

	if( !sArgs.genes_arg ) {
		cerr << "No genes given" << endl;
		return 1; }
	{
		CPCL	PCLGenes( false );

		if( !PCLGenes.Open( sArgs.genes_arg, 1 ) ) {
			cerr << "Could not open: " << sArgs.genes_arg << endl;
			return 1; }
		for( i = 0; i < PCLGenes.GetGenes( ); ++i )
			Genome.AddGene( PCLGenes.GetFeature( i, 1 ) );
	}

	vecpGenes.resize( sArgs.inputs_num ? sArgs.inputs_num : 1 );
	vecpYes.resize( vecpGenes.size( ) );
	vecpNo.resize( vecpGenes.size( ) );
	vecstrTmps.resize( vecpNo.size( ) );
	for( i = 0; i < vecpGenes.size( ); ++i ) {
		char	acTemp[ L_tmpnam + 1 ];

		vecpGenes[ i ]  = new CGenes( Genome );
		if( sArgs.inputs_num ) {
			ifstream	ifsm;

			ifsm.open( sArgs.inputs[ i ] );
			if( !vecpGenes[ i ]->Open( ifsm, false ) ) {
				cerr << "Couldn't open: " << sArgs.inputs[ i ] << endl;
				return 1; } }
		else
			vecpGenes[ i ]->Open( Genome.GetGeneNames( ), false );
		vecpYes[ i ] = new CDat( );
		vecpYes[ i ]->Open( Genome.GetGeneNames( ), false, ( (string)sArgs.output_arg + '/' +
			( sArgs.inputs_num ? CMeta::Basename( sArgs.inputs[ i ] ) : "global" ) + c_acDab ).c_str( ) );
		vecpNo[ i ] = new CDat( );
#pragma warning( disable : 4996 )
		vecstrTmps[ i ] = tmpnam( acTemp );
#pragma warning( default : 4996 )
		vecpNo[ i ]->Open( Genome.GetGeneNames( ), false, vecstrTmps[ i ].c_str( ) ); }

	veciGenes.resize( vecpYes[ 0 ]->GetGenes( ) );
	vecsData.resize( vecpGenes.size( ) );
	vecpthdThreads.resize( vecsData.size( ) );
	for( fFirst = true,iterDataset = mapstriDatasets.begin( ); iterDataset != mapstriDatasets.end( );
		fFirst = false,++iterDataset ) {
		CDataPair	Dat;
		string		strFile;

		if( !Dat.Open( ( strFile = ( (string)sArgs.directory_arg + '/' + iterDataset->first +
			c_acDab ) ).c_str( ), false, !!sArgs.memmap_flag ) ) {
			cerr << "Couldn't open: " << strFile << endl;
			return 1; }
		cerr << "Processing: " << strFile << endl;
		for( i = 0; i < veciGenes.size( ); ++i )
			veciGenes[ i ] = Dat.GetGene( vecpYes[ 0 ]->GetGene( i ) );
		for( iTerm = 0; iTerm < vecpGenes.size( ); iTerm += iThread ) {
			for( iThread = 0; ( ( sArgs.threads_arg == -1 ) || ( iThread < (size_t)sArgs.threads_arg ) ) &&
				( ( iTerm + iThread ) < vecpGenes.size( ) ); ++iThread ) {
				i = iTerm + iThread;
				vecsData[ i ].m_pBN = vecpBNs[ i ];
				vecsData[ i ].m_pDat = &Dat;
				vecsData[ i ].m_pGenes = vecpGenes[ i ];
				vecsData[ i ].m_pYes = vecpYes[ i ];
				vecsData[ i ].m_pNo = vecpNo[ i ];
				vecsData[ i ].m_iZero = ( ( iterZero = mapstriZeros.find( iterDataset->first ) ) ==
					mapstriZeros.end( ) ) ? -1 : iterZero->second;
				vecsData[ i ].m_iNode = iterDataset->second;
				vecsData[ i ].m_pveciGenes = &veciGenes;
				vecsData[ i ].m_fFirst = fFirst;
				vecsData[ i ].m_strName = sArgs.inputs_num ? sArgs.inputs[ i ] : "global";
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
			vecsData[ i ].m_strName = sArgs.inputs_num ? sArgs.inputs[ i ] : "global";
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

	return 0; }

void* evaluate( void* pData ) {
	SEvaluate*	psData;
	size_t		i, j, iOne, iTwo, iBin, iIndex;
	float*		adYes;
	float*		adNo;
	float		dNo, dYes;

	psData = (SEvaluate*)pData;
	if( psData->m_fFirst ) {
		float*	adBuffer;

		adBuffer = new float[ psData->m_pYes->GetGenes( ) ];
		for( i = 0; i < psData->m_pNo->GetGenes( ); ++i )
			adBuffer[ i ] = CMeta::GetNaN( );
		for( i = 0; i < psData->m_pNo->GetGenes( ); ++i ) {
			if( !( i % 1000 ) )
				cerr << "IN: " << psData->m_strName << ", " << i << endl;
			psData->m_pNo->Set( i, adBuffer ); }
		for( i = 0; i < psData->m_pYes->GetGenes( ); ++i ) {
			if( !( i % 1000 ) )
				cerr << "IY: " << psData->m_strName << ", " << i << endl;
			psData->m_pYes->Set( i, adBuffer ); }
		delete[] adBuffer; }

	dNo = log( psData->m_pBN->GetCPT( 0 ).Get( 0, 0 ) );
	dYes = log( psData->m_pBN->GetCPT( 0 ).Get( 1, 0 ) );
	adYes = new float[ psData->m_pYes->GetGenes( ) ];
	adNo = new float[ psData->m_pNo->GetGenes( ) ];
	for( i = 0; i < psData->m_pYes->GetGenes( ); ++i ) {
		if( !( i % 1000 ) )
			cerr << "C: " << psData->m_strName << ", " << i << endl;
		if( ( ( iOne = (*psData->m_pveciGenes)[ i ] ) == -1 ) && ( psData->m_iZero == -1 ) )
			continue;
		memcpy( adYes, psData->m_pYes->Get( i ), ( psData->m_pYes->GetGenes( ) - i - 1 ) * sizeof(*adYes) );
		memcpy( adNo, psData->m_pNo->Get( i ), ( psData->m_pNo->GetGenes( ) - i - 1 ) * sizeof(*adNo) );
		for( j = ( i + 1 ); j < psData->m_pYes->GetGenes( ); ++j ) {
			if( ( ( iTwo = (*psData->m_pveciGenes)[ j ] ) == -1 ) && ( psData->m_iZero == -1 ) )
				continue;
			if( ( iOne == -1 ) || ( iTwo == -1 ) ||
				( ( iBin = psData->m_pDat->Quantize( psData->m_pDat->Get( iOne, iTwo ) ) ) == -1 ) )
				iBin = psData->m_iZero;
			if( iBin == -1 )
				continue;
			if( CMeta::IsNaN( adYes[ iIndex = ( j - i - 1 ) ] ) ) {
				adYes[ iIndex ] = dYes;
				adNo[ iIndex ] = dNo; }
			adNo[ iIndex ] += log( psData->m_pBN->GetCPT( psData->m_iNode ).Get( 0, iBin ) );
			adYes[ iIndex ] += log( psData->m_pBN->GetCPT( psData->m_iNode ).Get( 1, iBin ) ); }
		psData->m_pNo->Set( i, adNo );
		psData->m_pYes->Set( i, adYes ); }
	delete[] adYes;
	delete[] adNo;

	pthread_exit( NULL );
	return NULL; }

void* finalize( void* pData ) {
	SEvaluate*	psData;
	size_t		i, j;
	float*		adYes;
	float*		adNo;
	float		dPrior;

	psData = (SEvaluate*)pData;
	dPrior = psData->m_pBN->GetCPT( 0 ).Get( 1, 0 );
	adYes = new float[ psData->m_pYes->GetGenes( ) ];
	adNo = new float[ psData->m_pNo->GetGenes( ) ];
	for( i = 0; i < psData->m_pYes->GetGenes( ); ++i ) {
		if( !( i % 1000 ) )
			cerr << "F: " << psData->m_strName << ", " << i << endl;
		memcpy( adYes, psData->m_pYes->Get( i ), ( psData->m_pYes->GetGenes( ) - i - 1 ) * sizeof(*adYes) );
		memcpy( adNo, psData->m_pNo->Get( i ), ( psData->m_pNo->GetGenes( ) - i - 1 ) * sizeof(*adNo) );
		for( j = 0; j < ( psData->m_pYes->GetGenes( ) - i - 1 ); ++j )
			adYes[ j ] = CMeta::IsNaN( adYes[ j ] ) ? dPrior :
				(float)( 1 / ( 1 + exp( (double)adNo[ j ] - (double)adYes[ j ] ) ) );
		psData->m_pYes->Set( i, adYes ); }
	delete[] adNo;
	delete[] adYes;

	pthread_exit( NULL );
	return NULL; }

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

struct SLearn {
	CBayesNetSmile*			m_pBN;
	const IDataset*			m_pData;
	const CGenes*			m_pGenes;
	const CGenes*			m_pNegatives;
	const CDat*				m_pAnswers;
	const CBayesNetSmile*	m_pBNDefault;
	const char*				m_szName;
	const vector<string>*	m_pvecstrNames;
	const vector<size_t>*	m_pveciZeros;
	bool					m_fZero;
};

void* learn( void* );

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info					sArgs;
	size_t								i, j, iTerm, iThread;
	map<string,size_t>					mapZeros;
	CDataPair							Answers;
	vector<CBayesNetSmile*>				vecpBNRoots, vecpBNs;
	vector<vector<CBayesNetSmile*>* >	vecpvecpBNData;
	CBayesNetSmile						BNDefault;
	vector<size_t>						veciZeros;
	CGenome								Genome;
	CDatasetCompact						Data;
	vector<string>						vecstrDummy;
	vector<CGenes*>						vecpGenes;
	string								strFile;
	vector<pthread_t>					vecpthdThreads;
	vector<SLearn>						vecsData;
	CGenes								GenesNeg( Genome );

#ifdef WIN32
	pthread_win32_process_attach_np( );
#endif // WIN32

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta = CMeta( sArgs.verbosity_arg, sArgs.random_arg );
#ifdef SMILEXML_LIB
	EnableXdslFormat( );
#endif

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

	if( sArgs.default_arg && !BNDefault.Open( sArgs.default_arg ) ) {
		cerr << "Couldn't open: " << sArgs.default_arg << endl;
		return 1; }

	if( !Answers.Open( sArgs.answers_arg, false, sArgs.memmap_flag && !sArgs.genex_arg ) ) {
		cerr << "Couldn't open: " << sArgs.answers_arg << endl;
		return 1; }
	if( sArgs.genex_arg && !Answers.FilterGenes( sArgs.genex_arg, CDat::EFilterExclude ) ) {
		cerr << "Couldn't open: " << sArgs.genex_arg << endl;
		return 1; }
	if( sArgs.negatives_arg && !GenesNeg.Open( sArgs.negatives_arg ) ) {
		cerr << "Couldn't open: " << sArgs.negatives_arg << endl;
		return 1; }

	vecpGenes.resize( sArgs.inputs_num );
	for( i = 0; i < vecpGenes.size( ); ++i ) {
		ifstream	ifsm;

		vecpGenes[ i ]  = new CGenes( Genome );
		ifsm.open( sArgs.inputs[ i ] );
		if( !vecpGenes[ i ]->Open( ifsm ) ) {
			cerr << "Couldn't open: " << sArgs.inputs[ i ] << endl;
			return 1; } }
	if( !Data.Open( Answers, vecstrDummy, true ) ) {
		cerr << "Couldn't open answer set" << endl;
		return 1; }
	vecstrDummy.push_back( "FR" );
	vecpBNRoots.resize( vecpGenes.size( ) );
	vecpthdThreads.resize( vecpBNRoots.size( ) );
	vecsData.resize( vecpthdThreads.size( ) );
	for( iTerm = 0; iTerm < vecpBNRoots.size( ); iTerm += iThread ) {
		cerr << "Learning root " << iTerm << '/' << vecpBNRoots.size( ) << endl;
		for( iThread = 0; ( ( sArgs.threads_arg == -1 ) || ( iThread < (size_t)sArgs.threads_arg ) ) &&
			( ( iTerm + iThread ) < vecpBNRoots.size( ) ); ++iThread ) {
			i = iTerm + iThread;
			vecsData[ i ].m_pBN = vecpBNRoots[ i ] = new CBayesNetSmile( );
			vecsData[ i ].m_pData = &Data;
			vecsData[ i ].m_pGenes = vecpGenes[ i ];
			vecsData[ i ].m_pNegatives = &GenesNeg;
			vecsData[ i ].m_pAnswers = &Answers;
			vecsData[ i ].m_pBNDefault = sArgs.default_arg ? &BNDefault : NULL;
			vecsData[ i ].m_szName = sArgs.inputs[ i ];
			vecsData[ i ].m_pvecstrNames = &vecstrDummy;
			vecsData[ i ].m_pveciZeros = &veciZeros;
			vecsData[ i ].m_fZero = false;
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
		vector<string>				vecstrNames;
		vector<CBayesNetSmile*>*	pvecpBNData;

		if( ( strFile.length( ) < ARRAYSIZE(c_acQuant) ) || ( ( i = strFile.rfind( c_acQuant ) ) !=
			( strFile.length( ) - ARRAYSIZE(c_acQuant) + 1 ) ) )
			continue;

		vecstrNames.push_back( (string)sArgs.directory_arg + "/" + strFile.substr( 0, i ) + c_acDab );
		if( !Data.Open( Answers, vecstrNames, true, sArgs.memmap_flag && !sArgs.randomize_flag ) ) {
			cerr << "Couldn't open: " << strFile << endl;
			return 1; }
		if( sArgs.randomize_flag )
			Data.Randomize( );
		vecstrNames.insert( vecstrNames.begin( ), sArgs.answers_arg );
		for( i = 0; i < vecstrNames.size( ); ++i )
			vecstrNames[ i ] = CMeta::Filename( CMeta::Deextension( CMeta::Basename(
				vecstrNames[ i ].c_str( ) ) ) );
		veciZeros.resize( vecstrNames.size( ) );
		for( i = 0; i < veciZeros.size( ); ++i ) {
			map<string,size_t>::const_iterator	iterZero;

			veciZeros[ i ] = ( ( iterZero = mapZeros.find( vecstrNames[ i ] ) ) == mapZeros.end( ) ) ? -1 :
				iterZero->second; }
		pvecpBNData = new vector<CBayesNetSmile*>( );
		pvecpBNData->resize( vecpGenes.size( ) );
		for( iTerm = 0; iTerm < vecpBNRoots.size( ); iTerm += iThread ) {
			cerr << "Learning term " << iTerm << '/' << vecpBNRoots.size( ) << endl;
			for( iThread = 0; ( ( sArgs.threads_arg == -1 ) || ( iThread < (size_t)sArgs.threads_arg ) ) &&
				( ( iTerm + iThread ) < vecpBNRoots.size( ) ); ++iThread ) {
				i = iTerm + iThread;
				vecsData[ i ].m_pBN = (*pvecpBNData)[ i ] = new CBayesNetSmile( );
				vecsData[ i ].m_pData = &Data;
				vecsData[ i ].m_pGenes = vecpGenes[ i ];
				vecsData[ i ].m_pNegatives = &GenesNeg;
				vecsData[ i ].m_pAnswers = &Answers;
				vecsData[ i ].m_pBNDefault = sArgs.default_arg ? &BNDefault : NULL;
				vecsData[ i ].m_szName = sArgs.inputs[ i ];
				vecsData[ i ].m_pvecstrNames = &vecstrNames;
				vecsData[ i ].m_pveciZeros = &veciZeros;
				vecsData[ i ].m_fZero = !!sArgs.zero_flag;
				if( pthread_create( &vecpthdThreads[ i ], NULL, learn, &vecsData[ i ] ) ) {
					cerr << "Couldn't create root thread: " << sArgs.inputs[ i ] << endl;
					return 1; } }
			for( i = 0; i < iThread; ++i )
				pthread_join( vecpthdThreads[ iTerm + i ], NULL ); }
		vecpvecpBNData.push_back( pvecpBNData ); }

#ifdef _MSC_VER
	FindClose( hSearch );
#else // _MSC_VER
	closedir( pDir );
#endif // _MSC_VER

	vecpBNs.resize( vecpvecpBNData.size( ) );
	for( i = 0; i < vecpGenes.size( ); ++i ) {
		CBayesNetSmile			BNOut;

		for( j = 0; j < vecpBNs.size( ); ++j )
			vecpBNs[ j ] = (*vecpvecpBNData[ j ])[ i ];
		if( !BNOut.Open( *vecpBNRoots[ i ], vecpBNs ) ) {
			cerr << "Couldn't merge networks: " << sArgs.inputs[ i ] << endl;
			return 1; }
		BNOut.Save( ( (string)sArgs.output_arg + "/" + CMeta::Deextension( CMeta::Basename(
			sArgs.inputs[ i ] ) ) + "." + ( sArgs.xdsl_flag ? "x" : "" ) + "dsl" ).c_str( ) ); }

	for( i = 0; i < vecpvecpBNData.size( ); ++i ) {
		for( j = 0; j < vecpvecpBNData[ i ]->size( ); ++j )
			delete (*vecpvecpBNData[ i ])[ j ];
		delete vecpvecpBNData[ i ]; }
	for( i = 0; i < vecpBNRoots.size( ); ++i ) {
		delete vecpBNRoots[ i ];
		delete vecpGenes[ i ]; }

	pthread_exit( NULL );
#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; }

void* learn( void* pData ) {
	CDataFilter	DataFilter, DataNegatives;
	SLearn*		psData;

	psData = (SLearn*)pData;
	if( psData->m_pNegatives->GetGenes( ) ) {
		DataNegatives.Attach( psData->m_pData, *psData->m_pNegatives, CDat::EFilterEdge, psData->m_pAnswers );
		DataFilter.Attach( &DataNegatives, *psData->m_pGenes, CDat::EFilterTerm, psData->m_pAnswers ); }
	else
		DataFilter.Attach( psData->m_pData, *psData->m_pGenes, CDat::EFilterTerm, psData->m_pAnswers );
	if( !psData->m_pBN->Open( &DataFilter, *psData->m_pvecstrNames, *psData->m_pveciZeros ) ) {
		cerr << "Couldn't create base network: " << psData->m_szName << endl;
		pthread_exit( NULL );
		return NULL; }
	if( psData->m_pBNDefault )
		psData->m_pBN->SetDefault( *psData->m_pBNDefault );
	if( !psData->m_pBN->Learn( &DataFilter, 1, psData->m_fZero ) )
		cerr << "Couldn't learn base network: " << psData->m_szName << endl;

	pthread_exit( NULL );
	return NULL; }

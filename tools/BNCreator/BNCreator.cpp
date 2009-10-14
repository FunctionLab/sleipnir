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

struct STerm {
	string					m_strInput;
	string					m_strOutput;
	CGenes*					m_pGenes;
	CBayesNetSmile			m_BNRoot;
	vector<CBayesNetSmile*>	m_vecpBNs;

	STerm( size_t iNodes, const string& strInput, const string& strOutput ) : m_strInput(strInput),
		m_strOutput(strOutput), m_pGenes(NULL) {

		m_vecpBNs.resize( iNodes ); }

	~STerm( ) {
		size_t	i;

		if( m_pGenes )
			delete m_pGenes;
		for( i = 0; i < m_vecpBNs.size( ); ++i )
			delete m_vecpBNs[ i ]; }

	bool Open( CGenome& Genome ) {
		ifstream	ifsm;

		m_pGenes = new CGenes( Genome );
		ifsm.open( m_strInput.c_str( ) );
		return m_pGenes->Open( ifsm ); }

	bool LearnRoot( const CDataPair& Answers, const IDataset* pData,
		const CBayesNetSmile* pBNDefault ) {
		CDataFilter		Data;
		vector<string>	vecstrDummy;
		vector<size_t>	veciZeros;

		if( m_pGenes ) {
			Data.Attach( pData, *m_pGenes, CDat::EFilterTerm, &Answers );
			pData = &Data; }
		vecstrDummy.push_back( "FR" );
		if( !m_BNRoot.Open( pData, vecstrDummy, veciZeros ) ) {
			cerr << "Couldn't create base network (" << m_strInput << ')' << endl;
			return false; }
		if( pBNDefault )
			m_BNRoot.SetDefault( *pBNDefault );
		if( !m_BNRoot.Learn( pData, 1 ) ) {
			cerr << "Couldn't learn base network (" << m_strInput << ')' << endl;
			return false; }

		return true; }

	bool LearnNode( size_t iNode, const CDataPair& Answers, const IDataset* pData, const vector<string>& vecstrNames,
		const vector<size_t>& veciZeros, const CBayesNetSmile* pBNDefault, bool fZero ) {
		CDataFilter		Data;
		CBayesNetSmile*	pBN;

		if( m_pGenes ) {
			Data.Attach( pData, *m_pGenes, CDat::EFilterTerm, &Answers );
			pData = &Data; }
		m_vecpBNs[iNode] = pBN = new CBayesNetSmile( );
		if( !pBN->Open( pData, vecstrNames, veciZeros ) ) {
			cerr << "Couldn't create network for (" << m_strInput << "): " << vecstrNames[ 1 ] << endl;
			return false; }
		if( pBNDefault )
			pBN->SetDefault( *pBNDefault );
		if( !pBN->Learn( pData, 1, fZero ) ) {
			cerr << "Couldn't learn network for (" << m_strInput << "): " << vecstrNames[ 1 ] << endl;
			return false; }

		return true; }

	bool Save( ) const {
		CBayesNetSmile	BNOut;

		if( !BNOut.Open( m_BNRoot, m_vecpBNs ) ) {
			cerr << "Couldn't merge networks (" << m_strInput << ')' << endl;
			return false; }
		BNOut.Save( m_strOutput.c_str( ) );
		return true; }
};

struct SLearn {
	size_t						m_iNode;
	string						m_strInput;
	const CDataPair*			m_pAnswers;
	const gengetopt_args_info*	m_psArgs;
	const map<string,size_t>*	m_pmapZeros;
	const CBayesNetSmile*		m_pBNDefault;
	vector<STerm*>*				m_pvecpsOutputs;
};

void* learn( void* );

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	size_t				i;
	map<string,size_t>	mapZeros;
	CBayesNetSmile		BNIn;
	vector<string>		vecstrNames;

#ifdef WIN32
	pthread_win32_process_attach_np( );
#endif // WIN32
	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg );

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

	if( sArgs.input_arg ) {
		vector<string>	vecstrFiles, vecstrGenes;
		CDat			DatYes, DatNo;
		CDatasetCompact	Data;
		vector<size_t>	veciGenes;
		size_t			j, k, iOne, iTwo, iBin, iZero, iIndex;
		double			d;
		CDataMatrix		MatCPT;
		ofstream		ofsm;
		char			acTemp[ L_tmpnam + 1 ];
		const char*		szTemp;
		float*			adYes;
		float*			adNo;
		CGenome			Genome;
		CGenes			GenesIn( Genome ), GenesEx( Genome );

		if( !BNIn.Open( sArgs.input_arg ) ) {
			cerr << "Couldn't open: " << sArgs.input_arg << endl;
			return 1; }

		BNIn.GetNodes( vecstrFiles );
		vecstrFiles.erase( vecstrFiles.begin( ) );
		for( i = 0; i < vecstrFiles.size( ); ++i )
			vecstrFiles[ i ] = (string)sArgs.directory_arg + '/' + vecstrFiles[ i ] + c_acDab;
		if( !Data.OpenGenes( vecstrFiles ) ) {
			cerr << "Couldn't open: " << sArgs.directory_arg << endl;
			return 1; }
		if( sArgs.genes_arg ) {
			ifstream	ifsm;

			ifsm.open( sArgs.genes_arg );
			if( !GenesIn.Open( ifsm ) ) {
				cerr << "Couldn't open: " << sArgs.genes_arg << endl;
				return 1; }
			ifsm.close( );
			for( i = 0; i < GenesIn.GetGenes( ); ++i )
				vecstrGenes.push_back( GenesIn.GetGene( i ).GetName( ) ); }

#pragma warning( disable : 4996 )
		if( !( szTemp = tmpnam( acTemp ) ) ) {
			cerr << "Couldn't create temp file: " << acTemp << endl;
			return 1; }
#pragma warning( default : 4996 )
		DatYes.Open( vecstrGenes.size( ) ? vecstrGenes : Data.GetGeneNames( ), false, sArgs.output_arg );
		DatNo.Open( DatYes.GetGeneNames( ), false ); // , szTemp );
		adYes = new float[ DatYes.GetGenes( ) ];
		adNo = new float[ DatNo.GetGenes( ) ];
		BNIn.GetCPT( 0, MatCPT );
		d = log( MatCPT.Get( 0, 0 ) );
		for( i = 0; i < DatNo.GetGenes( ); ++i )
			adNo[ i ] = (float)d;
		for( i = 0; i < DatNo.GetGenes( ); ++i )
			DatNo.Set( i, adNo );
		d = log( MatCPT.Get( 1, 0 ) );
		for( i = 0; i < DatYes.GetGenes( ); ++i )
			adYes[ i ] = (float)d;
		for( i = 0; i < DatYes.GetGenes( ); ++i )
			DatYes.Set( i, adYes );

		BNIn.GetNodes( vecstrNames );
		veciGenes.resize( DatYes.GetGenes( ) );
		for( i = 0; i < vecstrFiles.size( ); ++i ) {
			vector<string>						vecstrDatum;
			map<string,size_t>::const_iterator	iterZero;

			vecstrDatum.push_back( vecstrFiles[ i ] );
			if( !Data.Open( vecstrDatum, !!sArgs.memmap_flag ) ) {
				cerr << "Couldn't open: " << vecstrFiles[ i ] << endl;
				return 1; }
			iZero = ( ( iterZero = mapZeros.find( vecstrNames[ i + 1 ] ) ) == mapZeros.end( ) ) ? -1 :
					iterZero->second;
			BNIn.GetCPT( i + 1, MatCPT );
			for( j = 0; j < veciGenes.size( ); ++j )
				veciGenes[ j ] = Data.GetGene( DatYes.GetGene( j ) );
			for( j = 0; j < DatYes.GetGenes( ); ++j ) {
				iBin = -1;
				if( ( ( iOne = veciGenes[ j ] ) == -1 ) && ( iZero == -1 ) )
					continue;
				memcpy( adYes, DatYes.Get( j ), ( DatYes.GetGenes( ) - j - 1 ) * sizeof(*adYes) );
				memcpy( adNo, DatNo.Get( j ), ( DatNo.GetGenes( ) - j - 1 ) * sizeof(*adNo) );
				for( k = ( j + 1 ); k < DatYes.GetGenes( ); ++k ) {
					if( ( iOne == -1 ) || ( ( iTwo = veciGenes[ k ] ) == -1 ) ||
						!Data.IsExample( iOne, iTwo ) ||
						( ( iBin = Data.GetDiscrete( iOne, iTwo, 0 ) ) == -1 ) )
						iBin = iZero;
					if( iBin == -1 )
						continue;
					adNo[ iIndex = ( k - j - 1 ) ] += log( MatCPT.Get( iBin, 0 ) );
					adYes[ iIndex ] += log( MatCPT.Get( iBin, 1 ) ); }
				DatYes.Set( j, adYes );
				DatNo.Set( j, adNo ); } }
		for( i = 0; i < DatYes.GetGenes( ); ++i ) {
			memcpy( adYes, DatYes.Get( i ), ( DatYes.GetGenes( ) - i - 1 ) * sizeof(*adYes) );
			memcpy( adNo, DatNo.Get( i ), ( DatNo.GetGenes( ) - i - 1 ) * sizeof(*adNo) );
			for( j = 0; j < ( DatYes.GetGenes( ) - i - 1 ); ++j )
				adYes[ j ] = (float)( 1 / ( 1 + exp( (double)adNo[ j ] - (double)adYes[ j ] ) ) );
			DatYes.Set( i, adYes ); }
		_unlink( szTemp ); }
	else {
		size_t				iArg, iThread;
		CDataPair			Answers;
		CBayesNetSmile		BNDefault;
		CGenome				Genome;
		vector<STerm*>		vecpsOutputs;
		vector<string>		vecstrDummy;
		CDatasetCompact		Data;
		vector<pthread_t>	vecpthdThreads;
		vector<SLearn>		vecsData;

		if( sArgs.default_arg && !BNDefault.Open( sArgs.default_arg ) ) {
			cerr << "Couldn't open: " << sArgs.default_arg << endl;
			return 1; }

		if( !Answers.Open( sArgs.answers_arg, false ) ) {
			cerr << "Couldn't open: " << sArgs.answers_arg << endl;
			return 1; }
		if( sArgs.genes_arg && !Answers.FilterGenes( sArgs.genes_arg, CDat::EFilterInclude ) ) {
			cerr << "Couldn't open: " << sArgs.genes_arg << endl;
			return 1; }
		if( sArgs.genet_arg && !Answers.FilterGenes( sArgs.genet_arg, CDat::EFilterTerm ) ) {
			cerr << "Couldn't open: " << sArgs.genet_arg << endl;
			return 1; }
		if( sArgs.genee_arg && !Answers.FilterGenes( sArgs.genee_arg, CDat::EFilterEdge ) ) {
			cerr << "Couldn't open: " << sArgs.genee_arg << endl;
			return 1; }
		if( sArgs.genex_arg && !Answers.FilterGenes( sArgs.genex_arg, CDat::EFilterExclude ) ) {
			cerr << "Couldn't open: " << sArgs.genex_arg << endl;
			return 1; }

		if( sArgs.terms_arg ) {
			string			strFile;

			FOR_EACH_DIRECTORY_FILE((string)sArgs.terms_arg, strFile)
				if( strFile[ 0 ] == '.' )
					continue;

				vecpsOutputs.push_back( new STerm( sArgs.inputs_num, (string)sArgs.terms_arg + '/' + strFile,
					(string)sArgs.output_arg + '/' + strFile + ".xdsl" ) );
				if( !vecpsOutputs[ vecpsOutputs.size( ) - 1 ]->Open( Genome ) ) {
					cerr << "Could not open: " << strFile << endl;
					return 1; } } }
		else
			vecpsOutputs.push_back( new STerm( sArgs.inputs_num, "", sArgs.output_arg ) );

		if( !Data.Open( Answers, vecstrDummy ) ) {
			cerr << "Couldn't open answer set" << endl;
			return 1; }
		for( i = 0; i < vecpsOutputs.size( ); ++i ) {
			if( !( i % 50 ) )
				cerr << "Term " << i << '/' << vecpsOutputs.size( ) << endl;
			if( !vecpsOutputs[ i ]->LearnRoot( Answers, &Data, sArgs.default_arg ? &BNDefault : NULL ) )
				return 1; }

		vecpthdThreads.resize( sArgs.inputs_num );
		vecsData.resize( vecpthdThreads.size( ) );
		for( iArg = 0; iArg < sArgs.inputs_num; iArg += iThread ) {
			for( iThread = 0; ( ( sArgs.threads_arg == -1 ) || ( iThread < (size_t)sArgs.threads_arg ) ) &&
				( ( iArg + iThread ) < sArgs.inputs_num ); ++iThread ) {
				i = iArg + iThread;
				vecsData[ i ].m_iNode = i;
				vecsData[ i ].m_strInput = sArgs.inputs[ i ];
				vecsData[ i ].m_pAnswers = &Answers;
				vecsData[ i ].m_psArgs = &sArgs;
				vecsData[ i ].m_pmapZeros = &mapZeros;
				vecsData[ i ].m_pBNDefault = &BNDefault;
				vecsData[ i ].m_pvecpsOutputs = &vecpsOutputs;
				if( pthread_create( &vecpthdThreads[ i ], NULL, learn, &vecsData[ i ] ) ) {
					cerr << "Couldn't create thread: " << sArgs.inputs[ i ] << endl;
					return 1; } }
			for( i = 0; i < iThread; ++i )
				pthread_join( vecpthdThreads[ iArg + i ], NULL ); }

		for( i = 0; i < vecpsOutputs.size( ); ++i ) {
			vecpsOutputs[ i ]->Save( );
			delete vecpsOutputs[ i ]; } }

#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; }

void* learn( void* pData ) {
	CDatasetCompact	Data;
	SLearn*			psData;
	vector<string>	vecstrNames;
	vector<size_t>	veciZeros;
	size_t			i;

	psData = (SLearn*)pData;

	vecstrNames.push_back( psData->m_strInput );
	if( !Data.Open( *psData->m_pAnswers, vecstrNames, psData->m_psArgs->zero_flag || psData->m_psArgs->zeros_arg,
		!!psData->m_psArgs->memmap_flag, psData->m_psArgs->skip_arg, !!psData->m_psArgs->zscore_flag ) ) {
		cerr << "Couldn't open: " << psData->m_strInput << endl;
		return NULL; }
	vecstrNames.insert( vecstrNames.begin( ), psData->m_psArgs->answers_arg );
	for( i = 0; i < vecstrNames.size( ); ++i )
		vecstrNames[ i ] = CMeta::Filename( CMeta::Deextension( CMeta::Basename(
			vecstrNames[ i ].c_str( ) ) ) );
	veciZeros.resize( vecstrNames.size( ) );
	for( i = 0; i < veciZeros.size( ); ++i ) {
		map<string,size_t>::const_iterator	iterZero;

		veciZeros[ i ] = ( ( iterZero = psData->m_pmapZeros->find( vecstrNames[ i ] ) ) ==
			psData->m_pmapZeros->end( ) ) ? -1 : iterZero->second; }
	for( i = 0; i < psData->m_pvecpsOutputs->size( ); ++i ) {
		if( !( i % 50 ) )
			cerr << "Term " << i << '/' << psData->m_pvecpsOutputs->size( ) << endl;
		if( !(*psData->m_pvecpsOutputs)[ i ]->LearnNode( psData->m_iNode, *psData->m_pAnswers, &Data, vecstrNames, veciZeros,
			psData->m_psArgs->default_arg ? psData->m_pBNDefault : NULL, !!psData->m_psArgs->zero_flag ) )
			return NULL; }

	return NULL; }

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
#include "BNServer.h"
#include "dot.h"

static const char				c_szXDSL[]						= ".xdsl";
static const char				c_szDSL[]						= ".dsl";
static const char				c_szDOT[]						= ".dot";
static const char				c_szSVG[]						= ".svg";
const float						CBNServer::c_dCutoff			= 0.1f;
const float						CBNServer::c_adColorMin[]		= {0, 1, 0};
const float						CBNServer::c_adColorMax[]		= {1, 0, 0};
const CBNServer::TPFNProcessor	CBNServer::c_apfnProcessors[]	=
	{&CBNServer::ProcessInference, &CBNServer::ProcessData, &CBNServer::ProcessGraph,
	&CBNServer::ProcessContexts, &CBNServer::ProcessTermFinder, &CBNServer::ProcessDiseases,
	&CBNServer::ProcessGenes, &CBNServer::ProcessAssociation, &CBNServer::ProcessAssociations};

struct SPixie {
	size_t	m_iNode;
	float	m_dScore;

	SPixie( size_t iNode, float dScore ) : m_iNode(iNode), m_dScore(dScore) { }

	bool operator<( const SPixie& sPixie ) const {

		return ( m_dScore < sPixie.m_dScore ); }
};

int main( int iArgs, char** aszArgs ) {
	static const size_t					c_iBuffer	= 1024;
	gengetopt_args_info					sArgs;
	CServer								Server;
	CBayesNetSmile						BNSmile;
	CBayesNetMinimal					BNDefault;
	vector<CBayesNetMinimal>			vecBNs;
	ifstream							ifsm, ifsmGenes;
	istream*							pistm;
	char								acBuffer[ c_iBuffer ];
	vector<string>						vecstrLine;
	map<size_t, string>					mapistrBNs;
	map<size_t, string>::const_iterator	iterBN;
	size_t								i, j, iMax;
	CDatabase							Database;
	uint32_t							iSize;
	CDataMatrix							MatBackgrounds;
	ofstream							ofsm;
	int									iRet;
	COntologyKEGG						KEGG;
	COntologyGO							GOBP, GOMF, GOCC;
	CGenome								Genome;
	const IOntology*					apOntologies[]	= {&GOBP, &GOMF, &GOCC, &KEGG, NULL };
	vector<vector<size_t> >				vecveciDiseases, vecveciContexts;
	vector<size_t>						veciDiseases, veciContexts;
	set<size_t>							setiContexts;

	iRet = cmdline_parser2( iArgs, aszArgs, &sArgs, 0, 1, 0 );
	if( sArgs.config_arg )
		iRet = cmdline_parser_configfile( sArgs.config_arg, &sArgs, 0, 0, 1 ) && iRet;
	if( iRet ) {
		cmdline_parser_print_help( );
		return iRet; }
	CMeta Meta = CMeta( sArgs.verbosity_arg );
#ifdef SMILEXML_LIB
	EnableXdslFormat( );
#endif

	if( sArgs.go_onto_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.go_onto_arg );
		if( sArgs.go_anno_arg ) {
			ifsmGenes.clear( );
			ifsmGenes.open( sArgs.go_anno_arg ); }
		if( !COntologyGO::Open( ifsm, ifsmGenes, Genome, GOBP, GOMF, GOCC, false, true ) ) {
			cerr << "Could not open: " << sArgs.go_onto_arg << ", " << sArgs.go_anno_arg << endl;
			return 1; }
		ifsm.close( );
		if( sArgs.go_anno_arg )
			ifsmGenes.close( ); }

	if( sArgs.kegg_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.kegg_arg );
		if( !KEGG.Open( ifsm, Genome, sArgs.kegg_org_arg, true ) ) {
			cerr << "Could not open: " << sArgs.kegg_arg << endl;
			return 1; }
		ifsm.close( ); }

	if( !Database.Open( sArgs.database_arg ) ) {
		cerr << "Could not open: " << sArgs.database_arg << endl;
		return 1; }
	if( sArgs.minimal_in_flag ) {
		ifsm.clear( );
		ifsm.open( sArgs.networks_arg, ios_base::binary );
		if( !BNDefault.Open( ifsm ) ) {
			cerr << "Could not read: " << sArgs.networks_arg << endl;
			return 1; }
		ifsm.read( (char*)&iSize, sizeof(iSize) );
		vecBNs.resize( iSize );
		for( i = 0; i < vecBNs.size( ); ++i )
			if( !vecBNs[ i ].Open( ifsm ) ) {
				cerr << "Could not read: " << sArgs.networks_arg << " (" << i << ")" << endl;
				return 1; }
		ifsm.close( );
		cerr << "Loaded minimal networks: " << sArgs.networks_arg << endl; }
	else {
		if( sArgs.default_arg && !( BNSmile.Open( sArgs.default_arg ) && BNDefault.Open( BNSmile ) ) ) {
			cerr << "Could not open: " << sArgs.default_arg << endl;
			return 1; }
		BNDefault.SetID( sArgs.default_arg );
		if( sArgs.input_arg ) {
			ifsm.open( sArgs.input_arg );
			pistm = &ifsm; }
		else
			pistm = &cin;
		iMax = 0;
		while( !pistm->eof( ) ) {
			pistm->getline( acBuffer, c_iBuffer - 1 );
			acBuffer[ c_iBuffer - 1 ] = 0;
			if( !acBuffer[ 0 ] )
				continue;
			vecstrLine.clear( );
			CMeta::Tokenize( acBuffer, vecstrLine );
			if( vecstrLine.size( ) < 2 ) {
				cerr << "Ignoring line: " << acBuffer << endl;
				continue; }
			if( ( i = atoi( vecstrLine[ 0 ].c_str( ) ) ) > iMax )
				iMax = i;
			mapistrBNs[ i ] = vecstrLine[ 1 ]; }
		if( sArgs.input_arg )
			ifsm.close( );
		vecBNs.resize( iMax );
		for( iterBN = mapistrBNs.begin( ); iterBN != mapistrBNs.end( ); ++iterBN ) {
			if( !( BNSmile.Open( ( (string)sArgs.networks_arg + '/' + CMeta::Filename( iterBN->second ) +
				( sArgs.xdsl_flag ? c_szXDSL : c_szDSL ) ).c_str( ) ) &&
				vecBNs[ iterBN->first - 1 ].Open( BNSmile ) ) ) {
				cerr << "Could not open: " << iterBN->second << endl;
				return 1; }
			vecBNs[ iterBN->first - 1 ].SetID( iterBN->second ); } }

	if( sArgs.minimal_out_arg ) {
		ofsm.open( sArgs.minimal_out_arg, ios_base::binary );
		BNDefault.Save( ofsm );
		iSize = (uint32_t)vecBNs.size( );
		ofsm.write( (const char*)&iSize, sizeof(iSize) );
		for( i = 0; i < vecBNs.size( ); ++i )
			vecBNs[ i ].Save( ofsm );
		ofsm.close( ); }

	ifsm.clear( );
	ifsm.open( sArgs.contexts_arg );
	if( !ifsm.is_open( ) ) {
		cerr << "Could not open: " << sArgs.contexts_arg << endl;
		return 1; }
	vecveciContexts.resize( vecBNs.size( ) );
	while( !ifsm.eof( ) ) {
		size_t	iContext, iGene;

		ifsm.getline( acBuffer, c_iBuffer - 1 );
		acBuffer[ c_iBuffer - 1 ] = 0;
		if( !acBuffer[ 0 ] )
			continue;
		vecstrLine.clear( );
		CMeta::Tokenize( acBuffer, vecstrLine );
		if( vecstrLine.size( ) != 2 ) {
			cerr << "Invalid line in " << sArgs.contexts_arg << ": " << acBuffer << endl;
			return 1; }
		iContext = atoi( vecstrLine[ 0 ].c_str( ) ) - 1;
		iGene = atoi( vecstrLine[ 1 ].c_str( ) ) - 1;
		if( iContext >= vecveciContexts.size( ) ) {
			cerr << "Invalid context on line: " << acBuffer << endl;
			return 1; }
		setiContexts.insert( iGene );
		vecveciContexts[ iContext ].push_back( iGene ); }
	ifsm.close( );
	veciContexts.resize( setiContexts.size( ) );
	copy( setiContexts.begin( ), setiContexts.end( ), veciContexts.begin( ) );

	if( sArgs.diseases_arg ) {
		set<size_t>	setiDiseases;

		ifsm.clear( );
		ifsm.open( sArgs.diseases_arg );
		if( !ifsm.is_open( ) ) {
			cerr << "Could not open: " << sArgs.diseases_arg << endl;
			return 1; }
		while( !ifsm.eof( ) ) {
			size_t	iDisease, iGene;

			ifsm.getline( acBuffer, c_iBuffer - 1 );
			acBuffer[ c_iBuffer - 1 ] = 0;
			if( !acBuffer[ 0 ] )
				continue;
			vecstrLine.clear( );
			CMeta::Tokenize( acBuffer, vecstrLine );
			if( vecstrLine.size( ) != 2 ) {
				cerr << "Invalid line in " << sArgs.diseases_arg << ": " << acBuffer << endl;
				return 1; }
			iDisease = atoi( vecstrLine[ 0 ].c_str( ) ) - 1;
			iGene = atoi( vecstrLine[ 1 ].c_str( ) ) - 1;
			if( vecveciDiseases.size( ) <= iDisease )
				vecveciDiseases.resize( iDisease + 1 );
			setiDiseases.insert( iGene );
			vecveciDiseases[ iDisease ].push_back( iGene ); }
		ifsm.close( );
		veciDiseases.resize( setiDiseases.size( ) );
		copy( setiDiseases.begin( ), setiDiseases.end( ), veciDiseases.begin( ) ); }

	if( sArgs.backgrounds_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.backgrounds_arg, ios_base::binary );
		if( !MatBackgrounds.Open( ifsm, true ) ) {
			float*	adValues;

			cerr << "Calculating gene backgrounds" << endl;
			MatBackgrounds.Initialize( vecveciContexts.size( ) + 1, Database.GetGenes( ) );
			adValues = new float[ MatBackgrounds.GetColumns( ) ];
			for( i = 0; i < MatBackgrounds.GetColumns( ); ++i ) {
				if( !( i % 100 ) )
					cerr << "Gene " << i << '/' << MatBackgrounds.GetColumns( ) << endl;
				CBNServer::Get( i, 0, adValues, Database, vecBNs, BNDefault );
				MatBackgrounds.Set( 0, i, (float)CStatistics::Average( adValues, adValues +
					MatBackgrounds.GetColumns( ) ) ); }
			for( i = 0; i < vecveciContexts.size( ); ++i ) {
				cerr << "Context " << i << '/' << vecveciContexts.size( ) << endl;
				for( j = 0; j < MatBackgrounds.GetColumns( ); ++j ) {
					CBNServer::Get( j, i + 1, adValues, Database, vecBNs, BNDefault );
					MatBackgrounds.Set( i + 1, j, (float)CStatistics::Average( adValues, adValues +
						MatBackgrounds.GetColumns( ) ) ); } }
			delete[] adValues;

			ofsm.open( sArgs.backgrounds_arg, ios_base::binary );
			MatBackgrounds.Save( ofsm, true );
			ofsm.close( ); } }

	CBNServer	BNServer( BNDefault, vecBNs, vecveciContexts, 0, Database, "", sArgs.files_arg,
		sArgs.graphviz_arg, MatBackgrounds, Genome, apOntologies, vecveciDiseases, veciDiseases,
		veciContexts, sArgs.limit_arg );
	if( sArgs.networklets_flag && !BNServer.GenerateNetworkIcons( ) )
		return -1;

	Server.Initialize( sArgs.port_arg, sArgs.timeout_arg, &BNServer );
#ifdef WIN32
	pthread_win32_process_attach_np( );
#endif // WIN32
	Server.Start( );
#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32

	return 0; }

bool CBNServer::Get( size_t iGene, size_t iContext, float* adValues, const CDatabase& Database,
	const vector<CBayesNetMinimal>& vecBNs, const CBayesNetMinimal& BNDefault ) {
	const CBayesNetMinimal&	BNet		= ( iContext && vecBNs.size( ) ) ?
											vecBNs[ ( iContext - 1 ) % vecBNs.size( ) ] : BNDefault;
	vector<unsigned char>	vecbData;

	if( ( iGene < 1 ) || !Database.Get( iGene - 1, vecbData ) )
		return false;
	if( !BNet.Evaluate( vecbData, adValues, Database.GetGenes( ) ) )
		return false;
	adValues[ ( iGene - 1 ) % Database.GetGenes( ) ] = CMeta::GetNaN( );

	return true; }

CBNServer::CBNServer( const CBayesNetMinimal& BNDefault, const vector<CBayesNetMinimal>& vecBNs,
	const vector<vector<size_t> >& vecveciContexts, SOCKET iSocket, const CDatabase& Database,
	const string& strConnection, const char* szFiles, const char* szGraphviz,
	const CDataMatrix& MatBackgrounds, const CGenome& Genome, const IOntology** apOntologies,
	const vector<vector<size_t> >& vecveciDiseases, const vector<size_t>& veciDiseases,
	const vector<size_t>& veciContexts, size_t iLimit ) : m_BNDefault(BNDefault), m_vecBNs(vecBNs),
	m_vecveciContexts(vecveciContexts), m_iSocket(iSocket), m_Database(Database),
	m_strConnection(strConnection), m_adGenes(NULL), m_adContexts(NULL), m_adDiseases(NULL),
	m_strGraphviz(szGraphviz), m_strFiles(szFiles), m_MatBackgrounds(MatBackgrounds), m_Genome(Genome),
	m_apOntologies(apOntologies), m_vecveciDiseases(vecveciDiseases), m_veciDiseases(veciDiseases),
	m_veciContexts(veciContexts), m_iLimit(iLimit) {

	for( m_iOntologies = 0; m_apOntologies && m_apOntologies[ m_iOntologies ]; ++m_iOntologies );
	if( m_strConnection.length( ) > 0 )
		cerr << "New connection from: " << m_strConnection << endl; }

CBNServer::~CBNServer( ) {

	if( m_adGenes )
		delete[] m_adGenes;
	if( m_adContexts )
		delete[] m_adContexts;
	if( m_adDiseases )
		delete[] m_adDiseases; }

IServerClient* CBNServer::NewInstance( SOCKET iSocket, uint32_t iHost, uint16_t sPort ) {
	string	strConnection;
	char	acBuffer[ 16 ];
	in_addr	sAddr;

#pragma warning(disable : 4996)
	sprintf( acBuffer, "%hu", sPort );
#pragma warning(default : 4996)
	sAddr.s_addr = htonl( iHost );
	strConnection = (string)inet_ntoa( sAddr ) + ":" + acBuffer;
	return new CBNServer( m_BNDefault, m_vecBNs, m_vecveciContexts, iSocket, m_Database, strConnection,
		m_strFiles.c_str( ), m_strGraphviz.c_str( ), m_MatBackgrounds, m_Genome, m_apOntologies,
		m_vecveciDiseases, m_veciDiseases, m_veciContexts, m_iLimit ); }

void CBNServer::Destroy( ) {

	cerr << "Disconnected: " << m_strConnection << endl;

	delete this; }

bool CBNServer::ProcessMessage( const vector<unsigned char>& vecbMessage ) {
	size_t	iProcessed, iOffset;

	for( iOffset = 0; iOffset < vecbMessage.size( ); iOffset += ( iProcessed + 1 ) ) {
		if( vecbMessage[ iOffset ] >= ARRAYSIZE(c_apfnProcessors) ) {
			cerr << m_strConnection << " unknown opcode: " << (int)vecbMessage[ iOffset ] << endl;
			return false; }
		if( ( iProcessed = (this->*c_apfnProcessors[ vecbMessage[ iOffset ] ])( vecbMessage,
			iOffset + 1 ) ) == -1 )
			return false; }

	return true; }
	
size_t CBNServer::ProcessInference( const vector<unsigned char>& vecbMessage, size_t iOffset ) {
	size_t		iStart;
	uint32_t	iGene, iContext;

	if( ( iOffset + sizeof(iContext) ) > vecbMessage.size( ) )
		return -1;
	iStart = iOffset;
	iContext = *(uint32_t*)&vecbMessage[ iOffset ];
	for( iOffset += sizeof(iContext); ( iOffset + sizeof(iGene) ) <= vecbMessage.size( );
		iOffset += sizeof(iGene) ) {
		iGene = *(uint32_t*)&vecbMessage[ iOffset ];
		if( !Get( iGene, iContext ) )
			return -1; }

	return ( iOffset - iStart ); }

bool CBNServer::Get( size_t iGene, size_t iContext, float* adValues ) {
	float*		adTarget;
	uint32_t	iSize;

	cerr << m_strConnection << " inferring " << iGene  << " (" << m_Database.GetGene( iGene - 1 ) << ") in " <<
		iContext << endl;
	if( !( adTarget = adValues ) ) {
		InitializeGenes( );
		adTarget = m_adGenes; }
	if( !CBNServer::Get( iGene, iContext, adTarget, m_Database, m_vecBNs, m_BNDefault ) )
		return false;
	if( !adValues ) {
		iSize = (uint32_t)( GetGenes( ) * sizeof(*m_adGenes) );
		send( m_iSocket, (char*)&iSize, sizeof(iSize), 0 );
		send( m_iSocket, (char*)m_adGenes, iSize, 0 ); }

	return true; }

bool CBNServer::Get( size_t iGene, const vector<size_t>& veciGenes, size_t iContext, float* adValues ) {
	const CBayesNetMinimal&	BNet		= ( iContext && m_vecBNs.size( ) ) ?
											m_vecBNs[ ( iContext - 1 ) % m_vecBNs.size( ) ] : m_BNDefault;
	vector<unsigned char>	vecbData;

	return ( ( iGene >= 1 ) && m_Database.Get( iGene - 1, veciGenes, vecbData ) &&
		BNet.Evaluate( vecbData, adValues, veciGenes.size( ) ) ); }

size_t CBNServer::ProcessData( const vector<unsigned char>& vecbMessage, size_t iOffset ) {
	uint32_t				iSize, iOne, iTwo;
	vector<unsigned char>	vecbData, vecbSend;
	size_t					i;
	unsigned char			b, bValue;

	iSize = sizeof(iOne) + sizeof(iTwo);
	if( ( iOffset + iSize ) > vecbMessage.size( ) )
		return -1;
	iOne = *(uint32_t*)&vecbMessage[ iOffset ];
	iTwo = *(uint32_t*)&vecbMessage[ iOffset + sizeof(iOne) ];
	cerr << m_strConnection << " requested " << iOne << ',' << iTwo << endl;
	if( !m_Database.Get( iOne - 1, iTwo - 1, vecbData ) )
		return -1;

	vecbSend.resize( vecbData.size( ) * 2 );
	for( i = 0; i < vecbData.size( ); ++i ) {
		b = vecbData[ i ];
		if( ( bValue = ( b & 0xF ) ) == 0xF )
			bValue = -1;
		vecbSend[ 2 * i ] = bValue;
		if( ( bValue = ( ( b >> 4 ) & 0xF ) ) == 0xF )
			bValue = -1;
		vecbSend[ ( 2 * i ) + 1 ] = bValue; }
	iOne = (uint32_t)vecbSend.size( );
	send( m_iSocket, (char*)&iOne, sizeof(iOne), 0 );
	send( m_iSocket, (char*)&vecbSend[ 0 ], (int)vecbSend.size( ), 0 );

	return iSize; }

size_t CBNServer::ProcessGraph( const vector<unsigned char>& vecbMessage, size_t iOffset ) {
	size_t						iRet, i, iSize;
	float						dEdgeAggressiveness;
	uint32_t					iGene, iContext, iLimit;
	vector<size_t>				veciQuery, veciNeighbors;
	set<size_t>					setiQuery;
	CDat						DatGraph;
	vector<bool>				vecfQuery;
	set<size_t>::const_iterator	iterQuery;
	unsigned char				bFile;

	iSize = sizeof(bFile) + sizeof(iContext) + sizeof(iLimit) + sizeof(dEdgeAggressiveness);
	if( ( iOffset + iSize ) > vecbMessage.size( ) )
		return -1;
	bFile = vecbMessage[ iOffset ];
	iContext = *(uint32_t*)&vecbMessage[ iOffset + sizeof(bFile) ];
	iLimit = *(uint32_t*)&vecbMessage[ iOffset + sizeof(bFile) + sizeof(iContext) ];
	dEdgeAggressiveness = *(float*)&vecbMessage[ iOffset + sizeof(bFile) + sizeof(iContext) + sizeof(iLimit) ];
	for( i = iOffset + iSize; ( i + sizeof(iGene) ) <= vecbMessage.size( ); i += sizeof(iGene) ) {
		iGene = *(uint32_t*)&vecbMessage[ i ];
		setiQuery.insert( iGene ); }
	iRet = i - iOffset;

	veciQuery.reserve( setiQuery.size( ) );
	for( iterQuery = setiQuery.begin( ); iterQuery != setiQuery.end( ); ++iterQuery )
		veciQuery.push_back( *iterQuery );
	if( !( GraphCreate( veciQuery, iContext, iLimit, dEdgeAggressiveness, vecfQuery, veciNeighbors,
		DatGraph ) && GraphWrite( DatGraph, veciQuery, veciNeighbors, vecfQuery, iContext,
		bFile ? EGraphOutputFile : EGraphOutputSocket ) ) )
		return -1;

	return iRet; }

bool CBNServer::SelectNeighborsPixie( const vector<size_t>& veciQuery, const vector<bool>& vecfQuery,
	size_t iContext, size_t iLimit, const CDataMatrix& MatQuery, vector<size_t>& veciNeighbors ) const {
	vector<float>			vecdNeighbors;
	priority_queue<SPixie>	pqueNeighbors;
	float					d, dAve, dStd;
	size_t					i, j, iN;

	cerr << m_strConnection << " PIXIE query " << iContext << ':';
	for( i = 0; i < veciQuery.size( ); ++i )
		cerr << ' ' << veciQuery[ i ];
	cerr << endl;

	vecdNeighbors.resize( m_Database.GetGenes( ) );
	fill( vecdNeighbors.begin( ), vecdNeighbors.end( ), 0.0f );
	dAve = dStd = 0;
	for( iN = i = 0; i < veciQuery.size( ); ++i )
		for( j = 0; j < m_Database.GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = MatQuery.Get( i, j ) ) ) {
				iN++;
				dAve += d;
				dStd += d * d;
				if( !vecfQuery[ j ] )
					vecdNeighbors[ j ] += d; }
	for( i = 0; i < vecdNeighbors.size( ); ++i )
		if( ( d = vecdNeighbors[ i ] ) > 0 )
			pqueNeighbors.push( SPixie( i, d ) );
	while( !pqueNeighbors.empty( ) && ( ( veciQuery.size( ) + veciNeighbors.size( ) ) < iLimit ) ) {
		veciNeighbors.push_back( pqueNeighbors.top( ).m_iNode );
		pqueNeighbors.pop( ); }

	return true; }

bool CBNServer::SelectNeighborsRatio( const vector<size_t>& veciQuery, const vector<bool>& vecfQuery,
	size_t iContext, size_t iLimit, const CDataMatrix& MatQuery, vector<size_t>& veciNeighbors ) const {
	size_t					i, j;
	priority_queue<SPixie>	pqueNeighbors;
	float					d;
	vector<float>			vecdValues;

	cerr << m_strConnection << " RATIO query " << iContext << " (" << ( iContext ?
		m_vecBNs[ iContext - 1 ].GetID( ) : "total" ) << "):";
	for( i = 0; i < veciQuery.size( ); ++i )
		cerr << ' ' << veciQuery[ i ];
	cerr << endl;

	for( i = 0; i < MatQuery.GetColumns( ); ++i ) {
		if( vecfQuery[ i ] )
			continue;
		vecdValues.clear( );
		vecdValues.reserve( MatQuery.GetRows( ) );
		for( j = 0; j < MatQuery.GetRows( ); ++j )
			if( !CMeta::IsNaN( d = MatQuery.Get( j, i ) ) )
				vecdValues.push_back( d );
		if( !vecdValues.empty( ) ) {
			CStatistics::Winsorize( vecdValues );
			pqueNeighbors.push( SPixie( i, (float)CStatistics::Average( vecdValues ) /
				GetBackground( iContext, i ) ) ); } }
	while( !pqueNeighbors.empty( ) && ( ( veciQuery.size( ) + veciNeighbors.size( ) ) < iLimit ) ) {
		veciNeighbors.push_back( pqueNeighbors.top( ).m_iNode );
		pqueNeighbors.pop( ); }

	return true; }

bool CBNServer::GraphCreate( const vector<size_t>& veciQuery, size_t iContext, size_t iLimit,
	float dEdgeAggressiveness, vector<bool>& vecfQuery, vector<size_t>& veciNeighbors, CDat& DatGraph ) const {
	vector<float>			vecdNeighbors;
	float					d, dAve, dStd, dMin, dCutoff;
	size_t					i, j, iN, iMinOne, iMinTwo;
	vector<size_t>			veciDegree;
	bool					fDone;
	CDataMatrix				MatQuery, MatNeighbors;
	vector<string>			vecstrGenes;

	MatQuery.Initialize( veciQuery.size( ), GetGenes( ) );
	for( i = 0; i < veciQuery.size( ); ++i )
		if( !((CBNServer*)this)->Get( veciQuery[ i ], iContext, MatQuery.Get( i ) ) )
			return false;
	vecfQuery.resize( GetGenes( ) );
	fill( vecfQuery.begin( ), vecfQuery.end( ), false );
	for( i = 0; i < veciQuery.size( ); ++i )
		vecfQuery[ veciQuery[ i ] - 1 ] = true;
	veciNeighbors.clear( );
	fDone = m_MatBackgrounds.GetColumns( ) ?
		SelectNeighborsRatio( veciQuery, vecfQuery, iContext, iLimit, MatQuery, veciNeighbors ) :
		SelectNeighborsPixie( veciQuery, vecfQuery, iContext, iLimit, MatQuery, veciNeighbors );
	if( !fDone )
		return false;

	vecstrGenes.resize( veciQuery.size( ) + veciNeighbors.size( ) );
	vecfQuery.resize( vecstrGenes.size( ) );
	fill( vecfQuery.begin( ), vecfQuery.end( ), false );
	for( i = 0; i < veciQuery.size( ); ++i ) {
		vecstrGenes[ i ] = m_Database.GetGene( veciQuery[ i ] - 1 );
		vecfQuery[ i ] = true; }
	for( i = 0; i < veciNeighbors.size( ); ++i )
		vecstrGenes[ veciQuery.size( ) + i ] = m_Database.GetGene( veciNeighbors[ i ] );
	MatNeighbors.Initialize( veciNeighbors.size( ), veciNeighbors.size( ) );
	for( i = 0; i < veciNeighbors.size( ); ++i )
		if( !((CBNServer*)this)->Get( veciNeighbors[ i ] + 1, veciNeighbors, iContext,
			MatNeighbors.Get( i ) ) )
			return false;
	DatGraph.Open( vecstrGenes );
	for( i = 0; i < veciQuery.size( ); ++i ) {
		for( j = ( i + 1 ); j < veciQuery.size( ); ++j )
			DatGraph.Set( i, j, MatQuery.Get( i, ( veciQuery[ j ] - 1 ) % MatQuery.GetColumns( ) ) );
		for( j = 0; j < veciNeighbors.size( ); ++j )
			DatGraph.Set( i, veciQuery.size( ) + j, MatQuery.Get( i, veciNeighbors[ j ] ) ); }
	dAve = dStd = 0;
	for( iN = i = 0; i < veciNeighbors.size( ); ++i ) {
		const float*	adOne	= MatNeighbors.Get( i );
		size_t			iOne	= veciQuery.size( ) + i;

		for( j = ( i + 1 ); j < veciNeighbors.size( ); ++j )
			if( !CMeta::IsNaN( d = adOne[ j ] ) ) {
				iN++;
				dAve += d;
				dStd += d * d;
				DatGraph.Set( iOne, veciQuery.size( ) + j, adOne[ j ] ); } }

	dAve /= iN;
	dStd = sqrt( ( dStd / iN ) - ( dAve * dAve ) );
	dCutoff = (float)( dAve + ( dEdgeAggressiveness * dStd ) );
	veciDegree.resize( DatGraph.GetGenes( ) );
	fill( veciDegree.begin( ), veciDegree.end( ), veciDegree.size( ) - 1 );
	for( fDone = false; !fDone; ) {
		fDone = true;
		dMin = FLT_MAX;
		for( i = 0; i < DatGraph.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < DatGraph.GetGenes( ); ++j )
				if( !CMeta::IsNaN( d = DatGraph.Get( i, j ) ) && ( d < dCutoff ) && ( d < dMin ) &&
					( veciDegree[ i ] > c_iDegree ) && ( veciDegree[ j ] > c_iDegree ) ) {
					fDone = false;
					dMin = d;
					iMinOne = i;
					iMinTwo = j; }
		if( !fDone ) {
			veciDegree[ iMinOne ]--;
			veciDegree[ iMinTwo ]--;
			DatGraph.Set( iMinOne, iMinTwo, CMeta::GetNaN( ) ); } }

	return true; }

bool CBNServer::GraphWrite( const CDat& DatGraph, const vector<size_t>& veciQuery,
	const vector<size_t>& veciNeighbors, const vector<bool>& vecfQuery, size_t iContext,
	EGraphOutput eOutput ) const {
	static const size_t	c_iBuffer	= 1024;
	char		acBuffer[ c_iBuffer ];
	string		strCmd, strDotIn, strDotOut, strSvg;
	ofstream	ofsm;
	CDot		DotOut( DatGraph );
	CGenome		Genome;
	uint32_t	iSize;

	if( eOutput == EGraphOutputNamed )
		strDotIn = m_strFiles + "/" + CMeta::Filename( m_Database.GetGene( veciQuery[ 0 ] - 1 ) );
	else {
		sprintf_s( acBuffer, ( m_strFiles + "/inXXXXXX" ).c_str( ) );
		if( _mktemp_s( acBuffer ) )
			return false;
		strDotIn = acBuffer; }
	strDotIn += c_szDOT;
	ofsm.open( strDotIn.c_str( ) );
	if( !ofsm.is_open( ) )
		return false;
	Genome.Open( DatGraph.GetGeneNames( ) );
	DatGraph.SaveDOT( ofsm, CMeta::GetNaN( ), &Genome, false, false );
	ofsm.close( );
	if( eOutput == EGraphOutputNamed )
		return true;

	sprintf_s( acBuffer, ( m_strFiles + "/outXXXXXX" ).c_str( ) );
	if( _mktemp_s( acBuffer ) )
		return false;
	strDotOut = acBuffer;
	strDotOut += c_szDOT;
	strCmd = m_strGraphviz + " -Tdot -o" + strDotOut + ' ' + strDotIn;
	system( strCmd.c_str( ) );
// THE PAIN, IT BURNS!
// The Boost DOT parser doesn't handle continuation lines.
// sed doesn't handle newlines unless you beat it with a stick.
// Backslashes have to be triple escaped to get from here to sed.
// The combination makes me cry, but it works.
	system( ( "sed -c -i -n '1h;2,$H;${g;s/\\\\\\n//g;p}' " + strDotOut ).c_str( ) );

	if( !DotOut.Open( strDotOut.c_str( ) ) )
		return false;
	switch( eOutput ) {
		case EGraphOutputFile:
			sprintf_s( acBuffer, "svgXXXXXX" );
			if( _mktemp_s( acBuffer ) )
				return false;
			strSvg = m_strFiles + '/' + acBuffer + c_szSVG;
			ofsm.clear( );
			ofsm.open( strSvg.c_str( ) );
			if( !( ofsm.is_open( ) && DotOut.Save( ofsm, vecfQuery, iContext ) ) )
				return false;
			ofsm.close( );

			strSvg = acBuffer;
			strSvg += c_szSVG;
			iSize = (uint32_t)strSvg.length( ) + sizeof(iSize) + ( sizeof(iSize) * ( veciQuery.size( ) +
				veciNeighbors.size( ) ) );
			send( m_iSocket, (const char*)&iSize, sizeof(iSize), 0 );
			SendGenes( veciQuery, veciNeighbors );
			send( m_iSocket, strSvg.c_str( ), (int)strSvg.length( ), 0 );
			break;

		case EGraphOutputSocket:
			stringstream	sssm;

			if( !DotOut.Save( sssm, vecfQuery, iContext ) )
				return false;

			iSize = (uint32_t)( sssm.str( ).length( ) + sizeof(iSize) + ( sizeof(iSize) * ( veciQuery.size( ) +
				veciNeighbors.size( ) ) ) );
			send( m_iSocket, (const char*)&iSize, sizeof(iSize), 0 );
			SendGenes( veciQuery, veciNeighbors );
			send( m_iSocket, sssm.str( ).c_str( ), (int)sssm.str( ).length( ), 0 );
			break; }

	return true; }

bool CBNServer::SendGenes( const vector<size_t>& veciQuery, const vector<size_t>& veciNeighbors ) const {
	size_t		i;
	uint32_t	iSize;
	uint32_t*	aiGenes;

	iSize = (uint32_t)( veciQuery.size( ) + veciNeighbors.size( ) );
	send( m_iSocket, (const char*)&iSize, sizeof(iSize), 0 );
	aiGenes = new uint32_t[ iSize ];
	for( i = 0; i < veciQuery.size( ); ++i )
		aiGenes[ i ] = (uint32_t)veciQuery[ i ];
	for( i = 0; i < veciNeighbors.size( ); ++i )
		aiGenes[ veciQuery.size( ) + i ] = (uint32_t)veciNeighbors[ i ] + 1;
	send( m_iSocket, (const char*)aiGenes, iSize * sizeof(*aiGenes), 0 );
	delete[] aiGenes;

	return true; }

size_t CBNServer::ProcessContexts( const vector<unsigned char>& vecbMessage, size_t iOffset ) {
	uint32_t				iGene;
	size_t					iCount, iPlace, iContext;
	vector<unsigned char>	vecbData;

	iCount = InitializeContexts( );
	iGene = (uint32_t)( ( ( vecbMessage.size( ) - iOffset ) / sizeof(iGene) ) *
		( iCount * sizeof(*m_adContexts) ) );
	send( m_iSocket, (const char*)&iGene, sizeof(iGene), 0 );
	for( iPlace = iOffset; ( iPlace + sizeof(iGene) ) <= vecbMessage.size( ); iPlace += sizeof(iGene) ) {
		iGene = *(uint32_t*)&vecbMessage[ iPlace ];
		if( ( iGene < 1 ) || !m_Database.Get( iGene - 1, vecbData, true ) )
			return -1;
		cerr << m_strConnection << " contexting " << iGene  << " (" << m_Database.GetGene( iGene - 1 ) <<
			")" << endl;
		for( iContext = 0; iContext < m_vecBNs.size( ); ++iContext )
			if( !GetContext( iGene, vecbData, iContext ) )
				return -1;
		send( m_iSocket, (const char*)m_adContexts, (int)( iCount * sizeof(*m_adContexts) ), 0 ); }

	return ( iPlace - iOffset ); }

bool CBNServer::GetContext( size_t iGene, const vector<unsigned char>& vecbData, size_t iContext ) {
	const vector<size_t>&	veciContext	= m_vecveciContexts[ iContext ];
	size_t			i;
	float			d, dFrac;
	vector<float>	vecdValues;

	dFrac = ( veciContext.size( ) > m_iLimit ) ? ( (float)m_iLimit / veciContext.size( ) ) : 1;
	vecdValues.reserve( veciContext.size( ) );
	for( i = 0; i < veciContext.size( ); ++i ) {
		if( ( dFrac < 1 ) && ( ( (float)rand( ) / RAND_MAX ) > dFrac ) )
			continue;
		if( !CMeta::IsNaN( d = m_vecBNs[ iContext ].Evaluate( vecbData, veciContext[ i ] *
				( ( m_Database.GetDatasets( ) + 1 ) / 2 ) ) ) )
			vecdValues.push_back( d ); }
	CStatistics::Winsorize( vecdValues );
	m_adContexts[ iContext ] = (float)CStatistics::Average( vecdValues );
	m_adContexts[ m_vecveciContexts.size( ) + iContext ] = GetBackground( iContext + 1, iGene );

	return true; }

struct SSortFind {
	bool operator()( const STermFound& sOne, const STermFound& sTwo ) {

		return ( sOne.m_dP < sTwo.m_dP ); }
};

size_t CBNServer::ProcessTermFinder( const vector<unsigned char>& vecbMessage, size_t iOffset ) {
	uint32_t			iOntology, iGene, iSize;
	float				d, dP;
	size_t				i, j, iPlace;
	vector<string>		vecstrGenes;
	vector<STermFound>	vecsTerms;
	CGenes				Genes( *(CGenome*)&m_Genome );
	const IOntology*	pOnto;
	vector<size_t>		veciMapping;

	iSize = sizeof(iOntology) + sizeof(dP);
	if( !m_iOntologies || ( ( iOffset + iSize ) > vecbMessage.size( ) ) )
		return -1;
	iOntology = *(uint32_t*)&vecbMessage[ iOffset ];
	dP = *(float*)&vecbMessage[ iOffset + sizeof(iOntology) ];
	for( iPlace = iOffset + iSize; ( iPlace + sizeof(iGene) ) <= vecbMessage.size( );
		iPlace += sizeof(iGene) ) {
		iGene = *(uint32_t*)&vecbMessage[ iPlace ];
		vecstrGenes.push_back( m_Database.GetGene( iGene - 1 ) ); }

	pOnto = m_apOntologies[ iOntology % m_iOntologies ];
	cerr << m_strConnection << " TermFinder " << pOnto->GetID( ) << '@' << dP << ':';
	for( i = 0; i < vecstrGenes.size( ); ++i )
		cerr << ' ' << vecstrGenes[ i ];
	cerr << endl;

	Genes.Open( vecstrGenes, false );
	veciMapping.resize( Genes.GetGenes( ) );
	for( i = 0; i < Genes.GetGenes( ); ++i ) {
		const CGene&	Gene	= Genes.GetGene( i );

		veciMapping[ i ] = m_Database.GetGene( Gene.GetName( ) );
		for( j = 0; ( j < Gene.GetSynonyms( ) ) && ( veciMapping[ i ] == -1 ); ++j )
			veciMapping[ i ] = m_Database.GetGene( Gene.GetSynonym( j ) ); }
			
	pOnto->TermFinder( Genes, vecsTerms );
	sort( vecsTerms.begin( ), vecsTerms.end( ), SSortFind( ) );
	iSize = sizeof(iSize);
	for( i = 0; ( i < vecsTerms.size( ) ) && ( vecsTerms[ i ].m_dP <= dP ); ++i ) {
		iSize += (uint32_t)( sizeof(float) + ( 5 * sizeof(uint32_t) ) + 2 +
			pOnto->GetGloss( vecsTerms[ i ].m_iID ).length( ) +
			pOnto->GetID( vecsTerms[ i ].m_iID ).length( ) );
		for( j = 0; j < Genes.GetGenes( ); ++j )
			if( pOnto->IsAnnotated( vecsTerms[ i ].m_iID, Genes.GetGene( j ) ) && ( veciMapping[ j ] != -1 ) )
				iSize += sizeof(uint32_t); }

	send( m_iSocket, (const char*)&iSize, sizeof(iSize), 0 );
	iSize = (uint32_t)i;
	send( m_iSocket, (const char*)&iSize, sizeof(iSize), 0 );
	for( i = 0; ( i < vecsTerms.size( ) ) && ( vecsTerms[ i ].m_dP <= dP ); ++i ) {
		const string&		strID		= pOnto->GetID( vecsTerms[ i ].m_iID );
		const string&		strGloss	= pOnto->GetGloss( vecsTerms[ i ].m_iID );
		vector<uint32_t>	veciGenes;

		cerr << m_strConnection << " found " << strGloss << " @ " << vecsTerms[ i ].m_dP << endl;
		send( m_iSocket, strID.c_str( ), (int)strID.length( ) + 1, 0 );
		send( m_iSocket, strGloss.c_str( ), (int)strGloss.length( ) + 1, 0 );
		d = (float)vecsTerms[ i ].m_dP;
		send( m_iSocket, (const char*)&d, sizeof(d), 0 );
		iSize = (uint32_t)vecsTerms[ i ].m_iHitsTerm;
		send( m_iSocket, (const char*)&iSize, sizeof(iSize), 0 );
		iSize = (uint32_t)vecsTerms[ i ].m_iSizeTerm;
		send( m_iSocket, (const char*)&iSize, sizeof(iSize), 0 );
		iSize = (uint32_t)vecsTerms[ i ].m_iHitsTotal;
		send( m_iSocket, (const char*)&iSize, sizeof(iSize), 0 );
		iSize = (uint32_t)vecsTerms[ i ].m_iSizeTotal;
		send( m_iSocket, (const char*)&iSize, sizeof(iSize), 0 );

		for( j = 0; j < Genes.GetGenes( ); ++j )
			if( pOnto->IsAnnotated( vecsTerms[ i ].m_iID, Genes.GetGene( j ) ) && ( veciMapping[ j ] != -1 ) )
				veciGenes.push_back( veciMapping[ j ] );
		iSize = (uint32_t)veciGenes.size( );
		send( m_iSocket, (const char*)&iSize, sizeof(iSize), 0 );
		for( j = 0; j < veciGenes.size( ); ++j ) {
			iSize = (uint32_t)veciGenes[ j ] + 1;
			send( m_iSocket, (const char*)&iSize, sizeof(iSize), 0 ); } }

	return ( iPlace - iOffset ); }

size_t CBNServer::ProcessDiseases( const vector<unsigned char>& vecbMessage, size_t iOffset ) {
	uint32_t		iGene, iContext;
	size_t			iCount, iPlace, iDisease;
	vector<float>	vecdValues;

	if( ( iOffset + sizeof(iContext) ) > vecbMessage.size( ) )
		return -1;
	iContext = *(uint32_t*)&vecbMessage[ iOffset ];

	InitializeGenes( );
	iCount = InitializeDiseases( );
	iGene = (uint32_t)( ( ( vecbMessage.size( ) - iOffset - sizeof(iContext) ) / sizeof(iGene) ) *
		( iCount * sizeof(*m_adDiseases) ) );
	send( m_iSocket, (const char*)&iGene, sizeof(iGene), 0 );

	vecdValues.resize( m_veciDiseases.size( ) );
	for( iPlace = iOffset + sizeof(iContext); ( iPlace + sizeof(iGene) ) <= vecbMessage.size( );
		iPlace += sizeof(iGene) ) {
		iGene = *(uint32_t*)&vecbMessage[ iPlace ];
		cerr << m_strConnection << " diseasing " << iGene  << " (" << m_Database.GetGene( iGene - 1 ) <<
			") in " << iContext << endl;
		if( !Get( iGene, m_veciDiseases, iContext, &vecdValues.front( ) ) )
			return -1;
		for( iDisease = 0; iDisease < m_veciDiseases.size( ); ++iDisease )
			m_adGenes[ m_veciDiseases[ iDisease ] ] = vecdValues[ iDisease ];
		for( iDisease = 0; iDisease < m_vecveciDiseases.size( ); ++iDisease )
			if( !GetDisease( iGene, iContext, iDisease ) )
				return -1;
		send( m_iSocket, (const char*)m_adDiseases, (int)( iCount * sizeof(*m_adDiseases) ), 0 ); }

	return ( iPlace - iOffset ); }

bool CBNServer::GetDisease( size_t iGene, size_t iContext, size_t iDisease ) {
	size_t			i, iCur, iSize;
	vector<float>	vecdIn, vecdOut;
	float			d;

	vecdIn.reserve( iSize = m_vecveciDiseases[ iDisease ].size( ) );
	vecdOut.reserve( iSize );
	for( i = 0; i < iSize; ++i ) {
		iCur = m_vecveciDiseases[ iDisease ][ i ];
		if( !CMeta::IsNaN( d = m_adGenes[ iCur ] ) )
			vecdIn.push_back( d );
		vecdOut.push_back( GetBackground( iContext, iCur ) ); }

	CStatistics::Winsorize( vecdIn );
	CStatistics::Winsorize( vecdOut );
	m_adDiseases[ iDisease ] = (float)CStatistics::Average( vecdIn );
	m_adDiseases[ m_vecveciDiseases.size( ) + iDisease ] = (float)CStatistics::Average( vecdOut );

	return true; }

bool CBNServer::GenerateNetworkIcons( ) const {
	static const size_t	c_iSize					= 6;
	static const float	c_dEdgeAggressiveness	= 0;
	vector<size_t>	veciQuery, veciNeighbors;
	vector<bool>	vecfQuery;
	CDat			DatGraph;
	size_t			i;

	veciQuery.resize( 1 );
	for( i = 0; i < m_Database.GetGenes( ); ++i ) {
		veciQuery[ 0 ] = i + 1;
		if( !( GraphCreate( veciQuery, 0, c_iSize, c_dEdgeAggressiveness, vecfQuery, veciNeighbors,
			DatGraph ) && GraphWrite( DatGraph, veciQuery, veciNeighbors, vecfQuery, 0, EGraphOutputNamed ) ) )
			return false; }

	return true; }

size_t CBNServer::ProcessGenes( const vector<unsigned char>& vecbMessage, size_t iOffset ) {
	uint32_t				iContext, iSize, iGene;
	size_t					iPlace;
	vector<size_t>			veciGenes;

	if( ( iOffset + sizeof(iContext) ) > vecbMessage.size( ) )
		return -1;
	iContext = *(uint32_t*)&vecbMessage[ iOffset ];
	for( iPlace = iOffset + sizeof(iContext); ( iPlace + sizeof(iGene) ) <= vecbMessage.size( );
		iPlace += sizeof(iGene) ) {
		iGene = *(uint32_t*)&vecbMessage[ iPlace ];
		if( iGene < 1 )
			return -1;
		veciGenes.push_back( iGene ); }

	if( !GetGenes( veciGenes, iContext ) )
		return -1;
	iSize = (uint32_t)( GetGenes( ) * sizeof(*m_adGenes) );
	send( m_iSocket, (char*)&iSize, sizeof(iSize), 0 );
	send( m_iSocket, (char*)m_adGenes, iSize, 0 );

	return ( iPlace - iOffset ); }

bool CBNServer::GetGenes( const vector<size_t>& veciGenes, size_t iContext ) {
	CDataMatrix		MatGenes;
	size_t			i, j;
	vector<float>	vecdValues;
	float			d, dFrac;

	dFrac = ( veciGenes.size( ) > m_iLimit ) ? ( (float)m_iLimit / veciGenes.size( ) ) : 1;
	MatGenes.Initialize( veciGenes.size( ), GetGenes( ) );
	for( i = 0; i < MatGenes.GetRows( ); ++i ) {
		if( ( dFrac < 1 ) && ( ( (float)rand( ) / RAND_MAX ) > dFrac ) ) {
			fill( MatGenes.Get( i ), MatGenes.Get( i ) + MatGenes.GetColumns( ), CMeta::GetNaN( ) );
			continue; }
		if( !Get( veciGenes[ i ], iContext, MatGenes.Get( i ) ) )
			return false; }

	InitializeGenes( );
	for( i = 0; i < GetGenes( ); ++i ) {
		vecdValues.clear( );
		vecdValues.reserve( MatGenes.GetRows( ) );
		for( j = 0; j < MatGenes.GetRows( ); ++j )
			if( !CMeta::IsNaN( d = MatGenes.Get( j, i ) ) )
				vecdValues.push_back( d );
		CStatistics::Winsorize( vecdValues );
		m_adGenes[ i ] = (float)CStatistics::Average( vecdValues ) / GetBackground( iContext, i ); }

	return true; }

size_t CBNServer::ProcessAssociation( const vector<unsigned char>& vecbMessage, size_t iOffset ) {
	uint32_t				iGene, iContext;
	size_t					i, iPlace;
	vector<unsigned char>	vecbData;
	vector<size_t>			veciFore, veciBack;
	bool					fBack;
	vector<float>			vecdScores;
	float					dIn;

	if( ( iOffset + sizeof(iContext) ) > vecbMessage.size( ) )
		return -1;
	iContext = *(uint32_t*)&vecbMessage[ iOffset ];
	for( fBack = false,iPlace = ( iOffset + sizeof(iContext) );
		( iPlace + sizeof(iGene) ) <= vecbMessage.size( ); iPlace += sizeof(iGene) ) {
		iGene = *(uint32_t*)&vecbMessage[ iPlace ];
		if( iGene )
			( fBack ? veciBack : veciFore ).push_back( iGene - 1 );
		else if( fBack )
			return -1;
		else
			fBack = true; }
	if( veciFore.empty( ) || veciBack.empty( ) )
		return -1;

	cerr << m_strConnection << " association in " << iContext << ' ';
	for( i = 0; i < veciFore.size( ); ++i )
		cerr << ( i ? ' ' : '{' ) << veciFore[ i ];
	cerr << "} with ";
	for( i = 0; i < veciBack.size( ); ++i )
		cerr << ( i ? ' ' : '{' ) << veciBack[ i ];
	cerr << '}' << endl;

	vecdScores.resize( veciFore.size( ) );
	for( i = 0; i < veciFore.size( ); ++i ) {
		if( !( m_Database.Get( veciFore[ i ], vecbData, true ) &&
			GetAssociation( vecbData, veciBack, iContext, &dIn, NULL ) ) )
			return -1;
		vecdScores[ i ] = dIn / GetBackground( iContext, veciFore[ i ] ); }

	iGene = (uint32_t)( vecdScores.size( ) * sizeof(vecdScores.front( )) );
	send( m_iSocket, (char*)&iGene, sizeof(iGene), 0 );
	send( m_iSocket, (char*)&vecdScores.front( ), iGene, 0 );

	return ( iPlace - iOffset ); }

bool CBNServer::GetAssociation( const vector<unsigned char>& vecbData, const vector<size_t>& veciBack,
	size_t iContext, float* pdIn, float* pdOut ) {
	const CBayesNetMinimal&	BNet	= ( iContext && m_vecBNs.size( ) ) ?
										m_vecBNs[ ( iContext - 1 ) % m_vecBNs.size( ) ] : m_BNDefault;
	size_t			i;
	float			d, dFrac;
	vector<float>	vecdIn, vecdOut;

	if( pdIn )
		*pdIn = 0;
	if( pdOut )
		*pdOut = 0;
	dFrac = ( veciBack.size( ) > m_iLimit ) ? ( (float)m_iLimit / veciBack.size( ) ) : 1;
	for( i = 0; i < veciBack.size( ); ++i ) {
		if( ( dFrac < 1 ) && ( ( (float)rand( ) / RAND_MAX ) > dFrac ) )
			continue;
		if( pdIn && !CMeta::IsNaN( d = BNet.Evaluate( vecbData, veciBack[ i ] *
			( ( m_Database.GetDatasets( ) + 1 ) / 2 ) ) ) )
			vecdIn.push_back( d );
		if( pdOut && !CMeta::IsNaN( d = m_MatBackgrounds.Get( iContext, veciBack[ i ] ) ) )
			vecdOut.push_back( d ); }

	if( pdIn && !vecdIn.empty( ) ) {
		CStatistics::Winsorize( vecdIn );
		*pdIn = (float)CStatistics::Average( vecdIn ); }
	if( pdOut && !vecdOut.empty( ) ) {
		CStatistics::Winsorize( vecdOut );
		*pdOut = (float)CStatistics::Average( vecdOut ); }
	return true; }

size_t CBNServer::ProcessAssociations( const vector<unsigned char>& vecbMessage, size_t iOffset ) {
	uint32_t				iGene, iContext;
	size_t					i, j, k, iPlace;
	vector<size_t>			veciFore, veciBack;
	vector<unsigned char>	vecbData;
	unsigned char			bDiseases;
	float					dIn, dOut, dFrac;
	float*					ad;
	vector<float>			vecdValues;
	vector<vector<float> >	vecvecdIn, vecvecdOut;

	if( ( iOffset + sizeof(iContext) + sizeof(bDiseases) ) > vecbMessage.size( ) )
		return -1;
	bDiseases = vecbMessage[ iOffset ];
	iContext = *(uint32_t*)&vecbMessage[ iOffset + sizeof(bDiseases) ];
	for( iPlace = ( iOffset + sizeof(iContext) + sizeof(bDiseases) );
		( iPlace + sizeof(iGene) ) <= vecbMessage.size( ); iPlace += sizeof(iGene) ) {
		iGene = *(uint32_t*)&vecbMessage[ iPlace ];
		veciFore.push_back( iGene ); }

	cerr << m_strConnection << " associationsing ";
	for( i = 0; i < veciFore.size( ); ++i )
		cerr << ( i ? ' ' : '{' ) << veciFore[ i ];
	cerr << "} with " << ( bDiseases ? "diseases" : "contexts" ) << " in " << iContext << endl;

	bDiseases ?
		memset( m_adDiseases, 0, InitializeDiseases( ) * sizeof(*m_adDiseases) ) :
		memset( m_adContexts, 0, InitializeContexts( ) * sizeof(*m_adContexts) );
	vecvecdIn.resize( ( bDiseases ? m_vecveciDiseases : m_vecveciContexts ).size( ) );
	vecvecdOut.resize( vecvecdIn.size( ) );
	if( bDiseases )
		vecdValues.resize( m_veciDiseases.size( ) );
	dFrac = ( veciFore.size( ) > m_iLimit ) ? ( (float)m_iLimit / veciFore.size( ) ) : 1;
	for( i = 0; i < veciFore.size( ); ++i ) {
		if( ( dFrac < 1 ) && ( ( (float)rand( ) / RAND_MAX ) > dFrac ) )
			continue;
		if( bDiseases ) {
			cerr << m_strConnection << " associations " << veciFore[ i ] << " with diseases" << endl;
			if( !Get( veciFore[ i ], m_veciDiseases, iContext, &vecdValues.front( ) ) )
				return -1;
			memset( m_adGenes, 0, InitializeGenes( ) * sizeof(*m_adGenes) );
			for( j = 0; j < vecdValues.size( ); ++j )
				m_adGenes[ m_veciDiseases[ j ] ] = vecdValues[ j ];
			for( j = 0; j < m_vecveciDiseases.size( ); ++j )
				for( k = 0; k < m_vecveciDiseases[ j ].size( ); ++k ) {
					size_t	iCur	= m_vecveciDiseases[ j ][ k ];
					if( !( CMeta::IsNaN( dIn = m_adGenes[ iCur ] ) ||
						CMeta::IsNaN( dOut = m_MatBackgrounds.Get( iContext, iCur ) ) ) ) {
						vecvecdIn[ j ].push_back( dIn );
						vecvecdOut[ j ].push_back( dOut ); } } }
		else {
			cerr << m_strConnection << " associations " << veciFore[ i ] << " with contexts" << endl;
			if( !m_Database.Get( veciFore[ i ] - 1, vecbData, true ) )
				return -1;
			for( j = 0; j < m_vecveciContexts.size( ); ++j ) {
				if( !GetAssociation( vecbData, m_vecveciContexts[ j ],
					( iContext == -1 ) ? ( j + 1 ) : iContext, &dIn, &dOut ) )
					return -1;
				if( !( CMeta::IsNaN( dIn ) || CMeta::IsNaN( dOut ) ) ) {
					vecvecdIn[ j ].push_back( dIn );
					vecvecdOut[ j ].push_back( dOut ); } } } }
	ad = bDiseases ? m_adDiseases : m_adContexts;
	for( i = 0; i < vecvecdIn.size( ); ++i ) {
		if( vecvecdIn[ i ].empty( ) ) {
			ad[ i ] = ad[ vecvecdIn.size( ) + i ] = 0;
			continue; }
		CStatistics::Winsorize( vecvecdIn[ i ] );
		CStatistics::Winsorize( vecvecdOut[ i ] );
		ad[ i ] = (float)CStatistics::Average( vecvecdIn[ i ] );
		ad[ vecvecdIn.size( ) + i ] = (float)CStatistics::Average( vecvecdOut[ i ] ); }

	iGene = (uint32_t)( 2 * vecvecdIn.size( ) * sizeof(*ad) );
	send( m_iSocket, (char*)&iGene, sizeof(iGene), 0 );
	send( m_iSocket, (char*)ad, iGene, 0 );

	return ( iPlace - iOffset ); }

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
static const size_t				c_iContextOffset				= 2;
static const size_t				c_iParamOffset					= 4;
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
	gengetopt_args_info			sArgs;
	CServer						Server;
	CBayesNetSmile				BNSmile;
	ifstream					ifsm, ifsmGenes;
	vector<string>				vecstrLine;
	size_t						i, j, k;
	uint32_t					iSize;
	ofstream					ofsm;
	int							iRet;
	COntologyKEGG				KEGG;
	COntologyGO					GOBP, GOMF, GOCC;
	set<size_t>					setiContexts;
// Server data
	vector<float>				vecdPriors;
	CBayesNetMinimal			BNDefault;
	vector<CBayesNetMinimal>	vecBNs;
	CDatabase					Database;
	CDataMatrix					MatBackgrounds, MatParameters, MatWithinC, MatWithinD;
	CDataMatrix					MatBetweenCC, MatBetweenDD, MatBetweenDC;
	vector<vector<size_t> >		vecveciDiseases, vecveciContexts;
	vector<size_t>				veciDiseases, veciContexts;
	CGenome						Genome;
	const IOntology*			apOntologies[]	= {&GOBP, &GOMF, &GOCC, &KEGG, NULL };

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

// Open Gene Ontology
	if( sArgs.go_onto_arg ) {
		cerr << "Loading the Gene Ontology..." << endl;
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

// Open KEGG
	if( sArgs.kegg_arg ) {
		cerr << "Loading KEGG..." << endl;
		ifsm.clear( );
		ifsm.open( sArgs.kegg_arg );
		if( !KEGG.Open( ifsm, Genome, sArgs.kegg_org_arg, true ) ) {
			cerr << "Could not open: " << sArgs.kegg_arg << endl;
			return 1; }
		ifsm.close( ); }

// Open database
	cerr << "Loading the database..." << endl;
	if( !Database.Open( sArgs.database_arg ) ) {
		cerr << "Could not open: " << sArgs.database_arg << endl;
		return 1; }
// Open Bayes networks
	if( sArgs.minimal_in_flag ) {
		cerr << "Loading minimal Bayesian classifiers..." << endl;
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
		ifsm.close( ); }
	if( !sArgs.minimal_in_flag ) {
		if( !sArgs.default_arg ) {
			cmdline_parser_print_help( );
			return 1; }
		if( !( BNSmile.Open( sArgs.default_arg ) && BNDefault.Open( BNSmile ) ) ) {
			cerr << "Could not open: " << sArgs.default_arg << endl;
			return 1; }
		BNDefault.SetID( sArgs.default_arg ); }
// Open context IDs and p-value parameters
	{
		CPCL	PCLContexts( false );

		cerr << "Loading contexts..." << endl;
		if( !PCLContexts.Open( sArgs.input_arg, c_iParamOffset - 1 ) ) {
			cerr << "Could not open: " << ( sArgs.input_arg ? sArgs.input_arg : "standard input" ) << endl;
			return 1; }
		MatParameters.Initialize( PCLContexts.GetGenes( ) + 1, PCLContexts.GetExperiments( ) );
		for( i = 0; i < PCLContexts.GetGenes( ); ++i )
			MatParameters.Set( i + 1, PCLContexts.Get( i ) );

		if( !sArgs.minimal_in_flag ) {
			vecBNs.resize( PCLContexts.GetGenes( ) );
			for( i = 0; i < vecBNs.size( ); ++i ) {
				const string&	strContext	= PCLContexts.GetFeature( i, c_iContextOffset );

				if( !( BNSmile.Open( ( (string)sArgs.networks_arg + '/' + CMeta::Filename( strContext ) +
					( sArgs.xdsl_flag ? c_szXDSL : c_szDSL ) ).c_str( ) ) &&
					vecBNs[ i ].Open( BNSmile ) ) ) {
					cerr << "Could not open: " << strContext << endl;
					return 1; }
				vecBNs[ i ].SetID( strContext ); } }
	}
// Open global p-value parameters
	{
		CPCL	PCLParams( false );

		cerr << "Loading global p-value parameters..." << endl;
		if( !PCLParams.Open( sArgs.global_arg, 0 ) ) {
			cerr << "Could not open: " << sArgs.global_arg << endl;
			return 1; }
		MatParameters.Set( 0, 0, (float)atof( PCLParams.GetGene( 0 ).c_str( ) ) );
		for( i = 0; i < PCLParams.GetExperiments( ); ++i )
			MatParameters.Set( 0, i + 1, PCLParams.Get( 0, i ) );
	}
// Save minimal Bayes nets
	if( sArgs.minimal_out_arg ) {
		cerr << "Saving minimal Bayesian classifiers..." << endl;
		ofsm.open( sArgs.minimal_out_arg, ios_base::binary );
		BNDefault.Save( ofsm );
		iSize = (uint32_t)vecBNs.size( );
		ofsm.write( (const char*)&iSize, sizeof(iSize) );
		for( i = 0; i < vecBNs.size( ); ++i )
			vecBNs[ i ].Save( ofsm );
		ofsm.close( ); }
// Open context<->gene mapping
	{
		CPCL	PCLContextGenes( false );

		cerr << "Loading context contents..." << endl;
		if( !PCLContextGenes.Open( sArgs.contexts_arg, 1 ) ) {
			cerr << "Could not open: " << sArgs.contexts_arg << endl;
			return 1; }
		vecveciContexts.resize( vecBNs.size( ) );
		for( i = 0; i < PCLContextGenes.GetGenes( ); ++i ) {
			j = atoi( PCLContextGenes.GetGene( i ).c_str( ) ) - 1;
			k = atoi( PCLContextGenes.GetFeature( i, 1 ).c_str( ) ) - 1;
			if( j >= vecveciContexts.size( ) ) {
				cerr << "Unknown context: " << j << " (" << PCLContextGenes.GetGene( i ) << ')' << endl;
				return 1; }
			if( k >= Database.GetGenes( ) ) {
				cerr << "Unknown gene: " << k << " (" << PCLContextGenes.GetFeature( i, 1 ) << ')' << endl;
				return 1; }
			setiContexts.insert( k );
			vecveciContexts[ j ].push_back( k ); }
		veciContexts.resize( setiContexts.size( ) );
		copy( setiContexts.begin( ), setiContexts.end( ), veciContexts.begin( ) );
	}
// Open disease<->gene mapping
	{
		CPCL		PCLDiseaseGenes( false );

		cerr << "Loading disease contents..." << endl;
		if( !PCLDiseaseGenes.Open( sArgs.diseases_arg, 1 ) ) {
			cerr << "Could not open: " << sArgs.diseases_arg << endl;
			return 1; }
		for( i = j = 0; i < PCLDiseaseGenes.GetGenes( ); ++i )
			if( ( k = atoi( PCLDiseaseGenes.GetGene( i ).c_str( ) ) ) > j )
				j = k;
		vecveciDiseases.resize( j );
		setiContexts.clear( );
		for( i = 0; i < PCLDiseaseGenes.GetGenes( ); ++i ) {
			j = atoi( PCLDiseaseGenes.GetGene( i ).c_str( ) ) - 1;
			k = atoi( PCLDiseaseGenes.GetFeature( i, 1 ).c_str( ) ) - 1;
			if( j >= vecveciContexts.size( ) ) {
				cerr << "Unknown disease: " << j << " (" << PCLDiseaseGenes.GetGene( i ) << ')' << endl;
				return 1; }
			if( k >= Database.GetGenes( ) ) {
				cerr << "Unknown gene: " << k << " (" << PCLDiseaseGenes.GetFeature( i, 1 ) << ')' << endl;
				return 1; }
			setiContexts.insert( k );
			vecveciDiseases[ j ].push_back( k ); }
		veciDiseases.resize( setiContexts.size( ) );
		copy( setiContexts.begin( ), setiContexts.end( ), veciDiseases.begin( ) );
	}
// Open gene backgrounds
	cerr << "Loading gene backgrounds..." << endl;
	if( sArgs.backgrounds_arg && !MatBackgrounds.Open( sArgs.backgrounds_arg ) ) {
		cerr << "Could not open: " << sArgs.backgrounds_arg << endl;
		return 1; }
// Initialize context priors
	vecdPriors.resize( MatBackgrounds.GetRows( ) );
	fill( vecdPriors.begin( ), vecdPriors.end( ), 0.0f );
	for( i = 0; i < vecdPriors.size( ); ++i ) {
		for( j = 0; j < MatBackgrounds.GetColumns( ); ++j )
			vecdPriors[ i ] += MatBackgrounds.Get( i, j );
		vecdPriors[ i ] /= MatBackgrounds.GetColumns( ); }
// Open context withins
	if( sArgs.within_c_arg ) {
		cerr << "Loading within context scores..." << endl;
		if( !MatWithinC.Open( sArgs.within_c_arg ) ) {
			cerr << "Could not open: " << sArgs.within_c_arg << endl;
			return 1; } }
// Open disease withins
	if( sArgs.within_d_arg ) {
		cerr << "Loading within disease scores..." << endl;
		if( !MatWithinD.Open( sArgs.within_d_arg ) ) {
			cerr << "Could not open: " << sArgs.within_d_arg << endl;
			return 1; } }
// Open context betweens
	if( sArgs.between_cc_arg ) {
		cerr << "Loading between context scores..." << endl;
		if( !MatBetweenCC.Open( sArgs.between_cc_arg ) ) {
			cerr << "Could not open: " << sArgs.between_cc_arg << endl;
			return 1; } }
// Open disease betweens
	if( sArgs.between_dd_arg ) {
		cerr << "Loading between disease scores..." << endl;
		if( !MatBetweenDD.Open( sArgs.between_dd_arg ) ) {
			cerr << "Could not open: " << sArgs.between_dd_arg << endl;
			return 1; } }
// Open context-disease betweens
	if( sArgs.between_dc_arg ) {
		cerr << "Loading between disease-context scores..." << endl;
		if( !MatBetweenDC.Open( sArgs.between_dc_arg ) ) {
			cerr << "Could not open: " << sArgs.between_dc_arg << endl;
			return 1; } }

	SBNServerData	sData(
		apOntologies,
		BNDefault,
		Database,
		Genome,
		sArgs.limit_arg,
		MatBackgrounds,
		MatBetweenCC,
		MatBetweenDC,
		MatBetweenDD,
		MatParameters,
		MatWithinC,
		MatWithinD,
		sArgs.files_arg,
		sArgs.graphviz_arg,
		vecBNs,
		vecdPriors,
		veciContexts,
		veciDiseases,
		vecveciContexts,
		vecveciDiseases
		);
	CBNServer	BNServer( 0, "", sData );
	if( sArgs.networklets_flag )
		return ( BNServer.GenerateNetworkIcons( ) ? 0 : 1 );
	if( sArgs.assoc_diseases_arg )
		return ( BNServer.GenerateAssociations( sArgs.assoc_diseases_arg, sArgs.assoc_context_arg ) ? 0 : 1 );

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

	if( !( Database.Get( iGene, vecbData ) && BNet.Evaluate( vecbData, adValues, Database.GetGenes( ) ) ) )
		return false;
	adValues[ iGene % Database.GetGenes( ) ] = CMeta::GetNaN( );

	return true; }

CBNServer::CBNServer( SOCKET iSocket, const string& strConnection, const SBNServerData& sData ) :
	m_iSocket(iSocket), m_strConnection(strConnection), m_sData(sData), m_adContexts(NULL),
	m_adDiseases(NULL), m_adGenes(NULL) {

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
	return new CBNServer( iSocket, strConnection, m_sData ); }

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
		iOffset += sizeof(iGene) )
		if( !( ( iGene = *(uint32_t*)&vecbMessage[ iOffset ] ) && Get( iGene - 1, iContext ) ) )
			return -1;

	return ( iOffset - iStart ); }

bool CBNServer::Get( size_t iGene, size_t iContext, float* adValues ) {
	float*		adTarget;
	uint32_t	iSize;

	cerr << m_strConnection << " inferring " << iGene  << " (" << GetGene( iGene ) << ") in " <<
		iContext << endl;
	if( !( adTarget = adValues ) ) {
		InitializeGenes( );
		adTarget = m_adGenes; }
	if( !CBNServer::Get( iGene, iContext, adTarget, GetDatabase( ), m_sData.m_vecBNs, m_sData.m_BNDefault ) )
		return false;
	if( !adValues ) {
		iSize = (uint32_t)( GetGenes( ) * sizeof(*m_adGenes) );
		send( m_iSocket, (char*)&iSize, sizeof(iSize), 0 );
		send( m_iSocket, (char*)m_adGenes, iSize, 0 ); }

	return true; }

bool CBNServer::Get( size_t iGene, const vector<size_t>& veciGenes, size_t iContext, float* adValues ) {
	const CBayesNetMinimal&	BNet	= GetBN( iContext );
	vector<unsigned char>	vecbData;

	return ( GetDatabase( ).Get( iGene, veciGenes, vecbData ) &&
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
	if( !( iOne-- && iTwo-- ) )
		return -1;
	cerr << m_strConnection << " requested " << GetGene( iOne ) << ',' << GetGene( iTwo ) << endl;
	if( !GetDatabase( ).Get( iOne, iTwo, vecbData ) )
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
		if( !iGene-- )
			return -1;
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
		cerr << ' ' << GetGene( veciQuery[ i ] );
	cerr << endl;

	vecdNeighbors.resize( GetGenes( ) );
	fill( vecdNeighbors.begin( ), vecdNeighbors.end( ), 0.0f );
	dAve = dStd = 0;
	for( iN = i = 0; i < veciQuery.size( ); ++i )
		for( j = 0; j < GetGenes( ); ++j )
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
	size_t					i, j, iOut;
	priority_queue<SPixie>	pqueNeighbors;
	float					d, dOut;
	vector<float>			vecdValues;

	cerr << m_strConnection << " RATIO query " << iContext << " (" << ( iContext ?
		GetBN( iContext ).GetID( ) : "total" ) << "):";
	for( i = 0; i < veciQuery.size( ); ++i )
		cerr << ' ' << GetGene( veciQuery[ i ] );
	cerr << endl;

	for( dOut = 0, iOut = i = 0; i < veciQuery.size( ); ++i )
		if( !CMeta::IsNaN( d = GetBackground( iContext, veciQuery[ i ] ) ) ) {
			iOut++;
			dOut += d; }
	for( i = 0; i < MatQuery.GetColumns( ); ++i ) {
		if( vecfQuery[ i ] )
			continue;
		vecdValues.clear( );
		vecdValues.reserve( MatQuery.GetRows( ) );
		for( j = 0; j < MatQuery.GetRows( ); ++j )
			if( !CMeta::IsNaN( d = MatQuery.Get( j, i ) ) )
				vecdValues.push_back( d );
		if( !vecdValues.empty( ) ) {
			Winsorize( vecdValues );
			pqueNeighbors.push( SPixie( i, (float)CStatistics::Average( vecdValues ) * ( iOut + 1 ) /
				( dOut + GetBackground( iContext, i ) ) ) ); } }
	while( !pqueNeighbors.empty( ) && ( ( veciQuery.size( ) + veciNeighbors.size( ) ) < iLimit ) ) {
		veciNeighbors.push_back( pqueNeighbors.top( ).m_iNode );
		pqueNeighbors.pop( ); }

	return true; }

bool CBNServer::GraphCreate( const vector<size_t>& veciQuery, size_t iContext, size_t iLimit,
	float dEdgeAggressiveness, vector<bool>& vecfQuery, vector<size_t>& veciNeighbors, CDat& DatGraph ) const {
	vector<float>			vecdNeighbors, vecdEdges;
	float					d, dAve, dStd, dMin, dCutoff;
	size_t					i, j, iMinOne, iMinTwo;
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
		vecfQuery[ veciQuery[ i ] ] = true;
	veciNeighbors.clear( );
	fDone = SelectNeighborsRatio( veciQuery, vecfQuery, iContext, iLimit, MatQuery, veciNeighbors );
//		SelectNeighborsPixie( veciQuery, vecfQuery, iContext, iLimit, MatQuery, veciNeighbors );
	if( !fDone )
		return false;

	vecstrGenes.resize( veciQuery.size( ) + veciNeighbors.size( ) );
	vecfQuery.resize( vecstrGenes.size( ) );
	fill( vecfQuery.begin( ), vecfQuery.end( ), false );
	for( i = 0; i < veciQuery.size( ); ++i ) {
		vecstrGenes[ i ] = GetGene( veciQuery[ i ] );
		vecfQuery[ i ] = true; }
	for( i = 0; i < veciNeighbors.size( ); ++i )
		vecstrGenes[ veciQuery.size( ) + i ] = GetGene( veciNeighbors[ i ] );
	MatNeighbors.Initialize( veciNeighbors.size( ), veciNeighbors.size( ) );
	for( i = 0; i < veciNeighbors.size( ); ++i )
		if( !((CBNServer*)this)->Get( veciNeighbors[ i ], veciNeighbors, iContext,
			MatNeighbors.Get( i ) ) )
			return false;
	DatGraph.Open( vecstrGenes );
	for( i = 0; i < veciQuery.size( ); ++i ) {
		for( j = ( i + 1 ); j < veciQuery.size( ); ++j )
			DatGraph.Set( i, j, MatQuery.Get( i, veciQuery[ j ] % MatQuery.GetColumns( ) ) );
		for( j = 0; j < veciNeighbors.size( ); ++j )
			DatGraph.Set( i, veciQuery.size( ) + j, MatQuery.Get( i, veciNeighbors[ j ] ) ); }

	vecdEdges.reserve( veciNeighbors.size( ) * ( veciNeighbors.size( ) - 1 ) / 2 );
	for( i = 0; i < veciNeighbors.size( ); ++i ) {
		const float*	adOne	= MatNeighbors.Get( i );
		size_t			iOne	= veciQuery.size( ) + i;

		for( j = ( i + 1 ); j < veciNeighbors.size( ); ++j )
			if( !CMeta::IsNaN( d = adOne[ j ] ) ) {
				vecdEdges.push_back( d );
				DatGraph.Set( iOne, veciQuery.size( ) + j, adOne[ j ] ); } }

	Winsorize( vecdEdges );
	dAve = (float)CStatistics::Average( vecdEdges );
	dStd = (float)sqrt( CStatistics::Variance( vecdEdges, dAve ) );
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
		strDotIn = GetFiles( ) + "/" + CMeta::Filename( GetGene( veciQuery[ 0 ] ) );
	else {
		sprintf_s( acBuffer, ( GetFiles( ) + "/inXXXXXX" ).c_str( ) );
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

	sprintf_s( acBuffer, ( GetFiles( ) + "/outXXXXXX" ).c_str( ) );
	if( _mktemp_s( acBuffer ) )
		return false;
	strDotOut = acBuffer;
	strDotOut += c_szDOT;
	strCmd = m_sData.m_strGraphviz + " -Tdot -o" + strDotOut + ' ' + strDotIn;
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
			strSvg = GetFiles( ) + '/' + acBuffer + c_szSVG;
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
		aiGenes[ i ] = (uint32_t)veciQuery[ i ] + 1;
	for( i = 0; i < veciNeighbors.size( ); ++i )
		aiGenes[ veciQuery.size( ) + i ] = (uint32_t)veciNeighbors[ i ] + 1;
	send( m_iSocket, (const char*)aiGenes, iSize * sizeof(*aiGenes), 0 );
	delete[] aiGenes;

	return true; }

size_t CBNServer::ProcessContexts( const vector<unsigned char>& vecbMessage, size_t iOffset ) {
	uint32_t				iGene;
	size_t					iCount, iPlace, iContext;
	vector<unsigned char>	vecbData;
	float					dBetween, dBackground, dWithin;

	iCount = InitializeContexts( );
	iGene = (uint32_t)( ( ( vecbMessage.size( ) - iOffset ) / sizeof(iGene) ) *
		( iCount * sizeof(*m_adContexts) ) );
	send( m_iSocket, (const char*)&iGene, sizeof(iGene), 0 );
	for( iPlace = iOffset; ( iPlace + sizeof(iGene) ) <= vecbMessage.size( ); iPlace += sizeof(iGene) ) {
		iGene = *(uint32_t*)&vecbMessage[ iPlace ];
		if( !( iGene-- && GetDatabase( ).Get( iGene, vecbData, true ) ) )
			return -1;
		cerr << m_strConnection << " contexting " << iGene  << " (" << GetGene( iGene ) <<
			")" << endl;
		for( iContext = 0; iContext < GetContexts( ); ++iContext ) {
			if( !GetAssociation( iGene, vecbData, GetContext( iContext ), iContext + 1, true,
				&dBetween, &dBackground, &dWithin, NULL, NULL, GetWithinContext( iContext + 1, iContext ) ) )
				return -1;
			SetArray( m_adContexts, GetContexts( ), iContext, GetPValue( dBetween, dBackground,
				dWithin, iContext + 1, 1, GetContext( iContext ).size( ), GetContexts( ) ), dBetween,
				dBackground, dWithin ); }
		send( m_iSocket, (const char*)m_adContexts, (int)( iCount * sizeof(*m_adContexts) ), 0 ); }

	return ( iPlace - iOffset ); }

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
	CGenes				Genes( GetGenome( ) );
	const IOntology*	pOnto;
	vector<size_t>		veciMapping;

	iSize = sizeof(iOntology) + sizeof(dP);
	if( ( iOffset + iSize ) > vecbMessage.size( ) )
		return -1;
	iOntology = *(uint32_t*)&vecbMessage[ iOffset ];
	dP = *(float*)&vecbMessage[ iOffset + sizeof(iOntology) ];
	for( iPlace = iOffset + iSize; ( iPlace + sizeof(iGene) ) <= vecbMessage.size( );
		iPlace += sizeof(iGene) ) {
		iGene = *(uint32_t*)&vecbMessage[ iPlace ];
		vecstrGenes.push_back( GetGene( iGene - 1 ) ); }

	if( !( pOnto = GetOntology( iOntology ) ) )
		return -1;
	cerr << m_strConnection << " TermFinder " << pOnto->GetID( ) << '@' << dP << ':';
	for( i = 0; i < vecstrGenes.size( ); ++i )
		cerr << ' ' << vecstrGenes[ i ];
	cerr << endl;

	Genes.Open( vecstrGenes, false );
	veciMapping.resize( Genes.GetGenes( ) );
	for( i = 0; i < Genes.GetGenes( ); ++i ) {
		const CGene&	Gene	= Genes.GetGene( i );

		veciMapping[ i ] = GetGene( Gene.GetName( ) );
		for( j = 0; ( j < Gene.GetSynonyms( ) ) && ( veciMapping[ i ] == -1 ); ++j )
			veciMapping[ i ] = GetGene( Gene.GetSynonym( j ) ); }
			
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
	uint32_t				iGene, iContext;
	size_t					i, iCount, iPlace, iDisease, iChunk;
	vector<unsigned char>	vecbData, vecbTmp;
	float					dBackground, dWithin, dBetween;

	if( ( iOffset + sizeof(iContext) ) > vecbMessage.size( ) )
		return -1;
	iContext = *(uint32_t*)&vecbMessage[ iOffset ];

	iCount = InitializeDiseases( );
	iGene = (uint32_t)( ( ( vecbMessage.size( ) - iOffset - sizeof(iContext) ) / sizeof(iGene) ) *
		( iCount * sizeof(*m_adDiseases) ) );
	send( m_iSocket, (const char*)&iGene, sizeof(iGene), 0 );

	iChunk = GetDatachunk( );
	vecbData.resize( GetGenes( ) * iChunk );
	for( iPlace = iOffset + sizeof(iContext); ( iPlace + sizeof(iGene) ) <= vecbMessage.size( );
		iPlace += sizeof(iGene) ) {
		iGene = *(uint32_t*)&vecbMessage[ iPlace ];
		if( !iGene-- )
			return -1;
		cerr << m_strConnection << " diseasing " << iGene  << " (" << GetGene( iGene ) <<
			") in " << iContext << endl;

		if( !GetDatabase( ).Get( iGene, GetDiseaseGenes( ), vecbTmp, true ) )
			return -1;
		for( i = 0; i < GetDiseaseGenes( ).size( ); ++i )
// Hackedy hack hack hack
			memcpy( &vecbData[ GetDiseaseGenes( )[ i ] * iChunk ], &vecbTmp[ i * iChunk ], iChunk );
		for( iDisease = 0; iDisease < GetDiseases( ); ++iDisease ) {
			if( !GetAssociation( iGene, vecbData, GetDisease( iDisease ), iContext, true, &dBetween,
				&dBackground, &dWithin, NULL, NULL, GetWithinDisease( iContext, iDisease ) ) )
				return -1;
			SetArray( m_adDiseases, GetDiseases( ), iDisease, GetPValue( dBetween, dBackground,
				dWithin, iContext, 1, GetDisease( iDisease ).size( ), GetDiseases( ) ), dBetween,
				dBackground, dWithin ); }
		send( m_iSocket, (const char*)m_adDiseases, (int)( iCount * sizeof(*m_adDiseases) ), 0 ); }

	return ( iPlace - iOffset ); }

bool CBNServer::GenerateNetworkIcons( ) const {
	static const size_t	c_iSize					= 6;
	static const float	c_dEdgeAggressiveness	= 0;
	vector<size_t>	veciQuery, veciNeighbors;
	vector<bool>	vecfQuery;
	CDat			DatGraph;
	size_t			i;

	veciQuery.resize( 1 );
	for( i = 0; i < GetGenes( ); ++i ) {
		veciQuery[ 0 ] = i;
		if( !( GraphCreate( veciQuery, 0, c_iSize, c_dEdgeAggressiveness, vecfQuery, veciNeighbors,
			DatGraph ) && GraphWrite( DatGraph, veciQuery, veciNeighbors, vecfQuery, 0, EGraphOutputNamed ) ) )
			return false; }

	return true; }

bool CBNServer::GenerateAssociations( const char* szDiseases, size_t iContext ) {
	CPCL			PCLDiseases( false ), PCLAssociations;
	size_t			i, j;
	ofstream		ofsm;
	vector<string>	vecstrDCs;

	if( !PCLDiseases.Open( szDiseases, 0 ) ) {
		cerr << "Could not open: " << szDiseases << endl;
		return false; }
	if( PCLDiseases.GetGenes( ) != GetDiseases( ) ) {
		cerr << "Unexpected disease count " << PCLDiseases.GetGenes( ) << ", wanted " << GetDiseases( ) <<
			endl;
		return false; }
	vecstrDCs.resize( GetDiseases( ) + GetContexts( ) );
	copy( PCLDiseases.GetGeneNames( ).begin( ), PCLDiseases.GetGeneNames( ).end( ), vecstrDCs.begin( ) );
	for( i = 0; i < GetContexts( ); ++i )
		vecstrDCs[ GetDiseases( ) + i ] = GetBN( i + 1 ).GetID( );
	PCLAssociations.Open( vecstrDCs, vecstrDCs, vector<string>( ) );

	// Disease/disease
	InitializeDiseases( );
	for( i = 0; i < GetDiseases( ); ++i ) {
		memset( m_adDiseases, 0, GetDiseases( ) * c_iValues * sizeof(*m_adDiseases) );
		if( !GetAssociationsDC( 1, ESetDisease, i, iContext, true ) )
			return false;
		for( j = i; j < GetDiseases( ); ++j ) {
			PCLAssociations.Set( i, j, m_adDiseases[ j ] );
			PCLAssociations.Set( j, i, m_adDiseases[ j ] ); } }
	// Disease/context
	InitializeContexts( );
	for( i = 0; i < GetDiseases( ); ++i ) {
		memset( m_adContexts, 0, GetContexts( ) * c_iValues * sizeof(*m_adContexts) );
		if( !GetAssociationsDC( 0, ESetDisease, i, -1, true ) )
			return false;
		for( j = 0; j < GetContexts( ); ++j )
			PCLAssociations.Set( i, GetDiseases( ) + j, m_adContexts[ j ] ); }
	for( i = 0; i < GetContexts( ); ++i ) {
		memset( m_adDiseases, 0, GetDiseases( ) * c_iValues * sizeof(*m_adDiseases) );
		if( !GetAssociationsDC( 1, ESetContext, i, iContext, true ) )
			return false;
		for( j = 0; j < GetDiseases( ); ++j )
			PCLAssociations.Set( GetDiseases( ) + i, j, m_adDiseases[ j ] ); }
	// Context/context
	for( i = 0; i < GetContexts( ); ++i ) {
		memset( m_adContexts, 0, GetContexts( ) * c_iValues * sizeof(*m_adContexts) );
		if( !GetAssociationsDC( 0, ESetContext, i, -1, true ) )
			return false;
		for( j = 0; j < GetContexts( ); ++j )
			PCLAssociations.Set( GetDiseases( ) + i, GetDiseases( ) + j, m_adContexts[ j ] ); }

	ofsm.open( "associations.pcl" );
	PCLAssociations.Save( ofsm );
	return true; }

size_t CBNServer::ProcessGenes( const vector<unsigned char>& vecbMessage, size_t iOffset ) {
	uint32_t		iContext, iSize, iGene;
	size_t			i, iPlace;
	vector<size_t>	veciGenes;
	unsigned char	bType;

	if( ( iOffset + sizeof(iContext) + sizeof(bType) ) > vecbMessage.size( ) )
		return -1;
	iContext = *(uint32_t*)&vecbMessage[ iOffset ];
	bType = vecbMessage[ iOffset + sizeof(iContext) ];
	for( iPlace = iOffset + sizeof(iContext) + sizeof(bType);
		( iPlace + sizeof(iGene) ) <= vecbMessage.size( ); iPlace += sizeof(iGene) ) {
		if( !( iGene = *(uint32_t*)&vecbMessage[ iPlace ] ) )
			return -1;
		iGene--;
		if( bType == ESetGenes )
			veciGenes.push_back( iGene ); }

	cerr << m_strConnection << " genesing in " << iContext << ' ';
	if( bType == ESetGenes ) {
		for( i = 0; i < veciGenes.size( ); ++i )
			cerr << ( i ? ' ' : '{' ) << GetGene( veciGenes[ i ] );
		cerr << '}' << endl; }
	else
		cerr << ( ( bType == ESetContext ) ? "context" : "disease" ) << ' ' << ( iGene + 1 ) << endl;
	{
		const vector<size_t>&	veciSet = ( bType == ESetGenes ) ? veciGenes :
											GetGeneSets( bType - 1 )[ iGene ];
		float					dWithin;

		dWithin = ( bType == ESetGenes ) ? CMeta::GetNaN( ) : GetWithin( bType - 1, iContext, iGene );
		if( !GetGenes( veciSet, iContext, dWithin ) )
			return -1;
	}
	iSize = (uint32_t)( c_iValues * GetGenes( ) * sizeof(*m_adGenes) );
	send( m_iSocket, (char*)&iSize, sizeof(iSize), 0 );
	send( m_iSocket, (char*)m_adGenes, iSize, 0 );

	return ( iPlace - iOffset ); }

bool CBNServer::GetGenes( const vector<size_t>& veciGenes, size_t iContext, float dWithin ) {
	CDataMatrix		MatGenes;
	size_t			i, j;
	vector<float>	vecdBetween, vecdBackground;
	float			d, dFrac, dBackground;

	dFrac = GetFraction( veciGenes.size( ) );
	MatGenes.Initialize( veciGenes.size( ), GetGenes( ) );
	vecdBackground.resize( veciGenes.size( ) );
	for( i = 0; i < MatGenes.GetRows( ); ++i ) {
		vecdBackground[ i ] = GetBackground( iContext, veciGenes[ i ] );
		if( IsFraction( dFrac ) ) {
			fill( MatGenes.Get( i ), MatGenes.Get( i ) + MatGenes.GetColumns( ), CMeta::GetNaN( ) );
			continue; }
		if( !Get( veciGenes[ i ], iContext, MatGenes.Get( i ) ) )
			return false; }
	Winsorize( vecdBackground );
	dBackground = (float)CStatistics::Average( vecdBackground );
	i = veciGenes.size( ) * ( veciGenes.size( ) + 1 ) / 2;
	if( CMeta::IsNaN( dWithin ) ) {
		vector<float>	vecdWithin;

		vecdWithin.reserve( 1 + i );
		if( !GetWithin( veciGenes, iContext, NULL, &vecdWithin ) )
			return false;
// Alternative: calculate an average over all edges (weighted) rather than between sets (unweighted)
//		vecdWithin.push_back( GetPrior( iContext ) );
		Winsorize( vecdWithin );
//		dWithin = (float)CStatistics::Average( vecdWithin ); }
		dWithin = ( (float)CStatistics::Average( vecdWithin ) + GetPrior( iContext ) ) / 2; }
	else
//		dWithin = ( ( dWithin * i ) + GetPrior( iContext )  ) / ( i + 1 );
		dWithin = ( dWithin + GetPrior( iContext ) ) / 2;

	InitializeGenes( );
	for( i = 0; i < GetGenes( ); ++i ) {
		float	dBetween, dTmp;

		vecdBetween.clear( );
		vecdBetween.reserve( MatGenes.GetRows( ) );
		for( j = 0; j < MatGenes.GetRows( ); ++j )
			if( !CMeta::IsNaN( d = MatGenes.Get( j, i ) ) )
				vecdBetween.push_back( d );
		Winsorize( vecdBetween );
		dBetween = (float)CStatistics::Average( vecdBetween );
		dTmp = ( ( dBackground * veciGenes.size( ) ) + GetBackground( iContext, i ) ) /
			( veciGenes.size( ) + 1 );
		SetArray( m_adGenes, GetGenes( ), i, GetPValue( dBetween, dTmp, dWithin, iContext,
// Alternative: multiple hypothesis correction on all genes (gets too big)
//			1, veciGenes.size( ), GetGenes( ) ), dBetween, dTmp, dWithin ); }
			1, veciGenes.size( ), 0 ), dBetween, dTmp, dWithin ); }

	return true; }

size_t CBNServer::ProcessAssociation( const vector<unsigned char>& vecbMessage, size_t iOffset ) {
	uint32_t				iGene, iContext;
	size_t					i, iPlace, iBack;
	vector<unsigned char>	vecbData;
	vector<size_t>			veciFore, veciBack;
	float					dBetween, dBackground, dWithin, dWith;
	const vector<size_t>*	pveciBack;
	unsigned char			bType;

	if( ( iOffset + sizeof(iContext) ) > vecbMessage.size( ) )
		return -1;
	iContext = *(uint32_t*)&vecbMessage[ iOffset ];
	for( iPlace = iOffset + sizeof(iContext); ( iPlace + sizeof(iGene) ) <= vecbMessage.size( );
		iPlace += sizeof(iGene) ) {
		if( !( iGene = *(uint32_t*)&vecbMessage[ iPlace ] ) )
			break;
		veciFore.push_back( --iGene ); }
	if( ( ( iPlace += sizeof(iGene) ) + sizeof(bType) ) > vecbMessage.size( ) )
		return -1;
	bType = vecbMessage[ iPlace++ ];
	if( bType == ESetGenes ) {
		for( ; ( iPlace + sizeof(iGene) ) <= vecbMessage.size( ); iPlace += sizeof(iGene) ) {
			if( !( iGene = *(uint32_t*)&vecbMessage[ iPlace ] ) )
				return -1;
			veciBack.push_back( --iGene ); }
		iBack = -1;
		pveciBack = &veciBack; }
	else {
		if( ( ( iPlace + sizeof(iGene) ) > vecbMessage.size( ) ) ||
			!( iGene = *(uint32_t*)&vecbMessage[ iPlace ] ) )
			return -1;
		iPlace += sizeof(iGene);
		iBack = iGene - 1;
		pveciBack = &GetGeneSets( bType - 1 )[ iBack ]; }
	if( veciFore.empty( ) || pveciBack->empty( ) )
		return -1;

	cerr << m_strConnection << " association in " << iContext << ' ';
	for( i = 0; i < veciFore.size( ); ++i )
		cerr << ( i ? ' ' : '{' ) << GetGene( veciFore[ i ] );
	cerr << "} with ";
	if( iBack == -1 ) {
		for( i = 0; i < pveciBack->size( ); ++i )
			cerr << ( i ? ' ' : '{' ) << GetGene( pveciBack->at( i ) );
		cerr << '}'; }
	else
		cerr << ( ( bType == ESetContext ) ? "context" : "disease" ) << ' ' << ( iBack + 1 );
	cerr << endl;

	dWith = ( bType == ESetGenes ) ? CMeta::GetNaN( ) : GetWithin( bType - 1, iContext, iBack );
	InitializeGenes( );
	for( i = 0; i < veciFore.size( ); ++i ) {
		if( !( i % 100 ) )
			cerr << m_strConnection << " association " << i << '/' << veciFore.size( ) << ": " <<
				GetGene( veciFore[ i ] ) << endl;
		if( !( GetDatabase( ).Get( veciFore[ i ], *pveciBack, vecbData, true ) &&
			GetAssociation( veciFore[ i ], vecbData, *pveciBack, iContext, false, &dBetween, &dBackground,
				&dWithin, NULL, NULL, dWith ) ) )
			return -1;
		SetArray( m_adGenes, veciFore.size( ), i, GetPValue( dBetween, dBackground, dWithin,
			iContext, 1, pveciBack->size( ), veciFore.size( ) ), dBetween, dBackground, dWithin ); }

	iGene = (uint32_t)( c_iValues * veciFore.size( ) * sizeof(*m_adGenes) );
	send( m_iSocket, (const char*)&iGene, sizeof(iGene), 0 );
	send( m_iSocket, (const char*)m_adGenes, iGene, 0 );

	return ( iPlace - iOffset ); }

size_t CBNServer::ProcessAssociations( const vector<unsigned char>& vecbMessage, size_t iOffset ) {
	uint32_t		iGene, iCount;
	vector<size_t>	veciGenes;
	unsigned char	bType, bDiseases;
	size_t			i, iPlace, iContext;
	float*			ad;

	if( ( iOffset + sizeof(iGene) + sizeof(bType) + sizeof(bDiseases) ) > vecbMessage.size( ) )
		return -1;
	bDiseases = vecbMessage[ iOffset ];
	iGene = *(uint32_t*)&vecbMessage[ iOffset + sizeof(bDiseases) ];
	iContext = ( iGene == -1 ) ? -1 : (size_t)iGene;
	bType = vecbMessage[ iOffset + sizeof(bDiseases) + sizeof(iGene) ];
	for( iPlace = ( iOffset + sizeof(iGene) + sizeof(bDiseases) + sizeof(bType) );
		( iPlace + sizeof(iGene) ) <= vecbMessage.size( ); iPlace += sizeof(iGene) ) {
		if( !( iGene = *(uint32_t*)&vecbMessage[ iPlace ] ) )
			return -1;
		iGene--;
		if( bType == ESetGenes )
			veciGenes.push_back( iGene ); }

	cerr << m_strConnection << " associationsing ";
	if( bType == ESetGenes ) {
		for( i = 0; i < veciGenes.size( ); ++i )
			cerr << ( i ? ' ' : '{' ) << GetGene( veciGenes[ i ] );
		cerr << '}'; }
	else
		cerr << ( ( bType == ESetContext ) ? "context" : "disease" ) << ' ' << ( iGene + 1 );
	cerr << " with " << ( bDiseases ? "diseases" : "contexts" ) << " in " << iContext << endl;

	iCount = (uint32_t)( ( bDiseases ? InitializeDiseases( ) : InitializeContexts( ) ) * sizeof(*ad) );
	ad = bDiseases ? m_adDiseases : m_adContexts;
	memset( ad, 0, iCount );
	if( !( ( bType == ESetGenes ) ?
		GetAssociationsSet( bDiseases, veciGenes, iContext ) :
		GetAssociationsDC( bDiseases, bType, iGene, iContext ) ) )
		return -1;
	send( m_iSocket, (char*)&iCount, sizeof(iCount), 0 );
	send( m_iSocket, (char*)ad, iCount, 0 );

	return ( iPlace - iOffset ); }

bool CBNServer::GetAssociationsSet( unsigned char bDiseases, const vector<size_t>& veciGenes,
	size_t iContext ) const {
	const vector<vector<size_t> >&	vecveciSets	= GetGeneSets( bDiseases );
	size_t					i, j, k;
	vector<unsigned char>	vecbData;
	float					d, dFrac, dWithinGenes, dBetween;
	float*					ad;
	vector<vector<float> >	vecvecdBackground, vecvecdBetween;
	vector<float>			vecdCur;

	vecvecdBetween.resize( vecveciSets.size( ) );
	for( i = 0; i < vecvecdBetween.size( ); ++i ) {
		vecvecdBetween[ i ].resize( veciGenes.size( ) + vecveciSets[ i ].size( ) );
		fill( vecvecdBetween[ i ].begin( ), vecvecdBetween[ i ].end( ), 0.0f ); }
	vecvecdBackground.resize( vecveciSets.size( ) );
	dFrac = GetFraction( veciGenes.size( ) );
	for( i = 0; i < veciGenes.size( ); ++i ) {
		if( IsFraction( dFrac ) ) {
			for( j = 0; j < vecveciSets.size( ); ++j ) {
				k = GetContext( bDiseases, iContext, j );
				vecvecdBackground[ j ].push_back( GetBackground( k, veciGenes[ i ] ) ); }
			continue; }
		if( !GetDatabase( ).Get( veciGenes[ i ], vecbData, true ) )
			return false;
		cerr << m_strConnection << " associations " << GetGene( veciGenes[ i ] ) << endl;
		for( j = 0; j < vecveciSets.size( ); ++j ) {
			vecdCur.clear( );
			vecdCur.reserve( vecveciSets[ j ].size( ) );
			k = GetContext( bDiseases, iContext, j );
			if( !GetAssociation( veciGenes[ i ], vecbData, vecveciSets[ j ], k, true, &dBetween, NULL, NULL,
				&vecdCur, &vecvecdBackground[ j ], GetWithin( bDiseases, k, j ) ) )
				return false;
			vecvecdBetween[ j ][ i ] += dBetween;
			for( k = 0; k < vecveciSets[ j ].size( ); ++k )
				vecvecdBetween[ j ][ veciGenes.size( ) + k ] += vecdCur[ k ]; } }
	for( i = 0; i < vecvecdBetween.size( ); ++i ) {
		for( j = 0; j < veciGenes.size( ); ++j )
			vecvecdBetween[ i ][ j ] /= vecveciSets[ i ].size( );
		for( j = 0; j < vecveciSets[ i ].size( ); ++j )
			vecvecdBetween[ i ][ veciGenes.size( ) + j ] /= veciGenes.size( ); }
	for( i = 0; i < vecveciSets.size( ); ++i ) {
		k = GetContext( bDiseases, iContext, i );
		for( j = 0; j < vecveciSets[ i ].size( ); ++j )
			if( !CMeta::IsNaN( d = GetBackground( k, vecveciSets[ i ][ j ] ) ) )
				vecvecdBackground[ i ].push_back( d ); }
	if( ( iContext != -1 ) && !GetWithin( veciGenes, iContext, &dWithinGenes, NULL ) )
		return false;

	ad = bDiseases ? m_adDiseases : m_adContexts;
	for( i = 0; i < vecveciSets.size( ); ++i ) {
		float	dBackground, dWithin;
		size_t	iTmp;

		if( vecvecdBetween[ i ].empty( ) )
			continue;
		Winsorize( vecvecdBetween[ i ] );
		Winsorize( vecvecdBackground[ i ] );
		dBetween = (float)CStatistics::Average( vecvecdBetween[ i ] );
		dBackground = (float)CStatistics::Average( vecvecdBackground[ i ] );
		iTmp = GetContext( bDiseases, iContext, i );
		if( ( iContext == -1 ) && !GetWithin( veciGenes, iTmp, &dWithinGenes, NULL ) )
			return false;
		dWithin = ( dWithinGenes + GetWithin( bDiseases, iTmp, i ) ) / 2;
		SetArray( ad, vecveciSets.size( ), i, GetPValue( dBetween, dBackground, dWithin, iTmp,
			veciGenes.size( ), vecveciSets[ i ].size( ), vecveciSets.size( ) ), dBetween, dBackground,
			dWithin ); }

	return true; }

bool CBNServer::GetAssociationsDC( unsigned char bDiseases, unsigned char bType, size_t iDC,
	size_t iContext, bool fZ ) const {
	const vector<vector<size_t> >&	vecveciSets	= GetGeneSets( bDiseases );
	const vector<size_t>&			veciGenes	= GetGeneSets( bType - 1 )[ iDC ];
	size_t	i, j, k, iEdgesOne, iEdgesTwo;
	float	d, dBetween, dWithin, dBackground, dBackgroundOne;
	float*	ad;

	if( iContext == -1 )
		dBackgroundOne = CMeta::GetNaN( );
	else
		for( dBackgroundOne = 0,i = 0; i < veciGenes.size( ); ++i )
			dBackgroundOne += GetBackground( iContext, veciGenes[ i ] );
	iEdgesOne = veciGenes.size( ) * ( veciGenes.size( ) + 1 ) / 2;
	ad = bDiseases ? m_adDiseases : m_adContexts;
	for( i = 0; i < vecveciSets.size( ); ++i ) {
		iEdgesTwo = vecveciSets[ i ].size( ) * ( vecveciSets[ i ].size( ) + 1 ) / 2;
		k = GetContext( bDiseases, iContext, i );
		if( CMeta::IsNaN( dBackgroundOne ) )
			for( dBackground = 0,j = 0; j < veciGenes.size( ); ++j )
				dBackground += GetBackground( k, veciGenes[ j ] );
		else
			dBackground = dBackgroundOne;
		for( j = 0; j < vecveciSets[ i ].size( ); ++j )
			dBackground += GetBackground( k, vecveciSets[ i ][ j ] );
		dBackground /= veciGenes.size( ) + vecveciSets[ i ].size( );

		d = GetWithin( bType - 1, k, iDC );
		dWithin = GetWithin( bDiseases, k, i );
// Alternative: calculate an average over all edges (weighted) rather than between sets (unweighted)
//		dWithin = ( ( d * iEdgesOne ) + ( dWithin * iEdgesTwo ) ) / ( iEdgesOne + iEdgesTwo );
		dWithin = ( d + dWithin ) / 2;
		dBetween = GetBetween( k, bType - 1, iDC, bDiseases, i );

		SetArray( ad, vecveciSets.size( ), i, GetPValue( dBetween, dBackground, dWithin, k, veciGenes.size( ),
			vecveciSets[ i ].size( ), vecveciSets.size( ), fZ ), dBetween, dBackground, dWithin ); }

	return true; }

bool CBNServer::GetAssociation( size_t iGene, const vector<unsigned char>& vecbData,
	const vector<size_t>& veciGenes, size_t iContext, bool fAll, float* pdBetween, float* pdBackground,
	float* pdWithin, vector<float>* pvecdBetween, vector<float>* pvecdBackground, float dWithin ) const {
	const CBayesNetMinimal&	BNet	= GetBN( iContext );
	float					d, dFrac, dBackground;
	size_t					i;
	vector<float>			vecdBetween, vecdBackground;

	if( veciGenes.empty( ) ) {
		SetPointer( pdBetween, CMeta::GetNaN( ) );
		SetPointer( pdBackground, CMeta::GetNaN( ) );
		SetPointer( pdWithin, CMeta::GetNaN( ) );
		return true; }

	dFrac = GetFraction( veciGenes.size( ) );
	if( pdBetween )
		vecdBetween.reserve( veciGenes.size( ) );
	if( pdBackground )
		vecdBackground.reserve( veciGenes.size( ) + 1 );
	dBackground = GetBackground( iContext, iGene );
	PushPointer( pdBackground, vecdBackground, dBackground );
	PushPointer( pvecdBackground, dBackground );
	for( i = 0; i < veciGenes.size( ); ++i ) {
		PushPointer( pdBackground, vecdBackground, GetBackground( iContext, veciGenes[ i ] ) );
		if( veciGenes[ i ] == iGene )
			d = 1;
		else if( IsFraction( dFrac ) )
			continue;
		else
			d = BNet.Evaluate( vecbData, ( fAll ? veciGenes[ i ] : i ) * GetDatachunk( ) );
		if( !CMeta::IsNaN( d ) ) {
			PushPointer( pdBetween, vecdBetween, d );
			PushPointer( pvecdBetween, d ); } }
	Winsorize( vecdBetween );
	SetPointer( pdBetween, (float)CStatistics::Average( vecdBetween ) );
	Winsorize( vecdBackground );
	SetPointer( pdBackground, (float)CStatistics::Average( vecdBackground ) );

	if( pdWithin ) {
		i = veciGenes.size( ) * ( veciGenes.size( ) + 1 ) / 2;
		if( CMeta::IsNaN( dWithin ) ) {
			vector<float>	vecdWithin;

			vecdWithin.reserve( 1 + i );
// Alternative: calculate an average over all edges (weighted) rather than between sets (unweighted)
//			vecdWithin.push_back( GetPrior( iContext ) );
			if( !GetWithin( veciGenes, iContext, NULL, &vecdWithin ) )
				return false;
			Winsorize( vecdWithin );
//			SetPointer( pdWithin, (float)CStatistics::Average( vecdWithin ) ); }
			SetPointer( pdWithin, ( (float)CStatistics::Average( vecdWithin ) + GetPrior( iContext ) ) / 2 ); }
		else
//			SetPointer( pdWithin, ( GetPrior( iContext ) + ( i * dWithin ) ) / ( 1 + i ) ); }
			SetPointer( pdWithin, ( GetPrior( iContext ) + dWithin ) / 2 ); }

	return true; }

bool CBNServer::GetAssociation( const vector<size_t>& veciOne, const vector<size_t>& veciTwo,
	size_t iContext, float& dBetween, float& dBackground, float& dWithin ) const {
	const vector<size_t>&	veciSmall	= ( veciOne.size( ) < veciTwo.size( ) ) ? veciOne : veciTwo;
	const vector<size_t>&	veciBig		= ( veciOne.size( ) < veciTwo.size( ) ) ? veciTwo : veciOne;
	size_t					i, j;
	float					dFrac;
	vector<unsigned char>	vecbData;
	vector<float>			vecdBetween, vecdBackground, vecdWithin;

	dFrac = GetFraction( veciSmall.size( ) );
	vecdBetween.reserve( (size_t)( dFrac * veciSmall.size( ) * veciBig.size( ) ) + c_iOverestimate );
	vecdBackground.reserve( veciSmall.size( ) + veciBig.size( ) );
	for( i = 0; i < veciSmall.size( ); ++i ) {
		if( IsFraction( dFrac ) ) {
			vecdBackground.push_back( GetBackground( iContext, veciSmall[ i ] ) );
			continue; }
		if( !( GetDatabase( ).Get( veciSmall[ i ], veciBig, vecbData, true ) &&
			GetAssociation( veciSmall[ i ], vecbData, veciBig, iContext, false, NULL, NULL, NULL,
			&vecdBetween, &vecdBackground, CMeta::GetNaN( ) ) ) )
			return false; }
	for( i = 0; i < veciBig.size( ); ++i )
		vecdBackground.push_back( GetBackground( iContext, veciBig[ i ] ) );

	i = veciSmall.size( ) + 1;
	j = veciBig.size( ) + 1;
	vecdWithin.reserve( (size_t)( ( dFrac * i * i / 2 ) + ( dFrac * j * j / 2 ) ) + c_iOverestimate );
	if( !( GetWithin( veciSmall, iContext, NULL, &vecdWithin ) &&
		GetWithin( veciBig, iContext, NULL, &vecdWithin ) ) )
		return false;

	Winsorize( vecdBetween );
	dBetween = (float)CStatistics::Average( vecdBetween );
	Winsorize( vecdBackground );
	dBackground = (float)CStatistics::Average( vecdBackground );
	Winsorize( vecdWithin );
	dWithin = (float)CStatistics::Average( vecdWithin );

	return true; }

bool CBNServer::GetWithin( const vector<size_t>& veciGenes, size_t iContext, float* pdWithin,
	vector<float>* pvecdWithin ) const {
// Alternative: calculate priors in the current context rather than the global context
	iContext = 0;
	const CBayesNetMinimal&	BNet	= GetBN( iContext );
	size_t					i, j;
	vector<unsigned char>	vecbData;
	vector<float>			vecdWithin;
	float					d, dFrac;

	dFrac = GetFraction( veciGenes.size( ) );
	if( pdWithin ) {
		i = veciGenes.size( ) + 1;
		vecdWithin.reserve( (size_t)( dFrac * i * i / 2 ) + c_iOverestimate ); }
	for( i = 0; i < veciGenes.size( ); ++i ) {
		if( IsFraction( dFrac ) )
			continue;
		if( !GetDatabase( ).Get( veciGenes[ i ], veciGenes, vecbData, true ) )
			return false;
		d = GetPrior( iContext );
		PushPointer( pdWithin, vecdWithin, d );
		PushPointer( pvecdWithin, d );
		for( j = ( i + 1 ); j < veciGenes.size( ); ++j ) {
			if( IsFraction( dFrac ) )
				continue;
			if( !CMeta::IsNaN( d = BNet.Evaluate( vecbData, j * GetDatachunk( ) ) ) ) {
				PushPointer( pdWithin, vecdWithin, d );
				PushPointer( pvecdWithin, d ); } } }
	if( pdWithin ) {
		Winsorize( vecdWithin );
		SetPointer( pdWithin, (float)CStatistics::Average( vecdWithin ) ); }

	return true; }

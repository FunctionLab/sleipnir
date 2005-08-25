#include "stdafx.h"

#ifdef BAYESIAN_NETWORKS

#include "bayesnet.h"
#include "dat.h"
#include "dataset.h"
#include "meta.h"

namespace libBioUtils {

const char	CBayesNetSmileImpl::c_szGaussian[]	= "gaussian";

string CBayesNetImpl::EncodeDatum( const IDataset* pData, size_t iOne, size_t iTwo ) {
	string	strRet;
	size_t	i, iCur;

	for( i = 0; i < pData->GetExperiments( ); ++i )
		strRet += ( ( iCur = pData->GetDiscrete( iOne, iTwo, i ) ) == -1 ) ? '_' :
			(char)( 'A' + ( iCur & 0xFF ) );

	return strRet; }

void CBayesNetImpl::DecodeDatum( const string& strDatum, vector<size_t>& veciDatum ) {
	size_t	i;

	for( i = 0; i < strDatum.length( ); ++i )
		veciDatum[ i ] = ( strDatum[ i ] == '_' ) ? -1 : ( strDatum[ i ] - 'A' ); }

CBayesNetImpl::CBayesNetImpl( bool fGroup ) : m_fGroup(fGroup) { }

bool CBayesNetSmileImpl::IsGaussian( const DSL_network& BayesNet ) {
	int	i;

	if( ( i = ((DSL_network&)BayesNet).UserProperties( ).FindProperty( c_szGaussian ) ) < 0 )
		return false;

	return !!atoi( ((DSL_network&)BayesNet).UserProperties( ).GetPropertyValue( i ) ); }

CBayesNetSmileImpl::CBayesNetSmileImpl( bool fGroup ) : CBayesNetImpl(fGroup),
	m_fSmileNet(false) { }

CBayesNetSmile::CBayesNetSmile( bool fGroup ) : CBayesNetSmileImpl( fGroup ) { }

bool CBayesNetSmile::Open( const char* szDSL ) {

	return ( m_fSmileNet = !m_SmileNet.ReadFile( szDSL ) ); }

bool CBayesNetSmile::Save( const char* szDSL ) const {

	if( !m_fSmileNet )
		return false;

	return !((CBayesNetSmile*)this)->m_SmileNet.WriteFile( szDSL ); }

bool CBayesNetSmileImpl::LearnGrouped( const IDataset* pData, size_t iIterations,
	bool fZero ) {
	size_t							i, j, iVal, iIter, iDatum;
	string							strCur;
	map<string,size_t>				mapData;
	map<string,size_t>::iterator	iterDatum;
	DSL_Dmatrix*					pMat;
	vector<DSL_Dmatrix*>			vecpExpected;
	DSL_intArray					veciCoords;
	vector<size_t>					veciDatum;

	for( i = 0; i < pData->GetGenes( ); ++i )
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j ) {
			if( !pData->IsExample( i, j ) )
				continue;
			strCur = EncodeDatum( pData, i, j );
			if( ( iterDatum = mapData.find( strCur ) ) == mapData.end( ) )
				mapData[ strCur ] = 1;
			else
				iterDatum->second++; }

	veciDatum.resize( m_SmileNet.GetNumberOfNodes( ) );
	vecpExpected.resize( m_SmileNet.GetNumberOfNodes( ) );
	for( i = 0; i < vecpExpected.size( ); ++i )
		vecpExpected[ i ] = new DSL_Dmatrix( *m_SmileNet.GetNode( (int)i )->Definition(
			)->GetMatrix( ) );
	for( iIter = 0; iIter < iIterations; ++iIter ) {
		for( iDatum = i = 0; i < vecpExpected.size( ); ++i )
			vecpExpected[ i ]->FillWith( 0 );
		for( iterDatum = mapData.begin( ); iterDatum != mapData.end( ); ++iterDatum ) {
			if( !( iDatum++ % 50 ) )
				g_CatBioUtils.notice( "CBayesNetSmile::Learn( %d, %d ) iteration %d, datum %d/%d",
					iIterations, fZero, iIter, ( iDatum - 1 ), mapData.size( ) );
			DecodeDatum( iterDatum->first, veciDatum );
			m_SmileNet.ClearAllEvidence( );
			for( i = 0; i < m_SmileNet.GetNumberOfNodes( ); ++i ) {
				if( pData->IsHidden( i ) )
					continue;
				if( ( iVal = veciDatum[ i ] ) == -1 ) {
					if( fZero )
						iVal = 0;
					else
						continue; }
				m_SmileNet.GetNode( (int)i )->Value( )->SetEvidence( (int)iVal ); }
			m_SmileNet.UpdateBeliefs( );

			for( i = 0; i < m_SmileNet.GetNumberOfNodes( ); ++i )
				LearnExpected( m_SmileNet.GetNode( (int)i ), vecpExpected[ i ],
					iterDatum->second ); }
		for( i = 0; i < m_SmileNet.GetNumberOfNodes( ); ++i ) {
			pMat = m_SmileNet.GetNode( (int)i )->Definition( )->GetMatrix( );
			for( pMat->IndexToCoordinates( (int)( j = 0 ), veciCoords );
				j != DSL_OUT_OF_RANGE; j = pMat->NextCoordinates( veciCoords ) )
				pMat->Subscript( veciCoords ) = vecpExpected[ i ]->Subscript( veciCoords );
			pMat->Normalize( ); } }
	for( i = 0; i < vecpExpected.size( ); ++i )
		delete vecpExpected[ i ];

	return true; }

bool CBayesNetSmileImpl::LearnUngrouped( const IDataset* pData, size_t iIterations,
	bool fZero ) {
	size_t					iIter, i, j, k, iVal;
	DSL_Dmatrix*			pMat;
	vector<DSL_Dmatrix*>	vecpExpected;
	DSL_intArray			veciCoords;

	if( !m_fSmileNet || IsContinuous( ) )
		return false;

	vecpExpected.resize( m_SmileNet.GetNumberOfNodes( ) );
	for( i = 0; i < vecpExpected.size( ); ++i )
		vecpExpected[ i ] = new DSL_Dmatrix( *m_SmileNet.GetNode( (int)i )->Definition(
			)->GetMatrix( ) );
	for( iIter = 0; iIter < iIterations; ++iIter ) {
		for( i = 0; i < vecpExpected.size( ); ++i )
			vecpExpected[ i ]->FillWith( 0 );
		for( i = 0; i < pData->GetGenes( ); ++i ) {
			if( !( i % 50 ) )
				g_CatBioUtils.notice( "CBayesNetSmile::Learn( %d, %d ) iteration %d, gene %d/%d",
					iIterations, fZero, iIter, i, pData->GetGenes( ) );
			for( j = ( i + 1 ); j < pData->GetGenes( ); ++j ) {
				if( !pData->IsExample( i, j ) )
					continue;

				m_SmileNet.ClearAllEvidence( );
				for( k = 0; k < m_SmileNet.GetNumberOfNodes( ); ++k ) {
					if( pData->IsHidden( k ) )
						continue;
					if( ( iVal = pData->GetDiscrete( i, j, k ) ) == -1 ) {
						if( fZero )
							iVal = 0;
						else
							continue; }
					m_SmileNet.GetNode( (int)k )->Value( )->SetEvidence( (int)iVal ); }
				m_SmileNet.UpdateBeliefs( );

				for( k = 0; k < m_SmileNet.GetNumberOfNodes( ); ++k )
					LearnExpected( m_SmileNet.GetNode( (int)k ), vecpExpected[ k ] ); }
		for( i = 0; i < m_SmileNet.GetNumberOfNodes( ); ++i ) {
			pMat = m_SmileNet.GetNode( (int)i )->Definition( )->GetMatrix( );
			for( pMat->IndexToCoordinates( (int)( j = 0 ), veciCoords );
				j != DSL_OUT_OF_RANGE; j = pMat->NextCoordinates( veciCoords ) )
				pMat->Subscript( veciCoords ) = vecpExpected[ i ]->Subscript( veciCoords );
			pMat->Normalize( ); } } }
	for( i = 0; i < vecpExpected.size( ); ++i )
		delete vecpExpected[ i ];

	return true; }

bool CBayesNetSmile::Learn( const IDataset* pData, size_t iIterations, bool fZero ) {

	return ( m_fGroup ? LearnGrouped( pData, iIterations, fZero ) :
		LearnUngrouped( pData, iIterations, fZero ) ); }

void CBayesNetSmileImpl::LearnExpected( DSL_node* pNode, DSL_Dmatrix* pExpected,
	size_t iWeight ) {
	int				iEvid, iLast, i, j;
	DSL_intArray	veciParents, veciCoords;
	DSL_Dmatrix*	pDef;
	DSL_nodeValue*	pVal;
	double			dProd;

	veciParents = pNode->Parents( );
	pDef = pNode->Definition( )->GetMatrix( );
	pVal = pNode->Value( );
	iEvid = pVal->GetEvidence( );
	for( pDef->IndexToCoordinates( i = 0, veciCoords ); i != DSL_OUT_OF_RANGE;
		i = pDef->NextCoordinates( veciCoords ) ) {
		iLast = veciCoords[ veciCoords.GetSize( ) - 1 ];
		if( veciParents.NumItems( ) ) {
			if( iEvid == DSL_OUT_OF_RANGE ) {
				dProd = pVal->GetMatrix( )->Subscript( iLast );
				pVal->SetEvidence( iLast );
				m_SmileNet.UpdateBeliefs( ); }
			else if( iLast == iEvid )
				dProd = 1;
			else
				continue;

			for( j = 0; j < veciParents.NumItems( ); ++j )
				dProd *= m_SmileNet.GetNode( veciParents[ j ] )->Value( )->GetMatrix(
					)->Subscript( veciCoords[ j ] );
			if( iEvid == DSL_OUT_OF_RANGE ) {
				pVal->ClearEvidence( );
				m_SmileNet.UpdateBeliefs( ); } }
		else
			dProd = pVal->GetMatrix( )->Subscript( veciCoords[ 0 ] );

		pExpected->Subscript( veciCoords ) += dProd * iWeight; } }

bool CBayesNetSmile::Convert( CBayesNetPNL& BNPNL ) const {

	if( !m_fSmileNet )
		return false;

	return( ConvertGraph( BNPNL ) && ConvertCPTs( BNPNL ) ); }

bool CBayesNetSmileImpl::ConvertGraph( CBayesNetPNL& BNPNL ) const {
	CNodeType*						aNodeTypes;
	DSL_intArray					IntArrayNodes, IntArrayAdj;
	int								i, j, iSize;
	int**							aaiNeighborCounts;
	int*							aiNeighborCounts;
	int*							aiNodeAssociation;
	map<int,int>					mapSmilePNL, mapNodeTypes;
	ENeighborType**					aaeNeighborTypes;
	vector<int>						veciNeighborCounts;
	vector<ENeighborType>			veceNeighborTypes;
	map<int,int>::const_iterator	iterNodeTypes;
	DSL_node*						pNode;
	CGraph*							pGraph;

	m_SmileNet.GetAllNodes( IntArrayNodes );
	for( i = 0; i < IntArrayNodes.NumItems( ); ++i )
		mapSmilePNL[ IntArrayNodes[ i ] ] = i;

	aiNodeAssociation = new int[ IntArrayNodes.NumItems( ) ];
	for( i = 0; i < IntArrayNodes.NumItems( ); ++i ) {
		pNode = m_SmileNet.GetNode( IntArrayNodes[ i ] );
		if( IsGaussian( m_SmileNet ) )
			iSize = -1;
		else
			iSize = pNode->Definition( )->GetNumberOfOutcomes( );
		if( ( iterNodeTypes = mapNodeTypes.find( iSize ) ) == mapNodeTypes.end( ) ) {
			mapNodeTypes[ iSize ] = (int)mapNodeTypes.size( );
			iterNodeTypes = mapNodeTypes.find( iSize ); }
		aiNodeAssociation[ i ] = iterNodeTypes->second; }
	aNodeTypes = new CNodeType[ mapNodeTypes.size( ) ];
	for( iterNodeTypes = mapNodeTypes.begin( ); iterNodeTypes != mapNodeTypes.end( );
		++iterNodeTypes )
		aNodeTypes[ iterNodeTypes->second ].SetType( ( iterNodeTypes->first > 0 ),
			abs( iterNodeTypes->first ) );

	aiNeighborCounts = new int[ IntArrayNodes.NumItems( ) ];
	aaiNeighborCounts = new int*[ IntArrayNodes.NumItems( ) ];
	aaeNeighborTypes = new ENeighborType*[ IntArrayNodes.NumItems( ) ];
	for( i = 0; i < IntArrayNodes.NumItems( ); ++i ) {
		IntArrayAdj = m_SmileNet.GetParents( IntArrayNodes[ i ] );
		veciNeighborCounts.resize( iSize = aiNeighborCounts[ i ] = IntArrayAdj.NumItems( ) );
		veceNeighborTypes.resize( iSize );
		for( j = 0; j < IntArrayAdj.NumItems( ); ++j ) {
			veciNeighborCounts[ j ] = mapSmilePNL[ IntArrayAdj[ j ] ];
			veceNeighborTypes[ j ] = ntParent; }

		IntArrayAdj = m_SmileNet.GetChildren( IntArrayNodes[ i ] );
		veciNeighborCounts.resize( aiNeighborCounts[ i ] += IntArrayAdj.NumItems( ) );
		veceNeighborTypes.resize( aiNeighborCounts[ i ] );
		for( j = 0; j < IntArrayAdj.NumItems( ); ++j ) {
			veciNeighborCounts[ iSize + j ] = mapSmilePNL[ IntArrayAdj[ j ] ];
			veceNeighborTypes[ iSize + j ] = ntChild; }

		aaiNeighborCounts[ i ] = new int[ veciNeighborCounts.size( ) ];
		aaeNeighborTypes[ i ] = new ENeighborType[ veceNeighborTypes.size( ) ];
		for( j = 0; j < veciNeighborCounts.size( ); ++j ) {
			aaiNeighborCounts[ i ][ j ] = veciNeighborCounts[ j ];
			aaeNeighborTypes[ i ][ j ] = veceNeighborTypes[ j ]; } }

	if( BNPNL.m_pPNLNet )
		delete BNPNL.m_pPNLNet;
	pGraph = CGraph::Create( IntArrayNodes.NumItems( ), aiNeighborCounts,
		aaiNeighborCounts, aaeNeighborTypes );
	BNPNL.m_pPNLNet = CBNet::Create( IntArrayNodes.NumItems( ), (int)mapNodeTypes.size( ),
		aNodeTypes, aiNodeAssociation, pGraph );

	for( i = 0; i < IntArrayNodes.NumItems( ); ++i ) {
		delete[] aaiNeighborCounts[ i ];
		delete[] aaeNeighborTypes[ i ]; }
	delete[] aaiNeighborCounts;
	delete[] aaeNeighborTypes;
	delete[] aiNeighborCounts;
	delete[] aiNodeAssociation;
	delete[] aNodeTypes;

	return true; }

bool CBayesNetSmileImpl::ConvertCPTs( CBayesNetPNL& BNPNL ) const {
	int							i, j, k, iCPT;
	floatVector					vecdCPT;
	DSL_Dmatrix*				pCPT;
	DSL_intArray				IntArrayNodes, IntArrayCoords;
	DSL_node*					pNode;
	CNumericDenseMatrix<float>*	pMatrix;
	int							aiDims[ 2 ]	= { 1, 1 };
	float						dCPT;

	BNPNL.m_pPNLNet->AllocFactors( );
	m_SmileNet.GetAllNodes( IntArrayNodes );
	for( i = 0; i < IntArrayNodes.NumItems( ); ++i ) {
		pNode = m_SmileNet.GetNode( IntArrayNodes[ i ] );
		pCPT = pNode->Definition( )->GetMatrix( );
		if( IsGaussian( m_SmileNet ) ) {
			BNPNL.m_pPNLNet->AllocFactor( i );

			dCPT = (float)(*pCPT)[ iCPT = 0 ];
			pMatrix = CNumericDenseMatrix<float>::Create( 2, aiDims, &dCPT );
			BNPNL.m_pPNLNet->GetFactor( i )->AttachMatrix( pMatrix, matMean );

			dCPT = (float)(*pCPT)[ ++iCPT ];
			pMatrix = CNumericDenseMatrix<float>::Create( 2, aiDims, &dCPT );
			BNPNL.m_pPNLNet->GetFactor( i )->AttachMatrix( pMatrix, matCovariance );

			j = m_SmileNet.GetParents( IntArrayNodes[ i ] ).NumItems( );
			for( k = 0; k < j; ++k ) {
				dCPT = (float)(*pCPT)[ ++iCPT ];
				pMatrix = CNumericDenseMatrix<float>::Create( 2, aiDims, &dCPT );
				BNPNL.m_pPNLNet->GetFactor( i )->AttachMatrix( pMatrix, matWeights, k ); } }
		else {
			vecdCPT.resize( pCPT->GetSize( ) );
			for( j = 0; j < pCPT->GetSize( ); ++j )
				vecdCPT[ j ] = (float)(*pCPT)[ j ];
			BNPNL.m_pPNLNet->CreateTabularCPD( i, vecdCPT ); } }

	return true; }

vector<string> CBayesNetSmile::GetNodes( ) const {
	vector<string>	vecRet;
	size_t			i;

	if( m_fSmileNet )
		for( i = 0; i < m_SmileNet.GetNumberOfNodes( ); ++i )
			vecRet.push_back( m_SmileNet.GetNode( (int)i )->Info( ).Header( ).GetId( ) );

	return vecRet; }

bool CBayesNetSmile::IsContinuous( ) const {

	return CBayesNetSmileImpl::IsContinuous( ); }

bool CBayesNetSmileImpl::IsContinuous( ) const {

	return ( m_fSmileNet ? IsGaussian( m_SmileNet ) : false ); }

bool CBayesNetSmile::Evaluate( const IDataset* pData,
	vector<vector<float> >& vecvecdResults, bool fZero ) const {

	return CBayesNetSmileImpl::Evaluate( pData, NULL, &vecvecdResults, fZero ); }

bool CBayesNetSmile::Evaluate( const IDataset* pData,  CDat& DatOut, bool fZero ) const {

	return CBayesNetSmileImpl::Evaluate( pData, &DatOut, NULL, fZero ); }

bool CBayesNetSmileImpl::Evaluate( const IDataset* pData, CDat* pDatOut,
	vector<vector<float> >* pvecvecdOut, bool fZero ) const {
	size_t						i, j, k, iVal;
	DSL_intArray				IntArrayNodes;
	DSL_nodeValue*				pValue;
	string						strCur;
	map<string,float>			mapData;
	map<string,float>::iterator	iterDatum;

	if( !m_fSmileNet || IsContinuous( ) )
		return false;

	((CBayesNetSmileImpl*)this)->m_SmileNet.SetDefaultBNAlgorithm( DSL_ALG_BN_LAURITZEN );
	m_SmileNet.GetAllNodes( IntArrayNodes );
	for( i = 0; i < pData->GetGenes( ); ++i ) {
		if( !( i % 250 ) )
			g_CatBioUtils.notice( "CBayesNetSmile::Evaluate( %d ) %d/%d", fZero, i,
				pData->GetGenes( ) );
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j ) {
			if( !pData->IsExample( i, j ) )
				continue;
			if( m_fGroup ) {
				strCur = EncodeDatum( pData, i, j );
//cerr << (int)i << '\t' << (int)j << '\t' << strCur << endl;
				if( ( iterDatum = mapData.find( strCur ) ) != mapData.end( ) ) {
					if( pDatOut )
						pDatOut->Set( i, j, iterDatum->second );
					if( pvecvecdOut ) {
						pvecvecdOut->resize( pvecvecdOut->size( ) + 1 );
						(*pvecvecdOut)[ pvecvecdOut->size( ) - 1 ].push_back(
							iterDatum->second ); }
					continue; } }

			((CBayesNetSmileImpl*)this)->m_SmileNet.ClearAllEvidence( );
			for( k = 1; k < IntArrayNodes.NumItems( ); ++k ) {
				if( pData->IsHidden( k ) )
					continue;
				if( ( iVal = pData->GetDiscrete( i, j, k ) ) == -1 ) {
					if( fZero )
						iVal = 0;
					else
						continue; }
				m_SmileNet.GetNode( IntArrayNodes[ (int)k ] )->Value( )->SetEvidence(
					(int)iVal ); }

			((CBayesNetSmileImpl*)this)->m_SmileNet.UpdateBeliefs( );
			pValue = m_SmileNet.GetNode( IntArrayNodes[ 0 ] )->Value( );
			if( m_fGroup )
				mapData[ strCur ] = (float)(*pValue->GetMatrix( ))[ 0 ];
			if( pvecvecdOut ) {
				pvecvecdOut->resize( pvecvecdOut->size( ) + 1 );
				{
					vector<float>&	vecdCur	= (*pvecvecdOut)[ pvecvecdOut->size( ) - 1 ];

					for( k = 0; ( k + 1 ) < pValue->GetSize( ); ++k )
						vecdCur.push_back( (float)(*pValue->GetMatrix( ))[ (int)k ] );
				} }
			if( pDatOut )
				pDatOut->Set( i, j, (float)(*pValue->GetMatrix( ))[ 0 ] ); } }

	return true; }

const char	CBayesNetPNLImpl::c_szBN[]	= "bn";

CBayesNetPNL::CBayesNetPNL( bool fGroup ) : CBayesNetPNLImpl(fGroup) { }

CBayesNetPNLImpl::CBayesNetPNLImpl( bool fGroup ) : CBayesNetImpl(fGroup),
	m_pPNLNet(NULL) { }

CBayesNetPNLImpl::~CBayesNetPNLImpl( ) {

	if( m_pPNLNet )
		delete m_pPNLNet; }

bool CBayesNetPNL::Open( const char* szIn ) {
	CContextPersistence	ConPer;

	if( !ConPer.LoadXML( szIn ) )
		return false;
	if( m_pPNLNet )
		delete m_pPNLNet;
	return !!( m_pPNLNet = (CBNet*)ConPer.Get( c_szBN ) ); }

bool CBayesNetPNL::Save( const char* szOut ) const {
	CContextPersistence	ConPer;

	ConPer.Put( m_pPNLNet, c_szBN );
	return ConPer.SaveAsXML( szOut ); }

bool CBayesNetPNL::Learn( const IDataset* pData, size_t iIterations, bool fZero ) {
	CEMLearningEngineDumb*	pLearner;

	if( !m_pPNLNet )
		return false;

/*
	for( i = 0; i < pData->GetGenes( ); ++i )
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j ) {
			if( !pData->IsExample( i, j ) )
				continue;
			veciObserved.clear( );
			vecValues.clear( );
			for( k = 0; k < m_pPNLNet->GetNumberOfNodes( ); ++k ) {
				if( pData->IsHidden( k ) )
					continue;
				if( IsContinuous( ) ) {
					if( CMeta::IsNaN( d = pData->GetContinuous( i, j, k ) ) )
						continue;
					vecValues.resize( vecValues.size( ) + 1 );
					vecValues[ vecValues.size( ) - 1 ].SetFlt( d ); }
				else {
					if( ( iVal = pData->GetDiscrete( i, j, k ) ) == -1 )
						continue;
					vecValues.resize( vecValues.size( ) + 1 );
					vecValues[ vecValues.size( ) - 1 ].SetInt( (int)iVal ); }
				veciObserved.push_back( (int)k ); }
			if( veciObserved.empty( ) ) {
				g_CatBioUtils.error( "CBayesNetPNL::Learn( %g, %d ) found no evidence for %s, %s",
					dTolerance, iIterations, pData->GetGene( i ).c_str( ),
					pData->GetGene( j ).c_str( ) );
				return false; }
			vecEvidence.push_back( CEvidence::Create( m_pPNLNet, veciObserved, vecValues ) ); }
*/

	pLearner = CEMLearningEngineDumb::Create( m_pPNLNet );
	pLearner->SetMaxIterEM( (int)iIterations );
//	pLearner->SetData( (int)vecEvidence.size( ), &vecEvidence.front( ) );
	pLearner->Learn( pData, fZero );

	delete pLearner;
//	for( i = 0; i < vecEvidence.size( ); ++i )
//		delete vecEvidence[ i ];

	return true; }

vector<string> CBayesNetPNL::GetNodes( ) const {
	vector<string>	vecstrRet;

	return vecstrRet; }

bool CBayesNetPNL::IsContinuous( ) const {

	return CBayesNetPNLImpl::IsContinuous( ); }

bool CBayesNetPNLImpl::IsContinuous( ) const {

	return ( m_pPNLNet ? !m_pPNLNet->GetNodeType( 0 )->IsDiscrete( ) : false ); }

bool CBayesNetPNL::Evaluate( const IDataset* pData,
	vector<vector<float> >& vecvecdResults, bool fZero ) const {

	return CBayesNetPNLImpl::Evaluate( pData, NULL, &vecvecdResults, fZero ); }

bool CBayesNetPNL::Evaluate( const IDataset* pData, CDat& DatOut, bool fZero ) const {

	return CBayesNetPNLImpl::Evaluate( pData, &DatOut, NULL, fZero ); }

bool CBayesNetPNLImpl::Evaluate( const IDataset* pData, CDat* pDatOut,
	vector<vector<float> >* pvecvecdOut, bool fZero ) const {
	CInfEngine*							pInferrer;
	size_t								i, j, k, l, iVal;
	CEvidence*							pEvidence;
	intVector							veciObserved;
	valueVector							vecValues;
	int									iNode;
	const CFactor*						pFactor;
	const CMatrix<float>*				pMatrix;
	CMatrixIterator<float>*				pIter;
	float								d;
	const float*						pd;
	vector<float>*						pvecdCur;
	map<string,float>					mapData;
	map<string,float>::const_iterator	iterDatum;
	string								strCur;

	if( !m_pPNLNet )
		return false;

	pvecdCur = NULL;
	pInferrer = CJtreeInfEngine::Create( m_pPNLNet );
	iNode = 0;
	for( i = 0; i < pData->GetGenes( ); ++i ) {
		if( !( i % 250 ) )
			g_CatBioUtils.notice( "CBayesNetPNL::Evaluate( %d ) %d/%d", fZero, i,
				pData->GetGenes( ) );
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j ) {
			if( !pData->IsExample( i, j ) )
				continue;
			if( m_fGroup ) {
				strCur = EncodeDatum( pData, i, j );
				if( ( iterDatum = mapData.find( strCur ) ) != mapData.end( ) ) {
					if( pDatOut )
						pDatOut->Set( i, j, iterDatum->second );
					if( pvecvecdOut ) {
						pvecvecdOut->resize( pvecvecdOut->size( ) + 1 );
						(*pvecvecdOut)[ pvecvecdOut->size( ) - 1 ].push_back(
							iterDatum->second ); }
					continue; } }

			veciObserved.clear( );
			vecValues.clear( );
			for( k = 1; k < m_pPNLNet->GetNumberOfNodes( ); ++k ) {
				if( pData->IsHidden( k ) )
					continue;
				if( IsContinuous( ) ) {
					if( CMeta::IsNaN( d = pData->GetContinuous( i, j, k ) ) ) {
						if( fZero )
							d = 0;
						else
							continue; }
					vecValues.resize( vecValues.size( ) + 1 );
					vecValues[ vecValues.size( ) - 1 ].SetFlt( d ); }
				else {
					if( ( iVal = pData->GetDiscrete( i, j, k ) ) == -1 ) {
						if( fZero )
							iVal = 0;
						else
							continue; }
					vecValues.resize( vecValues.size( ) + 1 );
					vecValues[ vecValues.size( ) - 1 ].SetInt( (int)iVal ); }
				veciObserved.push_back( (int)k ); }

			pEvidence = CEvidence::Create( m_pPNLNet, veciObserved, vecValues );
			pInferrer->EnterEvidence( pEvidence );
			pInferrer->MarginalNodes( &iNode, 1 );
			pFactor = pInferrer->GetQueryJPD( );
			delete pEvidence;

			if( pvecvecdOut ) {
				pvecvecdOut->resize( pvecvecdOut->size( ) + 1 );
				pvecdCur = &(*pvecvecdOut)[ pvecvecdOut->size( ) - 1 ]; }
			if( pFactor->GetDistributionType( ) == dtTabular ) {
				pMatrix = pFactor->GetMatrix( matTable );
				pIter = pMatrix->InitIterator( );
				while( true ) {
					pd = pMatrix->Value( pIter );
					pMatrix->Next( pIter );
					if( !pMatrix->IsValueHere( pIter ) )
						break;
					mapData[ strCur ] = *pd;
					if( pvecdCur )
						pvecdCur->push_back( *pd );
					if( pDatOut ) {
						pDatOut->Set( i, j, *pd );
						break; } }
				delete pIter; }
			else {
				pMatrix = pFactor->GetMatrix( matMean );
				for( pIter = pMatrix->InitIterator( ); pMatrix->IsValueHere( pIter );
					pMatrix->Next( pIter ) ) {
					mapData[ strCur ] = *pMatrix->Value( pIter );
					if( pvecdCur )
						pvecdCur->push_back( *pMatrix->Value( pIter ) );
					if( pDatOut ) {
						pDatOut->Set( i, j, *pMatrix->Value( pIter ) );
						break; } }
				delete pIter;
				if( !pvecdCur )
					break;
				pMatrix = pFactor->GetMatrix( matCovariance );
				for( pIter = pMatrix->InitIterator( ); pMatrix->IsValueHere( pIter );
					pMatrix->Next( pIter ) )
					pvecdCur->push_back( *pMatrix->Value( pIter ) );
				delete pIter;

				veciObserved.clear( );
				pFactor->GetDomain( &veciObserved );
				for( l = k = 0; k < veciObserved.size( ); ++k )
					l += m_pPNLNet->GetGraph( )->GetNumberOfParents( veciObserved[ k ] );
				for( k = 0; k < l; ++k ) {
					pMatrix = pFactor->GetMatrix( matWeights, (int)k );
					for( pIter = pMatrix->InitIterator( ); pMatrix->IsValueHere( pIter );
						pMatrix->Next( pIter ) ) 
						pvecdCur->push_back( *pMatrix->Value( pIter ) );
					delete pIter; } } } }

	delete pInferrer;
	return true; }

}

#endif // BAYESIAN_NETWORKS

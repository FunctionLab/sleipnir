#include "stdafx.h"
#include "bayesnet.h"
#include "dat.h"
#include "dataset.h"
#include "meta.h"

namespace libBioUtils {

const char	CBayesNetSmileImpl::c_szGaussian[]	= "gaussian";

bool CBayesNetSmileImpl::IsGaussian( const DSL_network& BayesNet ) {
	int	i;

	if( ( i = ((DSL_network&)BayesNet).UserProperties( ).FindProperty( c_szGaussian ) ) < 0 )
		return false;

	return !!atoi( ((DSL_network&)BayesNet).UserProperties( ).GetPropertyValue( i ) ); }

bool CBayesNetSmileImpl::IsNaive( const DSL_network& BayesNet ) {
	size_t	i;

	{
		const DSL_intArray&	veciParents	= ((DSL_network&)BayesNet).GetNode( (int)0 )->Parents( );

		if( veciParents.NumItems( ) != 0 )
			return false;
	}
	for( i = 1; i < BayesNet.GetNumberOfNodes( ); ++i ) {
		const DSL_intArray&	veciParents	= ((DSL_network&)BayesNet).GetNode( (int)i )->Parents( );

		if( ( veciParents.NumItems( ) > 1 ) || ( veciParents[ 0 ] != 0 ) )
			return false; }

	return true; }

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
	size_t					i, j, iIter, iDatum;
	string					strCur;
	TMapData				mapData;
	TMapData::iterator		iterDatum;
	DSL_Dmatrix*			pMat;
	vector<DSL_Dmatrix*>	vecpExpected;
	DSL_intArray			veciCoords;
	vector<bool>			vecfHidden;

	vecfHidden.resize( pData->GetExperiments( ) );
	for( i = 0; i < vecfHidden.size( ); ++i )
		vecfHidden[ i ] = pData->IsHidden( i );
	EncodeData( pData, mapData );
	vecpExpected.resize( m_SmileNet.GetNumberOfNodes( ) );
	for( i = 0; i < vecpExpected.size( ); ++i )
		vecpExpected[ i ] = new DSL_Dmatrix( *m_SmileNet.GetNode( (int)i )->Definition(
			)->GetMatrix( ) );
	for( iIter = 0; iIter < iIterations; ++iIter ) {
		for( iDatum = i = 0; i < vecpExpected.size( ); ++i )
			vecpExpected[ i ]->FillWith( 0 );
		for( iterDatum = mapData.begin( ); iterDatum != mapData.end( ); ++iterDatum ) {
			if( !( iDatum++ % 50 ) )
				g_CatBioUtils.notice( "CBayesNetSmile::LearnGrouped( %d, %d ) iteration %d, datum %d/%d",
					iIterations, fZero, iIter, ( iDatum - 1 ), mapData.size( ) );
			FillCPTs( vecfHidden, iterDatum->first, fZero, true );
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

bool CBayesNetSmileImpl::FillCPTs( const IDataset* pData, size_t iOne, size_t iTwo, bool fZero, bool fLearn ) {
	size_t	i, iVal;

	if( !pData->IsExample( iOne, iTwo ) || ( fLearn && ( pData->GetDiscrete( iOne, iTwo, 0 ) == -1 ) ) )
		return false;

	m_SmileNet.ClearAllEvidence( );
	for( i = fLearn ? 0 : 1; i < m_SmileNet.GetNumberOfNodes( ); ++i ) {
		if( pData->IsHidden( i ) )
			continue;
		if( ( iVal = pData->GetDiscrete( iOne, iTwo, i ) ) == -1 ) {
			if( fZero )
				iVal = 0;
			else
				continue; }
		m_SmileNet.GetNode( (int)i )->Value( )->SetEvidence( (int)iVal ); }

	return true; }

bool CBayesNetSmileImpl::FillCPTs( const vector<bool>& vecfHidden, const string& strDatum, bool fZero,
	bool fLearn ) {
	size_t	i, iVal;

	if( fLearn && !IsAnswer( strDatum ) )
		return false;

	m_SmileNet.ClearAllEvidence( );
	for( i = fLearn ? 0 : 1; i < m_SmileNet.GetNumberOfNodes( ); ++i ) {
		if( vecfHidden[ i ] )
			continue;
		if( strDatum[ i ] == c_cMissing ) {
			if( !fZero )
				continue;
			iVal = 0; }
		else
			iVal = strDatum[ i ] - c_cBase;
		m_SmileNet.GetNode( (int)i )->Value( )->SetEvidence( (int)iVal ); }

	return true; }

bool CBayesNetSmileImpl::FillCPTs( const vector<bool>& vecfHidden, const vector<unsigned char>& vecbDatum,
	bool fZero, bool fLearn ) {
	size_t	i, iVal;

	if( fLearn && !vecbDatum[ 0 ] )
		return false;

	m_SmileNet.ClearAllEvidence( );
	for( i = fLearn ? 0 : 1; i < m_SmileNet.GetNumberOfNodes( ); ++i ) {
		if( vecfHidden[ i ] )
			continue;
		if( !vecbDatum[ i ] ) {
			if( !fZero )
				continue;
			iVal = 0; }
		else
			iVal = vecbDatum[ i ] - 1;
		m_SmileNet.GetNode( (int)i )->Value( )->SetEvidence( (int)iVal ); }

	return true; }

bool CBayesNetSmileImpl::LearnUngrouped( const IDataset* pData, size_t iIterations,
	bool fZero ) {
	size_t					iIter, i, j, k;
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
				g_CatBioUtils.notice( "CBayesNetSmile::LearnUngrouped( %d, %d ) iteration %d, gene %d/%d",
					iIterations, fZero, iIter, i, pData->GetGenes( ) );
			for( j = ( i + 1 ); j < pData->GetGenes( ); ++j ) {
				if( !FillCPTs( pData, i, j, fZero, true ) )
					continue;
				m_SmileNet.UpdateBeliefs( );

				for( k = 0; k < m_SmileNet.GetNumberOfNodes( ); ++k )
					LearnExpected( m_SmileNet.GetNode( (int)k ), vecpExpected[ k ] ); } }
		for( i = 0; i < m_SmileNet.GetNumberOfNodes( ); ++i ) {
			pMat = m_SmileNet.GetNode( (int)i )->Definition( )->GetMatrix( );
			for( pMat->IndexToCoordinates( (int)( j = 0 ), veciCoords );
				j != DSL_OUT_OF_RANGE; j = pMat->NextCoordinates( veciCoords ) )
				pMat->Subscript( veciCoords ) = vecpExpected[ i ]->Subscript( veciCoords );
			pMat->Normalize( ); } }
	for( i = 0; i < vecpExpected.size( ); ++i )
		delete vecpExpected[ i ];

	return true; }

bool CBayesNetSmileImpl::IsNaive( ) const {

	return ( m_fSmileNet ? CBayesNetSmileImpl::IsNaive( m_SmileNet ) : false ); }

bool CBayesNetSmile::Learn( const IDataset* pData, size_t iIterations, bool fZero, bool fELR ) {

	if( fELR )
		return LearnELR( pData, iIterations, fZero );
	if( IsNaive( ) )
		return LearnNaive( pData, fZero );

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
			j = (int)mapNodeTypes.size( );
			mapNodeTypes[ iSize ] = j;
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

unsigned char CBayesNetSmile::GetValues( size_t iNode ) const {

	return m_SmileNet.GetNode( (int)iNode )->Definition( )->GetNumberOfOutcomes( ); }

bool CBayesNetSmile::IsContinuous( size_t ) const {

	return IsContinuous( ); }

bool CBayesNetSmile::IsContinuous( ) const {

	return CBayesNetSmileImpl::IsContinuous( ); }

bool CBayesNetSmileImpl::IsContinuous( ) const {

	return ( m_fSmileNet ? IsGaussian( m_SmileNet ) : false ); }

bool CBayesNetSmile::Evaluate( const IDataset* pData,
	vector<vector<float> >& vecvecdResults, bool fZero ) const {

	return CBayesNetSmileImpl::Evaluate( pData, NULL, &vecvecdResults, fZero ); }

bool CBayesNetSmile::Evaluate( const IDataset* pData, CDat& DatOut, bool fZero ) const {

	return CBayesNetSmileImpl::Evaluate( pData, &DatOut, NULL, fZero ); }

bool CBayesNetSmileImpl::Evaluate( const IDataset* pData, CDat* pDatOut,
	vector<vector<float> >* pvecvecdOut, bool fZero ) const {
	size_t						i, j, k;
	DSL_nodeValue*				pValue;
	string						strCur;
	map<string,float>			mapData;
	map<string,float>::iterator	iterDatum;
	vector<bool>				vecfHidden;

	if( !m_fSmileNet || IsContinuous( ) )
		return false;

	vecfHidden.resize( pData->GetExperiments( ) );
	for( i = 0; i < vecfHidden.size( ); ++i )
		vecfHidden[ i ] = pData->IsHidden( i );
	for( i = 0; i < pData->GetGenes( ); ++i ) {
		if( !( i % 250 ) )
			g_CatBioUtils.notice( "CBayesNetSmile::Evaluate( %d ) %d/%d", fZero, i,
				pData->GetGenes( ) );
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j ) {
			if( !pData->IsExample( i, j ) )
				continue;
			strCur = EncodeDatum( pData, i, j );
			if( m_fGroup && ( ( iterDatum = mapData.find( strCur ) ) != mapData.end( ) ) ) {
				if( pDatOut )
					pDatOut->Set( i, j, iterDatum->second );
				if( pvecvecdOut ) {
					pvecvecdOut->resize( pvecvecdOut->size( ) + 1 );
					(*pvecvecdOut)[ pvecvecdOut->size( ) - 1 ].push_back(
						iterDatum->second ); }
				continue; }

			((CBayesNetSmileImpl*)this)->FillCPTs( vecfHidden, strCur, fZero, false );
			((CBayesNetSmileImpl*)this)->m_SmileNet.UpdateBeliefs( );
			pValue = m_SmileNet.GetNode( 0 )->Value( );
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

bool CBayesNetSmile::Evaluate( const vector<unsigned char>& vecbDatum, vector<float>& vecdOut, bool fZero ) const {
	vector<bool>	vecfHidden;
	DSL_nodeValue*	pValue;
	size_t			i;

	if( !m_fSmileNet || IsContinuous( ) )
		return false;

	vecfHidden.resize( vecbDatum.size( ) );
	for( i = 0; i < vecfHidden.size( ); ++i )
		vecfHidden[ i ] = false;
	((CBayesNetSmile*)this)->FillCPTs( vecfHidden, vecbDatum, fZero, false );
	((CBayesNetSmile*)this)->m_SmileNet.UpdateBeliefs( );
	pValue = m_SmileNet.GetNode( 0 )->Value( );
	for( i = 0; ( i + 1 ) < pValue->GetSize( ); ++i )
		vecdOut.push_back( (float)(*pValue->GetMatrix( ))[ (int)i ] );

	return true; }

void CBayesNetSmile::Randomize( ) {
	int	i;

	if( !m_fSmileNet )
		return;

	for( i = m_SmileNet.GetFirstNode( ); i != DSL_OUT_OF_RANGE;
		i = m_SmileNet.GetNextNode( i ) )
		Randomize( i ); }

void CBayesNetSmile::Randomize( size_t iNode ) {
	DSL_Dmatrix*	pMat;

	if( !m_fSmileNet )
		return;

	pMat = m_SmileNet.GetNode( (int)iNode )->Definition( )->GetMatrix( );

	{
		DSL_sysCoordinates	Coords( *pMat );

		Coords.GoFirst( );
		do
			Coords.CheckedValue( ) = (float)rand( ) / RAND_MAX;
		while( Coords.Next( ) != DSL_OUT_OF_RANGE );
	}

	pMat->Normalize( ); }

void CBayesNetSmile::Reverse( size_t iNode ) {
	int				iCoords;
	DSL_Dmatrix*	pMat;

	if( !m_fSmileNet )
		return;

	pMat = m_SmileNet.GetNode( (int)iNode )->Definition( )->GetMatrix( );
	{
		DSL_sysCoordinates	Coords( *pMat );

		iCoords = pMat->GetSizeOfDimension( pMat->GetLastDimension( ) );
		Coords.GoFirst( );
		do {
			DSL_intArray	veciCoords	= Coords.Coordinates( );
			int				iCoord;
			double			d;

			iCoord = veciCoords[ veciCoords.GetSize( ) - 1 ];
			if( iCoord >= ( iCoords / 2 ) )
				continue;
			d = Coords.CheckedValue( );
			veciCoords[ veciCoords.GetSize( ) - 1 ] = iCoords - iCoord - 1;
			Coords.CheckedValue( ) = (*pMat)[ veciCoords ];
			(*pMat)[ veciCoords ] = d; }
		while( Coords.Next( ) != DSL_OUT_OF_RANGE );
	} }

bool CBayesNetSmileImpl::LearnNaive( const IDataset* pData, bool fZero ) {
	vector<vector<size_t> >	vecveciCounts;
	size_t					i, j, k, iAnswer, iAnswers, iVal, iCount;
	DSL_nodeDefinition*		pDef;
	DSL_Dmatrix*			pMat;
	DSL_intArray			veciCoords;

	vecveciCounts.resize( m_SmileNet.GetNumberOfNodes( ) );
	iAnswers = m_SmileNet.GetNode( 0 )->Definition( )->GetNumberOfOutcomes( );
	vecveciCounts[ 0 ].resize( iAnswers );
	for( i = 1; i < vecveciCounts.size( ); ++i )
		vecveciCounts[ i ].resize( iAnswers *
			m_SmileNet.GetNode( (int)i )->Definition( )->GetNumberOfOutcomes( ) );
	for( i = 0; i < pData->GetGenes( ); ++i )
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j )
			if( pData->IsExample( i, j ) && ( ( iAnswer = pData->GetDiscrete( i, j, 0 ) ) != -1 ) ) {
				vecveciCounts[ 0 ][ iAnswer ]++;
				for( k = 1; k < pData->GetExperiments( ); ++k ) {
					if( ( iVal = pData->GetDiscrete( i, j, k ) ) == -1 ) {
						if( fZero )
							iVal = 0;
						else
							continue; }
					vecveciCounts[ k ][ ( iVal * iAnswers ) + iAnswer ]++; } }

	pMat = m_SmileNet.GetNode( 0 )->Definition( )->GetMatrix( );
	for( i = 0; i < iAnswers; ++i )
		(*pMat)[ (int)i ] = vecveciCounts[ 0 ][ (int)i ] + 1;
	pMat->Normalize( );
	for( i = 1; i < vecveciCounts.size( ); ++i ) {
		pDef = m_SmileNet.GetNode( (int)i )->Definition( );
		pMat = pDef->GetMatrix( );
		pMat->IndexToCoordinates( 0, veciCoords );
		for( j = 0; j < iAnswers; ++j ) {
			veciCoords[ 0 ] = (int)j;
			for( k = 0; k < pDef->GetNumberOfOutcomes( ); ++k ) {
				veciCoords[ 1 ] = (int)k;
				if( !( iCount = vecveciCounts[ i ][ ( k * iAnswers ) + j ] ) ) {
//* Indicates that unseen examples provide no information
					for( size_t m = 0; m < iAnswers; ++m )
						if( m != j )
							iCount += vecveciCounts[ i ][ ( k * iAnswers ) + m ];
					iCount /= iAnswers - 1; }
				if( !iCount ) {
//*/
					iCount = 1; }
				(*pMat)[ veciCoords ] = iCount; } }
		pMat->Normalize( ); }

	return true; }

}

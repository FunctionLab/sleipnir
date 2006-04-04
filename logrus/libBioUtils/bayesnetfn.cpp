#include "stdafx.h"
#include "bayesnet.h"
#include "statistics.h"

namespace libBioUtils {

const char	CBayesNetFNNode::c_szType[]	= "type";

CBayesNetFNNode* CBayesNetFNNode::Open( DSL_node* pNode ) {
	static const CBayesNetFNNode*	c_pDefault	= new CBayesNetFNNodeDiscrete( );
	static const CBayesNetFNNode*	c_apTypes[]	= {
		new CBayesNetFNNodeDiscrete( ),
		new CBayesNetFNNodeGaussian( ),
		NULL };
	CBayesNetFNNode*	pRet;
	int					i, j, iDim;
	const char*			szType;
	DSL_Dmatrix*		pMatrix;

	if( ( i = pNode->Info( ).UserProperties( ).FindProperty( c_szType ) ) < 0 )
		pRet = c_pDefault->New( pNode );
	else {
		pRet = NULL;
		szType = pNode->Info( ).UserProperties( ).GetPropertyValue( i );
		for( i = 0; c_apTypes[ i ]; ++i )
			if( !strcmp( szType, c_apTypes[ i ]->GetType( ) ) ) {
				pRet = c_apTypes[ i ]->New( pNode );
				break; } }
	pMatrix = pNode->Definition( )->GetMatrix( );
	if( pRet ) {
		DSL_intArray&	veciDims	= pMatrix->GetDimensions( );

		iDim = ( veciDims.GetSize( ) > 1 ) ? 1 : 0;
		pRet->m_strName = pNode->Info( ).Header( ).GetId( );
		pRet->m_Params.Initialize( iDim ? veciDims[ 1 ] : 1, veciDims[ iDim ] );
		for( i = 0; i < pRet->m_Params.GetRows( ); ++i ) {
			if( iDim )
				veciDims[ 0 ] = i;
			for( j = 0; j < pRet->m_Params.GetColumns( ); ++j ) {
				veciDims[ iDim ] = j;
				pRet->m_Params.Set( i, j, (float)(*pMatrix)[ veciDims ] ); } } }

	return pRet; }

const string& CBayesNetFNNode::GetName( ) const {

	return m_strName; }

unsigned char CBayesNetFNNode::GetParameters( ) const {

	return (unsigned char)m_Params.GetColumns( ); }

void CBayesNetFNNode::Reverse( ) {
	size_t			iCol, iRow;
	vector<float>	vecdRow;

	vecdRow.resize( m_Params.GetColumns( ) );
	for( iRow = 0; iRow < m_Params.GetRows( ); ++iRow ) {
		for( iCol = 0; iCol < vecdRow.size( ); ++iCol )
			vecdRow[ iCol ] = m_Params.Get( iRow, vecdRow.size( ) - iCol - 1 );
		for( iCol = 0; iCol < vecdRow.size( ); ++iCol )
			m_Params.Set( iRow, iCol, vecdRow[ iCol ] ); } }

bool CBayesNetFNNode::Save( DSL_node* pNode ) const {
	int				i;
	size_t			iCol, iRow;
	DSL_Dmatrix*	pMatrix		= pNode->Definition( )->GetMatrix( );
	DSL_intArray&	veciDims	= pMatrix->GetDimensions( );

	if( m_strName != pNode->Info( ).Header( ).GetId( ) )
		return false;
	if( ( i = pNode->Info( ).UserProperties( ).FindProperty( c_szType ) ) < 0 )
		pNode->Info( ).UserProperties( ).AddProperty( c_szType, GetType( ) );
	else
		pNode->Info( ).UserProperties( ).ChangePropertyValue( i, GetType( ) );

	i = ( veciDims.GetSize( ) > 1 ) ? 1 : 0;
	for( iRow = 0; iRow < m_Params.GetRows( ); ++iRow ) {
		if( i )
			veciDims[ 0 ] = (int)iRow;
		for( iCol = 0; iCol < m_Params.GetColumns( ); ++iCol ) {
			veciDims[ i ] = (int)iCol;
			(*pMatrix)[ veciDims ] = m_Params.Get( iRow, iCol ); } }

	return true; }

bool CBayesNetFNNode::Learn( const vector<size_t>& veciCounts ) {
	vector<float>	vecdRow;
	size_t			iRow, iCol, iSum;

	if( veciCounts.size( ) != m_Params.GetColumns( ) )
		return false;

	iSum = veciCounts.size( );
	for( iCol = 0; iCol < veciCounts.size( ); ++iCol )
		iSum += veciCounts[ iCol ];
	vecdRow.resize( veciCounts.size( ) );
	for( iCol = 0; iCol < vecdRow.size( ); ++iCol )
		vecdRow[ iCol ] = ( veciCounts[ iCol ] + 1.0f ) / iSum;
	for( iRow = 0; iRow < m_Params.GetRows( ); ++iRow )
		for( iCol = 0; iCol < m_Params.GetColumns( ); ++iCol )
			m_Params.Set( iRow, iCol, vecdRow[ iCol ] );

	return true; }

void CBayesNetFNNodeDiscrete::Randomize( ) {
	size_t		iCol, iRow, iSum;
	vector<int>	veciRow;

	veciRow.resize( m_Params.GetColumns( ) );
	for( iRow = 0; iRow < m_Params.GetRows( ); ++iRow ) {
		for( iSum = iCol = 0; iCol < veciRow.size( ); ++iCol )
			iSum += ( veciRow[ iCol ] = rand( ) );
		for( iCol = 0; iCol < veciRow.size( ); ++iCol )
			m_Params.Set( iRow, iCol, (float)veciRow[ iCol ] / iSum ); } }

bool CBayesNetFNNodeDiscrete::Learn( const IDataset* pData, size_t iNode, bool fZero ) {
	CFullMatrix<size_t>	MatCounts;
	vector<size_t>		veciSums;
	size_t				i, j, iAnswer, iVal;

	veciSums.resize( m_Params.GetRows( ) );
	MatCounts.Initialize( m_Params.GetRows( ), m_Params.GetColumns( ) );
	for( i = 0; i < MatCounts.GetRows( ); ++i ) {
		veciSums[ i ] = (unsigned int)MatCounts.GetColumns( );
		for( j = 0; j < MatCounts.GetColumns( ); ++j )
			MatCounts.Set( i, j, 1 ); }
	for( i = 0; i < pData->GetGenes( ); ++i )
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j )
			if( pData->IsExample( i, j ) && ( ( iAnswer = pData->GetDiscrete( i, j, 0 ) ) != -1 ) ) {
				if( ( iVal = pData->GetDiscrete( i, j, iNode ) ) == -1 ) {
					if( !fZero )
						continue;
					iVal = 0; }
				if( ( iAnswer >= MatCounts.GetRows( ) ) || ( iVal >= MatCounts.GetColumns( ) ) )
					return false;
				MatCounts.Get( iAnswer, iVal )++;
				veciSums[ iAnswer ]++; }

	for( i = 0; i < m_Params.GetRows( ); ++i )
		for( j = 0; j < m_Params.GetColumns( ); ++j )
			m_Params.Set( i, j, (float)MatCounts.Get( i, j ) / veciSums[ i ] );

	return true; }

bool CBayesNetFNNodeDiscrete::Evaluate( float dValue, vector<float>& vecdOut ) const {
	size_t	i;

	if( CMeta::IsNaN( dValue ) ) {
		if( m_Params.GetRows( ) != 1 )
			return false;
		vecdOut.resize( m_Params.GetColumns( ) );
		for( i = 0; i < vecdOut.size( ); ++i )
			vecdOut[ i ] = m_Params.Get( 0, i );
		return true; }

	vecdOut.resize( m_Params.GetRows( ) );
	for( i = 0; i < vecdOut.size( ); ++i )
		vecdOut[ i ] = m_Params.Get( i, (unsigned char)dValue );

	return true; }

void CBayesNetFNNodeGaussian::Randomize( ) {
	size_t	i, j;

	for( i = 0; i < m_Params.GetRows( ); ++i )
		for( j = 0; j < m_Params.GetColumns( ); ++j )
			m_Params.Set( i, c_iMu, (float)( rand( ) - ( RAND_MAX / 2 ) ) ); }

bool CBayesNetFNNodeGaussian::Learn( const IDataset* pData, size_t iNode, bool fZero ) {
	float	d, dVal;
	size_t	i, j, iAnswer, iCount;

	for( i = 0; i < m_Params.GetRows( ); ++i )
		for( j = 0; j < m_Params.GetColumns( ); ++j )
			m_Params.Set( i, j, 0 );
	for( iCount = i = 0; i < pData->GetGenes( ); ++i )
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j )
			if( pData->IsExample( i, j ) && ( ( iAnswer = pData->GetDiscrete( i, j, 0 ) ) != -1 ) ) {
				if( CMeta::IsNaN( dVal = pData->GetContinuous( i, j, iNode ) ) ) {
					if( !fZero )
						continue;
					dVal = 0; }
				if( iAnswer >= m_Params.GetRows( ) )
					return false;
				iCount++;
				m_Params.Get( iAnswer, c_iMu ) += dVal;
				m_Params.Get( iAnswer, c_iSigma ) += dVal * dVal; }
	if( !iCount )
		iCount = 1;
	for( i = 0; i < m_Params.GetRows( ); ++i ) {
		d = ( m_Params.Get( i, c_iMu ) /= iCount );
		d *= d;
		dVal = m_Params.Get( i, c_iSigma );
		dVal = ( dVal == d ) ? 1 : sqrt( ( dVal - d ) / ( iCount - 1 ) );
		m_Params.Set( i, c_iSigma, dVal ); }

	return true; }

bool CBayesNetFNNodeGaussian::Evaluate( float dValue, vector<float>& vecdOut ) const {
	float	dSum;
	size_t	i;

	vecdOut.resize( m_Params.GetRows( ) );
	dSum = 0;
	for( i = 0; i < vecdOut.size( ); ++i )
		dSum += ( vecdOut[ i ] = (float)CStatistics::NormalPDF( dValue, m_Params.Get( i, c_iMu ),
			m_Params.Get( i, c_iSigma ) ) );
	for( i = 0; i < vecdOut.size( ); ++i )
		vecdOut[ i ] /= dSum;

	return true; }

CBayesNetFNImpl::CBayesNetFNImpl( ) : CBayesNetImpl( true ), m_iNodes(0), m_apNodes(NULL), m_fSmileNet(false) { }

CBayesNetFNImpl::~CBayesNetFNImpl( ) {

	Reset( ); }

void CBayesNetFNImpl::Reset( ) {
	size_t	i;

	m_fSmileNet = false;
	if( m_apNodes ) {
		for( i = 0; i < m_iNodes; ++i )
			if( m_apNodes[ i ] )
				delete m_apNodes[ i ];
		delete[] m_apNodes; } }

bool CBayesNetFN::Open( const char* szFile ) {
	size_t	i;

	Reset( );
	if( !( ( m_fSmileNet = !m_SmileNet.ReadFile( szFile ) ) && CBayesNetSmileImpl::IsNaive( m_SmileNet ) ) )
		return false;

	m_iNodes = m_SmileNet.GetNumberOfNodes( );
	m_apNodes = new CBayesNetFNNode*[ m_iNodes ];
	for( i = 0; i < m_iNodes; ++i )
		if( !( m_apNodes[ i ] = CBayesNetFNNode::Open( m_SmileNet.GetNode( (int)i ) ) ) )
			return false;

	return true; }

bool CBayesNetFN::Save( const char* szFile ) const {
	size_t	i;

	if( !m_fSmileNet )
		return false;

	for( i = 0; i < m_iNodes; ++i )
		if( !m_apNodes[ i ]->Save( m_SmileNet.GetNode( (int)i ) ) )
			return false;

	return !((CBayesNetFN*)this)->m_SmileNet.WriteFile( szFile ); }

bool CBayesNetFN::Learn( const IDataset* pData, size_t iIterations, bool fZero, bool fELR ) {
	size_t			i, j, iAnswer;
	vector<size_t>	veciCounts;

	if( !m_iNodes )
		return false;

	veciCounts.resize( m_apNodes[ 0 ]->GetParameters( ) );
	for( i = 0; i < pData->GetGenes( ); ++i )
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j )
			if( pData->IsExample( i, j ) && ( ( iAnswer = pData->GetDiscrete( i, j, 0 ) ) != -1 ) )
				veciCounts[ iAnswer ]++;

	if( !m_apNodes[ 0 ]->Learn( veciCounts ) )
		return false;
	for( i = 1; i < m_iNodes; ++i )
		if( !m_apNodes[ i ]->Learn( pData, i, fZero ) )
			return false;

	return true; }

bool CBayesNetFN::Evaluate( const IDataset* pData, vector<vector<float> >& vecvecdOut, bool fZero ) const {

	return CBayesNetFNImpl::Evaluate( pData, NULL, &vecvecdOut, fZero ); }

bool CBayesNetFN::Evaluate( const IDataset* pData, CDat& DatOut, bool fZero ) const {

	return CBayesNetFNImpl::Evaluate( pData, &DatOut, NULL, fZero ); }

bool CBayesNetFNImpl::Evaluate( const IDataset* pData, CDat* pDatOut, vector<vector<float> >* pvecvecdOut,
	bool fZero ) const {
	size_t						i, j, k;
	string						strCur;
	vector<float>				vecdCur;
	map<string,float>			mapData;
	map<string,float>::iterator	iterDatum;
	bool						fContinuous;

	if( !m_iNodes )
		return false;

	for( i = 0; i < pData->GetGenes( ); ++i )
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j ) {
			if( !pData->IsExample( i, j ) )
				continue;
			fContinuous = false;
			for( k = 1; k < m_iNodes; ++k )
				if( m_apNodes[ k ]->IsContinuous( ) && !CMeta::IsNaN( pData->GetContinuous( i, j, k ) ) ) {
					fContinuous = true;
					break; }
			if( !fContinuous ) {
				strCur = EncodeDatum( pData, i, j );
				if( ( iterDatum = mapData.find( strCur ) ) != mapData.end( ) ) {
					if( pDatOut )
						pDatOut->Set( i, j, iterDatum->second );
					if( pvecvecdOut ) {
						pvecvecdOut->resize( pvecvecdOut->size( ) + 1 );
						(*pvecvecdOut)[ pvecvecdOut->size( ) - 1 ].push_back(
							iterDatum->second ); }
					continue; } }

			if( !Evaluate( pData, i, j, fZero, vecdCur ) )
				return false;
			if( !fContinuous )
				mapData[ strCur ] = vecdCur[ 0 ];
			if( pvecvecdOut ) {
				pvecvecdOut->resize( pvecvecdOut->size( ) + 1 );
				{
					vector<float>&	vecdOut	= (*pvecvecdOut)[ pvecvecdOut->size( ) - 1 ];

					for( k = 0; ( k + 1 ) < vecdCur.size( ); ++k )
						vecdOut.push_back( vecdCur[ k ] );
				} }
			if( pDatOut )
				pDatOut->Set( i, j, vecdCur[ 0 ] ); }

	return true; }

bool CBayesNetFNImpl::Evaluate( const IDataset* pData, size_t iOne, size_t iTwo, bool fZero,
	vector<float>& vecdOut ) const {
	vector<float>	vecdCur;
	float			dValue;
	size_t			i, j;

	if( !m_apNodes[ 0 ]->Evaluate( CMeta::GetNaN( ), vecdOut ) )
		return false;
	for( i = 0; i < vecdOut.size( ); ++i )
		vecdOut[ i ] = log( vecdOut[ i ] );
	for( i = 1; i < m_iNodes; ++i ) {
		dValue = 0;
		if( m_apNodes[ i ]->IsContinuous( ) ) {
			if( CMeta::IsNaN( dValue = pData->GetContinuous( iOne, iTwo, i ) ) && !fZero )
				continue; }
		else {
			if( ( ( j = pData->GetDiscrete( iOne, iTwo, i ) ) == -1 ) && !fZero )
				continue;
			dValue = (float)j; }
		if( !m_apNodes[ i ]->Evaluate( dValue, vecdCur ) || ( vecdCur.size( ) != vecdOut.size( ) ) )
			return false;
		for( j = 0; j < vecdCur.size( ); ++j )
			vecdOut[ j ] += log( vecdCur[ j ] ); }

	dValue = 0;
	for( i = 0; i < vecdOut.size( ); ++i )
		dValue += ( vecdOut[ i ] = exp( vecdOut[ i ] ) );
	for( i = 0; i < vecdOut.size( ); ++i )
		vecdOut[ i ] /= dValue;

	return true; }

bool CBayesNetFN::Evaluate( const vector<unsigned char>& vecbDatum, vector<float>& vecdOut, bool fZero ) const {
	vector<float>	vecdProd, vecdCur;
	float			dValue;
	size_t			i, j;

	if( !m_iNodes )
		return false;

	if( !m_apNodes[ 0 ]->Evaluate( CMeta::GetNaN( ), vecdProd ) )
		return false;
	for( i = 0; i < vecdProd.size( ); ++i )
		vecdProd[ i ] = log( vecdProd[ i ] );
	for( i = 1; i < m_iNodes; ++i ) {
		if( !( dValue = (float)vecbDatum[ i ] ) ) {
			if( !fZero )
				continue;
			dValue = 0; }
		else
			dValue--;
		if( !m_apNodes[ i ]->Evaluate( dValue, vecdCur ) )
			return false;
		for( j = 0; j < vecdCur.size( ); ++j )
			vecdProd[ j ] += log( vecdCur[ j ] ); }
	for( i = 0; ( i + 1 ) < vecdProd.size( ); ++i )
		vecdOut.push_back( exp( vecdProd[ i ] ) );

	return true; }

vector<string> CBayesNetFN::GetNodes( ) const {
	vector<string>	vecstrRet;
	size_t			i;

	vecstrRet.resize( m_iNodes );
	for( i = 0; i < vecstrRet.size( ); ++i )
		vecstrRet[ i ] = m_apNodes[ i ]->GetName( );

	return vecstrRet; }

unsigned char CBayesNetFN::GetValues( size_t iNode ) const {
	const CBayesNetFNNode*	pNode;

	pNode = m_apNodes[ iNode ];
	return ( pNode->IsContinuous( ) ? -1 : pNode->GetParameters( ) ); }

bool CBayesNetFN::IsContinuous( size_t iNode ) const {

	return m_apNodes[ iNode ]->IsContinuous( ); }

bool CBayesNetFN::IsContinuous( ) const {
	size_t	i;

	for( i = 0; i < m_iNodes; ++i )
		if( IsContinuous( i ) )
			return true;

	return false; }

void CBayesNetFN::Randomize( ) {
	size_t	i;

	for( i = 0; i < m_iNodes; ++i )
		Randomize( i ); }

void CBayesNetFN::Randomize( size_t iNode ) {

	m_apNodes[ iNode ]->Randomize( ); }

void CBayesNetFN::Reverse( size_t iNode ) {

	m_apNodes[ iNode ]->Reverse( ); }

}

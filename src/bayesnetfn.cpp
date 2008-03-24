#include "stdafx.h"
#include "bayesnet.h"
#include "statistics.h"

namespace Sleipnir {

#ifndef NO_SMILE

// CBayesNetFNNode //////////////////////////////////////////////////////////

const char	CBayesNetFNNode::c_szType[]	= "type";

CBayesNetFNNode* CBayesNetFNNode::Open( DSL_node* pNode ) {
	static const CBayesNetFNNode*	c_pDefault	= new CBayesNetFNNodeDiscrete( );
	static const CBayesNetFNNode*	c_apTypes[]	= {
		new CBayesNetFNNodeDiscrete( ),
		new CBayesNetFNNodeGaussian( ),
		new CBayesNetFNNodeBeta( ),
		new CBayesNetFNNodeExponential( ),
		new CBayesNetFNNodeMOG( ),
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
		pRet->m_Params.Initialize( iDim ? veciDims[ 0 ] : 1, veciDims[ iDim ] );
		for( i = 0; (size_t)i < pRet->m_Params.GetRows( ); ++i ) {
			if( iDim )
				veciDims[ 0 ] = i;
			for( j = 0; (size_t)j < pRet->m_Params.GetColumns( ); ++j ) {
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

// CBayesNetFNNodeDiscrete //////////////////////////////////////////////////

void CBayesNetFNNodeDiscrete::Randomize( ) {
	size_t		iCol, iRow, iSum;
	vector<int>	veciRow;

	veciRow.resize( m_Params.GetColumns( ) );
	for( iRow = 0; iRow < m_Params.GetRows( ); ++iRow ) {
		for( iSum = iCol = 0; iCol < veciRow.size( ); ++iCol )
			iSum += ( veciRow[ iCol ] = rand( ) );
		for( iCol = 0; iCol < veciRow.size( ); ++iCol )
			m_Params.Set( iRow, iCol, (float)veciRow[ iCol ] / iSum ); } }

bool CBayesNetFNNodeDiscrete::Learn( const IDataset* pData, size_t iNode, size_t iZero ) {
	CFullMatrix<size_t>	MatCounts;
	vector<size_t>		veciSums;
	size_t				i, j, iAnswer, iVal;

	veciSums.resize( m_Params.GetRows( ) );
	MatCounts.Initialize( m_Params.GetRows( ), m_Params.GetColumns( ) );
	for( i = 0; i < MatCounts.GetRows( ); ++i ) {
#pragma warning( disable : 4267 )
		veciSums[ i ] = MatCounts.GetColumns( );
#pragma warning( default : 4267 )
		for( j = 0; j < MatCounts.GetColumns( ); ++j )
			MatCounts.Set( i, j, 1 ); }
	for( i = 0; i < pData->GetGenes( ); ++i )
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j )
			if( pData->IsExample( i, j ) && ( ( iAnswer = pData->GetDiscrete( i, j, 0 ) ) != -1 ) ) {
				if( ( iVal = pData->GetDiscrete( i, j, iNode ) ) == -1 ) {
					if( iZero == -1 )
						continue;
					iVal = iZero; }
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

// CBayesNetFNNodeGaussian //////////////////////////////////////////////////

void CBayesNetFNNodeGaussian::Randomize( ) {
	size_t	i;

	for( i = 0; i < m_Params.GetRows( ); ++i ) {
		m_Params.Set( i, c_iMu, (float)( rand( ) - ( RAND_MAX / 2 ) ) );
		m_Params.Set( i, c_iSigma, (float)rand( ) / RAND_MAX ); } }

bool CBayesNetFNNodeGaussian::Learn( const IDataset* pData, size_t iNode, size_t iZero ) {
	float			d, dVal;
	size_t			i, j, iAnswer;
	vector<size_t>	veciCounts;

	for( i = 0; i < m_Params.GetRows( ); ++i )
		for( j = 0; j < m_Params.GetColumns( ); ++j )
			m_Params.Set( i, j, 0 );
	veciCounts.resize( m_Params.GetRows( ) );
	for( i = 0; i < pData->GetGenes( ); ++i )
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j )
			if( pData->IsExample( i, j ) && ( ( iAnswer = pData->GetDiscrete( i, j, 0 ) ) != -1 ) ) {
				if( CMeta::IsNaN( dVal = pData->GetContinuous( i, j, iNode ) ) ) {
					if( iZero == -1 )
						continue;
					dVal = (float)iZero; }
				if( iAnswer >= m_Params.GetRows( ) )
					return false;
				veciCounts[ iAnswer ]++;
				m_Params.Get( iAnswer, c_iMu ) += dVal;
				m_Params.Get( iAnswer, c_iSigma ) += dVal * dVal; }
	for( i = 0; i < m_Params.GetRows( ); ++i ) {
		if( !veciCounts[ i ] )
			veciCounts[ i ] = 1;
		d = ( m_Params.Get( i, c_iMu ) /= veciCounts[ i ] );
		d *= d;
		dVal = m_Params.Get( i, c_iSigma );
		dVal = ( dVal == d ) ? 1 : sqrt( ( dVal / ( veciCounts[ i ] - 1 ) ) - d );
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

// CBayesNetFNNodeBeta //////////////////////////////////////////////////////

void CBayesNetFNNodeBeta::Randomize( ) {
	size_t	i;

	for( i = 0; i < m_Params.GetRows( ); ++i ) {
		m_Params.Set( i, c_iMin, -( (float)rand( ) / RAND_MAX ) );
		m_Params.Set( i, c_iMax, (float)rand( ) / RAND_MAX );
		m_Params.Set( i, c_iAlpha, (float)rand( ) / RAND_MAX );
		m_Params.Set( i, c_iBeta, (float)rand( ) / RAND_MAX ); } }

bool CBayesNetFNNodeBeta::Learn( const IDataset* pData, size_t iNode, size_t iZero ) {
	float			d, dVal, dMean, dVariance;
	size_t			i, j, iAnswer;
	vector<size_t>	veciCounts;

	for( i = 0; i < m_Params.GetRows( ); ++i ) {
		m_Params.Set( i, c_iMin, FLT_MAX );
		m_Params.Set( i, c_iMax, -FLT_MAX );
		m_Params.Set( i, c_iAlpha, 0 );
		m_Params.Set( i, c_iBeta, 0 ); }
	veciCounts.resize( m_Params.GetRows( ) );
	for( i = 0; i < pData->GetGenes( ); ++i )
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j )
			if( pData->IsExample( i, j ) && ( ( iAnswer = pData->GetDiscrete( i, j, 0 ) ) != -1 ) ) {
				if( CMeta::IsNaN( dVal = pData->GetContinuous( i, j, iNode ) ) ) {
					if( iZero == -1 )
						continue;
					dVal = (float)iZero; }
				if( iAnswer >= m_Params.GetRows( ) )
					return false;
				veciCounts[ iAnswer ]++;
				m_Params.Get( iAnswer, c_iAlpha ) += dVal;
				m_Params.Get( iAnswer, c_iBeta ) += dVal * dVal;
				if( dVal < m_Params.Get( iAnswer, c_iMin ) )
					m_Params.Set( iAnswer, c_iMin, dVal );
				if( dVal > m_Params.Get( iAnswer, c_iMax ) )
					m_Params.Set( iAnswer, c_iMax, dVal ); }
	for( i = 0; i < m_Params.GetRows( ); ++i ) {
		if( !veciCounts[ i ] )
			veciCounts[ i ] = 1;
		d = m_Params.Get( i, c_iMax ) - m_Params.Get( i, c_iMin );
		dMean = m_Params.Get( i, c_iAlpha ) / veciCounts[ i ];
		dVariance = ( m_Params.Get( i, c_iBeta ) / ( veciCounts[ i ] - 1 ) ) - ( dMean * dMean );
		dMean = ( dMean - m_Params.Get( i, c_iMin ) ) / d;
		dVariance /= d * d;

		m_Params.Set( i, c_iAlpha, dMean * ( ( dMean * ( 1 - dMean ) / dVariance ) - 1 ) );
		m_Params.Set( i, c_iBeta, ( 1 - dMean ) * ( ( dMean * ( 1 - dMean ) / dVariance ) - 1 ) ); }

	return true; }

bool CBayesNetFNNodeBeta::Evaluate( float dValue, vector<float>& vecdOut ) const {
	float	d, dSum;
	size_t	i;

	vecdOut.resize( m_Params.GetRows( ) );
	dSum = 0;
	for( i = 0; i < vecdOut.size( ); ++i ) {
		if( ( d = dValue ) < m_Params.Get( i, c_iMin ) )
			d = m_Params.Get( i, c_iMin );
		else if( d > m_Params.Get( i, c_iMax ) )
			d = m_Params.Get( i, c_iMax );
		dSum += ( vecdOut[ i ] = (float)CStatistics::BetaPDF( d, m_Params.Get( i, c_iMin ),
			m_Params.Get( i, c_iMax ), m_Params.Get( i, c_iAlpha ), m_Params.Get( i, c_iBeta ) ) ); }
	for( i = 0; i < vecdOut.size( ); ++i )
		vecdOut[ i ] /= dSum;

	return true; }

// CBayesNetFNNodeExponential ///////////////////////////////////////////////

void CBayesNetFNNodeExponential::Randomize( ) {
	size_t	i;

	for( i = 0; i < m_Params.GetRows( ); ++i ) {
		m_Params.Set( i, c_iMin, (float)( rand( ) - ( RAND_MAX / 2 ) ) );
		m_Params.Set( i, c_iBeta, (float)rand( ) / RAND_MAX ); } }

bool CBayesNetFNNodeExponential::Learn( const IDataset* pData, size_t iNode, size_t iZero ) {
	float			dVal;
	size_t			i, j, iAnswer;
	vector<size_t>	veciCounts;

	for( i = 0; i < m_Params.GetRows( ); ++i ) {
		m_Params.Set( i, c_iMin, FLT_MAX );
		m_Params.Set( i, c_iBeta, 0 ); }
	veciCounts.resize( m_Params.GetRows( ) );
	for( i = 0; i < pData->GetGenes( ); ++i )
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j )
			if( pData->IsExample( i, j ) && ( ( iAnswer = pData->GetDiscrete( i, j, 0 ) ) != -1 ) ) {
				if( CMeta::IsNaN( dVal = pData->GetContinuous( i, j, iNode ) ) ) {
					if( iZero == -1 )
						continue;
					dVal = (float)iZero; }
				if( iAnswer >= m_Params.GetRows( ) )
					return false;
				veciCounts[ iAnswer ]++;
				m_Params.Get( iAnswer, c_iBeta ) += dVal;
				if( dVal < m_Params.Get( iAnswer, c_iMin ) )
					m_Params.Set( iAnswer, c_iMin, dVal ); }
	for( i = 0; i < m_Params.GetRows( ); ++i ) {
		if( !veciCounts[ i ] )
			veciCounts[ i ] = 1;
		m_Params.Set( i, c_iBeta, ( m_Params.Get( i, c_iBeta ) / veciCounts[ i ] ) - m_Params.Get( i, c_iMin ) ); }

	return true; }

bool CBayesNetFNNodeExponential::Evaluate( float dValue, vector<float>& vecdOut ) const {
	float	d, dSum;
	size_t	i;

	vecdOut.resize( m_Params.GetRows( ) );
	dSum = 0;
	for( i = 0; i < vecdOut.size( ); ++i ) {
		d = ( dValue < m_Params.Get( i, c_iMin ) ) ? c_iMin : dValue;
		dSum += ( vecdOut[ i ] = exp( ( m_Params.Get( i, c_iMin ) - d ) / m_Params.Get( i, c_iBeta ) ) /
			m_Params.Get( i, c_iBeta ) ); }
	for( i = 0; i < vecdOut.size( ); ++i )
		vecdOut[ i ] /= dSum;

	return true; }

// CBayesNetFNNodeMOG ///////////////////////////////////////////////////////

void CBayesNetFNNodeMOG::Randomize( ) {
	size_t	i, j;

	for( i = 0; i < m_Params.GetRows( ); ++i )
		for( j = 0; j < m_Params.GetColumns( ); j += 2 ) {
			m_Params.Set( i, j + c_iMu, (float)( rand( ) - ( RAND_MAX / 2 ) ) );
			m_Params.Set( i, j + c_iSigma, (float)rand( ) / RAND_MAX ); } }

bool CBayesNetFNNodeMOG::Learn( const IDataset* pData, size_t iNode, size_t iZero ) {

	return false; }
/*
	float			d, dVal;
	size_t			i, j, iAnswer;
	vector<size_t>	veciCounts;

	for( i = 0; i < m_Params.GetRows( ); ++i )
		for( j = 0; j < m_Params.GetColumns( ); ++j )
			m_Params.Set( i, j, 0 );
	veciCounts.resize( m_Params.GetRows( ) );
	for( i = 0; i < pData->GetGenes( ); ++i )
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j )
			if( pData->IsExample( i, j ) && ( ( iAnswer = pData->GetDiscrete( i, j, 0 ) ) != -1 ) ) {
				if( CMeta::IsNaN( dVal = pData->GetContinuous( i, j, iNode ) ) ) {
					if( iZero == -1 )
						continue;
					dVal = iZero; }
				if( iAnswer >= m_Params.GetRows( ) )
					return false;
				veciCounts[ iAnswer ]++;
				m_Params.Get( iAnswer, c_iMu ) += dVal;
				m_Params.Get( iAnswer, c_iSigma ) += dVal * dVal; }
	for( i = 0; i < m_Params.GetRows( ); ++i ) {
		if( !veciCounts[ i ] )
			veciCounts[ i ] = 1;
		d = ( m_Params.Get( i, c_iMu ) /= veciCounts[ i ] );
		d *= d;
		dVal = m_Params.Get( i, c_iSigma );
		dVal = ( dVal == d ) ? 1 : sqrt( ( dVal / ( veciCounts[ i ] - 1 ) ) - d );
		m_Params.Set( i, c_iSigma, dVal ); }

	return true; }
*/

bool CBayesNetFNNodeMOG::Evaluate( float dValue, vector<float>& vecdOut ) const {
	float	dCur, dSum;
	size_t	i, j;

	vecdOut.resize( m_Params.GetRows( ) );
	dSum = 0;
	for( i = 0; i < vecdOut.size( ); ++i ) {
		dCur = 0;
		for( j = 0; j < m_Params.GetColumns( ); j += 2 )
			dCur += (float)CStatistics::NormalPDF( dValue, m_Params.Get( i, j + c_iMu ),
				m_Params.Get( i, j + c_iSigma ) );
		dSum += ( vecdOut[ i ] = dCur ); }
	for( i = 0; i < vecdOut.size( ); ++i )
		vecdOut[ i ] /= dSum;

	return true; }

// CBayesNetFN //////////////////////////////////////////////////////////////

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
	size_t			i, j, iAnswer, iZero;
	vector<size_t>	veciCounts;
	int				iProp;

	if( !m_iNodes )
		return false;

	veciCounts.resize( m_apNodes[ 0 ]->GetParameters( ) );
	for( i = 0; i < pData->GetGenes( ); ++i )
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j )
			if( pData->IsExample( i, j ) && ( ( iAnswer = pData->GetDiscrete( i, j, 0 ) ) != -1 ) )
				veciCounts[ iAnswer ]++;

	if( !m_apNodes[ 0 ]->Learn( veciCounts ) )
		return false;
	for( i = 1; i < m_iNodes; ++i ) {
		DSL_userProperties&	Props	= m_SmileNet.GetNode( (int)i )->Info( ).UserProperties( );

		if( ( iProp = Props.FindProperty( c_szZero ) ) < 0 )
			iZero = fZero ? 0 : -1;
		else
			iZero = atoi( Props.GetPropertyValue( iProp ) );
		if( !m_apNodes[ i ]->Learn( pData, i, iZero ) )
			return false; }

	return true; }

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

	for( i = 0; i < pData->GetGenes( ); ++i ) {
		if( !( i % 250 ) )
			g_CatSleipnir.notice( "CBayesNetFN::Evaluate( %d ) %d/%d", fZero, i, pData->GetGenes( ) );
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
				pDatOut->Set( i, j, vecdCur[ 0 ] ); } }

	return true; }

bool CBayesNetFNImpl::Evaluate( const IDataset* pData, size_t iOne, size_t iTwo, bool fZero,
	vector<float>& vecdOut ) const {
	vector<float>	vecdCur;
	float			dValue;
	size_t			i, j, iZero;
	int				iProp;

	if( !m_apNodes[ 0 ]->Evaluate( CMeta::GetNaN( ), vecdOut ) )
		return false;
	for( i = 0; i < vecdOut.size( ); ++i )
		vecdOut[ i ] = log( vecdOut[ i ] );
	for( i = 1; i < m_iNodes; ++i ) {
		DSL_userProperties&	Props	= m_SmileNet.GetNode( (int)i )->Info( ).UserProperties( );

		if( ( iProp = Props.FindProperty( c_szZero ) ) < 0 )
			iZero = fZero ? 0 : -1;
		else
			iZero = atoi( Props.GetPropertyValue( iProp ) );

		dValue = 0;
		if( m_apNodes[ i ]->IsContinuous( ) ) {
			if( CMeta::IsNaN( dValue = pData->GetContinuous( iOne, iTwo, i ) ) ) {
				if( iZero == -1 )
					continue;
				dValue = (float)iZero; } }
		else {
			if( ( j = pData->GetDiscrete( iOne, iTwo, i ) ) == -1 ) {
				if( iZero == -1 )
					continue;
				j = iZero; }
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

bool CBayesNetFN::Evaluate( const vector<unsigned char>& vecbDatum, vector<float>& vecdOut, bool fZero,
	size_t iNode, bool fNoData ) const {
	vector<float>	vecdProd, vecdCur;
	float			dValue;
	size_t			i, j, iZero;
	int				iProp;

	if( !m_iNodes )
		return false;

	if( !m_apNodes[ iNode ]->Evaluate( CMeta::GetNaN( ), vecdProd ) )
		return false;
	for( i = 0; i < vecdProd.size( ); ++i )
		vecdProd[ i ] = log( vecdProd[ i ] );
	for( i = 1; i < m_iNodes; ++i ) {
		DSL_userProperties&	Props	= m_SmileNet.GetNode( (int)i )->Info( ).UserProperties( );

		if( ( iProp = Props.FindProperty( c_szZero ) ) < 0 )
			iZero = fZero ? 0 : -1;
		else
			iZero = atoi( Props.GetPropertyValue( iProp ) );

		if( !( dValue = (float)vecbDatum[ i ] ) ) {
			if( fNoData || ( iZero == -1 ) )
				continue;
			dValue = (float)iZero; }
		else
			dValue--;
		if( !m_apNodes[ i ]->Evaluate( dValue, vecdCur ) )
			return false;
		for( j = 0; j < vecdCur.size( ); ++j )
			vecdProd[ j ] += log( vecdCur[ j ] ); }
	for( i = 0; ( i + 1 ) < vecdProd.size( ); ++i )
		vecdOut.push_back( exp( vecdProd[ i ] ) );

	return true; }

void CBayesNetFN::GetNodes( std::vector<std::string>& vecstrNodes ) const {
	size_t	i;

	for( i = 0; i < m_iNodes; ++i )
		vecstrNodes.push_back( m_apNodes[ i ]->GetName( ) ); }

unsigned char CBayesNetFN::GetValues( size_t iNode ) const {
	const CBayesNetFNNode*	pNode;

	pNode = m_apNodes[ iNode ];
	return ( pNode->IsContinuous( ) ? -1 : pNode->GetParameters( ) ); }

bool CBayesNetFN::IsContinuous( ) const {
	size_t	i;

	for( i = 0; i < m_iNodes; ++i )
		if( IsContinuous( i ) )
			return true;

	return false; }

// CBayesNetMinimal //////////////////////////////////////////////////////////

/*!
 * \brief
 * Construct a new minimal Bayes net from the given SMILE-based network.
 * 
 * \param BNSmile
 * SMILE-based network from which to copy node parameters; must have naive structure.
 * 
 * \returns
 * True if Bayes net was successfully constructed.
 * 
 * \remarks
 * BNSmile must have only discrete nodes and naive structure.
 */
bool CBayesNetMinimal::Open( const CBayesNetSmile& BNSmile ) {
	CDataMatrix		Mat;
	vector<string>	vecstrNodes;
	size_t			i;

	if( m_adNY ) {
		delete[] m_adNY;
		m_adNY = NULL; }
	BNSmile.GetNodes( vecstrNodes );
	if( !vecstrNodes.size( ) )
		return false;
	BNSmile.GetCPT( 0, m_MatRoot );
	m_vecNodes.resize( vecstrNodes.size( ) - 1 );
	for( i = 0; i < m_vecNodes.size( ); ++i ) {
		BNSmile.GetCPT( i + 1, m_vecNodes[ i ].m_MatCPT );
		if( m_vecNodes[ i ].m_MatCPT.GetColumns( ) != m_MatRoot.GetRows( ) )
			return false;
		m_vecNodes[ i ].m_bDefault = BNSmile.GetDefault( i + 1 ); }
	m_adNY = new long double[ m_MatRoot.GetRows( ) ];

	return true; }

#endif // NO_SMILE

/*!
 * \brief
 * Load a minimal Bayes net from the given binary stream.
 * 
 * \param istm
 * Stream from which Bayes net is loaded.
 * 
 * \returns
 * True if Bayes net was successfully loaded.
 * 
 * \remarks
 * istm must be binary and contain a minimal Bayes net stored by CBayesNetMinimal::Save.
 */
bool CBayesNetMinimal::Open( std::istream& istm ) {
	uint32_t	iSize;
	size_t		i;

	if( m_adNY ) {
		delete[] m_adNY;
		m_adNY = NULL; }
	istm.read( (char*)&iSize, sizeof(iSize) );
	m_strID.resize( iSize );
	istm.read( &m_strID[ 0 ], iSize * sizeof(*m_strID.c_str( )) );
	if( !m_MatRoot.Open( istm, true ) )
		return false;
	m_adNY = new long double[ m_MatRoot.GetRows( ) ];
	istm.read( (char*)&iSize, sizeof(iSize) );
	m_vecNodes.resize( iSize );
	for( i = 0; i < m_vecNodes.size( ); ++i ) {
		CBayesNetMinimalNode&	BNNode	= m_vecNodes[ i ];

		istm.read( (char*)&BNNode.m_bDefault, sizeof(BNNode.m_bDefault) );
		if( !BNNode.m_MatCPT.Open( istm, true ) )
			return false; }

	return true; }

/*!
 * \brief
 * Save a minimal Bayes net to the given binary stream.
 * 
 * \param ostm
 * Stream to which Bayes net is saved.
 * 
 * \see
 * CBayesNetMinimal::Open
 */
void CBayesNetMinimal::Save( std::ostream& ostm ) const {
	uint32_t	iSize;
	size_t		i;

	iSize = (uint32_t)m_strID.length( );
	ostm.write( (const char*)&iSize, sizeof(iSize) );
	ostm.write( m_strID.c_str( ), iSize * sizeof(*m_strID.c_str( )) );
	m_MatRoot.Save( ostm, true );
	iSize = (uint32_t)m_vecNodes.size( );
	ostm.write( (const char*)&iSize, sizeof(iSize) );
	for( i = 0; i < m_vecNodes.size( ); ++i ) {
		const CBayesNetMinimalNode&	BNNode	= m_vecNodes[ i ];

		ostm.write( (const char*)&BNNode.m_bDefault, sizeof(BNNode.m_bDefault) );
		BNNode.m_MatCPT.Save( ostm, true ); } }

/*!
 * \brief
 * Perform Bayesian inference to obtain the class probability given evidence for some number of nodes.
 * 
 * \param vecbDatum
 * Values for each evidence node; 0xF indicates missing data (no evidence) for a particular node.  Note
 * that each evidence value is stored in <b>four bits</b>, not a full byte.
 * 
 * \param iOffset
 * Position of the first piece of evidence within vecbDatum; zero by default.  This can be used to
 * store multiple data in a single vector and rapidly perform inference for each subsequent data setting.
 * 
 * \returns
 * Posterior probability of the largest value of the class node given the evidence (generally the
 * probability of functional relationship).
 * 
 * \remarks
 * Evidence is stored in nibbles, not full bytes, so for a network containing N evidence (non-root) nodes,
 * vecbDatum must be of size at least iOffset + ceil(N/2).
 * 
 * \see
 * IBayesNet::Evaluate
 */
float CBayesNetMinimal::Evaluate( const vector<unsigned char>& vecbDatum, size_t iOffset ) const {
	long double		dNum, dDen;
	size_t			i, j;
	unsigned char	c;

	if( !m_adNY )
		return CMeta::GetNaN( );

	for( i = 0; i < m_MatRoot.GetRows( ); ++i )
		m_adNY[ i ] = m_MatRoot.Get( i, 0 );
	for( i = 0; i < m_vecNodes.size( ); ++i ) {
		c = vecbDatum[ ( i / 2 ) + iOffset ];
		if( i % 2 )
			c >>= 4;
		if( ( ( c &= 0xF ) == 0xF ) && ( ( c = m_vecNodes[ i ].m_bDefault ) == 0xFF ) )
			continue;

		const CDataMatrix&	MatCPT	= m_vecNodes[ i ].m_MatCPT;
		for( j = 0; j < MatCPT.GetColumns( ); ++j ) {
			if( c >= MatCPT.GetRows( ) ) {
				g_CatSleipnir.error( "CBayesNetMinimal::Evaluate( %d ) illegal value: %d/%d in %d", iOffset, c, MatCPT.GetRows( ), i );
				return CMeta::GetNaN( ); }
			m_adNY[ j ] *= MatCPT.Get( c, j ); } }

	dNum = dDen = m_adNY[ m_MatRoot.GetRows( ) - 1 ];
	for( i = 0; ( i + 1 ) < m_MatRoot.GetRows( ); ++i )
		dDen += m_adNY[ i ];

	return (float)( dNum / dDen ); }

/*!
 * \brief
 * Repeatedly perform Bayesian inference to obtain the class probability given evidence for some number of
 * nodes.
 * 
 * \param vecbData
 * Values for each evidence node; 0xF indicates missing data (no evidence) for a particular node.  Note
 * that each evidence value is stored in <b>four bits</b>, not a full byte.  Multiple sets of evidence can
 * be included in vecbData, e.g. for N nodes, entries 0 through floor(N/2) comprise one set of evidence,
 * floor(N/2)+1 through N the next, and so forth.
 * 
 * \param adResults
 * Array into which posterior probabilities of the largest value of the class node are inserted given the
 * evidence (generally probabilities of functional relationships).
 * 
 * \param iGenes
 * Number of inferences to perform and probabilities to generate.
 * 
 * \param iStart
 * First gene to process; this means that the first output probability is placed into the iStart element of
 * adResults, and the first element read from vecbDatum is at iStart * ceil(N/2).
 * 
 * \returns
 * True if evaluation was successful.
 * 
 * Perform Bayesian inference iGenes - iStart times using evidence from vecbData, which consists of zero or
 * more sets of evidence values for the N non-root nodes in the Bayes net.  In pseudocode:
 * \code
 * for( i = iStart; i < iGenes; ++i )
 *   adValues[ i ] = Evaluate( vecbData, i * floor((N+1)/2) );
 * \endcode
 * 
 * \remarks
 * Evidence is stored in nibbles, not full bytes, so for a network containing N evidence (non-root) nodes,
 * vecbData must be of size at least iGenes * ceil(N/2).
 * 
 * \see
 * IBayesNet::Evaluate
 */
bool CBayesNetMinimal::Evaluate( const vector<unsigned char>& vecbData, float* adResults,
	size_t iGenes, size_t iStart ) const {
	size_t	iGene, iOffset, iChunk;

	if( !adResults )
		return false;

	iChunk = ( m_vecNodes.size( ) + 1 ) / 2;
	iOffset = iChunk * iStart;
	for( iGene = iStart; ( iGene < iGenes ) && ( iOffset < vecbData.size( ) ); ++iGene,iOffset += iChunk )
		adResults[ iGene ] = Evaluate( vecbData, iOffset );

	return true; }

}

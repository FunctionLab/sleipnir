#include "stdafx.h"
#include "bayesnet.h"
#include "statistics.h"

namespace libBioUtils {

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

	for( i = 0; i < pData->GetGenes( ); ++i ) {
		if( !( i % 250 ) )
			g_CatBioUtils.notice( "CBayesNetFN::Evaluate( %d ) %d/%d", fZero, i, pData->GetGenes( ) );
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

bool CBayesNetFN::Evaluate( const vector<unsigned char>& vecbDatum, vector<float>& vecdOut, bool fZero ) const {
	vector<float>	vecdProd, vecdCur;
	float			dValue;
	size_t			i, j, iZero;
	int				iProp;

	if( !m_iNodes )
		return false;

	if( !m_apNodes[ 0 ]->Evaluate( CMeta::GetNaN( ), vecdProd ) )
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
			if( iZero == -1 )
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

void CBayesNetFN::GetNodes( vector<string>& vecstrNodes ) const {
	size_t	i;

	for( i = 0; i < m_iNodes; ++i )
		vecstrNodes.push_back( m_apNodes[ i ]->GetName( ) ); }

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

bool CBayesNetFN::GetCPT( size_t iNode, CDataMatrix& MatCPT ) const {

	return CBayesNetSmileImpl::GetCPT( m_SmileNet.GetNode( (int)iNode ), MatCPT ); }

bool CBayesNetFN::Evaluate( const CPCLPair&, CPCL&, bool, int ) const {

	return false; }

}

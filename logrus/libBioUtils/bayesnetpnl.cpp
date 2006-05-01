#include "stdafx.h"
#include "bayesnet.h"
#include "dat.h"
#include "dataset.h"
#include "meta.h"

namespace libBioUtils {

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

bool CBayesNetPNL::Learn( const IDataset* pData, size_t iIterations, bool fZero, bool fELR ) {
	CEMLearningEngineDumb*	pLearner;

	if( !m_pPNLNet || fELR )
		return false;

	pLearner = CEMLearningEngineDumb::Create( m_pPNLNet );
	pLearner->SetMaxIterEM( (int)iIterations );
	pLearner->Learn( pData, fZero );

	delete pLearner;
	return true; }

vector<string> CBayesNetPNL::GetNodes( ) const {
	vector<string>	vecstrRet;

	return vecstrRet; }

bool CBayesNetPNL::IsContinuous( size_t ) const {

	return IsContinuous( ); }

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

void CBayesNetPNL::Randomize( ) { }

void CBayesNetPNL::Randomize( size_t ) { }

void CBayesNetPNL::Reverse( size_t ) { }

bool CBayesNetPNL::Evaluate( const vector<unsigned char>&, vector<float>&, bool ) const {

	return false; }

unsigned char CBayesNetPNL::GetValues( size_t ) const {

	return 0; }

bool CBayesNetPNL::GetCPT( size_t, CDataMatrix& ) const {

	return false; }

}

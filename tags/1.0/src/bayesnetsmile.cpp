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
#include "bayesnet.h"
#include "dat.h"
#include "dataset.h"
#include "meta.h"

#ifndef NO_SMILE

#if ( defined(_MSC_VER) && defined(_DEBUG) )
extern "C" void __cdecl _invalid_parameter_noinfo( ) { }
#endif // ( defined(_MSC_VER) && defined(_DEBUG) )

namespace Sleipnir {

const char	CBayesNetSmileImpl::c_szGaussian[]	= "gaussian";

bool CBayesNetSmileImpl::GetCPT( DSL_node* pNode, CDataMatrix& MatCPT ) {
	DSL_Dmatrix*	pMat;
	DSL_intArray	veciCoord;

	pMat = pNode->Definition( )->GetMatrix( );
	const DSL_intArray&	veciDims	= pMat->GetDimensions( );

	if( veciDims.GetSize( ) > 2 )
		return false;
	pMat->IndexToCoordinates( 0, veciCoord );
	if( veciDims.GetSize( ) == 1 ) {
		MatCPT.Initialize( veciDims[ 0 ], 1 );
		for( veciCoord[ 0 ] = 0; veciCoord[ 0 ] < veciDims[ 0 ]; ++veciCoord[ 0 ] )
			MatCPT.Set( veciCoord[ 0 ], 0, (float)(*pMat)[ veciCoord ] );
		return true; }

	MatCPT.Initialize( veciDims[ 1 ], veciDims[ 0 ] );
	for( veciCoord[ 0 ] = 0; veciCoord[ 0 ] < veciDims[ 0 ]; ++veciCoord[ 0 ] )
		for( veciCoord[ 1 ] = 0; veciCoord[ 1 ] < veciDims[ 1 ]; ++veciCoord[ 1 ] )
			MatCPT.Set( veciCoord[ 1 ], veciCoord[ 0 ], (float)(*pMat)[ veciCoord ] );
	return true; }

bool CBayesNetSmileImpl::IsGaussian( const DSL_network& BayesNet ) {
	int	i;

	if( ( i = ((DSL_network&)BayesNet).UserProperties( ).FindProperty( c_szGaussian ) ) < 0 )
		return false;

	return !!atoi( ((DSL_network&)BayesNet).UserProperties( ).GetPropertyValue( i ) ); }

bool CBayesNetSmileImpl::IsNaive( const DSL_network& BayesNet ) {
	int	i;

	{
		const DSL_intArray&	veciParents	= ((DSL_network&)BayesNet).GetNode( 0 )->Parents( );

		if( veciParents.NumItems( ) != 0 )
			return false;
	}
	for( i = 1; i < BayesNet.GetNumberOfNodes( ); ++i ) {
		const DSL_intArray&	veciParents	= ((DSL_network&)BayesNet).GetNode( i )->Parents( );

		if( ( veciParents.NumItems( ) > 1 ) || ( veciParents[ 0 ] != 0 ) )
			return false; }

	return true; }

CBayesNetSmileImpl::CBayesNetSmileImpl( bool fGroup ) : CBayesNetImpl(fGroup),
	m_fSmileNet(false), m_pDefaults(NULL) { }

/*!
 * \brief
 * Construct a new SMILE-based Bayes net.
 * 
 * \param fGroup
 * If true, group identical learning/evaluation examples together into a single heavily weighted example.
 * 
 * \remarks
 * There's essentially never a reason to set fGroup to false.
 */
CBayesNetSmile::CBayesNetSmile( bool fGroup ) : CBayesNetSmileImpl( fGroup ) { }

bool CBayesNetSmileImpl::LearnGrouped( const IDataset* pData, size_t iIterations, bool fZero ) {
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
				g_CatSleipnir.notice( "CBayesNetSmile::LearnGrouped( %d, %d ) iteration %d, datum %d/%d",
					iIterations, fZero, iIter, ( iDatum - 1 ), mapData.size( ) );
			FillCPTs( vecfHidden, iterDatum->first, fZero, true );
			m_SmileNet.UpdateBeliefs( );

			for( i = 0; i < (size_t)m_SmileNet.GetNumberOfNodes( ); ++i )
				LearnExpected( m_SmileNet.GetNode( (int)i ), vecpExpected[ i ],
					iterDatum->second ); }
		for( i = 0; i < (size_t)m_SmileNet.GetNumberOfNodes( ); ++i ) {
			pMat = m_SmileNet.GetNode( (int)i )->Definition( )->GetMatrix( );
			for( pMat->IndexToCoordinates( (int)( j = 0 ), veciCoords );
				j != DSL_OUT_OF_RANGE; j = pMat->NextCoordinates( veciCoords ) )
				pMat->Subscript( veciCoords ) = vecpExpected[ i ]->Subscript( veciCoords );
			pMat->Normalize( ); } }
	for( i = 0; i < vecpExpected.size( ); ++i )
		delete vecpExpected[ i ];

	return true; }

bool CBayesNetSmileImpl::FillCPTs( const IDataset* pData, size_t iOne, size_t iTwo, bool fZero, bool fLearn ) {
	size_t	i, iVal, iZero;
	int		iProp;

	if( !pData->IsExample( iOne, iTwo ) || ( fLearn && ( pData->GetDiscrete( iOne, iTwo, 0 ) == -1 ) ) )
		return false;

	m_SmileNet.ClearAllEvidence( );
	for( i = fLearn ? 0 : 1; i < (size_t)m_SmileNet.GetNumberOfNodes( ); ++i ) {
		if( pData->IsHidden( i ) )
			continue;

		DSL_userProperties&	Props	= m_SmileNet.GetNode( (int)i )->Info( ).UserProperties( );

		if( ( iProp = Props.FindProperty( c_szZero ) ) < 0 )
			iZero = fZero ? 0 : -1;
		else
			iZero = atoi( Props.GetPropertyValue( iProp ) );

		if( ( iVal = pData->GetDiscrete( iOne, iTwo, i ) ) == -1 ) {
			if( iZero == -1 )
				continue;
			iVal = iZero; }
		m_SmileNet.GetNode( (int)i )->Value( )->SetEvidence( (int)iVal ); }

	return true; }

bool CBayesNetSmileImpl::FillCPTs( const std::vector<bool>& vecfHidden, const std::string& strDatum,
	bool fZero, bool fLearn, bool fAll ) {
	size_t	i, iVal, iZero;
	int		iProp;

	if( !fAll && fLearn && !IsAnswer( strDatum ) )
		return false;

	m_SmileNet.ClearAllEvidence( );
	for( i = ( fAll || fLearn ) ? 0 : 1; i < (size_t)m_SmileNet.GetNumberOfNodes( ); ++i ) {
		if( vecfHidden[ i ] )
			continue;

		DSL_userProperties&	Props	= m_SmileNet.GetNode( (int)i )->Info( ).UserProperties( );

		if( ( iProp = Props.FindProperty( c_szZero ) ) < 0 )
			iZero = fZero ? 0 : -1;
		else
			iZero = atoi( Props.GetPropertyValue( iProp ) );

		if( strDatum[ i ] == c_cMissing ) {
			if( iZero == -1 )
				continue;
			iVal = iZero; }
		else
			iVal = strDatum[ i ] - c_cBase;
		m_SmileNet.GetNode( (int)i )->Value( )->SetEvidence( (int)iVal ); }

	return true; }

bool CBayesNetSmileImpl::FillCPTs( const vector<bool>& vecfHidden, const vector<unsigned char>& vecbDatum,
	bool fZero, bool fLearn, bool fNoData ) {
	size_t	i, iVal, iZero;
	int		iProp;

	if( fLearn && !vecbDatum[ 0 ] )
		return false;

	m_SmileNet.ClearAllEvidence( );
	for( i = fLearn ? 0 : 1; i < (size_t)m_SmileNet.GetNumberOfNodes( ); ++i ) {
		if( vecfHidden[ i ] )
			continue;

		DSL_userProperties&	Props	= m_SmileNet.GetNode( (int)i )->Info( ).UserProperties( );

		if( ( iProp = Props.FindProperty( c_szZero ) ) < 0 )
			iZero = fZero ? 0 : -1;
		else
			iZero = atoi( Props.GetPropertyValue( iProp ) );

		if( !vecbDatum[ i ] ) {
			if( fNoData || ( iZero == -1 ) )
				continue;
			iVal = iZero; }
		else
			iVal = vecbDatum[ i ] - 1;
		m_SmileNet.GetNode( (int)i )->Value( )->SetEvidence( (int)iVal ); }

	return true; }

bool CBayesNetSmileImpl::LearnUngrouped( const IDataset* pData, size_t iIterations, bool fZero ) {
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
				g_CatSleipnir.notice( "CBayesNetSmile::LearnUngrouped( %d, %d ) iteration %d, gene %d/%d",
					iIterations, fZero, iIter, i, pData->GetGenes( ) );
			for( j = ( i + 1 ); j < pData->GetGenes( ); ++j ) {
				if( !FillCPTs( pData, i, j, fZero, true ) )
					continue;
				m_SmileNet.UpdateBeliefs( );

				for( k = 0; k < (size_t)m_SmileNet.GetNumberOfNodes( ); ++k )
					LearnExpected( m_SmileNet.GetNode( (int)k ), vecpExpected[ k ] ); } }
		for( i = 0; i < (size_t)m_SmileNet.GetNumberOfNodes( ); ++i ) {
			pMat = m_SmileNet.GetNode( (int)i )->Definition( )->GetMatrix( );
			for( pMat->IndexToCoordinates( (int)( j = 0 ), veciCoords );
				j != DSL_OUT_OF_RANGE; j = pMat->NextCoordinates( veciCoords ) )
				pMat->Subscript( veciCoords ) = vecpExpected[ i ]->Subscript( veciCoords );
			pMat->Normalize( ); } }
	for( i = 0; i < vecpExpected.size( ); ++i )
		delete vecpExpected[ i ];

	return true; }

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

#ifdef PNL_ENABLED

/*!
 * \brief
 * Generate a PNL-based Bayes net equivalent to the current SMILE-based network.
 * 
 * \param BNPNL
 * PNL-based network into which the current SMILE-based network is copied.
 * 
 * If PNL is enabled, generate a PNL-based network with equivalent structure and parameters to the current
 * SMILE-based network.
 */
bool CBayesNetSmile::Convert( CBayesNetPNL& BNPNL ) const {

	if( !m_fSmileNet )
		return false;

	return( ConvertGraph( BNPNL ) && ConvertCPTs( BNPNL ) ); }

#endif // PNL_ENABLED

void CBayesNetSmile::GetNodes( std::vector<std::string>& vecstrNodes ) const {
	int	i;

	if( m_fSmileNet )
		for( i = 0; i < m_SmileNet.GetNumberOfNodes( ); ++i )
			vecstrNodes.push_back( m_SmileNet.GetNode( i )->Info( ).Header( ).GetId( ) ); }

bool CBayesNetSmileImpl::Evaluate( const IDataset* pData, CDat* pDatOut, TVecVecD* pvecvecdOut,
	bool fZero ) const {
	size_t						i, j, k, iOne, iTwo;
	DSL_nodeValue*				pValue;
	string						strCur;
	map<string,float>			mapData;
	map<string,float>::iterator	iterDatum;
	vector<bool>				vecfHidden;
	bool						fZeroable;
	float						dPrior;
	vector<size_t>				veciGenes;

	if( !m_fSmileNet || IsContinuous( ) )
		return false;

	vecfHidden.resize( pData->GetExperiments( ) );
	for( i = 0; i < vecfHidden.size( ); ++i )
		vecfHidden[ i ] = pData->IsHidden( i );
	if( !( fZeroable = fZero ) )
		for( i = 1; i < (size_t)m_SmileNet.GetNumberOfNodes( ); ++i ) {
			DSL_userProperties&	Props	= m_SmileNet.GetNode( (int)i )->Info( ).UserProperties( );
			if( Props.FindProperty( c_szZero ) >= 0 ) {
				fZeroable = true;
				break; } }
	if( pDatOut ) {
		veciGenes.resize( pData->GetGenes( ) );
		for( i = 0; i < pData->GetGenes( ); ++i )
			veciGenes[ i ] = pDatOut->GetGene( pData->GetGene( i ) );
		((CBayesNetSmileImpl*)this)->m_SmileNet.UpdateBeliefs( );
		pValue = m_SmileNet.GetNode( 0 )->Value( );
		dPrior = (float)(*pValue->GetMatrix( ))[ 0 ]; }
	for( i = 0; i < pData->GetGenes( ); ++i ) {
		if( !( i % 250 ) )
			g_CatSleipnir.notice( "CBayesNetSmile::Evaluate( %d ) %d/%d", fZero, i,
				pData->GetGenes( ) );
		if( pDatOut && !pvecvecdOut && ( ( iOne = veciGenes[ i ] ) == -1 ) )
			continue;
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j ) {
			if( pDatOut && !pvecvecdOut && ( ( iTwo = veciGenes[ j ] ) == -1 ) )
				continue;
			if( !( fZeroable || pData->IsExample( i, j ) ) ) {
				if( pDatOut && ( iOne != -1 ) && ( iTwo != -1 ) )
					pDatOut->Set( iOne, iTwo, dPrior );
				continue; }
			strCur = EncodeDatum( pData, i, j );
			if( m_fGroup && ( ( iterDatum = mapData.find( strCur ) ) != mapData.end( ) ) ) {
				if( pDatOut && ( iOne != -1 ) && ( iTwo != -1 ) )
					pDatOut->Set( iOne, iTwo, iterDatum->second );
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

					for( k = 0; ( k + 1 ) < (size_t)pValue->GetSize( ); ++k )
						vecdCur.push_back( (float)(*pValue->GetMatrix( ))[ (int)k ] );
				} }
			if( pDatOut && ( iOne != -1 ) && ( iTwo != -1 ) )
				pDatOut->Set( iOne, iTwo, (float)(*pValue->GetMatrix( ))[ 0 ] ); } }

	return true; }

bool CBayesNetSmile::Evaluate( const vector<unsigned char>& vecbDatum, vector<float>& vecdResults, bool fZero,
	size_t iNode, bool fIgnoreMissing ) const {
	vector<bool>	vecfHidden;
	DSL_nodeValue*	pValue;
	size_t			i;

	if( !m_fSmileNet || IsContinuous( ) )
		return false;

	vecfHidden.resize( vecbDatum.size( ) );
	for( i = 0; i < vecfHidden.size( ); ++i )
		vecfHidden[ i ] = false;
	((CBayesNetSmile*)this)->FillCPTs( vecfHidden, vecbDatum, fZero, false, fIgnoreMissing );
	((CBayesNetSmile*)this)->m_SmileNet.UpdateBeliefs( );
	pValue = m_SmileNet.GetNode( iNode )->Value( );
	for( i = 0; ( i + 1 ) < (size_t)pValue->GetSize( ); ++i )
		vecdResults.push_back( (float)(*pValue->GetMatrix( ))[ (int)i ] );

	return true; }

/*!
 * \brief
 * Evaluate the output of a Bayesian classifier given only a single node's evidence value.
 * 
 * \param iNode
 * Node for which evidence is set.
 * 
 * \param bValue
 * Value of evidence to set.
 * 
 * \returns
 * Posterior probabillity of classifier node's first value given the data.
 * 
 * Evaluates the posterior probability of the Bayesian network's first node (i.e. the class node) given
 * only a single piece of evidence.  This can be used to calculate the impact of a single dataset on
 * predicted probabilities, for example.
 * 
 * \remarks
 * Unlike other evaluation methods, ignores default values for all nodes.
 */
float CBayesNetSmile::Evaluate( size_t iNode, unsigned char bValue ) const {
	vector<bool>			vecfHidden;
	vector<unsigned char>	vecbDatum;
	DSL_nodeValue*			pValue;
	size_t					i;

	if( !m_fSmileNet || IsContinuous( ) )
		return CMeta::GetNaN( );

	vecbDatum.resize( m_SmileNet.GetNumberOfNodes( ) );
	vecbDatum[ iNode ] = bValue + 1;
	vecfHidden.resize( vecbDatum.size( ) );
	for( i = 0; i < vecbDatum.size( ); ++i )
		vecfHidden[ i ] = ( i != iNode );
	((CBayesNetSmile*)this)->FillCPTs( vecfHidden, vecbDatum, false, false );
	((CBayesNetSmile*)this)->m_SmileNet.UpdateBeliefs( );
	pValue = m_SmileNet.GetNode( 0 )->Value( );

	return (float)(*pValue->GetMatrix( ))[ 0 ]; }

/*!
 * \brief
 * Returns the default value (if any) for the requested node.
 * 
 * \param iNode
 * Node whose default value should be retrieved.
 * 
 * \returns
 * Default value of the requested node, or -1 if none exists.
 */
unsigned char CBayesNetSmile::GetDefault( size_t iNode ) const {
	int	i;

	if( !m_fSmileNet ||
		( ( i = ((DSL_network&)m_SmileNet).GetNode(
		iNode )->Info( ).UserProperties( ).FindProperty( c_szZero ) ) < 0 ) )
		return -1;

	return atoi( ((DSL_network&)m_SmileNet).GetNode(
		iNode )->Info( ).UserProperties( ).GetPropertyValue( i ) ); }

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
	DSL_Dmatrix*			pDefault;
	DSL_intArray			veciCoords;
	vector<size_t>			veciZeros;
	int						iProp;
	bool					fZeroable, fFallback;
	float					dLambda;
	double					dCount;

	vecveciCounts.resize( m_SmileNet.GetNumberOfNodes( ) );
	iAnswers = m_SmileNet.GetNode( 0 )->Definition( )->GetNumberOfOutcomes( );
	vecveciCounts[ 0 ].resize( iAnswers );
	for( i = 1; i < vecveciCounts.size( ); ++i )
		vecveciCounts[ i ].resize( iAnswers *
			m_SmileNet.GetNode( (int)i )->Definition( )->GetNumberOfOutcomes( ) );
	veciZeros.resize( m_SmileNet.GetNumberOfNodes( ) );
	fZeroable = fZero;
	for( i = 0; i < veciZeros.size( ); ++i ) {
		DSL_userProperties&	Props	= m_SmileNet.GetNode( (int)i )->Info( ).UserProperties( );

		if( ( iProp = Props.FindProperty( c_szZero ) ) < 0 )
			veciZeros[ i ] = fZero ? 0 : -1;
		else {
			fZeroable = true;
			veciZeros[ i ] = atoi( Props.GetPropertyValue( iProp ) ); } }
	for( iCount = i = 0; i < pData->GetGenes( ); ++i )
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j )
			if( ( fZeroable || pData->IsExample( i, j ) ) &&
				( ( iAnswer = pData->GetDiscrete( i, j, 0 ) ) != -1 ) ) {
				vecveciCounts[ 0 ][ iAnswer ]++;
				iCount++;
				for( k = 1; k < pData->GetExperiments( ); ++k ) {
					if( ( iVal = pData->GetDiscrete( i, j, k ) ) == -1 ) {
						if( veciZeros[ k ] == -1 )
							continue;
						iVal = veciZeros[ k ]; }
//iVal = iVal % m_SmileNet.GetNode( k )->Definition( )->GetNumberOfOutcomes( );
					vecveciCounts[ k ][ ( iVal * iAnswers ) + iAnswer ]++; } }

	fFallback = m_pDefaults && ( iCount < c_iMinimum );
	pMat = m_SmileNet.GetNode( 0 )->Definition( )->GetMatrix( );
	for( i = 0; i < iAnswers; ++i )
		(*pMat)[ (int)i ] = ( j = vecveciCounts[ 0 ][ (int)i ] ) ? j : ( fFallback ? 0 : 1 );
	if( fFallback ) {
		g_CatSleipnir.warn( "CBayesNetSmile::LearnNaive( %d ) insufficient data for node %s",
			fZero, m_SmileNet.GetNode( 0 )->Info( ).Header( ).GetId( ) );
		dLambda = 1 - ( (float)iCount / c_iMinimum );
		pMat->Normalize( );
		pDefault = m_pDefaults->m_SmileNet.GetNode( 0 )->Definition( )->GetMatrix( );
		for( i = 0; i < iAnswers; ++i )
			(*pMat)[ (int)i ] = ( ( 1 - dLambda ) * (*pMat)[ (int)i ] ) +
				( dLambda * (*pDefault)[ (int)i ] ); }
	pMat->Normalize( );
	for( i = 1; i < vecveciCounts.size( ); ++i ) {
		pDef = m_SmileNet.GetNode( (int)i )->Definition( );
		pMat = pDef->GetMatrix( );
		pMat->IndexToCoordinates( 0, veciCoords );
		pDefault = m_pDefaults ? m_pDefaults->m_SmileNet.GetNode( (int)i )->Definition( )->GetMatrix( ) : NULL;
		for( j = 0; j < iAnswers; ++j ) {
			veciCoords[ 0 ] = (int)j;
			for( k = 0; k < (size_t)pDef->GetNumberOfOutcomes( ); ++k ) {
				veciCoords[ 1 ] = (int)k;
				(*pMat)[ veciCoords ] = vecveciCounts[ i ][ ( k * iAnswers ) + j ]; } }
		if( pDefault )
			for( j = 0; j < iAnswers; ++j ) {
				veciCoords[ 0 ] = (int)j;
				for( dCount = k = 0; k < (size_t)pDef->GetNumberOfOutcomes( ); ++k ) {
					veciCoords[ 1 ] = (int)k;
					dCount += (*pMat)[ veciCoords ]; }
				if( dCount < c_iMinimum ) {
					g_CatSleipnir.warn( "CBayesNetSmile::LearnNaive( %d ) insufficient data for node %s, column %d",
						fZero, m_SmileNet.GetNode( (int)i )->Info( ).Header( ).GetId( ), j );
					dLambda = 1 - ( (float)dCount / c_iMinimum );
					for( k = 0; k < (size_t)pDef->GetNumberOfOutcomes( ); ++k ) {
						veciCoords[ 1 ] = (int)k;
						(*pMat)[ veciCoords ] = ( dCount ? ( ( 1 - dLambda ) * (*pMat)[ veciCoords ] /
							dCount ) : 0 ) + ( dLambda * (*pDefault)[ veciCoords ] ); } }
				else
					for( k = 0; k < (size_t)pDef->GetNumberOfOutcomes( ); ++k ) {
						veciCoords[ 1 ] = (int)k;
						if( !(*pMat)[ veciCoords ] )
							(*pMat)[ veciCoords ] = 1; } }
		else
			for( j = 0; j < iAnswers; ++j ) {
				veciCoords[ 0 ] = (int)j;
				for( k = 0; k < (size_t)pDef->GetNumberOfOutcomes( ); ++k ) {
					veciCoords[ 1 ] = (int)k;
					if( !(*pMat)[ veciCoords ] )
						(*pMat)[ veciCoords ] = 1; } }
		pMat->Normalize( ); }

	return true; }

bool CBayesNetSmile::Evaluate( const CPCLPair& PCLData, CPCL& PCLResults, bool fZero, int iAlgorithm ) const {
	size_t									i, j, k, iExp;
	string									strCur;
	map<string, vector<float> >				mapData;
	map<string, vector<float> >::iterator	iterDatum;
	vector<size_t>							veciMap;
	vector<bool>							vecfHidden;
	int										iPrev;

	if( !m_fSmileNet || IsContinuous( ) )
		return false;

	iPrev = ((CBayesNetSmile*)this)->m_SmileNet.GetDefaultBNAlgorithm( );
	veciMap.resize( m_SmileNet.GetNumberOfNodes( ) );
	vecfHidden.resize( veciMap.size( ) );
	for( i = 0; i < veciMap.size( ); ++i ) {
		veciMap[ i ] = -1;
		vecfHidden[ i ] = true;
		for( j = 0; j < PCLData.GetExperiments( ); ++j )
			if( PCLData.GetExperiment( j ) == m_SmileNet.GetNode( (int)i )->Info( ).Header( ).GetId( ) ) {
				vecfHidden[ i ] = false;
				veciMap[ i ] = (unsigned int)j;
				break; } }
	((CBayesNetSmile*)this)->m_SmileNet.SetDefaultBNAlgorithm( iAlgorithm );
	for( i = 0; i < PCLResults.GetGenes( ); ++i ) {
		if( !( i % 1 ) )
			g_CatSleipnir.notice( "CBayesNetSmile::Evaluate( %d ) %d/%d", fZero, i,
				PCLResults.GetGenes( ) );
		strCur = EncodeDatum( PCLData, PCLData.GetGene( PCLResults.GetGene( i ) ), veciMap );
		if( m_fGroup && ( ( iterDatum = mapData.find( strCur ) ) != mapData.end( ) ) ) {
			for( j = 0; j < iterDatum->second.size( ); ++j )
				PCLResults.Set( i, j, iterDatum->second[ j ] );
			continue; }

		((CBayesNetSmile*)this)->FillCPTs( vecfHidden, strCur, fZero, false, true );
		((CBayesNetSmile*)this)->m_SmileNet.UpdateBeliefs( );
		for( iExp = j = 0; j < veciMap.size( ); ++j ) {
			DSL_Dmatrix*	pMatrix;

			if( veciMap[ j ] != -1 )
				continue;
			pMatrix = m_SmileNet.GetNode( (int)j )->Value( )->GetMatrix( );
			for( k = 0; k < GetValues( j ); ++k )
				PCLResults.Set( i, iExp++, (float)(*pMatrix)[ (int)k ] ); }
		if( m_fGroup ) {
			vector<float>	vecfCur;

			vecfCur.resize( PCLResults.GetExperiments( ) );
			for( j = 0; j < vecfCur.size( ); ++j )
				vecfCur[ j ] = PCLResults.Get( i, j );
			mapData[ strCur ] = vecfCur; } }
	((CBayesNetSmile*)this)->m_SmileNet.SetDefaultBNAlgorithm( iPrev );

	return true; }

/*!
 * \brief
 * Construct a new SMILE-based naive Bayes net with nodes corresponding to the given datasets.
 * 
 * \param vecstrFiles
 * Filenames of datasets, one per node.
 * 
 * \param iValues
 * Number of values into which each dataset will be quantized.
 * 
 * \returns
 * True if Bayes net was successfully constructed.
 * 
 * This version of Open can be used to quickly construct a uniform, naive Bayes net corresponding to a
 * particular set of data.  These data files are usually PCLs or DATs containing microarray data, since
 * large numbers of microarray datasets can be processed in this manner.  In addition to one node per
 * given file, one additional class node will be created at the top of the naive model with two possible
 * values (generally corresponding to functional unrelatedness or relatedness).
 * 
 * \see
 * CDat | CDataPair | CPCL | CPCLPair
 */
bool CBayesNetSmile::Open( const std::vector<std::string>& vecstrFiles, size_t iValues ) {
	size_t			i, j;
	DSL_stringArray	vecstrOutcomes;
	string			strCur;

	m_fSmileNet = true;
	m_SmileNet.DeleteAllNodes( );
	m_SmileNet.AddNode( DSL_CPT, (char*)c_szFR );
	vecstrOutcomes.Add( ( (string)c_szFR + "No" ).c_str( ) );
	vecstrOutcomes.Add( ( (string)c_szFR + "Yes" ).c_str( ) );
	m_SmileNet.GetNode( 0 )->Definition( )->SetNumberOfOutcomes( vecstrOutcomes );
	for( i = 0; i < vecstrFiles.size( ); ++i ) {
		m_SmileNet.AddNode( DSL_CPT, (char*)( strCur =
			CMeta::Filename( CMeta::Deextension( vecstrFiles[ i ] ) ) ).c_str( ) );
		vecstrOutcomes.Flush( );
		for( j = 0; j < iValues; ++j ) {
			char	acNum[ 8 ];

#pragma warning( disable : 4996 )
			sprintf( acNum, "%02d", j );
#pragma warning( default : 4996 )
			vecstrOutcomes.Add( ( strCur + acNum ).c_str( ) ); }
		m_SmileNet.GetNode( (int)i + 1 )->Definition( )->SetNumberOfOutcomes( vecstrOutcomes );
		m_SmileNet.AddArc( 0, (int)i + 1 ); }

	return true; }

/*!
 * \brief
 * Construct a new SMILE-based naive Bayes net with nodes corresponding to the given datasets.
 * 
 * \param pData
 * Datasets from which new Bayes net nodes should be constructed.
 * 
 * \param vecstrNames
 * String identifiers of the newly constructed nodes.
 * 
 * \param veciDefaults
 * Default values (if any) for missing data from each dataset.  -1 is ignored, any other value is used as a
 * default value when data is missing for the corresponding node.
 * 
 * \returns
 * True if Bayes net was successfully constructed.
 * 
 * Constructs a naive Bayes classifier from the given datasets, with one node per dataset plus one
 * additional class node at the top of the naive model.  This class node corresponds to the first dataset
 * in pData and will take two values, generally corresponding to function unrelatedness and relatedness.
 * Each other node is named as indicated and takes the number of discrete values indicated by the dataset.
 * Each value in veciDefaults not equal to -1 is used as a default when data is missing for the corresponding
 * node.
 * 
 * \remarks
 * The order and length of pData, vecstrNames, and veciDefaults must be identical.
 */
bool CBayesNetSmile::Open( const IDataset* pData, const std::vector<std::string>& vecstrNames,
	const vector<size_t>& veciDefaults ) {
	size_t			i, j;
	DSL_stringArray	vecstrOutcomes;
	char			acNum[ 8 ];

	if( pData->GetExperiments( ) != vecstrNames.size( ) )
		return false;

	m_fSmileNet = true;
	m_SmileNet.DeleteAllNodes( );
	m_SmileNet.AddNode( DSL_CPT, (char*)c_szFR );
	vecstrOutcomes.Add( ( (string)c_szFR + "No" ).c_str( ) );
	vecstrOutcomes.Add( ( (string)c_szFR + "Yes" ).c_str( ) );
	m_SmileNet.GetNode( 0 )->Definition( )->SetNumberOfOutcomes( vecstrOutcomes );
	for( i = 1; i < pData->GetExperiments( ); ++i ) {
		m_SmileNet.AddNode( DSL_CPT, (char*)vecstrNames[ i ].c_str( ) );
		vecstrOutcomes.Flush( );
		for( j = 0; j < pData->GetBins( i ); ++j ) {
#pragma warning( disable : 4996 )
			sprintf( acNum, "%02d", j );
#pragma warning( default : 4996 )
			vecstrOutcomes.Add( ( vecstrNames[ i ] + acNum ).c_str( ) ); }
		m_SmileNet.GetNode( (int)i )->Definition( )->SetNumberOfOutcomes( vecstrOutcomes );
		if( veciDefaults[ i ] != -1 ) {
#pragma warning( disable : 4996 )
			sprintf( acNum, "%d", veciDefaults[ i ] );
#pragma warning( default : 4996 )
			m_SmileNet.GetNode( (int)i )->Info( ).UserProperties( ).AddProperty( c_szZero, acNum ); }
		m_SmileNet.AddArc( 0, (int)i ); }

	return true; }

/*!
 * \brief
 * Construct a new SMILE-based naive Bayes net by merging the given class and data nodes.
 * 
 * \param BNPrior
 * Bayes net from which class (root) node is taken.
 * 
 * \param vecpBNs
 * Bayes nets from which data (child) nodes are taken.
 * 
 * \returns
 * True if Bayes net was successfully constructed.
 * 
 * Constructs a new SMILE-based Bayes net by merging the root (prior or class) node from one Bayes net with
 * the child (non-root) nodes from zero or more other networks.  In other words, suppose BNPrior was a naive
 * network with root P1 and children P2 and P3.  vecpBNs contains two networks, one with root A1 and
 * data nodes A2 and A3 and one with root B1 and child node B2.  The newly constructed Bayes net would have
 * a root node with P1's parameters and three children with A2, A3, and B2's parameters.  This can be used
 * to merge multiple naive classifiers created independently from the same answer set.
 * 
 * \remarks
 * In the prior (class) network, only the root (first) node is used.  In the data (child) networks, only
 * the root (first) node is ignored, and the rest are copied into the new network as child nodes.
 */
bool CBayesNetSmile::Open( const CBayesNetSmile& BNPrior, const vector<CBayesNetSmile*>& vecpBNs ) {
	DSL_node*	pFrom;
	size_t		iNet, iNode;
	int			iTo, iProp;

	if( !BNPrior.m_fSmileNet )
		return false;
	for( iNet = 0; iNet < vecpBNs.size( ); ++iNet )
		if( !vecpBNs[ iNet ]->m_fSmileNet )
			return false;

	m_fSmileNet = true;
	m_SmileNet.DeleteAllNodes( );
	pFrom = BNPrior.m_SmileNet.GetNode( 0 );
	m_SmileNet.AddNode( pFrom->Definition( )->GetType( ), pFrom->Info( ).Header( ).GetId( ) );
	m_SmileNet.GetNode( 0 )->Definition( )->SetNumberOfOutcomes( *pFrom->Definition( )->GetOutcomesNames( ) );
	m_SmileNet.GetNode( 0 )->Definition( )->SetDefinition( *pFrom->Definition( )->GetMatrix( ) );
	for( iNet = 0; iNet < vecpBNs.size( ); ++iNet )
		for( iNode = 1; iNode < (size_t)vecpBNs[ iNet ]->m_SmileNet.GetNumberOfNodes( ); ++iNode ) {
			pFrom = vecpBNs[ iNet ]->m_SmileNet.GetNode( iNode );
			m_SmileNet.AddNode( pFrom->Definition( )->GetType( ), pFrom->Info( ).Header( ).GetId( ) );
			m_SmileNet.AddArc( 0, iTo = ( m_SmileNet.GetNumberOfNodes( ) - 1 ) );
			for( iProp = 0; iProp < pFrom->Info( ).UserProperties( ).GetNumberOfProperties( ); ++iProp )
				m_SmileNet.GetNode( iTo )->Info( ).UserProperties( ).AddProperty(
					pFrom->Info( ).UserProperties( ).GetPropertyName( iProp ),
					pFrom->Info( ).UserProperties( ).GetPropertyValue( iProp ) );
			m_SmileNet.GetNode( iTo )->Definition( )->SetNumberOfOutcomes( *pFrom->Definition( )->GetOutcomesNames( ) );
			m_SmileNet.GetNode( iTo )->Definition( )->SetDefinition( *pFrom->Definition( )->GetMatrix( ) ); }

	return true; }

/*!
 * \brief
 * Creates a SMILE Bayes net equivalent to the given minimal naive Bayesian classifier.
 * 
 * \param BNMinimal
 * Minimal naive Bayesian classifier to be copied into the new SMILE network.
 * 
 * \param vecstrNames
 * Node IDs to be assigned to the SMILE network nodes.
 * 
 * \returns
 * True if Bayes net was successfully constructed.
 * 
 * \remarks
 * vecstrNames must contain the same number of strings as BNMinimal has non-root nodes.
 */
bool CBayesNetSmile::Open( const CBayesNetMinimal& BNMinimal, const std::vector<std::string>& vecstrNames ) {
	DSL_stringArray	vecstrOutcomes;
	char			acNum[ 8 ];
	size_t			i, j, k;
	string			strCur;
	DSL_Dmatrix*	pMat;

	m_fSmileNet = true;
	m_SmileNet.DeleteAllNodes( );
	m_SmileNet.AddNode( DSL_CPT, (char*)c_szFR );
	for( i = 0; i < BNMinimal.GetCPT( 0 ).GetRows( ); ++i ) {
#pragma warning( disable : 4996 )
		sprintf( acNum, "%02d", i );
		vecstrOutcomes.Add( ( (string)c_szFR + acNum ).c_str( ) ); }
	m_SmileNet.GetNode( 0 )->Definition( )->SetNumberOfOutcomes( vecstrOutcomes );
	pMat = m_SmileNet.GetNode( 0 )->Definition( )->GetMatrix( );
	for( i = 0; i < BNMinimal.GetCPT( 0 ).GetRows( ); ++i )
		(*pMat)[ i ] = BNMinimal.GetCPT( 0 ).Get( i, 0 );
	for( i = 1; i < BNMinimal.GetNodes( ); ++i ) {
		m_SmileNet.AddNode( DSL_CPT, (char*)( strCur = CMeta::Filename( vecstrNames[ i - 1 ] ) ).c_str( ) );
		vecstrOutcomes.Flush( );
		for( j = 0; j < BNMinimal.GetCPT( i ).GetRows( ); ++j ) {
			sprintf( acNum, "%02d", j );
#pragma warning( default : 4996 )
			vecstrOutcomes.Add( ( strCur + acNum ).c_str( ) ); }
		m_SmileNet.GetNode( (int)i )->Definition( )->SetNumberOfOutcomes( vecstrOutcomes );
		m_SmileNet.AddArc( 0, (int)i );
		pMat = m_SmileNet.GetNode( i )->Definition( )->GetMatrix( );
		for( j = 0; j < BNMinimal.GetCPT( i ).GetColumns( ); ++j )
			for( k = 0; k < BNMinimal.GetCPT( i ).GetRows( ); ++k )
				(*pMat)[ ( j * BNMinimal.GetCPT( i ).GetRows( ) ) + k ] =
					BNMinimal.GetCPT( i ).Get( k, j );
		if( BNMinimal.GetDefault( i ) != 0xFF ) {
			char	acNum[ 16 ];

#pragma warning( disable : 4996 )
			sprintf( acNum, "%d", BNMinimal.GetDefault( i ) );
#pragma warning( default : 4996 )
			m_SmileNet.GetNode( i )->Info( ).UserProperties( ).AddProperty( c_szZero, acNum ); } }
	
	return true; }

}

#endif // NO_SMILE

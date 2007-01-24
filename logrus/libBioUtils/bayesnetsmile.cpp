#include "stdafx.h"
#include "bayesnet.h"
#include "dat.h"
#include "dataset.h"
#include "meta.h"

#if ( defined(_MSC_VER) && defined(_DEBUG) )
extern "C" void __cdecl _invalid_parameter_noinfo( ) { }
#endif // ( defined(_MSC_VER) && defined(_DEBUG) )

namespace libBioUtils {

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

CBayesNetSmile::CBayesNetSmile( bool fGroup ) : CBayesNetSmileImpl( fGroup ) { }

bool CBayesNetSmile::Open( const char* szDSL ) {

	return ( m_fSmileNet = !m_SmileNet.ReadFile( szDSL ) ); }

bool CBayesNetSmile::Save( const char* szDSL ) const {

	if( !m_fSmileNet )
		return false;

	return !((CBayesNetSmile*)this)->m_SmileNet.WriteFile( szDSL ); }

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
				g_CatBioUtils.notice( "CBayesNetSmile::LearnGrouped( %d, %d ) iteration %d, datum %d/%d",
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

bool CBayesNetSmileImpl::FillCPTs( const vector<bool>& vecfHidden, const string& strDatum, bool fZero,
	bool fLearn, bool fAll ) {
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
	bool fZero, bool fLearn ) {
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
			if( iZero == -1 )
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
				g_CatBioUtils.notice( "CBayesNetSmile::LearnUngrouped( %d, %d ) iteration %d, gene %d/%d",
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

void CBayesNetSmile::GetNodes( vector<string>& vecstrNodes ) const {
	int	i;

	if( m_fSmileNet )
		for( i = 0; i < m_SmileNet.GetNumberOfNodes( ); ++i )
			vecstrNodes.push_back( m_SmileNet.GetNode( i )->Info( ).Header( ).GetId( ) ); }

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

					for( k = 0; ( k + 1 ) < (size_t)pValue->GetSize( ); ++k )
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
	for( i = 0; ( i + 1 ) < (size_t)pValue->GetSize( ); ++i )
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
					vecveciCounts[ k ][ ( iVal * iAnswers ) + iAnswer ]++; } }

	fFallback = m_pDefaults && ( iCount < c_iMinimum );
	pMat = m_SmileNet.GetNode( 0 )->Definition( )->GetMatrix( );
	for( i = 0; i < iAnswers; ++i )
		(*pMat)[ (int)i ] = ( j = vecveciCounts[ 0 ][ (int)i ] ) ? j : ( fFallback ? 0 : 1 );
	if( fFallback ) {
		g_CatBioUtils.warn( "CBayesNetSmile::LearnNaive( %d ) insufficient data for node %s",
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
					g_CatBioUtils.warn( "CBayesNetSmile::LearnNaive( %d ) insufficient data for node %s, column %d",
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

void CBayesNetSmile::SetDefault( const CBayesNetSmile& Defaults ) {

	m_pDefaults = &Defaults; }

bool CBayesNetSmile::GetCPT( size_t iNode, CDataMatrix& MatCPT ) const {

	return CBayesNetSmileImpl::GetCPT( m_SmileNet.GetNode( (int)iNode ), MatCPT ); }

bool CBayesNetSmile::Evaluate( const CPCLPair& PCLIn, CPCL& PCLOut, bool fZero, int iAlgorithm ) const {
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
		for( j = 0; j < PCLIn.GetExperiments( ); ++j )
			if( PCLIn.GetExperiment( j ) == m_SmileNet.GetNode( (int)i )->Info( ).Header( ).GetId( ) ) {
				vecfHidden[ i ] = false;
				veciMap[ i ] = (unsigned int)j;
				break; } }
	((CBayesNetSmile*)this)->m_SmileNet.SetDefaultBNAlgorithm( iAlgorithm );
	for( i = 0; i < PCLOut.GetGenes( ); ++i ) {
		if( !( i % 1 ) )
			g_CatBioUtils.notice( "CBayesNetSmile::Evaluate( %d ) %d/%d", fZero, i,
				PCLOut.GetGenes( ) );
		strCur = EncodeDatum( PCLIn, PCLIn.GetGene( PCLOut.GetGene( i ) ), veciMap );
		if( m_fGroup && ( ( iterDatum = mapData.find( strCur ) ) != mapData.end( ) ) ) {
			for( j = 0; j < iterDatum->second.size( ); ++j )
				PCLOut.Set( i, j, iterDatum->second[ j ] );
			continue; }

		((CBayesNetSmile*)this)->FillCPTs( vecfHidden, strCur, fZero, false, true );
		((CBayesNetSmile*)this)->m_SmileNet.UpdateBeliefs( );
		for( iExp = j = 0; j < veciMap.size( ); ++j ) {
			DSL_Dmatrix*	pMatrix;

			if( veciMap[ j ] != -1 )
				continue;
			pMatrix = m_SmileNet.GetNode( (int)j )->Value( )->GetMatrix( );
			for( k = 0; k < GetValues( j ); ++k )
				PCLOut.Set( i, iExp++, (float)(*pMatrix)[ (int)k ] ); }
		if( m_fGroup ) {
			vector<float>	vecfCur;

			vecfCur.resize( PCLOut.GetExperiments( ) );
			for( j = 0; j < vecfCur.size( ); ++j )
				vecfCur[ j ] = PCLOut.Get( i, j );
			mapData[ strCur ] = vecfCur; } }
	((CBayesNetSmile*)this)->m_SmileNet.SetDefaultBNAlgorithm( iPrev );

	return true; }

bool CBayesNetSmile::Open( const vector<string>& vecstrPCLs, size_t iBins ) {
	size_t			i, j;
	DSL_stringArray	vecstrOutcomes;
	string			strCur;

	m_fSmileNet = true;
	m_SmileNet.DeleteAllNodes( );
	m_SmileNet.AddNode( DSL_CPT, (char*)c_szFR );
	vecstrOutcomes.Add( ( (string)c_szFR + "No" ).c_str( ) );
	vecstrOutcomes.Add( ( (string)c_szFR + "Yes" ).c_str( ) );
	m_SmileNet.GetNode( 0 )->Definition( )->SetNumberOfOutcomes( vecstrOutcomes );
	for( i = 0; i < vecstrPCLs.size( ); ++i ) {
		m_SmileNet.AddNode( DSL_CPT, (char*)( strCur =
			CMeta::Filename( CMeta::Deextension( vecstrPCLs[ i ] ) ) ).c_str( ) );
		vecstrOutcomes.Flush( );
		for( j = 0; j < iBins; ++j ) {
			char	acNum[ 8 ];

#pragma warning( disable : 4996 )
			sprintf( acNum, "%02d", j );
#pragma warning( default : 4996 )
			vecstrOutcomes.Add( ( strCur + acNum ).c_str( ) ); }
		m_SmileNet.GetNode( (int)i + 1 )->Definition( )->SetNumberOfOutcomes( vecstrOutcomes );
		m_SmileNet.AddArc( 0, (int)i + 1 ); }

	return true; }

bool CBayesNetSmile::Open( const IDataset* pData, const vector<string>& vecstrNames,
	const vector<size_t>& veciZeros ) {
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
		if( veciZeros[ i ] != -1 ) {
#pragma warning( disable : 4996 )
			sprintf( acNum, "%d", veciZeros[ i ] );
#pragma warning( default : 4996 )
			m_SmileNet.GetNode( (int)i )->Info( ).UserProperties( ).AddProperty( c_szZero, acNum ); }
		m_SmileNet.AddArc( 0, (int)i ); }

	return true; }

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

}

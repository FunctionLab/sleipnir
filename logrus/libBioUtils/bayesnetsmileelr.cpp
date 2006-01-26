#include "stdafx.h"
#include "bayesnet.h"
#include "dataset.h"

namespace libBioUtils {

bool CBayesNetSmileImpl::LearnELR( const IDataset* pData, size_t iIterations, bool fZero ) {
	size_t			i, j, iParameters;
	TVecVecF		vecvecfBeta, vecvecfGradient, vecvecfPrev, vecvecfDirection, vecvecfOriginal;
	DSL_Dmatrix*	pMatrix;
	DSL_intArray	veciCoords;
	float			dAX, dBX, dDotNew, dDotOld, dDotON, dB;

	if( !m_fSmileNet )
		return false;

	iParameters = ELRCountParameters( );
	vecvecfBeta.resize( m_SmileNet.GetNumberOfNodes( ) );
	vecvecfGradient.resize( m_SmileNet.GetNumberOfNodes( ) );
	vecvecfPrev.resize( m_SmileNet.GetNumberOfNodes( ) );
	vecvecfDirection.resize( m_SmileNet.GetNumberOfNodes( ) );
	vecvecfOriginal.resize( m_SmileNet.GetNumberOfNodes( ) );
	for( i = 0; i < vecvecfBeta.size( ); ++i ) {
		pMatrix = m_SmileNet.GetNode( (int)i )->Definition( )->GetMatrix( );
		vecvecfOriginal[ i ].resize( pMatrix->GetSize( ) );
		vecvecfDirection[ i ].resize( pMatrix->GetSize( ) );
		vecvecfPrev[ i ].resize( pMatrix->GetSize( ) );
		vecvecfGradient[ i ].resize( pMatrix->GetSize( ) );
		vecvecfBeta[ i ].resize( pMatrix->GetSize( ) );
		for( j = 0; j < vecvecfBeta[ i ].size( ); ++j )
			vecvecfBeta[ i ][ j ] = (float)log( (*pMatrix)[ (int)j ] ); }
	dAX = 0;
	dBX = 0.1f;

	copy( vecvecfBeta.begin( ), vecvecfBeta.end( ), vecvecfOriginal.begin( ) );
	ELRComputeGradient( pData, fZero, vecvecfGradient );
	copy( vecvecfGradient.begin( ), vecvecfGradient.end( ), vecvecfPrev.begin( ) );
	copy( vecvecfGradient.begin( ), vecvecfGradient.end( ), vecvecfDirection.begin( ) );
	dDotNew = ELRDot( vecvecfGradient, vecvecfPrev );

	for( i = 0; i < iIterations; ++i ) {
		g_CatBioUtils.notice( "CBayesNetSmile::LearnELR( %d, %d ) iteration %d/%d",
			iIterations, fZero, i, iIterations );
		if( !dDotNew )
			continue;
		ELRNormalizeDirection( vecvecfDirection );
		ELRLineSearch( pData, vecvecfDirection, vecvecfOriginal, vecvecfBeta, dAX, dBX, fZero );
		copy( vecvecfBeta.begin( ), vecvecfBeta.end( ), vecvecfOriginal.begin( ) );
		ELRCopyParameters( vecvecfBeta );
		ELRComputeGradient( pData, fZero, vecvecfGradient );

		dDotOld = dDotNew;
		dDotON = ELRDot( vecvecfGradient, vecvecfPrev );
		copy( vecvecfGradient.begin( ), vecvecfGradient.end( ), vecvecfPrev.begin( ) );
		dDotNew = ELRDot( vecvecfGradient, vecvecfPrev );
		dB = ( dDotNew - dDotON ) / dDotOld;
		if( !( ( i + 1 ) % iParameters ) || ( dB <= 0 ) )
			copy( vecvecfGradient.begin( ), vecvecfGradient.end( ), vecvecfDirection.begin( ) );
		else
			ELRComputeDirection( dB, vecvecfGradient, vecvecfDirection ); }

	return true; }

size_t CBayesNetSmileImpl::ELRCountParameters( ) const {
	size_t	i, iRet;

	for( iRet = i = 0; i < m_SmileNet.GetNumberOfNodes( ); ++i )
		iRet += m_SmileNet.GetNode( (int)i )->Definition( )->GetNumberOfOutcomes( );

	return iRet; }

void CBayesNetSmileImpl::ELRCopyParameters( TVecVecF& vecvecfBeta ) {
	static const float	c_dMinimum	= log( FLT_MIN );
	size_t				i, j, k, iIndex, iDomains;
	float				dSum, dCur, dMax;
	DSL_nodeDefinition*	pDef;
	DSL_Dmatrix*		pDefs;

	for( i = 0; i < vecvecfBeta.size( ); ++i ) {
		pDef = m_SmileNet.GetNode( (int)i )->Definition( );
		pDefs = pDef->GetMatrix( );
		iDomains = pDef->GetNumberOfOutcomes( );
		for( j = 0; j < vecvecfBeta[ i ].size( ); j += iDomains ) {
			dSum = 0;
			dMax = vecvecfBeta[ i ][ j ];
			for( k = 1; k < iDomains; ++k )
				if( ( dCur = vecvecfBeta[ i ][ j + k ] ) > dMax )
					dMax = dCur;
			for( k = 0; k < iDomains; ++k ) {
				iIndex = j + k;
				if( ( vecvecfBeta[ i ][ iIndex ] -= dMax ) < c_dMinimum )
					vecvecfBeta[ i ][ iIndex ] = c_dMinimum;
				dSum += exp( vecvecfBeta[ i ][ iIndex ] ); }
			if( dSum )
				for( k = 0; k < iDomains; ++k ) {
					iIndex = j + k;
					(*pDefs)[ (int)iIndex ] = exp( vecvecfBeta[ i ][ iIndex ] ) / dSum; } } } }

// EXPENSIVE
void CBayesNetSmileImpl::ELRComputeGradient( const IDataset* pData, bool fZero, TVecVecF& vecvecfGradient ) {
	size_t	i, j, iVal;

	for( i = 0; i < vecvecfGradient.size( ); ++i )
		for( j = 0; j < vecvecfGradient[ i ].size( ); ++j )
			vecvecfGradient[ i ][ j ] = 0;

	for( i = 0; i < pData->GetGenes( ); ++i )
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j ) {
			if( ( ( iVal = pData->GetDiscrete( i, j, 0 ) ) == -1 ) || !FillCPTs( pData, i, j, fZero, false ) )
				continue;
			ELRUpdateGradient( -1, vecvecfGradient );
			m_SmileNet.GetNode( 0 )->Value( )->SetEvidence( (int)iVal );
			ELRUpdateGradient( 1, vecvecfGradient ); } }

void CBayesNetSmileImpl::ELRUpdateGradient( float dRate, TVecVecF& vecvecfGradient ) {
	size_t				i, j;
	double				dSum, dPosterior;
	int					k, iTmp, iDomain, iDomains, iEvidence;
	DSL_node*			pNode;
	DSL_nodeValue*		pValue;
	DSL_nodeDefinition*	pDef;
	DSL_Dmatrix*		pValues;
	DSL_Dmatrix*		pDefs;
	DSL_intArray		veciParents, veciCoords;

	m_SmileNet.UpdateBeliefs( );
	for( i = 0; i < vecvecfGradient.size( ); ++i ) {
		pNode = m_SmileNet.GetNode( (int)i );
		veciParents = pNode->Parents( );
		pValue = pNode->Value( );
		iEvidence = pValue->GetEvidence( );
		pValues = pValue->GetMatrix( );
		pDef = pNode->Definition( );
		iDomains = pDef->GetNumberOfOutcomes( );
		pDefs = pDef->GetMatrix( );

		dSum = 0;
		for( j = 0,pDefs->IndexToCoordinates( iTmp = 0, veciCoords ); iTmp != DSL_OUT_OF_RANGE;
			++j,iTmp = pDefs->NextCoordinates( veciCoords ) ) {
			iDomain = veciCoords[ veciCoords.GetSize( ) - 1 ];
			if( veciParents.NumItems( ) ) {
				if( iEvidence == DSL_OUT_OF_RANGE )
					dPosterior = (*pValues)[ iDomain ];
				else if( iDomain == iEvidence ) {
					dPosterior = 1;
					for( k = 0; k < veciParents.NumItems( ); ++k )
						dPosterior *= m_SmileNet.GetNode( veciParents[ k ] )->Value(
							)->GetMatrix( )->Subscript( veciCoords[ k ] ); }
				else
					dPosterior = 0; }
			else
				dPosterior = (*pValues)[ veciCoords[ 0 ] ];
			dSum += dPosterior;
			vecvecfGradient[ i ][ j ] += dRate * (float)dPosterior;

			if( ( iDomain + 1 ) == iDomains ) {
				for( k = 0; k < iDomains; ++k ) {
					veciCoords[ veciCoords.GetSize( ) - 1 ] = k;
					vecvecfGradient[ i ][ k + j - iDomains + 1 ] -=
						(float)( dRate * dSum * (*pDefs)[ veciCoords ] ); }
				dSum = 0; } } } }

float CBayesNetSmileImpl::ELRDot( const TVecVecF& vecvecfOne, const TVecVecF& vecvecfTwo ) {
	size_t	i, j;
	float	dRet;

	dRet = 0;
	for( i = 0; i < vecvecfOne.size( ); ++i )
		for( j = 0; j < vecvecfOne[ i ].size( ); ++j )
			dRet += vecvecfOne[ i ][ j ] * vecvecfTwo[ i ][ j ];

	return dRet; }

void CBayesNetSmileImpl::ELRNormalizeDirection( TVecVecF& vecvecfDirection ) const {
	size_t	i, j;
	float	d, dSum;
	int		k, iDomains;

	for( i = 0; i < vecvecfDirection.size( ); ++i ) {
		iDomains = m_SmileNet.GetNode( (int)i )->Definition( )->GetNumberOfOutcomes( );
		for( j = 0; j < vecvecfDirection[ i ].size( ); j += iDomains ) {
			dSum = 0;
			for( k = 0; k < iDomains; ++k ) {
				d = vecvecfDirection[ i ][ j + k ];
				dSum += d * d; }

			if( dSum = sqrt( dSum ) )
				for( k = 0; k < iDomains; ++k )
					vecvecfDirection[ i ][ j + k ] /= dSum; } } }

// EXPENSIVE
float CBayesNetSmileImpl::ELRLineSearch( const IDataset* pData, const TVecVecF& vecvecfDirection,
	const TVecVecF& vecvecfOriginal, TVecVecF& vecvecfBeta, float& dAX, float& dBX, bool fZero ) {
	float	dFA, dFB, dFC, dCX;

	ELRBracket( pData, vecvecfDirection, vecvecfOriginal, vecvecfBeta, dAX, dBX, dCX, dFA, dFB, dFC, fZero );
	return ELRBrent( pData, vecvecfDirection, vecvecfOriginal, vecvecfBeta, dAX, dBX, dCX, dFA, dFB, dFC, fZero ); }

float CBayesNetSmileImpl::ELREvalFunction( const IDataset* pData, float dX, const TVecVecF& vecvecfDirection,
	const TVecVecF& vecvecfOriginal, TVecVecF& vecvecfBeta, bool fZero  ) {
	size_t	i, j;

	for( i = 0; i < vecvecfBeta.size( ); ++i )
		for( j = 0; j < vecvecfBeta[ i ].size( ); ++j )
			vecvecfBeta[ i ][ j ] = ( dX * vecvecfDirection[ i ][ j ] ) + vecvecfOriginal[ i ][ j ];
	ELRCopyParameters( vecvecfBeta );

	return -ELRConditionalLikelihood( pData, fZero ); }

float CBayesNetSmileImpl::ELRAvoidZero( float d ) {
	static const float	c_dTiny	= 1e-10f;

	if( d >= 0 ) {
		if( d < c_dTiny )
			return c_dTiny; }
	else if( d > -c_dTiny )
		return -c_dTiny;

	return d; }

void CBayesNetSmileImpl::ELRBracket( const IDataset* pData, const TVecVecF& vecvecfDirection,
	const TVecVecF& vecvecfOriginal, TVecVecF& vecvecfBeta, float& dAX, float& dBX, float& dCX, float& dFA,
	float& dFB, float& dFC, bool fZero ) {
	static const float	c_dGolden		= 1.618034f;
	static const float	c_dLimit		= 100;
	static const size_t	c_iIterations	= 100;
	size_t	i;
	float	dR, dQ, dU, dFU, dULim;

	dFA = ELREvalFunction( pData, dAX, vecvecfDirection, vecvecfOriginal, vecvecfBeta, fZero );
	dFB = ELREvalFunction( pData, dBX, vecvecfDirection, vecvecfOriginal, vecvecfBeta, fZero );
	if( dFB > dFA ) {
		swap( dAX, dBX );
		swap( dFA, dFB ); }

	dCX = dBX + ( c_dGolden * ( dBX - dAX ) );
	dFC = ELREvalFunction( pData, dCX, vecvecfDirection, vecvecfOriginal, vecvecfBeta, fZero );
	for( i = 0; ( i < c_iIterations ) && ( dFB > dFC ); ++i ) {
		/* finding the minimum point dU of the extrapolated parabola*/
		dR = ( dBX - dAX ) * ( dFB - dFC );
		dQ = ( dBX - dCX ) * ( dFB - dFA );
		dU = dBX - ( ( ( ( dBX - dCX ) * dQ ) - ( ( dBX - dAX ) * dR ) ) / ( 2 * ELRAvoidZero( dQ - dR ) ) );
		dULim = dBX + ( c_dLimit * ( dCX - dBX ) );

		/* bx>u and u>cx, or u is between b and c */
		if( ( ( dBX - dU ) * ( dU - dCX ) ) > 0 ) {
			dFU = ELREvalFunction( pData, dU, vecvecfDirection, vecvecfOriginal, vecvecfBeta, fZero );

			/* the minimum is between b and c, because c is lower than b,
			 * otherwise we can not be in the loop. So u is lower than both b
			 * and c. We are done. */
			if( dFU < dFC ) {
				dAX = dBX;
				dFA = dFB;
				dBX = dU;
				dFB = dFU;
				break; }
			/* b is lower than both a and u. we are done*/
			else if( dFU > dFB ) {
				dCX = dU;
				dFC = dFU;
				break; }
			/* u is lower than b, but higher than c. useless.
			 * extend it beyound c to prepare for the next iteration
			 * u will be copied to c. */

			dU = dCX + ( c_dGolden * ( dCX - dBX ) );
			dFU = ELREvalFunction( pData, dU, vecvecfDirection, vecvecfOriginal, vecvecfBeta, fZero ); }
		/* u is between c and the far limit */
		else if( ( ( dCX - dU ) * ( dU - dULim ) ) > 0 ) {
			dFU = ELREvalFunction( pData, dU, vecvecfDirection, vecvecfOriginal, vecvecfBeta, fZero );
			/* if u is still lower than c, extend the points beyond u and shift the values */
			if( dFU < dFC ) {
				dBX = dCX;
				dCX = dU;
				dU = dCX + ( c_dGolden * ( dCX - dBX ) );
				dFB = dFC;
				dFC = dFU;
				dFU = ELREvalFunction( pData, dU, vecvecfDirection, vecvecfOriginal, vecvecfBeta, fZero ); } }
		/* u is beyond the far limit. ulim is between u and c. c and ulim defines the direction. */
		else if( ( ( dU - dULim ) * ( dULim - dCX ) ) > 0 )
			dFU = ELREvalFunction( pData, dU = dULim, vecvecfDirection, vecvecfOriginal, vecvecfBeta, fZero );
		/* other cases... what could be left? */
		else
			dFU = ELREvalFunction( pData, dU = dCX + ( c_dGolden * ( dCX - dBX ) ), vecvecfDirection,
				vecvecfOriginal, vecvecfBeta, fZero );

		/*eliminate oldest point and repeat the process*/
		dAX = dBX;
		dBX = dCX;
		dCX = dU;
		dFA = dFB;
		dFB = dFC;
		dFC = dFU; } }

float CBayesNetSmileImpl::ELRBrent( const IDataset* pData, const TVecVecF& vecvecfDirection,
	const TVecVecF& vecvecfOriginal, TVecVecF& vecvecfBeta, float& dAX, float& dBX, float dCX, float dFA,
	float dFB, float dFC, bool fZero ) {
	static const size_t	c_iIterations	= 100;
	static const float	c_dZEPS			= 1e-8f;
	static const float	c_dIGolden		= 0.3819660f;
	static const float	c_dTolerance	= 1e-4f;
	/* d, e, etemp are the step. e is the step before d, etemp is the step
	 * before e
	 * tol,tol1,tol2 are the error tolerance. for simple thinking, consider them 0. */
	float	a,b,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u=0,v,w,x,xm,e,d;
	size_t	iter;

	/* a, b brackets the minimum. ensure a<=b*/
	if( dAX < dCX ) {
		a = dAX;
		b = dCX; }
	else {
		a = dCX;
		b = dAX; }

	/* initialize all three lowest points to the middle variable of the bracket*/
	e = d = b - a;
	x = dBX;
	fx = dFB;
	if( dFA < dFC ) {
		w = dAX;
		fw = dFA;
		v = dCX;
		fv = dFC; }
	else {
		w = dCX;
		fw = dFC;
		v = dAX;
		fv = dFA; }

	for( iter = 0; iter < c_iIterations; ++iter ) {
		xm = 0.5f * ( a + b );
		/* error tolerance:2*tol*x.  note that tol, tol1, tol2 are all positive*/
		tol2 = 2 * ( tol1 = ( ( c_dTolerance * abs( x ) ) + c_dZEPS ) );
		/* the sum of x - xm and 0.5*(b-a) smaller than tolerance*/
		if( ( fw - fx ) <= ( c_dTolerance * fx ) ) {
			dAX = a;
			dBX = b;
			return fx; }

		/* e is used to check to see whether it is ok to take a parabolic step
		 * parabolic step will not be considered if the step before previous
		 * step was too small , made  no progress.
		 * I'm not sure when could this happen... if progress is too small,
		 * it should be done.... unless it's the first time */
		if( abs( e ) > tol1 ) {
			/* fit a parabola*/
			r = ( x - w ) * ( fx - fv );
			q = ( x - v ) * ( fx - fw );
			p = ( ( x - v ) * q ) - ( ( x - w ) * r );
			q = 2 * ( q - r );

			/*ensure q is positive, and reverse the ratio p/q
			 * these extra steps are used to avoid dividing by zero, making the
			 * denominator positive ensure the later inequality holds because it is
			 * multiplied by a positive number */
			if( q > 0 )
				p *= -1;
			else
				q = -q;

			etemp = e;
			/*the previouse step size*/
			e = d;

			/* after this iteration, e will be either the previous step size
			 * (if the current step is a parabolic step), or the size of the
			 * current larger segment, (if the current step is a GOLD step).
			 * therefore, on the next iteration, a parabolic step will be
			 * taken only if the parabolic step is smaller than the
			 * previous step, or smaller than the current larger segment */

			/* the first part test if the current step is greater than the step
			 * before the previous step. If yes, it may not be making progress
			 * and may not converge using parabola, so use GOLD step instead.
			 *
			 * the 2nd part test the new minimum of the parabola x+p/q (which
			 * is the old x-p/q) is outside of (a,b)
			 * e is set to be the larger of the 2 segments */
			if( ( abs( p ) >= abs( 0.5 * q * etemp ) ) ||
				( ( p <= ( q * ( a - x ) ) ) || ( p >= ( q * ( b - x ) ) ) ) )
				d = c_dIGolden * ( e = ( x >= xm ? ( a - x ) : ( b - x ) ) );
			/* parabola is good. update the minimum u */
			else {
				/* parabolic step*/
				d = p / q;
				/* notice that u's value is not fixed yet, it's only set for
				 * testing purpose now. it will be fixed later
				 * Bad coding. Here should really make the code cleaner by using another
				 * variable, say uTest, or boundaryTest. */
				u = x + d;

				/* if the minimum is too close to a or b, step away by tol1 to
				 * ensure there is progress made.*/
				if( ( u - a ) < tol2 )
					d += tol1;
				else if( ( b - u ) < tol2 )
					d -= tol1; } }
		/* take a GOLD step in case of problem*/
		else
			d = c_dIGolden * ( e = ( ( x > xm ) ? ( a - x ) : ( b - x ) ) );

		/* update the new minimum
		 * if d<tolerance, step by tolerance, because otherwise it would give
		 * the "same" result and make no progress */
		if( abs( d ) >=tol1 )
			u = x + d;
		else {
			if( d > 0 )
				u = x + tol1;
			else
				u = x - tol1; }
		fu = ELREvalFunction( pData, u, vecvecfDirection, vecvecfOriginal, vecvecfBeta, fZero );

		/* u is the lowest point of the function,
		 * lower than the previous lowest point x.
		 * so make x one of the boundary points.
		 * then old point v will be discarded.
		 * and one of old a,b will be discarded, and will be the 2nd lowest
		 * point w, which was the old lowest point x.
		 * so w will be one of the boundary points a,b. */
		if( fu < fx ) {
			if( u >= x )
				a = x;
			else
				b = x;
			/*update three lowest points x,v,w*/
			v = w;
			w = x;
			x = u;
			fv = fw;
			fw = fx;
			fx = fu; }
		/* x is still the lowest point*/
		else {
			/* we know u is between a and b from before,
			 * so we can reduce the bracket */
			if( u < x )
				a = u;
			else
				b = u;

			/* if u is lower than 2nd lowest point w, make u the 2nd lowest
			 * point, and w the 3rd
			 * if the lowest point is same as the 2nd lowest point */
			if( ( fu <= fw ) || ( w == x ) ) {
				v = w;
				fv = fw;
				w = u;
				fw = fu; }
			/* if u is lower than the 3rd lowest point v,
			 * or if the 3rd lowest point is same as the lowest or 2nd lowest point */
			else if( ( fu <= fv ) || ( v == x ) || ( v == w ) ) {
				v = u;
				fv = fu; } } }

	g_CatBioUtils.warn( "CBayesNetSmileImpl::ELRBrent( ) too many iterations" );
	dAX = a;
	dBX = b;
	return fx; }

float CBayesNetSmileImpl::ELRConditionalLikelihood( const IDataset* pData, bool fZero ) {
	size_t			i, j, iCount;
	float			dRet;
	DSL_Dmatrix*	pMatrix;

	dRet = 0;
	for( iCount = i = 0; i < pData->GetGenes( ); ++i )
		for( j = ( i + 1 ); j < pData->GetGenes( ); ++j ) {
			if( ( pData->GetDiscrete( i, j, 0 ) == -1 ) || !FillCPTs( pData, i, j, fZero, false ) )
				continue;
			iCount++;
			m_SmileNet.UpdateBeliefs( );
			pMatrix = m_SmileNet.GetNode( 0 )->Value( )->GetMatrix( );
			dRet += (float)log( (*pMatrix)[ (int)pData->GetDiscrete( i, j, 0 ) ] ); }

	return ( dRet / iCount ); }

void CBayesNetSmileImpl::ELRComputeDirection( float dB, const TVecVecF& vecvecfGradient,
	TVecVecF& vecvecfDirection ) {
	size_t	i, j;

	for( i = 0; i < vecvecfDirection.size( ); ++i )
		for( j = 0; j < vecvecfDirection[ i ].size( ); ++j )
			vecvecfDirection[ i ][ j ] = vecvecfGradient[ i ][ j ] + ( dB * vecvecfDirection[ i ][ j ] ); }

}

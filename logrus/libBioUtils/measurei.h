#ifndef MEASUREI_H
#define MEASUREI_H

#include "fullmatrix.h"

namespace libBioUtils {

class IMeasure;

class CMeasureImpl {
protected:
	CMeasureImpl( const IMeasure* );

	const IMeasure*	m_pMeasure;
};

class CMeasureSigmoidImpl : protected CMeasureImpl {
protected:
	CMeasureSigmoidImpl( const IMeasure*, float );

	float	m_dMult;
};

class CMeasureSpearmanImpl {
protected:
	CMeasureSpearmanImpl( bool );

	bool	m_fTransformed;
};

class CMeasureKendallsTauImpl {
protected:
	struct SKendallsFirst {
		const float*	m_adX;
		const float*	m_adY;

		SKendallsFirst( const float* adX, const float* adY ) : m_adX(adX), m_adY(adY) { }

		bool operator( )( size_t iX, size_t iY ) const {
			float	dX, dY;

			dX = m_adX[ iX ];
			dY = m_adX[ iY ];
			if( dX > dY )
				return true;
			if( dX < dY )
				return false;

			dX = m_adY[ iX ];
			dY = m_adY[ iY ];
			return ( dX > dY ); }
	};

	struct SKendallsSecond {
		const float*	m_adX;
		const float*	m_adY;

		SKendallsSecond( const float* adX, const float* adY ) : m_adX(adX), m_adY(adY) { }

		char operator( )( size_t iX, size_t iY ) const {
			float	dX, dY;

			dX = m_adY[ iX ];
			dY = m_adY[ iY ];
			return ( ( dX > dY ) ? -1 : ( ( dX < dY ) ? 1 : 0 ) ); }
	};

	static double MeasureWeighted( const float*, const float*, size_t, const float*,
		const float* );
	static double MeasureUnweighted( const float*, const float*, size_t );
	static size_t CountExchanges( size_t*, size_t, size_t*, const SKendallsSecond&,
		size_t = 0 );
};

}

#endif // MEASUREI_H

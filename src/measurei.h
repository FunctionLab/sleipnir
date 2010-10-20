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
#ifndef MEASUREI_H
#define MEASUREI_H

#include "fullmatrix.h"

namespace Sleipnir {

class CMeasureImpl {
protected:
	friend class CMeasureKendallsTau;
	friend class CMeasureKolmogorovSmirnov;
	friend class CMeasureSpearman;

	static double MeasureTrim( const IMeasure*, const float*, size_t, const float*, size_t, const IMeasure::EMap,
		const float*, const float*, bool );
	static bool IsNaN( const float*, size_t );

	CMeasureImpl( const IMeasure*, bool );
	virtual ~CMeasureImpl( );

	IMeasure*	m_pMeasure;
	bool		m_fMemory;
};

class CMeasureSigmoidImpl : protected CMeasureImpl {
protected:
	CMeasureSigmoidImpl( const IMeasure* pMeasure, bool fMemory, float dMult ) :
		m_dMult(dMult), CMeasureImpl( pMeasure, fMemory ) { }

	float	m_dMult;
};

class CMeasureSpearmanImpl {
protected:
	CMeasureSpearmanImpl( bool fTransformed ) : m_fTransformed(fTransformed) { }

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

class CMeasurePearNormImpl {
protected:
	CMeasurePearNormImpl( double dAve, double dStd ) : m_dAverage(dAve), m_dStdDev(dStd) { }

	double	m_dAverage;
	double	m_dStdDev;
};

}

#endif // MEASUREI_H

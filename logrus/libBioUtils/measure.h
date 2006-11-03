#ifndef MEASURE_H
#define MEASURE_H

namespace libBioUtils {

class IMeasure {
public:
	enum EMap {
		EMapNone	= 0,
		EMapCenter	= EMapNone + 1,
		EMapAbs		= EMapCenter + 1
	};

	virtual const char* GetName( ) const = 0;
	virtual bool IsRank( ) const = 0;
	virtual double Measure( const float*, size_t, const float*, size_t, EMap = EMapCenter,
		const float* = NULL, const float* = NULL ) const = 0;
	virtual IMeasure* Clone( ) const = 0;
};

}

#include "measurei.h"

namespace libBioUtils {

class CMeasureSigmoid : CMeasureSigmoidImpl, public IMeasure {
public:
	CMeasureSigmoid( const IMeasure*, bool, float dDiv = 1 );

	const char* GetName( ) const;
	bool IsRank( ) const;
	double Measure( const float*, size_t, const float*, size_t, EMap, const float*,
		const float* ) const;
	IMeasure* Clone( ) const;
};

class CMeasureNegate : CMeasureImpl, public IMeasure {
public:
	CMeasureNegate( const IMeasure*, bool );

	const char* GetName( ) const;
	bool IsRank( ) const;
	double Measure( const float*, size_t, const float*, size_t, EMap, const float*,
		const float* ) const;
	IMeasure* Clone( ) const;
};

class CMeasureAutocorrelate : CMeasureImpl, public IMeasure {
public:
	CMeasureAutocorrelate( const IMeasure*, bool );

	const char* GetName( ) const;
	bool IsRank( ) const;
	double Measure( const float*, size_t, const float*, size_t, EMap, const float*,
		const float* ) const;
	IMeasure* Clone( ) const;
};

class CMeasureEuclidean : public IMeasure {
public:
	const char* GetName( ) const;
	bool IsRank( ) const;
	double Measure( const float*, size_t, const float*, size_t, EMap, const float*,
		const float* ) const;
	IMeasure* Clone( ) const;
};

class CMeasurePearson : public IMeasure {
public:
	static double Pearson( const float*, size_t, const float*, size_t, EMap, const float*,
		const float* );

	const char* GetName( ) const;
	bool IsRank( ) const;
	double Measure( const float*, size_t, const float*, size_t, EMap, const float*,
		const float* ) const;
	IMeasure* Clone( ) const;
};

class CMeasureQuickPearson : public IMeasure {
public:
	const char* GetName( ) const;
	bool IsRank( ) const;
	double Measure( const float*, size_t, const float*, size_t, EMap, const float*,
		const float* ) const;
	IMeasure* Clone( ) const;
};

class CMeasureKolmogorovSmirnov : public IMeasure {
public:
	const char* GetName( ) const;
	bool IsRank( ) const;
	double Measure( const float*, size_t, const float*, size_t, EMap, const float*,
		const float* ) const;
	IMeasure* Clone( ) const;
};

class CMeasureKendallsTau : CMeasureKendallsTauImpl, public IMeasure {
public:
	const char* GetName( ) const;
	bool IsRank( ) const;
	double Measure( const float*, size_t, const float*, size_t, EMap, const float*,
		const float* ) const;
	IMeasure* Clone( ) const;
};

class CMeasureSpearman : CMeasureSpearmanImpl, public IMeasure {
public:
	CMeasureSpearman( bool = false );

	const char* GetName( ) const;
	bool IsRank( ) const;
	double Measure( const float*, size_t, const float*, size_t, EMap, const float*,
		const float* ) const;
	IMeasure* Clone( ) const;
};

class CMeasurePearNorm : CMeasurePearNormImpl, public IMeasure {
public:
	CMeasurePearNorm( );
	CMeasurePearNorm( double, double );

	const char* GetName( ) const;
	bool IsRank( ) const;
	double Measure( const float*, size_t, const float*, size_t, EMap, const float*,
		const float* ) const;
	IMeasure* Clone( ) const;
};

class CMeasureHypergeometric : public IMeasure {
public:
	const char* GetName( ) const;
	bool IsRank( ) const;
	double Measure( const float*, size_t, const float*, size_t, EMap, const float*,
		const float* ) const;
	IMeasure* Clone( ) const;
};

}

#endif // MEASURE_H

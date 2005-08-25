#ifndef EXAMPLEI_H
#define EXAMPLEI_H

#include <utility>

#include "datapair.h"

namespace libBioUtils {

class CExampleImpl {
public:
	CExampleImpl( );
	~CExampleImpl( );

	void Set( size_t, float, const CDataPair&, size_t );
	bool Equals( const CExampleImpl&, size_t ) const;
	size_t GetDiscrete( size_t ) const;
	float GetContinuous( size_t ) const;
	void Reset( );
	bool IsEvidence( bool, size_t ) const;
	bool IsSet( ) const;

protected:
	union UDatum {
		float	m_d;
		size_t	m_i;
	};

	UDatum*	m_auData;
};

}

#endif // EXAMPLEI_H

#ifndef EXAMPLEI_H
#define EXAMPLEI_H

#include <utility>

#include "datapair.h"

namespace Sleipnir {

class CExampleImpl {
public:
	CExampleImpl( );
	~CExampleImpl( );

	void Set( size_t, float, const CDataPair&, size_t );
	bool Equals( const CExampleImpl&, size_t ) const;
	void Reset( );
	bool IsEvidence( size_t ) const;

	size_t GetDiscrete( size_t iDatum ) const {

		return m_auData[ iDatum ].m_i; }

	float GetContinuous( size_t iDatum ) const {

		return m_auData[ iDatum ].m_d; }

	bool IsSet( ) const {

		return !!m_auData; }

protected:
	union UDatum {
		float	m_d;
		size_t	m_i;
	};

	UDatum*	m_auData;
};

}

#endif // EXAMPLEI_H

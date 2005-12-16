#ifndef BAYESNET_H
#define BAYESNET_H

#ifdef BAYESIAN_NETWORKS

#include "bayesneti.h"

namespace libBioUtils {

class IDataset;

class IBayesNet {
public:
	virtual bool Open( const char* ) = 0;
	virtual bool Save( const char* ) const = 0;
	virtual bool Learn( const IDataset*, size_t, bool = false ) = 0;
	virtual bool Evaluate( const IDataset*, std::vector<std::vector<float> >&,
		bool = false ) const = 0;
	virtual bool Evaluate( const IDataset*, CDat&, bool = false ) const = 0;
	virtual std::vector<std::string> GetNodes( ) const = 0;
	virtual bool IsContinuous( ) const = 0;
	virtual void Randomize( ) = 0;
	virtual void Randomize( size_t ) = 0;
	virtual void Reverse( size_t ) = 0;
};

class CBayesNetSmile : CBayesNetSmileImpl, public IBayesNet {
public:
	CBayesNetSmile( bool = true );

	bool Open( const char* );
	bool Save( const char* ) const;
	bool Learn( const IDataset*, size_t, bool = false );
	bool Evaluate( const IDataset*, std::vector<std::vector<float> >&, bool = false ) const;
	bool Evaluate( const IDataset*, CDat&, bool = false ) const;
	std::vector<std::string> GetNodes( ) const;
	bool IsContinuous( ) const;
	void Randomize( );
	void Randomize( size_t );
	void Reverse( size_t );

	bool Convert( CBayesNetPNL& ) const;
};

class CBayesNetPNL : public CBayesNetPNLImpl, public IBayesNet {
public:
	CBayesNetPNL( bool = true );

	bool Open( const char* );
	bool Save( const char* ) const;
	bool Learn( const IDataset*, size_t, bool = false );
	bool Evaluate( const IDataset*, std::vector<std::vector<float> >&, bool = false ) const;
	bool Evaluate( const IDataset*, CDat&, bool = false ) const;
	std::vector<std::string> GetNodes( ) const;
	bool IsContinuous( ) const;
	void Randomize( );
	void Randomize( size_t );
	void Reverse( size_t );

	bool Convert( CBayesNetSmile& ) const;
};

}

#endif // BAYESIAN_NETWORKS

#endif // BAYESNET_H

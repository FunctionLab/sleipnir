#ifndef BAYESNET_H
#define BAYESNET_H

#include "bayesneti.h"

namespace libBioUtils {

class IDataset;

class IBayesNet {
public:
	virtual bool Open( const char* ) = 0;
	virtual bool Save( const char* ) const = 0;
	virtual bool Learn( const IDataset*, size_t, bool = false, bool = false ) = 0;
	virtual bool Evaluate( const IDataset*, std::vector<std::vector<float> >&,
		bool = false ) const = 0;
	virtual bool Evaluate( const IDataset*, CDat&, bool = false ) const = 0;
	virtual bool Evaluate( const std::vector<unsigned char>&, std::vector<float>&, bool = false, size_t = 0 ) const = 0;
	virtual bool Evaluate( const CPCLPair&, CPCL&, bool = false, int = -1 ) const = 0;
	virtual void GetNodes( std::vector<std::string>& ) const = 0;
	virtual unsigned char GetValues( size_t ) const = 0;
	virtual bool IsContinuous( ) const = 0;
	virtual bool IsContinuous( size_t ) const = 0;
	virtual void Randomize( ) = 0;
	virtual void Randomize( size_t ) = 0;
	virtual void Reverse( size_t ) = 0;
	virtual bool GetCPT( size_t, CDataMatrix& ) const = 0;
};

class CBayesNetSmile : public CBayesNetSmileImpl, public IBayesNet {
public:
	CBayesNetSmile( bool = true );

	bool Open( const std::vector<std::string>&, size_t );
	bool Convert( CBayesNetPNL& ) const;
	void SetDefault( const CBayesNetSmile& );
	bool Open( const IDataset*, const std::vector<std::string>&, const std::vector<size_t>& );
	bool Open( const CBayesNetSmile&, const std::vector<CBayesNetSmile*>& );

	bool Open( const char* );
	bool Save( const char* ) const;
	bool Learn( const IDataset*, size_t, bool = false, bool = false );
	bool Evaluate( const IDataset*, std::vector<std::vector<float> >&, bool = false ) const;
	bool Evaluate( const IDataset*, CDat&, bool = false ) const;
	bool Evaluate( const std::vector<unsigned char>&, std::vector<float>&, bool = false, size_t = 0 ) const;
	bool Evaluate( const CPCLPair&, CPCL&, bool = false, int = DSL_ALG_BN_LAURITZEN ) const;
	void GetNodes( std::vector<std::string>& ) const;
	unsigned char GetValues( size_t ) const;
	bool IsContinuous( ) const;
	bool IsContinuous( size_t ) const;
	void Randomize( );
	void Randomize( size_t );
	void Reverse( size_t );
	bool GetCPT( size_t, CDataMatrix& ) const;
};

class CBayesNetPNL : public CBayesNetPNLImpl, public IBayesNet {
public:
	CBayesNetPNL( bool = true );

	bool Convert( CBayesNetSmile& ) const;

	bool Open( const char* );
	bool Save( const char* ) const;
	bool Learn( const IDataset*, size_t, bool = false, bool = false );
	bool Evaluate( const IDataset*, std::vector<std::vector<float> >&, bool = false ) const;
	bool Evaluate( const IDataset*, CDat&, bool = false ) const;
	bool Evaluate( const std::vector<unsigned char>&, std::vector<float>&, bool = false, size_t = 0 ) const;
	bool Evaluate( const CPCLPair&, CPCL&, bool = false, int = -1 ) const;
	void GetNodes( std::vector<std::string>& ) const;
	unsigned char GetValues( size_t ) const;
	bool IsContinuous( ) const;
	bool IsContinuous( size_t ) const;
	void Randomize( );
	void Randomize( size_t );
	void Reverse( size_t );
	bool GetCPT( size_t, CDataMatrix& ) const;
};

// MEFIT OFF

class CBayesNetFN : CBayesNetFNImpl, public IBayesNet {
public:
	bool Open( const char* );
	bool Save( const char* ) const;
	bool Learn( const IDataset*, size_t, bool = false, bool = false );
	bool Evaluate( const IDataset*, std::vector<std::vector<float> >&, bool = false ) const;
	bool Evaluate( const IDataset*, CDat&, bool = false ) const;
	bool Evaluate( const std::vector<unsigned char>&, std::vector<float>&, bool = false, size_t = 0 ) const;
	bool Evaluate( const CPCLPair&, CPCL&, bool = false, int = -1 ) const;
	void GetNodes( std::vector<std::string>& ) const;
	unsigned char GetValues( size_t ) const;
	bool IsContinuous( ) const;
	bool IsContinuous( size_t ) const;
	void Randomize( );
	void Randomize( size_t );
	void Reverse( size_t );
	bool GetCPT( size_t, CDataMatrix& ) const;
};

// MEFIT ON

}

#endif // BAYESNET_H

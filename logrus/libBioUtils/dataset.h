#ifndef DATASET_H
#define DATASET_H

#include "dataseti.h"

namespace libBioUtils {

class IBayesNet;

class IDataset {
public:
	virtual bool IsHidden( size_t ) const = 0;
	virtual size_t GetDiscrete( size_t, size_t, size_t ) const = 0;
	virtual float GetContinuous( size_t, size_t, size_t ) const = 0;
	virtual const std::string& GetGene( size_t ) const = 0;
	virtual size_t GetGenes( ) const = 0;
	virtual bool IsExample( size_t, size_t ) const = 0;
	virtual const std::vector<std::string>& GetGeneNames( ) const = 0;
	virtual size_t GetExperiments( ) const = 0;
	virtual size_t GetGene( const std::string& ) const = 0;
};

class CDataset : CDatasetImpl, public IDataset {
public:
	bool Open( const char*, const char*, const IBayesNet* );
	bool Open( const CDataPair&, const char*, const IBayesNet* );
	bool Open( const char*, const std::vector<std::string>& );
	bool Open( const std::vector<std::string>& );
//	bool Open( const char*, const IBayesNet* );

	bool IsHidden( size_t ) const;
	size_t GetDiscrete( size_t, size_t, size_t ) const;
	float GetContinuous( size_t, size_t, size_t ) const;
	const std::string& GetGene( size_t ) const;
	size_t GetGenes( ) const;
	bool IsExample( size_t, size_t ) const;
	const std::vector<std::string>& GetGeneNames( ) const;
	size_t GetExperiments( ) const;
	size_t GetGene( const std::string& ) const;
};

class CDatasetCompact : CDatasetCompactImpl, public IDataset {
public:
	bool Open( const CDataPair&, const char*, const IBayesNet* );
	bool Open( const char*, const IBayesNet* );
	bool Open( const char*, const IBayesNet*, const CGenes&, const CGenes& );
	bool Open( const std::vector<std::string>& );

	bool IsHidden( size_t ) const;
	size_t GetDiscrete( size_t, size_t, size_t ) const;
	float GetContinuous( size_t, size_t, size_t ) const;
	const std::string& GetGene( size_t ) const;
	size_t GetGenes( ) const;
	bool IsExample( size_t, size_t ) const;
	const std::vector<std::string>& GetGeneNames( ) const;
	size_t GetExperiments( ) const;
	size_t GetGene( const std::string& ) const;
};

class CDataMask : CDataMaskImpl, public IDataset {
public:
	void AttachRandom( const IDataset*, float );
	void AttachComplement( const CDataMask& );
	void AttachGenes( const IDataset*, std::istream& );

	bool IsHidden( size_t ) const;
	size_t GetDiscrete( size_t, size_t, size_t ) const;
	float GetContinuous( size_t, size_t, size_t ) const;
	const std::string& GetGene( size_t ) const;
	size_t GetGenes( ) const;
	bool IsExample( size_t, size_t ) const;
	const std::vector<std::string>& GetGeneNames( ) const;
	size_t GetExperiments( ) const;
	size_t GetGene( const std::string& ) const;
};

class CDataSubset : CDataSubsetImpl, public IDataset {
public:
	bool Initialize( const char*, const IBayesNet*, size_t );
	bool Initialize( const std::vector<std::string>&, size_t );
	bool Open( size_t );

	bool IsHidden( size_t ) const;
	size_t GetDiscrete( size_t, size_t, size_t ) const;
	float GetContinuous( size_t, size_t, size_t ) const;
	const std::string& GetGene( size_t ) const;
	size_t GetGenes( ) const;
	bool IsExample( size_t, size_t ) const;
	const std::vector<std::string>& GetGeneNames( ) const;
	size_t GetExperiments( ) const;
	size_t GetGene( const std::string& ) const;
};

}

#endif // DATASET_H

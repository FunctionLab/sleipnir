#ifndef SVM_H
#define SVM_H

#include "svmi.h"

namespace libBioUtils {

class CSVM : CSVMImpl {
public:
	enum EKernel {
		EKernelLinear		= 0,
		EKernelPolynomial	= EKernelLinear + 1,
		EKernelRBF			= EKernelPolynomial + 1
	};

	bool OpenAlphas( std::istream& );
	bool Open( std::istream& );
	bool Save( std::ostream& ) const;
	bool Learn( const CPCLSet&, const CDataPair& );
	bool Learn( const IDataset*, const CDataPair& );
	bool Learn( const char* );
	bool Evaluate( const char*, CDat& ) const;
	bool Evaluate( const CPCLSet&, CDat& ) const;
	bool Evaluate( const IDataset*, CDat& ) const;
	bool Evaluate( const CPCLSet&, const CGenes&, CDat& ) const;
	bool Evaluate( const IDataset*, const CGenes&, CDat& ) const;
	void SetIterations( size_t );
	void SetCache( size_t );
	void SetKernel( EKernel );
};

}

#endif // SVM_H

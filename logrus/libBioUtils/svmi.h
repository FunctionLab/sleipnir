#ifndef SVMI_H
#define SVMI_H

extern "C" {
#include "svm_common.h"
}

namespace libBioUtils {

class CDat;
class CDataPair;
class CGenes;
class CPCLSet;

class CSVMImpl {
protected:
	struct SLearn : LEARN_PARM {
		SLearn( );
	};

	struct SKernel : KERNEL_PARM {
		SKernel( );
	};

	static DOC* CreateDoc( const CPCLSet&, size_t, size_t, size_t );
	static size_t GetWords( const CPCLSet& );

	CSVMImpl( );
	~CSVMImpl( );

	void Reset( bool, bool, bool );
	bool Initialize( const CPCLSet&, const CDataPair& );
	bool Evaluate( const CPCLSet&, const CGenes*, CDat& ) const;

	MODEL*	m_pModel;
	DOC**	m_apDocs;
	size_t	m_iDocs;
	double*	m_adLabels;
	double*	m_adAlphas;
	size_t	m_iAlphas;
	SLearn	m_sLearn;
	SKernel	m_sKernel;
};

}

#endif // SVMI_H

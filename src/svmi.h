#ifndef SVMI_H
#define SVMI_H

extern "C" {
#include "svm_common.h"
}

#include "typesi.h"

namespace Sleipnir {

class CDat;
class CDataPair;
class CGenes;
class CPCL;
class CPCLSet;
class IDataset;

class CSVMImpl {
protected:
	static const size_t	c_iWords	= 512;

	struct SLearn : LEARN_PARM {
		SLearn( );
	};

	struct SKernel : KERNEL_PARM {
		SKernel( );
	};

	struct SData {
		enum {
			EPCL,
			EPCLs,
			EData,
			EFile
		}	m_eType;
		union {
			const CPCLSet*	m_pPCLs;
			const IDataset*	m_pData;
			const char*		m_szFile;
			const CPCL*		m_pPCL;
		}	m_uData;
		union {
			const CDataPair*	m_pAnswers;
			const CGenes*		m_pGenes;
		}	m_uAnswers;
		const CGenes*	m_pNegative;
	};

	static SWORD	s_asWords[ c_iWords ];

	CSVMImpl( );
	~CSVMImpl( );

	void Reset( bool, bool, bool );
	bool Initialize( const SData& );
	bool Evaluate( const SData&, const CGenes*, CDat& ) const;
	bool EvaluateFile( const char*, CDat& ) const;
	bool Learn( const SData& );
	size_t GetWords( const SData& ) const;
	DOC* CreateDoc( const SData&, size_t, size_t, size_t ) const;
	DOC* CreateDoc( const SData&, size_t ) const;

	MODEL*		m_pModel;
	DOC**		m_apDocs;
	uint32_t	m_iDocs;
	uint32_t	m_iWords;
	double*		m_adLabels;
	double*		m_adAlphas;
	size_t		m_iAlphas;
	SLearn		m_sLearn;
	SKernel		m_sKernel;
};

}

#endif // SVMI_H

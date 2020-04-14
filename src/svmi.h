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
#ifndef SVMI_H
#define SVMI_H

namespace SVMPerf {

    extern "C" {
#undef int64_t

#include <svm_light/svm_common.h>

    }

}
using namespace SVMPerf;

namespace Sleipnir {

    class CDat;

    class CDataPair;

    class CGenes;

    class CPCL;

    class CPCLSet;

    class IDataset;

    class CSVMImpl {
    protected:
        static const size_t c_iWords = 512;

        struct SLearn : LEARN_PARM {
            SLearn();
        };

        struct SKernel : KERNEL_PARM {
            SKernel();
        };

        struct SData {
            enum {
                EPCL,
                EPCLs,
                EData,
                EFile
            } m_eType;
            union {
                const CPCLSet *m_pPCLs;
                const IDataset *m_pData;
                const char *m_szFile;
                const CPCL *m_pPCL;
            } m_uData;
            union {
                const CDataPair *m_pAnswers;
                const CGenes *m_pGenes;
            } m_uAnswers;
            const CGenes *m_pNegative;
        };

        static SVMPerf::WORD s_asWords[c_iWords];

        CSVMImpl();

        ~CSVMImpl();

        void Reset(bool, bool, bool);

        bool Initialize(const SData &);

        bool Evaluate(const SData &, const CGenes *, CDat &) const;

        bool EvaluateFile(const char *, CDat &) const;

        bool Learn(const SData &);

        size_t GetWords(const SData &) const;

        DOC *CreateDoc(const SData &, size_t, size_t, size_t) const;

        DOC *CreateDoc(const SData &, size_t) const;

        MODEL *m_pModel;
        DOC **m_apDocs;
        uint32_t m_iDocs;
        uint32_t m_iWords;
        double *m_adLabels;
        double *m_adAlphas;
        size_t m_iAlphas;
        SLearn m_sLearn;
        SKernel m_sKernel;
    };

}

#endif // SVMI_H

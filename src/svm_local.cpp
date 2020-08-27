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
#include "stdafx.h"
#include "svm_local.h"
#include "pclset.h"
#include "dataset.h"
#include "meta.h"
#include "genome.h"

#ifndef NO_SVM_PERF

extern "C" {
KERNEL_CACHE *kernel_cache_init(long, long);
void kernel_cache_cleanup(KERNEL_CACHE *);
void svm_learn_classification(DOC **, double *, long, long, LEARN_PARM *, KERNEL_PARM *,
                              KERNEL_CACHE *, MODEL *, double *);
void svm_learn_regression(DOC **, double *, long, long, LEARN_PARM *, KERNEL_PARM *,
                          KERNEL_CACHE **, MODEL *);
void svm_learn_ranking(DOC **, double *, long, long, LEARN_PARM *, KERNEL_PARM *,
                       KERNEL_CACHE **, MODEL *);
void svm_learn_optimization(DOC **, double *, long, long, LEARN_PARM *, KERNEL_PARM *,
                            KERNEL_CACHE *, MODEL *, double *);
}

namespace Sleipnir {

    bool read_documents_bin(char *szFile, DOC ***papDocs, double **padLabels,
                            uint32_t *piWords, uint32_t *piDocs) {
        char szComment[1024];
        FILE *pfileDoc;
        SVMPerf::WORD *aWords;
        uint32_t i, iDoc, iWord;
        float d;
        float *ad;

        g_CatSleipnir().info("CSVM::read_documents_bin( ) Reading binary examples into memory...");

#pragma warning( disable : 4996 )
        if (!(pfileDoc = fopen(szFile, "rb"))) {
#pragma warning( default : 4996 )
            g_CatSleipnir().error("CSVM::read_documents_bin( ) Could not open: %s", szFile);
            return false;
        }

        fread(piWords, sizeof(*piWords), 1, pfileDoc);
        fread(piDocs, sizeof(*piDocs), 1, pfileDoc);

        (*papDocs) = (DOC **) my_malloc(sizeof(DOC * ) * (*piDocs));
        (*padLabels) = (double *) my_malloc(sizeof(double) * (*piDocs));
        ad = (float *) my_malloc(sizeof(*ad) * (*piWords));
        aWords = (SVMPerf::WORD *) my_malloc(sizeof(SVMPerf::WORD) * ((*piWords) + 1));
        for (iWord = 0; iWord < (*piWords); ++iWord)
            aWords[iWord].wnum = iWord + 1;
        aWords[iWord].wnum = 0;
        for (iDoc = 0; iDoc < (*piDocs); ++iDoc) {
            if (!(iDoc % 100000))
                g_CatSleipnir().info("CSVM::read_documents_bin( ) Read %d/%d", iDoc, *piDocs);
            fread(&d, sizeof(d), 1, pfileDoc);
            (*padLabels)[iDoc] = d;
            fread(ad, sizeof(*ad), (*piWords), pfileDoc);
            for (iWord = 0; iWord < (*piWords); ++iWord)
                aWords[iWord].weight = ad[iWord];
            fread(&i, sizeof(i), 1, pfileDoc);
            if (i)
                fread(szComment, sizeof(*szComment), i, pfileDoc);
            szComment[i] = 0;
            (*papDocs)[iDoc] = create_example(iDoc, 0, 0, 1, create_svector(aWords, szComment, 1));
        }
        free(aWords);
        free(ad);

        fclose(pfileDoc);
        return true;
    }

    SVMPerf::WORD    CSVMImpl::s_asWords[CSVMImpl::c_iWords];

    CSVMImpl::SLearn::SLearn() {

        predfile[0] = 0;
        alphafile[0] = 0;
        biased_hyperplane = 1;
        sharedslack = 0;
        remove_inconsistent = 0;
        skip_final_opt_check = 0;
        svm_maxqpsize = 10;
        svm_newvarsinqp = 0;
        svm_iter_to_shrink = -1;
        maxiter = 100000;
        kernel_cache_size = 40;
        svm_c = 0;
        eps = 0.1;
        transduction_posratio = -1.0;
        svm_costratio = 0;
        svm_costratio_unlab = 1;
        svm_unlabbound = 1e-5;
        epsilon_crit = 0.001;
        epsilon_a = 1e-15;
        compute_loo = 0;
        rho = 1;
        xa_depth = 0;
        type = CLASSIFICATION;
    }

    CSVMImpl::SKernel::SKernel() {

        kernel_type = 0;
        poly_degree = 3;
        rbf_gamma = 1;
        coef_lin = 1;
        coef_const = 1;
        custom[0] = 0;
    }

    CSVMImpl::CSVMImpl() : m_apDocs(NULL), m_iDocs(0), m_adAlphas(NULL), m_iAlphas(0),
                           m_pModel(NULL), m_adLabels(NULL) {

        verbosity = 2;
    }

    CSVMImpl::~CSVMImpl() {

        Reset(true, true, true);
    }

    void CSVMImpl::Reset(bool fData, bool fModel, bool fAlphas) {
        size_t i;

        if (fModel && m_pModel) {
            free_model(m_pModel, 0);
            m_pModel = NULL;
        }
        if (fAlphas && m_adAlphas) {
            free(m_adAlphas);
            m_adAlphas = NULL;
        }
        if (fData) {
            if (m_apDocs) {
                for (i = 0; i < m_iDocs; ++i)
                    free_example(m_apDocs[i], 1);
                delete[] m_apDocs;
                m_apDocs = NULL;
            }
            if (m_adLabels) {
                delete[] m_adLabels;
                m_adLabels = NULL;
            }
        }
    }

    size_t CSVMImpl::GetWords(const SData &sData) const {
        size_t i, iRet;

        switch (sData.m_eType) {
            case SData::EPCLs:
                for (iRet = i = 0; i < sData.m_uData.m_pPCLs->GetPCLs(); ++i)
                    iRet += sData.m_uData.m_pPCLs->Get(i).GetExperiments();
                return (iRet * 2);

            case SData::EData:
                return sData.m_uData.m_pData->GetExperiments();

            case SData::EFile:
                return m_iWords;

            case SData::EPCL:
                return sData.m_uData.m_pPCL->GetExperiments();
        }

        return -1;
    }

    DOC *CSVMImpl::CreateDoc(const SData &sData, size_t iOne, size_t iTwo, size_t iDoc) const {
        SVMPerf::WORD *asWords;
        size_t i, j, iWord, iWords;
        float d;
        DOC *pRet;

        iWords = GetWords(sData);
        asWords = (iWords >= c_iWords) ? new SVMPerf::WORD[iWords + 1] : s_asWords;
        for (i = 0; i < iWords; ++i)
            asWords[i].wnum = i + 1;
        asWords[i].wnum = 0;
        if (sData.m_eType == SData::EPCLs) {
            const CPCLSet &PCLs = *sData.m_uData.m_pPCLs;

            for (iWord = i = 0; i < PCLs.GetPCLs(); ++i) {
                for (j = 0; j < PCLs.Get(i).GetExperiments(); ++j) {
                    if (CMeta::IsNaN(d = PCLs.Get(i, iOne, j)))
                        d = 0;
                    assert((iWord + j) < iWords);
                    asWords[iWord + j].weight = d;
                    if (CMeta::IsNaN(d = PCLs.Get(i, iTwo, j)))
                        d = 0;
                    assert((iWord + (iWords / 2) + j) < iWords);
                    asWords[iWord + (iWords / 2) + j].weight = d;
                }
                iWord += PCLs.Get(i).GetExperiments();
            }
        } else {
            const IDataset *pData = sData.m_uData.m_pData;

            for (i = 0; i < pData->GetExperiments(); ++i) {
                if (CMeta::IsNaN(d = pData->GetContinuous(iOne, iTwo, i)))
                    d = 0;
                asWords[i].weight = d;
            }
        }

        char empty[] = "";

        pRet = create_example(iDoc, 0, 0, 1, create_svector(asWords, empty, 1));
        if (asWords != s_asWords)
            delete[] asWords;
        return pRet;
    }

    DOC *CSVMImpl::CreateDoc(const SData &sData, size_t iGene) const {
        SVMPerf::WORD *asWords;
        size_t i, iWords;
        DOC *pRet;
        float d;
        char empty[] = "";

        iWords = GetWords(sData);
        asWords = (iWords >= c_iWords) ? new SVMPerf::WORD[iWords + 1] : s_asWords;
        for (i = 0; i < iWords; ++i)
            asWords[i].wnum = i + 1;
        asWords[i].wnum = 0;

        for (i = 0; i < iWords; ++i)
            asWords[i].weight = CMeta::IsNaN(d = sData.m_uData.m_pPCL->Get(iGene, i)) ? 0 : d;
        pRet = create_example(i, 0, 0, 1, create_svector(asWords, empty, 1));

        if (asWords != s_asWords)
            delete[] asWords;
        return pRet;
    }

/*!
 * \brief
 * Open an initial file of alphas to be used as a starting point during learning.
 * 
 * \param istm
 * Stream containing a file of alphas, one number per line, to be used to initialize learning.
 * 
 * \returns
 * True if the file was opened successfully.
 * 
 * \remarks
 * Equivalent to the svm_learn -y parameter.
 */
    bool CSVM::OpenAlphas(std::istream &istm) {
        static const size_t c_iBuf = 1024;
        char szBuf[c_iBuf];
        vector<float> vecdAlphas;
        size_t i;

        Reset(false, false, true);
        while (istm.peek() != EOF) {
            istm.getline(szBuf, c_iBuf - 1);
            vecdAlphas.push_back((float) atof(szBuf));
        }
        m_adAlphas = new double[m_iAlphas = vecdAlphas.size()];
        for (i = 0; i < m_iAlphas; ++i)
            m_adAlphas[i] = vecdAlphas[i];

        return true;
    }

    bool CSVMImpl::Initialize(const SData &sData) {
        size_t i, j, iOne, iTwo, iDoc;
        vector <size_t> veciGenes;
        float d;

        Reset(true, false, false);
        if (sData.m_eType == SData::EFile) {
            read_documents_bin((char *) sData.m_uData.m_szFile, &m_apDocs, &m_adLabels,
                               &m_iWords, &m_iDocs);
            return true;
        }
        if (sData.m_eType == SData::EPCL) {
            for (m_iDocs = i = 0; i < sData.m_uData.m_pPCL->GetGenes(); ++i)
                if (!sData.m_uData.m_pPCL->IsMasked(i) && (!sData.m_pNegative ||
                                                           sData.m_uAnswers.m_pGenes->IsGene(
                                                                   sData.m_uData.m_pPCL->GetGene(i)) ||
                                                           sData.m_pNegative->IsGene(sData.m_uData.m_pPCL->GetGene(i))))
                    m_iDocs++;
            m_apDocs = new DOC *[m_iDocs];
            m_adLabels = new double[m_iDocs];
            for (i = j = 0; i < sData.m_uData.m_pPCL->GetGenes(); ++i)
                if (!sData.m_uData.m_pPCL->IsMasked(i)) {
                    const string &strGene = sData.m_uData.m_pPCL->GetGene(i);

                    d = 0;
                    if (!sData.m_pNegative)
                        d = sData.m_uAnswers.m_pGenes->IsGene(strGene) ? 1.0f : -1.0f;
                    else if (sData.m_uAnswers.m_pGenes->IsGene(strGene))
                        d = 1;
                    else if (sData.m_pNegative->IsGene(strGene))
                        d = -1;
                    if (d) {
                        m_apDocs[j] = CreateDoc(sData, i);
                        m_adLabels[j++] = d;
                    }
                }
            return true;
        }

        veciGenes.resize((sData.m_eType == SData::EPCLs) ?
                         sData.m_uData.m_pPCLs->GetGenes() : sData.m_uData.m_pData->GetGenes());
        for (i = 0; i < veciGenes.size(); ++i)
            veciGenes[i] = sData.m_uAnswers.m_pAnswers->GetGene((sData.m_eType == SData::EPCLs) ?
                                                                sData.m_uData.m_pPCLs->GetGene(i)
                                                                                                : sData.m_uData.m_pData->GetGene(
                            i));
        for (m_iDocs = i = 0; i < veciGenes.size(); ++i)
            if ((iOne = veciGenes[i]) != -1)
                for (j = (i + 1); j < veciGenes.size(); ++j)
                    if (((iTwo = veciGenes[j]) != -1) &&
                        !CMeta::IsNaN(sData.m_uAnswers.m_pAnswers->Get(iOne, iTwo)))
                        m_iDocs++;
        m_apDocs = new DOC *[m_iDocs];
        m_adLabels = new double[m_iDocs];

        for (iDoc = i = 0; i < veciGenes.size(); ++i)
            if ((iOne = veciGenes[i]) != -1)
                for (j = (i + 1); j < veciGenes.size(); ++j)
                    if (((iTwo = veciGenes[j]) != -1) &&
                        !CMeta::IsNaN(d = sData.m_uAnswers.m_pAnswers->Get(iOne, iTwo))) {
                        m_adLabels[iDoc] = d ? 1 : -1;
                        iDoc++;
                        m_apDocs[iDoc] = CreateDoc(sData, i, j, iDoc);
                    }
        assert(iDoc == m_iDocs);

        return true;
    }

/*!
 * \brief
 * Learn an SVM recognizing the given genes' values in the given PCL.
 * 
 * \param PCL
 * PCL from which features are read.
 * 
 * \param GenesPositive
 * Set of positive examples to be learned.
 * 
 * \returns
 * True if the model was learned successfully.
 * 
 * Learns an SVM model using the current settings, assuming each column of the given PCL is a feature and
 * each row a record.  Genes in the given positive set become positive examples, and all other genes
 * are assumed to be negative.
 */
    bool CSVM::Learn(const CPCL &PCL, const CGenes &GenesPositive) {
        CGenes GenesNeg(GenesPositive.GetGenome());

        return Learn(PCL, GenesPositive, GenesNeg);
    }

/*!
 * \brief
 * Learn an SVM recognizing the given genes' values in the given PCL.
 * 
 * \param PCL
 * PCL from which features are read.
 * 
 * \param GenesPositive
 * Set of positive examples to be learned.
 * 
 * \param GenesNegative
 * Set of negative examples to be learned.
 * 
 * \returns
 * True if the model was learned successfully.
 * 
 * Learns an SVM model using the current settings, assuming each column of the given PCL is a feature and
 * each row a record.  Genes in the given positive set become positive examples, genes in the given negative
 * set become negative examples, and all other rows are unlabeled.
 */
    bool CSVM::Learn(const CPCL &PCL, const CGenes &GenesPositive, const CGenes &GenesNegative) {
        SData sData;

        sData.m_eType = SData::EPCL;
        sData.m_uData.m_pPCL = &PCL;
        sData.m_uAnswers.m_pGenes = &GenesPositive;
        sData.m_pNegative = GenesNegative.GetGenes() ? &GenesNegative : NULL;

        return CSVMImpl::Learn(sData);
    }

    bool CSVMImpl::Learn(const SData &sData) {
        KERNEL_CACHE *pCache;
        size_t i, iNeg, iPos, iWords;

        Reset(false, true, false);
        m_pModel = (MODEL *) calloc(1, sizeof(*m_pModel));
        if (!Initialize(sData))
            return false;
        if (!m_sLearn.svm_costratio) {
            for (iNeg = iPos = i = 0; i < m_iDocs; ++i)
                if (m_adLabels[i] == 1)
                    iPos++;
                else
                    iNeg++;
            m_sLearn.svm_costratio = (float) iNeg / iPos;
        }
        if (m_sLearn.svm_iter_to_shrink < 0)
            m_sLearn.svm_iter_to_shrink = (m_sKernel.kernel_type == LINEAR) ? 2 : 100;
        iWords = GetWords(sData);

        pCache = (m_sKernel.kernel_type == LINEAR) ? NULL :
                 kernel_cache_init(m_iDocs, m_sLearn.kernel_cache_size);
        switch (m_sLearn.type) {
            case CLASSIFICATION:
                svm_learn_classification(m_apDocs, m_adLabels, m_iDocs, iWords,
                                         (LEARN_PARM * ) & m_sLearn, (KERNEL_PARM * ) & m_sKernel, pCache, m_pModel,
                                         m_adAlphas);
                break;

            case REGRESSION:
                svm_learn_regression(m_apDocs, m_adLabels, m_iDocs, iWords,
                                     (LEARN_PARM * ) & m_sLearn, (KERNEL_PARM * ) & m_sKernel, &pCache, m_pModel);
                break;

            case RANKING:
                svm_learn_ranking(m_apDocs, m_adLabels, m_iDocs, iWords, (LEARN_PARM * ) & m_sLearn,
                                  (KERNEL_PARM * ) & m_sKernel, &pCache, m_pModel);
                break;

            case OPTIMIZATION:
                svm_learn_optimization(m_apDocs, m_adLabels, m_iDocs, iWords,
                                       (LEARN_PARM * ) & m_sLearn, (KERNEL_PARM * ) & m_sKernel, pCache, m_pModel,
                                       m_adAlphas);
                break;
        }

        if (pCache)
            kernel_cache_cleanup(pCache);

        return true;
    }

/*!
 * \brief
 * Save an SVM model file to the given stream.
 * 
 * \param ostm
 * Stream to which model file is saved.
 * 
 * \returns
 * True if the model file was saved successfully.
 * 
 * \remarks
 * Equivalent to the model file generated by svm_learn or input to svm_classify.
 * 
 * \see
 * Open
 */
    bool CSVM::Save(std::ostream &ostm) const {
        size_t i, j;
        SVECTOR *pVec;

        if (!m_pModel)
            return false;

        ostm << "SVM-light Version " << VERSION << endl;
        ostm << m_pModel->kernel_parm.kernel_type << " # kernel type" << endl;
        ostm << m_pModel->kernel_parm.poly_degree << " # kernel parameter -d" << endl;
        ostm << m_pModel->kernel_parm.rbf_gamma << " # kernel parameter -g" << endl;
        ostm << m_pModel->kernel_parm.coef_lin << " # kernel parameter -s" << endl;
        ostm << m_pModel->kernel_parm.coef_const << " # kernel parameter -r" << endl;
        ostm << m_pModel->kernel_parm.custom << "# kernel parameter -u" << endl;
        ostm << m_pModel->totwords << " # highest feature index" << endl;
        ostm << m_pModel->totdoc << " # number of training documents" << endl;

        for (i = j = 1; i < (size_t) m_pModel->sv_num; ++i)
            for (pVec = m_pModel->supvec[i]->fvec; pVec; pVec = pVec->next)
                j++;
        ostm << j << " # number of support vectors plus 1" << endl;
        ostm << m_pModel->b <<
             " # threshold b, each following line is a SV (starting with alpha*y)" << endl;

        for (i = 1; i < (size_t) m_pModel->sv_num; ++i)
            for (pVec = m_pModel->supvec[i]->fvec; pVec; pVec = pVec->next) {
                ostm << (m_pModel->alpha[i] * pVec->factor) << ' ';
                for (j = 0; pVec->words[j].wnum; ++j)
                    ostm << pVec->words[j].wnum << ':' << pVec->words[j].weight << ' ';
                ostm << '#' << endl;
            }
//			ostm << '#' << pVec->userdefined << endl; }

        return true;
    }

    bool CSVMImpl::Evaluate(const SData &sData, const CGenes *pGenesIn, CDat &DatOut) const {
        size_t i, j, iGenes;
        DOC *pDoc;

        if (!m_pModel)
            return false;
        if (m_pModel->kernel_parm.kernel_type == 0)
            add_weight_vector_to_linear_model(m_pModel);

        if (sData.m_eType == SData::EFile)
            return EvaluateFile(sData.m_uData.m_szFile, DatOut);

        iGenes = (sData.m_eType == SData::EPCLs) ? sData.m_uData.m_pPCLs->GetGenes() :
                 sData.m_uData.m_pData->GetGenes();
        for (i = 0; i < iGenes; ++i) {
            const string &strGeneOne = (sData.m_eType == SData::EPCLs) ?
                                       sData.m_uData.m_pPCLs->GetGene(i) : sData.m_uData.m_pData->GetGene(i);

            if (!(i % 10))
                g_CatSleipnir().notice("CSVMImpl::Evaluate( ) gene %d/%d", i, iGenes);
            if (pGenesIn && !pGenesIn->IsGene(strGeneOne))
                continue;
            for (j = (i + 1); j < iGenes; ++j) {
                const string &strGeneTwo = (sData.m_eType == SData::EPCLs) ?
                                           sData.m_uData.m_pPCLs->GetGene(j) : sData.m_uData.m_pData->GetGene(j);

                if (pGenesIn && !pGenesIn->IsGene(strGeneTwo))
                    continue;
                if (!(pDoc = CreateDoc(sData, i, j, 0)))
                    return false;
                DatOut.Set(i, j, (float) (m_pModel->kernel_parm.kernel_type ?
                                          classify_example(m_pModel, pDoc) :
                                          classify_example_linear(m_pModel, pDoc)));
                free_example(pDoc, 1);
            }
        }

        return true;
    }

/*!
 * \brief
 * Evaluate the SVM's output for each row of the given PCL.
 * 
 * \param PCL
 * PCL from which features are read.
 * 
 * \param vecdResults
 * Values output by the SVM for each row of the given PCL.
 * 
 * \returns
 * True if the evaluation completed successfully.
 * 
 * Evaluates the current SVM model, assuming each column of the given PCL is a feature and each row a record.
 */
    bool CSVM::Evaluate(const CPCL &PCL, vector<float> &vecdResults) const {
        size_t i;
        DOC *pDoc;
        SData sData;

        if (!m_pModel)
            return false;
        if (m_pModel->kernel_parm.kernel_type == 0)
            add_weight_vector_to_linear_model(m_pModel);

        sData.m_eType = SData::EPCL;
        sData.m_uData.m_pPCL = &PCL;
        for (i = 0; i < PCL.GetGenes(); ++i) {
            if (!(i % 1000))
                g_CatSleipnir().notice("CSVMImpl::Evaluate( ) gene %d/%d", i, PCL.GetGenes());
            if (PCL.IsMasked(i))
                continue;

            if (!(pDoc = CreateDoc(sData, i)))
                return false;
            vecdResults.push_back((float) (m_pModel->kernel_parm.kernel_type ?
                                           classify_example(m_pModel, pDoc) :
                                           classify_example_linear(m_pModel, pDoc)));
            free_example(pDoc, 1);
        }

        return true;
    }

    bool CSVMImpl::EvaluateFile(const char *szFile, CDat &DatOut) const {
        static const size_t c_iSize = 512;
        char szGene[c_iSize];
        SVMPerf::WORD *asWords;
        char *pc;
        ifstream ifsm;
        vector <string> vecstrGenes;
        uint32_t i, j, k, iDocs, iWords, iGenes;
        float *ad;
        DOC *pDoc;
        char empty[] = "";

        ifsm.open(szFile, ios_base::binary);
        if (!ifsm.is_open())
            return false;
        ifsm.read((char *) &iWords, sizeof(iWords));
        ifsm.read((char *) &iDocs, sizeof(iDocs));
        ifsm.seekg(iDocs * (((iWords + 1) * sizeof(float)) +
                            (3 * sizeof(iDocs))), ios_base::cur);
        ifsm.read((char *) &iGenes, sizeof(iGenes));
        vecstrGenes.resize(iGenes);
        for (i = 0; i < iGenes; ++i) {
            for (pc = szGene;; ++pc) {
                ifsm.read(pc, 1);
                if (!*pc)
                    break;
            }
            vecstrGenes[i] = szGene;
        }
        DatOut.Open(vecstrGenes);

        asWords = (iWords >= c_iWords) ? new SVMPerf::WORD[iWords + 1] : s_asWords;
        for (i = 0; i < iWords; ++i)
            asWords[i].wnum = i + 1;
        asWords[i].wnum = 0;

        ad = new float[iWords + 1];
        ifsm.seekg(2 * sizeof(iDocs), ios_base::beg);
        for (i = 0; i < iDocs; ++i) {
            if (!(i % 1000))
                g_CatSleipnir().notice("CSVMImpl::EvaluateFile( %s ) pair %d/%d", szFile, i,
                                       iDocs);
            ifsm.read((char *) ad, (iWords + 1) * sizeof(*ad));
            for (j = 0; j < iWords; ++j)
                asWords[j].weight = ad[j + 1];
            pDoc = create_example(i, 0, 0, 1, create_svector(asWords, empty, 1));
            ifsm.read((char *) &j, sizeof(j));
            ifsm.read((char *) &j, sizeof(j));
            ifsm.read((char *) &k, sizeof(k));
            DatOut.Set(j, k, (float) (m_pModel->kernel_parm.kernel_type ?
                                      classify_example(m_pModel, pDoc) :
                                      classify_example_linear(m_pModel, pDoc)));
            free_example(pDoc, 1);
        }
        delete[] ad;

        if (asWords != s_asWords)
            delete[] asWords;

        return true;
    }

/*!
 * \brief
 * Open an SVM model file in the given stream.
 * 
 * \param istm
 * Stream from which model file is loaded.
 * 
 * \returns
 * True if the model file was loaded successfully.
 * 
 * \remarks
 * Equivalent to the model file generated by svm_learn or input to svm_classify.
 * 
 * \see
 * Save
 */
    bool CSVM::Open(std::istream &istm) {
        static const size_t c_iBuf = 131072;
        char szBuf[c_iBuf];
        vector <string> vecstrLine, vecstrToken;
        SVMPerf::WORD *asWords;
        size_t i, j;
        char empty[] = "";

        Reset(false, true, true);
        m_pModel = (MODEL *) calloc(1, sizeof(*m_pModel));

        istm.getline(szBuf, c_iBuf - 1);
        istm >> m_pModel->kernel_parm.kernel_type;
        istm.getline(szBuf, c_iBuf - 1);
        istm >> m_pModel->kernel_parm.poly_degree;
        istm.getline(szBuf, c_iBuf - 1);
        istm >> m_pModel->kernel_parm.rbf_gamma;
        istm.getline(szBuf, c_iBuf - 1);
        istm >> m_pModel->kernel_parm.coef_lin;
        istm.getline(szBuf, c_iBuf - 1);
        istm >> m_pModel->kernel_parm.coef_const;
        istm.getline(szBuf, c_iBuf - 1);
        istm.getline(szBuf, c_iBuf - 1);
        CMeta::Tokenize(szBuf, vecstrLine, "#", true);
        if (vecstrLine.size() > 1)
            strcpy_s(m_pModel->kernel_parm.custom, 49, vecstrLine[0].c_str());
        istm >> m_pModel->totwords;
        istm.getline(szBuf, c_iBuf - 1);
        istm >> m_pModel->totdoc;
        istm.getline(szBuf, c_iBuf - 1);
        istm >> m_pModel->sv_num;
        istm.getline(szBuf, c_iBuf - 1);
        istm >> m_pModel->b;
        istm.getline(szBuf, c_iBuf - 1);

        m_pModel->supvec = (DOC **) malloc(m_pModel->sv_num * sizeof(*m_pModel->supvec));
        m_pModel->alpha = (double *) malloc(m_pModel->sv_num * sizeof(*m_pModel->alpha));
        m_pModel->index = nullptr;
        m_pModel->lin_weights = nullptr;

        asWords = new SVMPerf::WORD[m_pModel->totwords + 1];
        asWords[m_pModel->totwords].wnum = 0;
        for (i = 1; i < (size_t) m_pModel->sv_num; ++i) {
            istm.getline(szBuf, c_iBuf - 1);
            szBuf[c_iBuf - 1] = 0;
            vecstrLine.clear();
            CMeta::Tokenize(szBuf, vecstrLine, CMeta::c_szWS, true);
            if (vecstrLine.size() != (m_pModel->totwords + 2)) {
                g_CatSleipnir().error("CSVM::Open( ) wanted %d words but only found %d on line: %s",
                                      (m_pModel->totwords + 2), vecstrLine.size(), szBuf);
                delete[] asWords;
                return false;
            }
            m_pModel->alpha[i] = strtod(vecstrLine[0].c_str(), nullptr);
            for (j = 1; (j + 1) < vecstrLine.size(); ++j) {
                vecstrToken.clear();
                CMeta::Tokenize(vecstrLine[j].c_str(), vecstrToken, ":", true);
                if (vecstrToken.size() != 2) {
                    g_CatSleipnir().error("CSVM::Open( ) found illegal token \"%s\" on line: %s",
                                          vecstrLine[j].c_str(), szBuf);
                    delete[] asWords;
                    return false;
                }
                asWords[j - 1].wnum = strtol(vecstrToken[0].c_str(), nullptr, 10);
                asWords[j - 1].weight = (float) strtof(vecstrToken[1].c_str(), nullptr);
            }
            m_pModel->supvec[i] = create_example(-1, 0, 0, 0, create_svector(asWords, empty,
                                                                             1));
        }

        delete[] asWords;
        return true;
    }

}

#endif // NO_SVM_PERF

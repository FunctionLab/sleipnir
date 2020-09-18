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
#include "cmdline.h"

struct SProcessor {
    const CDat &m_Dat;
    const vector <vector<size_t>> &m_vecveciSets1;
    const vector <vector<size_t>> &m_vecveciSets2;
    long double *m_adResults;
    size_t m_iThreads;

    SProcessor(const CDat &Dat, const vector <vector<size_t>> &vecveciSets1,
               const vector <vector<size_t>> &vecveciSets2, long double *adResults, size_t iThreads) : m_Dat(Dat),
                                                                                                       m_vecveciSets1(
                                                                                                               vecveciSets1),
                                                                                                       m_vecveciSets2(
                                                                                                               vecveciSets2),
                                                                                                       m_adResults(
                                                                                                               adResults),
                                                                                                       m_iThreads(
                                                                                                               iThreads) {}
};

typedef int (TFnProcessor)(const SProcessor &);

static const char c_szDab[] = ".dab";

struct SWithin {
    const CDat *m_pDat;
    const vector <vector<size_t>> *m_pvecveciSets;
    long double *m_adResults;
    long double m_dPrior;
    size_t m_iStart;
    size_t m_iLength;
};

struct SBetween {
    const CDat *m_pDat;
    const vector <vector<size_t>> *m_pvecveciSets1;
    const vector <vector<size_t>> *m_pvecveciSets2;
    long double *m_adResults;
    size_t m_iStart;
    size_t m_iLength;
};

struct SBackground {
    const CDat *m_pDat;
    const vector <vector<size_t>> *m_pvecveciSets;
    long double *m_adResults;
    size_t m_iStart;
    size_t m_iLength;
};

struct SSummary {
    size_t m_iCount;
    float m_dSum;
    float m_dSumSquares;

    SSummary() {

        Clear();
    }

    void Add(float d) {

        m_iCount++;
        m_dSum += d;
        m_dSumSquares += d * d;
    }

    void Add(const SSummary &sSummary) {

        m_iCount += sSummary.m_iCount;
        m_dSum += sSummary.m_dSum;
        m_dSumSquares += sSummary.m_dSumSquares;
    }

    void Multiply(float d) {

        m_iCount = (size_t)(m_iCount * d);
        m_dSum *= d;
        m_dSumSquares *= d;
    }

    void Clear() {

        m_iCount = 0;
        m_dSum = m_dSumSquares = 0;
    }

    float GetAverage() const {

        return (m_iCount ? (m_dSum / m_iCount) : CMeta::GetNaN());
    }

    float GetStdev() const {

        return (m_iCount ? sqrt((m_dSumSquares / (max((size_t) 2, m_iCount) - 1)) - pow(GetAverage(), 2))
                         : CMeta::GetNaN());
    }
};

struct SDatum {
    SSummary m_sHubbiness;
    SSummary m_sCliquiness;
    vector <pair<size_t, float>> m_vecprSpecific;

    void Clear() {

        m_sHubbiness.Clear();
        m_sCliquiness.Clear();
        m_vecprSpecific.clear();
    }
};

struct STerm {
    string m_strFile;
    size_t m_iGenes;
    size_t m_iTotal;
    const CDat *m_pDat;
    const CPCL *m_pPCLWeights;
    string m_strClip;
    const vector <size_t> *m_pveciDat2PCL;
    const vector <SSummary> *m_pvecsHubs;
    pthread_mutex_t *m_pmutx;
};

static size_t hubs(const CDat &, vector <SSummary> &, const CPCL &, const vector <size_t> &);

static size_t cliques(const CDat &, size_t, const vector <SSummary> &, size_t, SDatum &, const CGenes *, const CPCL &,
                      const vector <size_t> &);

static void enset(const CDat &, const vector <vector<string>> &, vector <vector<size_t>> &);

static int sets(const char *, const vector <string> &, vector <vector<string>> &);

static int process(const char *, bool, bool, const vector <vector<string>> &, const vector <vector<string>> &,
                   long double *, TFnProcessor *, size_t);

static TFnProcessor background;
static TFnProcessor within;
static TFnProcessor between;

static void *between_thread(void *);

static void *within_thread(void *);

static void *background_thread(void *);

static void *term_thread(void *);

int main(int iArgs, char **aszArgs) {
#ifdef WIN32
    pthread_win32_process_attach_np( );
#endif // WIN32
    gengetopt_args_info sArgs{};
    CDat Dat;
    size_t i, j, iGenes, iTotal, iThread;
    vector <SSummary> vecsHubs;
    SDatum sDatum;
    CPCL PCLWeights;
    vector <size_t> veciDat2PCL;
    vector <STerm> vecsTerms;
    vector <pthread_t> vecpthdThreads;
    pthread_mutex_t mutxTerms;

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }
    CMeta Meta(sArgs.verbosity_arg);

    if (sArgs.output_arg) {
        CDataMatrix MatBackgrounds;
        int iRet;
        vector <string> vecstrGenes, vecstrContexts;
        vector <vector<string>> vecvecstrSets1, vecvecstrSets2;
        TFnProcessor *pfnProcessor;
        long double *adResults;

        if (!(sArgs.contexts_arg && sArgs.genelist_arg)) {
            cmdline_parser_print_help();
            return 1;
        }
        {
            static const size_t c_iContext = 2;
            CPCL PCLContexts(false);

            if (!PCLContexts.Open(sArgs.contexts_arg, c_iContext)) {
                cerr << "Could not open: " << sArgs.contexts_arg << endl;
                return 1;
            }
            vecstrContexts.resize(PCLContexts.GetGenes());
            for (i = 0; i < vecstrContexts.size(); ++i)
                vecstrContexts[i] = PCLContexts.GetFeature(i, c_iContext);
        }
        {
            CPCL PCLGenes(false);

            if (!PCLGenes.Open(sArgs.genelist_arg, 1)) {
                cerr << "Could not open: " << sArgs.genelist_arg << endl;
                return 1;
            }
            vecstrGenes.resize(PCLGenes.GetGenes());
            for (i = 0; i < vecstrGenes.size(); ++i)
                vecstrGenes[i] = PCLGenes.GetFeature(i, 1);
        }
        if (sArgs.genesets_arg) {
            if ((iRet = sets(sArgs.genesets_arg, vecstrGenes, vecvecstrSets1)))
                return iRet;
            if (sArgs.between_arg) {
                if ((iRet = sets(sArgs.between_arg, vecstrGenes, vecvecstrSets2)))
                    return iRet;
                pfnProcessor = between;
            } else
                pfnProcessor = within;
        } else {
            vecvecstrSets1.resize(vecstrGenes.size());
            for (i = 0; i < vecvecstrSets1.size(); ++i) {
                vecvecstrSets1[i].resize(1);
                vecvecstrSets1[i][0] = vecstrGenes[i];
            }
            pfnProcessor = background;
        }

        MatBackgrounds.Initialize(vecstrContexts.size() + 1, vecvecstrSets1.size() * max((size_t) 1,
                                                                                         vecvecstrSets2.size()));
        adResults = new long double[MatBackgrounds.GetColumns()];
        if ((iRet = process(sArgs.input_arg, !!sArgs.memmap_flag, !!sArgs.normalize_flag, vecvecstrSets1,
                           vecvecstrSets2, adResults, pfnProcessor, sArgs.threads_arg)))
            return iRet;
        for (i = 0; i < MatBackgrounds.GetColumns(); ++i)
            MatBackgrounds.Set(0, i, (float) adResults[i]);
        for (i = 1; i < MatBackgrounds.GetRows(); ++i) {
            if ((iRet = process(((string) sArgs.directory_arg + '/' +
                                CMeta::Filename(vecstrContexts[i - 1]) + ".dab").c_str(), !!sArgs.memmap_flag,
                               !!sArgs.normalize_flag, vecvecstrSets1, vecvecstrSets2, adResults, pfnProcessor,
                               sArgs.threads_arg)))
                return iRet;
            for (j = 0; j < MatBackgrounds.GetColumns(); ++j)
                MatBackgrounds.Set(i, j, (float) adResults[j]);
        }

        MatBackgrounds.Save(sArgs.output_arg, true);
        return 0;
    }

    if (sArgs.input_arg) {
        if (!Dat.Open(sArgs.input_arg,
                      !(sArgs.memmap_flag || sArgs.normalize_flag || sArgs.genin_arg || sArgs.genex_arg))) {
            cerr << "Could not open: " << sArgs.input_arg << endl;
            return 1;
        }
    } else if (!Dat.Open(cin, CDat::EFormatText)) {
        cerr << "Could not open input" << endl;
        return 1;
    }
    if (sArgs.genin_arg && !Dat.FilterGenes(sArgs.genin_arg, CDat::EFilterInclude)) {
        cerr << "Could not open: " << sArgs.genin_arg << endl;
        return 1;
    }
    if (sArgs.genex_arg && !Dat.FilterGenes(sArgs.genex_arg, CDat::EFilterExclude)) {
        cerr << "Could not open: " << sArgs.genex_arg << endl;
        return 1;
    }
    if (sArgs.normalize_flag)
        Dat.Normalize(CDat::ENormalizeSigmoid);

    if (sArgs.weights_arg) {
        if (!PCLWeights.Open(sArgs.weights_arg, 0)) {
            cerr << "Could not open weights: " << sArgs.weights_arg << endl;
            return 1;
        }
        veciDat2PCL.resize(Dat.GetGenes());
        for (i = 0; i < veciDat2PCL.size(); ++i)
            veciDat2PCL[i] = PCLWeights.GetGene(Dat.GetGene(i));
    }

    iTotal = hubs(Dat, vecsHubs, PCLWeights, veciDat2PCL);
    if (sArgs.genes_arg == -1) {
        cout << "Function";
        for (i = 0; i < Dat.GetGenes(); ++i)
            cout << '\t' << Dat.GetGene(i);
    } else {
        cliques(Dat, iTotal, vecsHubs, 0, sDatum, nullptr, PCLWeights, veciDat2PCL);
        cout << "name	size	hubbiness	hubbiness std.	hubbiness n	cliquiness	cliquiness std.	cliquiness n"
             << endl;
        cout << "total	" << iTotal << '\t' << sDatum.m_sHubbiness.GetAverage() << '\t' <<
             sDatum.m_sHubbiness.GetStdev() << '\t' << sDatum.m_sHubbiness.m_iCount << '\t'
             << sDatum.m_sCliquiness.GetAverage() << '\t' <<
             sDatum.m_sCliquiness.GetStdev() << '\t' << sDatum.m_sCliquiness.m_iCount;
    }
    cout << endl;

    pthread_mutex_init(&mutxTerms, nullptr);
    vecsTerms.resize(sArgs.threads_arg);
    vecpthdThreads.resize(sArgs.threads_arg);
    for (iGenes = 0; iGenes < sArgs.inputs_num; iGenes += sArgs.threads_arg) {
        for (iThread = 0;
             (iThread < (size_t) sArgs.threads_arg) && ((iGenes + iThread) < sArgs.inputs_num); ++iThread) {
            if (!((iGenes + iThread) % 25))
                cerr << (iGenes + iThread) << '/' << sArgs.inputs_num << endl;
            vecsTerms[iThread].m_strFile = sArgs.inputs[iGenes + iThread];
            vecsTerms[iThread].m_iGenes = sArgs.genes_arg;
            vecsTerms[iThread].m_pDat = &Dat;
            vecsTerms[iThread].m_iTotal = iTotal;
            vecsTerms[iThread].m_pPCLWeights = &PCLWeights;
            vecsTerms[iThread].m_strClip = sArgs.clip_arg ? sArgs.clip_arg : "";
            vecsTerms[iThread].m_pveciDat2PCL = &veciDat2PCL;
            vecsTerms[iThread].m_pvecsHubs = &vecsHubs;
            vecsTerms[iThread].m_pmutx = &mutxTerms;
            if (pthread_create(&vecpthdThreads[iThread], nullptr, term_thread, &vecsTerms[iThread])) {
                cerr << "Couldn't create thread for term: " << sArgs.input_arg[iGenes + iThread] << endl;
                return 1;
            }
        }
        for (i = 0; (i < vecpthdThreads.size()) && ((i + iGenes) < sArgs.inputs_num); ++i)
            pthread_join(vecpthdThreads[i], nullptr);
    }
    pthread_mutex_destroy(&mutxTerms);

#ifdef WIN32
    pthread_win32_process_detach_np( );
#endif // WIN32
    return 0;
}

static inline float get_weight(size_t iGene, const vector <size_t> &veciDat2PCL, const CPCL &PCLWeights) {
    size_t iPCL;

    if (iGene >= veciDat2PCL.size())
        return 1;

    iPCL = veciDat2PCL[iGene];
    return (((iPCL == -1) || (iPCL >= PCLWeights.GetGenes())) ? 1 : PCLWeights.Get(iPCL, 0));
}

size_t hubs(const CDat &Dat, vector <SSummary> &vecsHubs, const CPCL &PCLWeights, const vector <size_t> &veciDat2PCL) {
    size_t i, j, iRet;
    float d;

    vecsHubs.resize(Dat.GetGenes());
    for (i = 0; i < Dat.GetGenes(); ++i)
        for (j = (i + 1); j < Dat.GetGenes(); ++j) {
            if (CMeta::IsNaN(d = Dat.Get(i, j)))
                continue;
            d *= get_weight(i, veciDat2PCL, PCLWeights) * get_weight(j, veciDat2PCL, PCLWeights);
            vecsHubs[i].Add(d);
            vecsHubs[j].Add(d);
        }
    for (iRet = i = 0; i < vecsHubs.size(); ++i)
        if (vecsHubs[i].m_iCount)
            iRet++;

    return iRet;
}

struct SSorter {

    bool operator()(const pair<size_t, float> &prOne, const pair<size_t, float> &prTwo) const {

        return (prOne.second > prTwo.second);
    }
};

size_t cliques(const CDat &Dat, size_t iGenes, const vector <SSummary> &vecsHubs, size_t iSort, SDatum &sDatum,
               const CGenes *pGenes, const CPCL &PCLWeights, const vector <size_t> &veciDat2PCL) {
    size_t i, j, iRet;
    float d;
    vector <SSummary> vecsClique;
    vector<bool> vecfOutside;

    vecfOutside.resize(Dat.GetGenes());
    if (pGenes) {
        for (i = 0; i < vecfOutside.size(); ++i)
            if (!pGenes->IsGene(Dat.GetGene(i)))
                vecfOutside[i] = true;
    }
    vecsClique.resize(Dat.GetGenes());
    sDatum.Clear();
    for (i = 0; i < Dat.GetGenes(); ++i)
        for (j = (i + 1); j < Dat.GetGenes(); ++j) {
            if (CMeta::IsNaN(d = Dat.Get(i, j)))
                continue;
            d *= get_weight(i, veciDat2PCL, PCLWeights) * get_weight(j, veciDat2PCL, PCLWeights);
            // Added 9/4/13 to fix Arjun Issue #17, ignore precomputed hubs -- avoids double counting with summary output
            if (!vecfOutside[i] || !vecfOutside[j]) {//if either is not outside, edge is incident to context
                sDatum.m_sHubbiness.Add(d);
            }
            if (!vecfOutside[i])
                vecsClique[j].Add(d);
            if (!vecfOutside[j])
                vecsClique[i].Add(d);
        }
    for (iRet = i = 0; i < vecsClique.size(); ++i)
        if (!vecfOutside[i]) {
            if (vecsClique[i].m_iCount)
                iRet++;
            sDatum.m_sCliquiness.Add(vecsClique[i]);
        }
    sDatum.m_sCliquiness.Multiply(0.5);

    /* Removed 9/4/13 to fix Arjun Issue #17
	for( i = 0; i < Dat.GetGenes( ); ++i )
		if( !vecfOutside[ i ] )
			sDatum.m_sHubbiness.Add( vecsHubs[ i ] );
    */

    if (iSort) {
        sDatum.m_vecprSpecific.resize(Dat.GetGenes());
        for (i = 0; i < sDatum.m_vecprSpecific.size(); ++i) {
            sDatum.m_vecprSpecific[i].first = i;
            //divide each genes cliquiness (within) by precomputed hubbiness (overall average) to calculate group membership
            sDatum.m_vecprSpecific[i].second = vecsClique[i].GetAverage() / vecsHubs[i].GetAverage();
        }
        if (iSort != -1)
            sort(sDatum.m_vecprSpecific.begin(), sDatum.m_vecprSpecific.end(), SSorter());
    }

    return iRet;
}

int process(const char *szFile, bool fMemmap, bool fNormalize, const vector <vector<string>> &vecvecstrSets1,
            const vector <vector<string>> &vecvecstrSets2, long double *adResults, TFnProcessor *pfnProcessor,
            size_t iThreads) {
    CDat Dat;
    vector <vector<size_t>> vecveciSets1, vecveciSets2;

    if (!Dat.Open(szFile, fMemmap && !fNormalize)) {
        cerr << "Could not open: " << szFile << endl;
        return 1;
    }
    cerr << "Calculating for: " << szFile << endl;
    if (fNormalize)
        Dat.Normalize(CDat::ENormalizeSigmoid);

    enset(Dat, vecvecstrSets1, vecveciSets1);
    enset(Dat, vecvecstrSets2, vecveciSets2);
    return pfnProcessor(SProcessor(Dat, vecveciSets1, vecveciSets2, adResults, iThreads));
}

int background(const SProcessor &sProcessor) {
    const CDat &Dat = sProcessor.m_Dat;
    const vector <vector<size_t>> &vecveciSets = sProcessor.m_vecveciSets1;
    long double *adResults = sProcessor.m_adResults;
    vector <pthread_t> vecpthdThreads;
    vector <SBackground> vecsBackground;
    size_t i, iChunk;

    vecpthdThreads.resize(min(sProcessor.m_iThreads, vecveciSets.size()));
    vecsBackground.resize(vecpthdThreads.size());
    iChunk = 1 + ((vecveciSets.size() - 1) / vecpthdThreads.size());
    for (i = 0; i < vecpthdThreads.size(); ++i) {
        vecsBackground[i].m_adResults = adResults;
        vecsBackground[i].m_iLength = iChunk;
        vecsBackground[i].m_iStart = i * iChunk;
        vecsBackground[i].m_pDat = &Dat;
        vecsBackground[i].m_pvecveciSets = &vecveciSets;
        if (pthread_create(&vecpthdThreads[i], nullptr, background_thread, &vecsBackground[i])) {
            cerr << "Couldn't create thread for set group: " << i << endl;
            return 1;
        }
    }
    for (i = 0; i < vecpthdThreads.size(); ++i)
        pthread_join(vecpthdThreads[i], nullptr);

    return 0;
}

void *background_thread(void *pData) {
    SBackground *psData;
    vector<float> vecdBackgrounds, vecdBackground;
    size_t i, j, k;
    float d;

    psData = (SBackground *) pData;
    for (i = 0; (i < psData->m_iLength) && ((psData->m_iStart + i) < psData->m_pvecveciSets->size());
         ++i) {
        size_t iSet = psData->m_iStart + i;
        const vector <size_t> &veciSet = psData->m_pvecveciSets->at(iSet);

        vecdBackgrounds.resize(veciSet.size());
        for (j = 0; j < veciSet.size(); ++j) {
            size_t iOne = veciSet[j];

            vecdBackground.clear();
            vecdBackground.reserve(psData->m_pDat->GetGenes());
            for (k = 0; k < psData->m_pDat->GetGenes(); ++k)
                if (!CMeta::IsNaN(d = psData->m_pDat->Get(iOne, k)))
                    vecdBackground.push_back(d);
            CStatistics::Winsorize(vecdBackground, (vecdBackground.size() / 10) + 1);
            vecdBackgrounds[j] = (float) CStatistics::Average(vecdBackground);
        }
//			vecdBackgrounds[ j ] = (float)CStatistics::Median( vecdBackground ); }
//		psData->m_adResults[ iSet ] = CStatistics::Median( vecdBackgrounds ); }
        CStatistics::Winsorize(vecdBackgrounds, (vecdBackgrounds.size() / 10) + 1);
        psData->m_adResults[iSet] = CStatistics::Average(vecdBackgrounds);
    }

    return nullptr;
}

int within(const SProcessor &sProcessor) {
    const CDat &Dat = sProcessor.m_Dat;
    const vector <vector<size_t>> &vecveciSets = sProcessor.m_vecveciSets1;
    long double *adResults = sProcessor.m_adResults;
    size_t i, j, iCount, iChunk;
    float d;
    long double dPrior;
    vector <pthread_t> vecpthdThreads;
    vector <SWithin> vecsWithin;

    for (dPrior = 0, iCount = i = 0; i < Dat.GetGenes(); ++i)
        for (j = (i + 1); j < Dat.GetGenes(); ++j)
            if (!CMeta::IsNaN(d = Dat.Get(i, j))) {
                iCount++;
                dPrior += d;
            }
    dPrior /= iCount ? iCount : 1;

    vecpthdThreads.resize(min(sProcessor.m_iThreads, vecveciSets.size()));
    vecsWithin.resize(vecpthdThreads.size());
    iChunk = 1 + ((vecveciSets.size() - 1) / vecpthdThreads.size());
    for (i = 0; i < vecpthdThreads.size(); ++i) {
        vecsWithin[i].m_dPrior = dPrior;
        vecsWithin[i].m_pDat = &Dat;
        vecsWithin[i].m_adResults = adResults;
        vecsWithin[i].m_pvecveciSets = &vecveciSets;
        vecsWithin[i].m_iStart = i * iChunk;
        vecsWithin[i].m_iLength = iChunk;
        if (pthread_create(&vecpthdThreads[i], nullptr, within_thread, &vecsWithin[i])) {
            cerr << "Couldn't create thread for set group: " << i << endl;
            return 1;
        }
    }
    for (i = 0; i < vecpthdThreads.size(); ++i)
        pthread_join(vecpthdThreads[i], nullptr);

    return 0;
}

void *within_thread(void *pData) {
    SWithin *psData;
    size_t i, j, k;
    float d;
    vector<float> vecdWithin;

    psData = (SWithin *) pData;
    for (i = 0; (i < psData->m_iLength) && ((psData->m_iStart + i) < psData->m_pvecveciSets->size());
         ++i) {
        size_t iSet = psData->m_iStart + i;
        const vector <size_t> &veciSet = psData->m_pvecveciSets->at(iSet);

        vecdWithin.clear();
        vecdWithin.reserve(veciSet.size() * (veciSet.size() + 1) / 2);
        for (j = 0; j < veciSet.size(); ++j) {
            vecdWithin.push_back((float) psData->m_dPrior);
            for (k = (j + 1); k < veciSet.size(); ++k)
                if (!CMeta::IsNaN(d = psData->m_pDat->Get(veciSet[j], veciSet[k])))
                    vecdWithin.push_back(d);
        }
        CStatistics::Winsorize(vecdWithin, (vecdWithin.size() / 10) + 1);
        psData->m_adResults[iSet] = CStatistics::Average(vecdWithin);
    }
//		psData->m_adResults[ iSet ] = CStatistics::Median( vecdWithin ); }

    return nullptr;
}

int between(const SProcessor &sProcessor) {
    const CDat &Dat = sProcessor.m_Dat;
    const vector <vector<size_t>> &vecveciSets1 = sProcessor.m_vecveciSets1;
    const vector <vector<size_t>> &vecveciSets2 = sProcessor.m_vecveciSets2;
    long double *adResults = sProcessor.m_adResults;
    vector <pthread_t> vecpthdThreads;
    vector <SBetween> vecsBetween;
    size_t i, iChunk;

    vecpthdThreads.resize(min(sProcessor.m_iThreads, vecveciSets1.size()));
    vecsBetween.resize(vecpthdThreads.size());
    iChunk = 1 + ((vecveciSets1.size() - 1) / vecpthdThreads.size());
    for (i = 0; i < vecpthdThreads.size(); ++i) {
        vecsBetween[i].m_adResults = adResults;
        vecsBetween[i].m_iLength = iChunk;
        vecsBetween[i].m_iStart = i * iChunk;
        vecsBetween[i].m_pDat = &Dat;
        vecsBetween[i].m_pvecveciSets1 = &vecveciSets1;
        vecsBetween[i].m_pvecveciSets2 = &vecveciSets2;
        if (pthread_create(&vecpthdThreads[i], nullptr, between_thread, &vecsBetween[i])) {
            cerr << "Couldn't create thread for set group: " << i << endl;
            return 1;
        }
    }
    for (i = 0; i < vecpthdThreads.size(); ++i)
        pthread_join(vecpthdThreads[i], nullptr);

    return 0;
}

void *between_thread(void *pData) {
    size_t i, iSetOne, iSetTwo, iGeneOne, iGeneTwo, iOne, iTwo;
    float d;
    SBetween *psData;
    vector<long double> vecdBetween;

    psData = (SBetween *) pData;
    for (iSetOne = 0; (iSetOne < psData->m_iLength) &&
                      ((psData->m_iStart + iSetOne) < psData->m_pvecveciSets1->size()); ++iSetOne) {
        const vector <size_t> &veciOne = (*psData->m_pvecveciSets1)[psData->m_iStart + iSetOne];

        cerr << iSetOne << '/' << psData->m_iLength << endl;
        for (iSetTwo = 0; iSetTwo < psData->m_pvecveciSets2->size(); ++iSetTwo) {
            const vector <size_t> &veciTwo = (*psData->m_pvecveciSets2)[iSetTwo];

            vecdBetween.resize(veciOne.size() + veciTwo.size());
            fill(vecdBetween.begin(), vecdBetween.end(), 0);
            for (iGeneOne = 0; iGeneOne < veciOne.size(); ++iGeneOne) {
                iOne = veciOne[iGeneOne];
                for (iGeneTwo = 0; iGeneTwo < veciTwo.size(); ++iGeneTwo) {
                    iTwo = veciTwo[iGeneTwo];
                    d = (iOne == iTwo) ? 1 : psData->m_pDat->Get(iOne, iTwo);
                    vecdBetween[iGeneOne] += d;
                    vecdBetween[veciOne.size() + iGeneTwo] += d;
                }
            }
            for (i = 0; i < veciOne.size(); ++i)
                vecdBetween[i] /= max((size_t) 1, veciTwo.size());
            for (i = 0; i < veciTwo.size(); ++i)
                vecdBetween[veciOne.size() + i] /= max((size_t) 1, veciOne.size());
            CStatistics::Winsorize(vecdBetween, (vecdBetween.size() / 10) + 1);
            psData->m_adResults[((psData->m_iStart + iSetOne) * psData->m_pvecveciSets2->size()) +
                                iSetTwo] = CStatistics::Average(vecdBetween);
        }
    }
//				iSetTwo ] = CStatistics::Median( vecdBetween ); } }

    return nullptr;
}

int sets(const char *szFile, const vector <string> &vecstrGenes, vector <vector<string>> &vecvecstrSets) {
    CPCL PCLSets(false);
    size_t i, iSet, iSets;

    if (!PCLSets.Open(szFile, 1)) {
        cerr << "Could not open: " << szFile << endl;
        return 1;
    }
    for (iSets = i = 0; i < PCLSets.GetGenes(); ++i) {
        if (!(iSet = strtol(PCLSets.GetGene(i).c_str(), nullptr, 10))) {
            cerr << "Invalid set: " << PCLSets.GetGene(i) << endl;
            return 1;
        }
        if (iSet > iSets)
            iSets = iSet;
    }
    vecvecstrSets.resize(iSets);
    for (i = 0; i < PCLSets.GetGenes(); ++i) {
        size_t iGene;

        iSet = strtol(PCLSets.GetGene(i).c_str(), nullptr, 10) - 1;
        iGene = strtol(PCLSets.GetFeature(i, 1).c_str(), nullptr, 10) - 1;
        if (iGene >= vecstrGenes.size()) {
            cerr << "Unknown gene: " << iGene << " (" << PCLSets.GetFeature(i, 1) << ") in set " <<
                 PCLSets.GetGene(i) << endl;
            return 1;
        }
        vecvecstrSets[iSet].push_back(vecstrGenes[iGene]);
    }

    return 0;
}

void enset(const CDat &Dat, const vector <vector<string>> &vecvecstrSets,
           vector <vector<size_t>> &vecveciSets) {
    size_t i, j, iGene;
    map <string, size_t> mapstriGenes;

    for (i = 0; i < vecvecstrSets.size(); ++i)
        for (j = 0; j < vecvecstrSets[i].size(); ++j) {
            const string &strGene = vecvecstrSets[i][j];

            if (mapstriGenes.find(strGene) == mapstriGenes.end())
                mapstriGenes[strGene] = Dat.GetGene(strGene);
        }

    vecveciSets.resize(vecvecstrSets.size());
    for (i = 0; i < vecveciSets.size(); ++i) {
        vector <size_t> &veciSet = vecveciSets[i];
        const vector <string> &vecstrSet = vecvecstrSets[i];

        veciSet.clear();
        veciSet.reserve(vecstrSet.size());
        for (j = 0; j < vecstrSet.size(); ++j)
            if ((iGene = mapstriGenes[vecstrSet[j]]) != -1)
                veciSet.push_back(iGene);
    }
}

void *term_thread(void *pData) {
    CGenome Genome;
    CGenes Genes(Genome);
    ifstream ifsm;
    size_t i, iCur;
    STerm *psData;
    SDatum sDatum;

    psData = (STerm *) pData;
    ifsm.open(psData->m_strFile.c_str());
    if (!Genes.Open(ifsm)) {
        cerr << "Could not open: " << psData->m_strFile << endl;
        return nullptr;
    }
    ifsm.close();
    iCur = cliques(*psData->m_pDat, psData->m_iTotal, *psData->m_pvecsHubs, psData->m_iGenes, sDatum, &Genes,
                   *psData->m_pPCLWeights, *psData->m_pveciDat2PCL);
    pthread_mutex_lock(psData->m_pmutx);
    cout << CMeta::Basename(psData->m_strFile.c_str());
    if (psData->m_iGenes == -1)
        for (i = 0; i < sDatum.m_vecprSpecific.size(); ++i)
            cout << '\t' << (sDatum.m_vecprSpecific[i].second *
                             (Genes.IsGene(psData->m_pDat->GetGene(sDatum.m_vecprSpecific[i].first)) ? -1 : 1));
    else {
        cout << '\t' << iCur << '\t' << sDatum.m_sHubbiness.GetAverage() << '\t' <<
             sDatum.m_sHubbiness.GetStdev() << '\t' << sDatum.m_sHubbiness.m_iCount << '\t'
             << sDatum.m_sCliquiness.GetAverage() << '\t' <<
             sDatum.m_sCliquiness.GetStdev() << '\t' << sDatum.m_sCliquiness.m_iCount;
        for (i = 0; i < min(psData->m_iGenes, sDatum.m_vecprSpecific.size()); ++i)
            cout << '\t' << psData->m_pDat->GetGene(sDatum.m_vecprSpecific[i].first) << '|' <<
                 sDatum.m_vecprSpecific[i].second << '|' <<
                 (Genes.IsGene(psData->m_pDat->GetGene(sDatum.m_vecprSpecific[i].first)) ? 1 : 0);
    }
    cout << endl;
    pthread_mutex_unlock(psData->m_pmutx);

    if (psData->m_strClip.length()) {
        CDat DatOut;
        size_t j;
        vector<bool> vecfGenes;

        vecfGenes.resize(psData->m_pDat->GetGenes());
        for (i = 0; i < vecfGenes.size(); ++i)
            vecfGenes[i] = Genes.IsGene(psData->m_pDat->GetGene(i));
        DatOut.Open(psData->m_pDat->GetGeneNames());
        for (i = 0; i < DatOut.GetGenes(); ++i)
            for (j = (i + 1); j < DatOut.GetGenes(); ++j)
                if (vecfGenes[i] || vecfGenes[j])
                    DatOut.Set(i, j, psData->m_pDat->Get(i, j));
        DatOut.Save((psData->m_strClip + '/' + CMeta::Deextension(CMeta::Basename(
                psData->m_strFile.c_str())) + c_szDab).c_str());
    }

    return nullptr;
}

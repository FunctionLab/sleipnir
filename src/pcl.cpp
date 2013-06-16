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
#include "pcl.h"
#include "meta.h"
#include "statistics.h"
#include "genome.h"
#include "measure.h"
#include "dat.h"

namespace Sleipnir {

const char CPCLImpl::c_szEWEIGHT[] = "EWEIGHT";
const char CPCLImpl::c_szGENE[] = "GENE";
const char CPCLImpl::c_szGID[] = "GID";
const char CPCLImpl::c_szGWEIGHT[] = "GWEIGHT";
const char CPCLImpl::c_szNAME[] = "NAME";
const char CPCLImpl::c_szOne[] = "1";
const char CPCLImpl::c_szExtension[] = ".pcl";
const char CPCLImpl::c_szBinExtension[] = ".bin";
const char CPCLImpl::c_szDabExtension[] = ".dab";

struct SNeighbors {
	CFullMatrix<pair<size_t, float> > m_MatNeighbors;
	vector<size_t> m_veciMin;
	vector<size_t> m_veciColumns;

	void Initialize(const float* adValues, size_t iValues, size_t iNeighbors) {
		size_t i, j;

		m_veciColumns.clear();
		for (i = 0; i < iValues; ++i)
			if (CMeta::IsNaN(adValues[i]))
				m_veciColumns.push_back(i);

		m_veciMin.resize(m_veciColumns.size());
		m_MatNeighbors.Initialize(iNeighbors, m_veciColumns.size());
		for (i = 0; i < m_MatNeighbors.GetRows(); ++i)
			for (j = 0; j < m_MatNeighbors.GetColumns(); ++j)
				m_MatNeighbors.Get(i, j).second = -FLT_MAX;
	}

	bool Add(size_t iNeighbor, float dSim, const float* adValues) {
		size_t i, j, iCol;
		bool fRet;

		if (!m_MatNeighbors.GetRows())
			return false;
		for (fRet = false, i = 0; i < m_veciColumns.size(); ++i) {
			iCol = m_veciColumns[i];
			if (!CMeta::IsNaN(adValues[iCol]) && (dSim > m_MatNeighbors.Get(
					m_veciMin[i], i).second)) {
				fRet = true;
				m_MatNeighbors.Get(m_veciMin[i], i).first = iNeighbor;
				m_MatNeighbors.Get(m_veciMin[i], i).second = dSim;

				for (m_veciMin[i] = 0, j = 1; j < m_MatNeighbors.GetRows(); ++j)
					if (m_MatNeighbors.Get(j, i).second < m_MatNeighbors.Get(
							m_veciMin[i], i).second)
						m_veciMin[i] = j;
			}
		}

		return fRet;
	}

	size_t GetColumn(size_t iColumn) const {
		size_t i;

		for (i = 0; i < m_veciColumns.size(); ++i)
			if (m_veciColumns[i] == iColumn)
				return i;

		return -1;
	}
};

/*!
 * \brief
 * Kitchen sink method for completely loading a PCL and calculating its pairwise similarity scores into a CDat.
 *
 * \param szFile
 * If non-null, file from which PCL is to be loaded; if null, standard input is used.
 *
 * \param iSkip
 * Number of columns to skip in the PCL file between gene IDs and experimental data.
 *
 * \param szSimilarityMeasure
 * String identifier of similarity measure to use for CDat generation.
 *
 * \param fNormalize
 * If true, normalize the generated CDat to the range [0, 1].
 *
 * \param fZScore
 * If true, normalize the generated CDat to z-scores (subtract mean, divide by standard deviation).
 *
 * \param fAutocorrelate
 * If true, autocorrelate the requested similarity measure.
 *
 * \param szGeneFile
 * If non-null, only convert genes in the given file to pairwise scores in the CDat.
 *
 * \param dCutoff
 * If finite, remove all pairwise scores less than the given cutoff.
 *
 * \param iLimit
 * If not equal to -1 and the PCL contains more genes than this limit, do not precalculate pairwise scores;
 * instead, configure the CDat to calculate scores on the fly as needed from the PCL.
 *
 * \param PCL
 * Output PCL with the loaded data.
 *
 * \param Dat
 * Output CDat with the calculated pairwise scores.
 *
 * \param eMap
 * Way in which returned value should be centered (implementation-specific).
 *
 * \param fFrequencyWeight
 * If true, weight each condition by the frequency with which it is nonzero over all genes.
 *
 * \returns
 * 0 on successes, a nonzero value on failure.
 *
 * The (many) steps performed by this method are as follows:
 * <ol>
 * <li>A PCL is loaded from szFile into PCL using Open.  If szFile is null, the PCL is loaded from
 * standard input.</li>
 * <li>An IMeasure object is constructed by iterating over the available implementations and finding
 * one whose name corresponds with szSimilarityMeasure.  Names which are distance measures (e.g. Euclidean)
 * are automatically inverted.  If requested, the measure is autocorrelated.</li>
 * <li>If given, szGeneFile is loaded into a CGenes object.</li>
 * <li>If iLimit is not -1 and the PCL contains more than the limiting number of genes, Dat is given a
 * reference to PCL and configured to calculate pairwise scores as needed.  Processing then stops.</li>
 * <li>Otherwise, an empty CDat is initialized to contain either the genes in szGeneFile (if given) or
 * all of the genes in the PCL.</li>
 * <li>Gene pairs are assigned scores in the CDat using the requested similarity measure.</li>
 * <li>If given, scores below dCutoff are replaced with missing values.</li>
 * <li>If requested, the remaining scores are normalized either to the range [0, 1] or to z-scores.</li>
 * </ol>
 *
 * \remarks
 * This method is written to make it easy to construct tools that must load a PCL and immediately convert
 * it to pairwise scores using some user-selected similarity measure.
 */
int CPCL::Distance(const char* szFile, size_t iSkip,
		const char* szSimilarityMeasure, bool fNormalize, bool fZScore,
		bool fAutocorrelate, const char* szGeneFile, float dCutoff,
		size_t iLimit, CPCL& PCL, CDat& Dat, IMeasure::EMap eMap,
		bool fFrequencyWeight, float dAlpha, int nThreads) {


	size_t i, j, iOne, iTwo;
	float d;
	ifstream ifsm;
	vector<string> vecstrGenes;
	CGenome Genome;
	CGenes GenesIn(Genome);
	vector<size_t> veciGenes;
	const float* adOne;
	const float* adWeights;
	vector<float> vecdWeights;
	IMeasure* pMeasure;
	CMeasurePearson Pearson;
	CMeasureEuclidean Euclidean;
	CMeasureKendallsTau KendallsTau;
	CMeasureKolmogorovSmirnov KolmSmir;
	CMeasureSpearman Spearman(true);
	CMeasurePearNorm PearNorm;
	CMeasureHypergeometric Hypergeom;
	CMeasureQuickPearson PearQuick;
	CMeasureInnerProduct InnerProd;
	CMeasureBinaryInnerProduct BinInnerProd;
	CMeasureMutualInformation MutualInfo;
	CMeasureRelativeAUC RelAuc;
	CMeasurePearsonSignificance PearSig;
	CMeasureDice Dice( dAlpha );
	CMeasureDistanceCorrelation DCor;
	CMeasureSignedDistanceCorrelation SDCor;
	if (szFile) {
        g_CatSleipnir().debug("Opening PCL for distance");
        if (!PCL.Open(szFile, iSkip, false, false)) {
			g_CatSleipnir().error(
					"CPCL::Distance( %s, %d, %s, %d, %d, %d, %s, %g ) failed to open PCL",
					szFile, iSkip, szSimilarityMeasure, fNormalize, fZScore,
					fAutocorrelate, szGeneFile ? szGeneFile : "", dCutoff);
			return 1;
		}
		ifsm.close();
	} else if (!PCL.Open(cin, iSkip)) {
		g_CatSleipnir().error(
				"CPCL::Distance( %s, %d, %s, %d, %d, %d, %s, %g ) failed to open PCL",
				"stdin", iSkip, szSimilarityMeasure, fNormalize, fZScore,
				fAutocorrelate, szGeneFile ? szGeneFile : "", dCutoff);
		return 1;
	}
	if (!szSimilarityMeasure)
		return 0;

	CMeasureSigmoid
			EuclideanSig(&Euclidean, false, 1.0f / PCL.GetExperiments());
	IMeasure* apMeasures[] = { &Pearson, &EuclideanSig, &KendallsTau,
			&KolmSmir, &Spearman, &PearNorm, &Hypergeom, &PearQuick,
			&InnerProd, &BinInnerProd, &MutualInfo, &RelAuc, &PearSig, &Dice,&DCor,&SDCor, NULL };

	pMeasure = NULL;
	for (i = 0; apMeasures[i]; ++i)
		if (!strcmp(apMeasures[i]->GetName(), szSimilarityMeasure)) {
			pMeasure = apMeasures[i];
			break;
		}
	if (!pMeasure)
		return 1;
	g_CatSleipnir().info("Method: %s ", pMeasure->GetName(), i);
	if (fFrequencyWeight) {
		vecdWeights.resize(PCL.GetExperiments());
		for (i = 0; i < vecdWeights.size(); ++i) {
			for (iOne = j = 0; j < PCL.GetGenes(); ++j)
				if (!CMeta::IsNaN(d = PCL.Get(j, i)) && (d > 0))
					iOne++;
			vecdWeights[i] = (float) ((PCL.GetGenes() + 1) - iOne)
					/ PCL.GetGenes();
		}
	}
	adWeights = vecdWeights.empty() ? NULL : &vecdWeights[0];

	CMeasureAutocorrelate Autocorrelate(pMeasure, false);
	if (fAutocorrelate)
		pMeasure = &Autocorrelate;

	if (szGeneFile) {
		ifsm.clear();
		ifsm.open(szGeneFile);
		if (!GenesIn.Open(ifsm)) {
			g_CatSleipnir().error(
					"CPCL::Distance( %s, %d, %s, %d, %d, %d, %s, %g ) failed to open genes",
					szFile ? szFile : "stdin", iSkip, szSimilarityMeasure,
					fNormalize, fZScore, fAutocorrelate, szGeneFile, dCutoff);
			return 1;
		}
		ifsm.close();
	} else
		GenesIn.Open(PCL.GetGeneNames());
	veciGenes.resize(GenesIn.GetGenes());
	for (i = 0; i < veciGenes.size(); ++i)
		veciGenes[i] = szGeneFile ? PCL.GetGene(GenesIn.GetGene(i).GetName())
				: i;

	if (pMeasure->IsRank())
		PCL.RankTransform();

	if ((iLimit != -1) && (PCL.GetGenes() > iLimit))
		Dat.Open(PCL, pMeasure->Clone(), true);
	else {
		Dat.Open(GenesIn.GetGeneNames());
		for (i = 0; i < Dat.GetGenes(); ++i)
			for (j = (i + 1); j < Dat.GetGenes(); ++j)
				Dat.Set(i, j, CMeta::GetNaN());
		int origNThreads = omp_get_num_threads();
		omp_set_num_threads(nThreads);
		for (i = 0; i < GenesIn.GetGenes(); ++i) {
			if (!(i % 100))
				g_CatSleipnir().info(
						"CPCL::Distance( %s, %d, %s, %d, %d, %d, %s, %g ) processing gene %d/%d",
						szFile ? szFile : "stdin", iSkip, szSimilarityMeasure,
						fNormalize, fZScore, fAutocorrelate,
						szGeneFile ? szGeneFile : "", dCutoff, i,
						GenesIn.GetGenes());
			if ((iOne = veciGenes[i]) == -1)
				continue;
			adOne = PCL.Get(iOne);
			for (j = (i + 1); j < GenesIn.GetGenes(); ++j){
				if ((iTwo = veciGenes[j]) != -1){
					Dat.Set(i, j, (float) pMeasure->Measure(adOne,
							PCL.GetExperiments(), PCL.Get(iTwo),
							PCL.GetExperiments(), eMap, adWeights, adWeights));
				}
			}
		}
		omp_set_num_threads(origNThreads);
		if (fNormalize || fZScore)
			Dat.Normalize(fZScore ? CDat::ENormalizeZScore
					: CDat::ENormalizeMinMax);
		if (!CMeta::IsNaN(dCutoff))
			for (i = 0; i < Dat.GetGenes(); ++i)
				for (j = (i + 1); j < Dat.GetGenes(); ++j)
					if (!CMeta::IsNaN(d = Dat.Get(i, j)) && (d < dCutoff))
						Dat.Set(i, j, CMeta::GetNaN());

		if(pMeasure->GetName()=="pearson" || pMeasure->GetName()=="pearnorm"){
			vector<unsigned long> bins;
			bins.resize(55);
			float upper = 0;
			float lower = 0;
			if(pMeasure->GetName()=="pearson"){
				upper = 1.0;
				lower = -1.0;
			}else if(pMeasure->GetName()=="pearnorm"){
				upper = 5.0;
				lower = -5.0;
			}
			float bin_size = (upper - lower) / 50;
			for(i=0; i<55; i++)
				bins[i] = 0;
			for(i=0; i<Dat.GetGenes(); i++){
				for(j=i+1; j<Dat.GetGenes(); j++){
					d = Dat.Get(i,j);
					if(CMeta::IsNaN(d)) continue;
					int b = (int) ((d - lower) / bin_size);
					if(b<0){
						bins[0]++;
						continue;
					}
					if(b>=55){
						bins[54]++;
						continue;
					}
					bins[b]++;
				}
			}
			g_CatSleipnir().info(
			"Distribution of distances: bin size: %.5f, number of bins: %d, min bin value: %.5f, max bin value: %.5f",
			bin_size, 55, lower, upper);
			for(i=0; i<55; i++){
				g_CatSleipnir().info("%lu\t%lu", i, bins[i]);
			}
		}

	}

	return 0;
}

/*!
 * \brief
 * Kitchen sink method for completely loading a PCL and calculating its pairwise similarity scores into a CDat.  This method supports weighted metrics as opposed to the above method which does not.  This allows distancer to take a weights vector (currently included in the documentation but it was not previously functioning).  The weights vector is used only if freqweight is false.
 *
 * \param szFile
 * If non-null, file from which PCL is to be loaded; if null, standard input is used.
 *
 * \param iSkip
 * Number of columns to skip in the PCL file between gene IDs and experimental data.
 *
 * \param szWeights
 * If non-null, file from which weights are to be loaded; if null, this is ignored.
 *
 * \param szSimilarityMeasure
 * String identifier of similarity measure to use for CDat generation.
 *
 * \param fNormalize
 * If true, normalize the generated CDat to the range [0, 1].
 *
 * \param fZScore
 * If true, normalize the generated CDat to z-scores (subtract mean, divide by standard deviation).
 *
 * \param fAutocorrelate
 * If true, autocorrelate the requested similarity measure.
 *
 * \param szGeneFile
 * If non-null, only convert genes in the given file to pairwise scores in the CDat.
 *
 * \param dCutoff
 * If finite, remove all pairwise scores less than the given cutoff.
 *
 * \param iLimit
 * If not equal to -1 and the PCL contains more genes than this limit, do not precalculate pairwise scores;
 * instead, configure the CDat to calculate scores on the fly as needed from the PCL.
 *
 * \param PCL
 * Output PCL with the loaded data.
 *
 * \param Dat
 * Output CDat with the calculated pairwise scores.
 *
 * \param eMap
 * Way in which returned value should be centered (implementation-specific).
 *
 * \param fFrequencyWeight
 * If true, weight each condition by the frequency with which it is nonzero over all genes.
 *
 * \returns
 * 0 on successes, a nonzero value on failure.
 *
 * The (many) steps performed by this method are as follows:
 * <ol>
 * <li>A PCL is loaded from szFile into PCL using Open.  If szFile is null, the PCL is loaded from
 * standard input.</li>
 * <li>An IMeasure object is constructed by iterating over the available implementations and finding
 * one whose name corresponds with szSimilarityMeasure.  Names which are distance measures (e.g. Euclidean)
 * are automatically inverted.  If requested, the measure is autocorrelated.</li>
 * <li>If given, szGeneFile is loaded into a CGenes object.</li>
 * <li>If iLimit is not -1 and the PCL contains more than the limiting number of genes, Dat is given a
 * reference to PCL and configured to calculate pairwise scores as needed.  Processing then stops.</li>
 * <li>Otherwise, an empty CDat is initialized to contain either the genes in szGeneFile (if given) or
 * all of the genes in the PCL.</li>
 * <li>Gene pairs are assigned scores in the CDat using the requested similarity measure.</li>
 * <li>If given, scores below dCutoff are replaced with missing values.</li>
 * <li>If requested, the remaining scores are normalized either to the range [0, 1] or to z-scores.</li>
 * </ol>
 *
 * \remarks
 * This method is written to make it easy to construct tools that must load a PCL and immediately convert
 * it to pairwise scores using some user-selected similarity measure.
 */
int CPCL::Distance(const char* szFile, size_t iSkip, const char* szWeights,
		const char* szSimilarityMeasure, bool fNormalize, bool fZScore,
		bool fAutocorrelate, const char* szGeneFile, float dCutoff,
		size_t iLimit, CPCL& PCL, CDat& Dat, IMeasure::EMap eMap,
		bool fFrequencyWeight, float dAlpha, int nThreads) {
	size_t i, j, iOne, iTwo;
	float d;
	ifstream ifsm;
	vector<string> vecstrGenes;
	CGenome Genome;
	CGenes GenesIn(Genome);
	vector<size_t> veciGenes;
	const float* adOne;
	const float* adWeights;
	vector<float> vecdWeights;
	vector<string> vecstrWeights;
    string acStr;
	IMeasure* pMeasure;
	CMeasurePearson Pearson;
	CMeasureEuclidean Euclidean;
	CMeasureKendallsTau KendallsTau;
	CMeasureKolmogorovSmirnov KolmSmir;
	CMeasureSpearman Spearman(true);
	CMeasurePearNorm PearNorm;
	CMeasureHypergeometric Hypergeom;
	CMeasureQuickPearson PearQuick;
	CMeasureInnerProduct InnerProd;
	CMeasureBinaryInnerProduct BinInnerProd;
	CMeasureMutualInformation MutualInfo;
	CMeasureRelativeAUC RelAuc;
	CMeasurePearsonSignificance PearSig;
	CMeasureDice Dice( dAlpha );
	CMeasureDistanceCorrelation DCor;
	CMeasureSignedDistanceCorrelation SDCor;
	if (szFile) {
        g_CatSleipnir().debug("Opening PCL for distance");
        if (!PCL.Open(szFile, iSkip, false, false)) {
			g_CatSleipnir().error(
					"CPCL::Distance( %s, %d, %s, %d, %d, %d, %s, %g ) failed to open PCL",
					szFile, iSkip, szSimilarityMeasure, fNormalize, fZScore,
					fAutocorrelate, szGeneFile ? szGeneFile : "", dCutoff);
			return 1;
		}
		ifsm.close();
	} else if (!PCL.Open(cin, iSkip)) {
		g_CatSleipnir().error(
				"CPCL::Distance( %s, %d, %s, %d, %d, %d, %s, %g ) failed to open PCL",
				"stdin", iSkip, szSimilarityMeasure, fNormalize, fZScore,
				fAutocorrelate, szGeneFile ? szGeneFile : "", dCutoff);
		return 1;
	}
	if (!szSimilarityMeasure)
		return 0;

	CMeasureSigmoid
			EuclideanSig(&Euclidean, false, 1.0f / PCL.GetExperiments());
	IMeasure* apMeasures[] = { &Pearson, &EuclideanSig, &KendallsTau,
			&KolmSmir, &Spearman, &PearNorm, &Hypergeom, &PearQuick,
			&InnerProd, &BinInnerProd, &MutualInfo, &RelAuc, &PearSig, &Dice,&DCor,&SDCor, NULL };

	pMeasure = NULL;
	for (i = 0; apMeasures[i]; ++i)
		if (!strcmp(apMeasures[i]->GetName(), szSimilarityMeasure)) {
			pMeasure = apMeasures[i];
			break;
		}
	if (!pMeasure)
		return 1;
	g_CatSleipnir().info("Method: %s ", pMeasure->GetName(), i);
	if (fFrequencyWeight) {
		vecdWeights.resize(PCL.GetExperiments());
		for (i = 0; i < vecdWeights.size(); ++i) {
			for (iOne = j = 0; j < PCL.GetGenes(); ++j)
				if (!CMeta::IsNaN(d = PCL.Get(j, i)) && (d > 0))
					iOne++;
			vecdWeights[i] = (float) ((PCL.GetGenes() + 1) - iOne)
					/ PCL.GetGenes();
		}
	}
	adWeights = vecdWeights.empty() ? NULL : &vecdWeights[0];

    if ( adWeights == NULL ) {
        if (szWeights != NULL) {
            ifsm.open(szWeights);
            getline(ifsm, acStr);
            ifsm.close();
            CMeta::Tokenize( acStr.c_str(), vecstrWeights, CMeta::c_szWS, true);
            if ( vecstrWeights.size() != PCL.GetExperiments()) {
                g_CatSleipnir().error(
                        "CPCL::Distance( %s, %d, %s, %d, %d, %d, %s, %g ) length of weights file (%s, length: %d) did not match number of experiments (%d).",
                        szFile, iSkip, szSimilarityMeasure, fNormalize, fZScore,
                        fAutocorrelate, szGeneFile ? szGeneFile : "", dCutoff, szWeights,
                        vecstrWeights.size(), PCL.GetExperiments());
                return 1;
            }
		    vecdWeights.resize(PCL.GetExperiments());
            for (i = 0; i < vecstrWeights.size(); ++i) {
                vecdWeights[i] = atof(vecstrWeights[i].c_str());
            }
	        adWeights = vecdWeights.empty() ? NULL : &vecdWeights[0];
        }
    }

	CMeasureAutocorrelate Autocorrelate(pMeasure, false);
	if (fAutocorrelate)
		pMeasure = &Autocorrelate;

	if (szGeneFile) {
		ifsm.clear();
		ifsm.open(szGeneFile);
		if (!GenesIn.Open(ifsm)) {
			g_CatSleipnir().error(
					"CPCL::Distance( %s, %d, %s, %d, %d, %d, %s, %g ) failed to open genes",
					szFile ? szFile : "stdin", iSkip, szSimilarityMeasure,
					fNormalize, fZScore, fAutocorrelate, szGeneFile, dCutoff);
			return 1;
		}
		ifsm.close();
	} else
		GenesIn.Open(PCL.GetGeneNames());
	veciGenes.resize(GenesIn.GetGenes());
	for (i = 0; i < veciGenes.size(); ++i)
		veciGenes[i] = szGeneFile ? PCL.GetGene(GenesIn.GetGene(i).GetName())
				: i;

	if (pMeasure->IsRank())
		PCL.RankTransform();

	if ((iLimit != -1) && (PCL.GetGenes() > iLimit))
		Dat.Open(PCL, pMeasure->Clone(), true);
	else {
		Dat.Open(GenesIn.GetGeneNames());
		for (i = 0; i < Dat.GetGenes(); ++i)
			for (j = (i + 1); j < Dat.GetGenes(); ++j)
				Dat.Set(i, j, CMeta::GetNaN());
		for (i = 0; i < GenesIn.GetGenes(); ++i) {
			if (!(i % 100))
				g_CatSleipnir().info(
						"CPCL::Distance( %s, %d, %s, %d, %d, %d, %s, %g ) processing gene %d/%d",
						szFile ? szFile : "stdin", iSkip, szSimilarityMeasure,
						fNormalize, fZScore, fAutocorrelate,
						szGeneFile ? szGeneFile : "", dCutoff, i,
						GenesIn.GetGenes());
			if ((iOne = veciGenes[i]) == -1)
				continue;
			adOne = PCL.Get(iOne);
			#pragma omp parallel for num_threads(nThreads)
			for (j = (i + 1); j < GenesIn.GetGenes(); ++j){
				if ((iTwo = veciGenes[j]) != -1){
					Dat.Set(i, j, (float) pMeasure->Measure(adOne,
							PCL.GetExperiments(), PCL.Get(iTwo),
							PCL.GetExperiments(), eMap, adWeights, adWeights));
				}
			}
		}

		if (fNormalize || fZScore)
			Dat.Normalize(fZScore ? CDat::ENormalizeZScore
					: CDat::ENormalizeMinMax);

		if (!CMeta::IsNaN(dCutoff))
			for (i = 0; i < Dat.GetGenes(); ++i)
				for (j = (i + 1); j < Dat.GetGenes(); ++j)
					if (!CMeta::IsNaN(d = Dat.Get(i, j)) && (d < dCutoff))
						Dat.Set(i, j, CMeta::GetNaN());

		if(pMeasure->GetName()=="pearson" || pMeasure->GetName()=="pearnorm"){
			vector<unsigned long> bins;
			bins.resize(55);
			float upper = 0;
			float lower = 0;
			if(pMeasure->GetName()=="pearson"){
				upper = 1.0;
				lower = -1.0;
			}else if(pMeasure->GetName()=="pearnorm"){
				upper = 5.0;
				lower = -5.0;
			}
			float bin_size = (upper - lower) / 50;
			for(i=0; i<55; i++)
				bins[i] = 0;
			for(i=0; i<Dat.GetGenes(); i++){
				for(j=i+1; j<Dat.GetGenes(); j++){
					d = Dat.Get(i,j);
					if(CMeta::IsNaN(d)) continue;
					int b = (int) ((d - lower) / bin_size);
					if(b<0){
						bins[0]++;
						continue;
					}
					if(b>=55){
						bins[54]++;
						continue;
					}
					bins[b]++;
				}
			}
			g_CatSleipnir().info(
			"Distribution of distances: bin size: %.5f, number of bins: %d, min bin value: %.5f, max bin value: %.5f",
			bin_size, 55, lower, upper);
			for(i=0; i<55; i++){
				g_CatSleipnir().info("%lu\t%lu", i, bins[i]);
			}
		}
	}
	return 0;
}

CPCLImpl::~CPCLImpl() {

	Reset();
}

void CPCLImpl::Reset() {

	m_Data.Reset();
	m_vecstrGenes.clear();
	m_vecstrExperiments.clear();
	m_vecstrFeatures.clear();
	m_vecvecstrFeatures.clear();
	m_setiGenes.clear();
	m_mapstriGenes.clear();
}

/*!
 * \brief
 * Create a new PCL by reading from text or binary
 *
 * \param szFile
 * File from which PCL file is loaded.
 *
 * \returns
 * True if the PCL was opened successfully.
 *
 */
//bool CPCL::Open(const char* szFile, size_t iSkip, bool Memmap) {
//	std::ifstream ifsm;
//	if (szFile)
//		ifsm.open(szFile);
//
//	// is this a binary binary file?
//	if (!strcmp(szFile + strlen(szFile) - strlen(c_szBinExtension),
//			c_szBinExtension))
//		return OpenBinary(szFile ? ifsm : cin, Memmap);
//	else
//		// is this a text based PCL file?
//		return Open(szFile ? ifsm : cin, iSkip);
//}


bool CPCL::Open(const char* szFile, size_t iSkip, bool Memmap, bool rTable) {
	bool isBinary = false;
	bool isDAB = false;
	ifstream ifsm;
	if (!strcmp(szFile + strlen(szFile) - strlen(c_szBinExtension), c_szBinExtension)) {
		isBinary = true;
	}
	if (!strcmp(szFile + strlen(szFile) - strlen(c_szDabExtension), c_szDabExtension)) {
		isDAB = true;
	}
	if (isBinary && Memmap) {
		g_CatSleipnir().notice("CPCL::Open, openning with memory mapping");
		Reset();
		if (!CMeta::MapRead(m_abData, m_hndlData, m_iData, szFile)) {
			g_CatSleipnir().error("CPCL::Open( %s, %d ) failed memory mapping",
					szFile, Memmap);

			return false;
		}
		bool fRet = OpenHelper();
		return fRet;
	}
	else if (isBinary) {
		ifsm.open(szFile, ios::binary);
		return OpenBinary(ifsm);
	} else if (isDAB) {
		CDat dat;
		if (!dat.Open(szFile,false, 0, false, false)) {
		    cerr << "Could not open dab to make pcl." << endl;
		    return false;
		}
		std::vector<std::string> features;
		std::vector<std::string> genes = dat.GetGeneNames();
		Open(genes, genes, features);
		for (size_t i = 0; i < dat.GetGenes(); i++) {
		    for (size_t j = i; j < dat.GetGenes(); j++) {
			if (i == j) {
			    Set(i, j, 1);
			}
			else {
			    float pairVal = dat.Get(i,j);
			    Set(i, j, pairVal);
			    Set(j, i, pairVal);
			}
		    }
		}
		return true;
	} else {
		ifstream ifsm;
		if (szFile)
		  ifsm.open(szFile);
		return Open(szFile ? ifsm : cin, iSkip, rTable);
	}
}



/*!
 * \brief
 * Save a PCL to the given text file.
 *
 * \param szFile
 * File into which PCL file is saved.
 *
 * \remarks
 * If null, output defaults to stdout.
 *
 * \see
 * Save
 */
void CPCL::Save(const char* szFile) {
	std::ofstream ofsm;

	if (szFile)
		ofsm.open(szFile);

	// Save as binary PCL file?
	if (!strcmp(szFile + strlen(szFile) - strlen(c_szExtension), c_szBinExtension))
		SaveBinary(szFile ? ofsm : cout);
	else
		Save(szFile ? ofsm : cout);
		// Save as text file

}

/*!
 * \brief
 * Create a new PCL by copying the given one.
 *
 * \param PCL
 * PCL to be copied into the current one.
 *
 * \remarks
 * All values are copied into newly allocated memory within the current PCL, so it's safe to destroy the
 * input PCL after the new one is opened.
 */
void CPCL::Open(const CPCL& PCL) {
	size_t i, j;
	TSetI::const_iterator iterGene;
	TMapStrI::const_iterator iterGene2;

	Reset();
	m_fHeader = PCL.m_fHeader;
	m_Data.Initialize(PCL.m_Data.GetRows(), PCL.m_Data.GetColumns());
	for (i = 0; i < m_Data.GetRows(); ++i)
		for (j = 0; j < m_Data.GetColumns(); ++j)
			m_Data.Set(i, j, PCL.m_Data.Get(i, j));

	for (iterGene = PCL.m_setiGenes.begin(); iterGene != PCL.m_setiGenes.end(); ++iterGene)
		m_setiGenes.insert(*iterGene);
	m_vecstrExperiments.resize(PCL.m_vecstrExperiments.size());
	copy(PCL.m_vecstrExperiments.begin(), PCL.m_vecstrExperiments.end(),
			m_vecstrExperiments.begin());
	m_vecstrFeatures.resize(PCL.m_vecstrFeatures.size());
	copy(PCL.m_vecstrFeatures.begin(), PCL.m_vecstrFeatures.end(),
			m_vecstrFeatures.begin());
	m_vecstrGenes.resize(PCL.m_vecstrGenes.size());
	copy(PCL.m_vecstrGenes.begin(), PCL.m_vecstrGenes.end(),
			m_vecstrGenes.begin());
	for (iterGene2 = PCL.m_mapstriGenes.begin(); iterGene2
			!= PCL.m_mapstriGenes.end(); ++iterGene2)
		m_mapstriGenes[iterGene2->first] = iterGene2->second;
	m_vecvecstrFeatures.resize(PCL.m_vecvecstrFeatures.size());
	for (i = 0; i < m_vecvecstrFeatures.size(); ++i) {
		m_vecvecstrFeatures[i].resize(PCL.m_vecvecstrFeatures[i].size());
		copy(PCL.m_vecvecstrFeatures[i].begin(),
				PCL.m_vecvecstrFeatures[i].end(),
				m_vecvecstrFeatures[i].begin());
	}
}

/*!
 * \brief
 * Load a PCL from the given text stream.
 *
 * \param istm
 * Stream from which PCL file is loaded.
 *
 * \param iSkip
 * Number of feature columns to skip between the gene IDs and first experimental column.
 *
 * \param rTable
 * Is this PCL file generated by R language (This means its missing the "YORF" label in the first line, thus has 1 less token compared to sleipnir PCL format)
 *
 * \returns
 * True if the PCL was opened successfully.
 *
 * \see
 * Save
 */

bool CPCL::Open(std::istream& istm, size_t iSkip, bool rTable) {
	vector<float> vecdData;
	size_t i;
	char* acBuf;
	bool fRet;
	string acStr;

	if (!istm.good())
		return false;

	acStr.reserve(c_iBufferSize);
	
	if (!OpenExperiments(istm, iSkip, acStr, rTable))
		fRet = false;
	else {
		m_vecvecstrFeatures.resize(m_vecstrFeatures.size() - 1);
		while (OpenGene(istm, vecdData, acStr))
			;
		for (fRet = !!GetGenes(), i = 0; i < GetGenes(); ++i)
			if (GetGene(i).empty() || !isprint(GetGene(i)[0])) {
				fRet = false;
				g_CatSleipnir().error(
						"CPCL::Open( %d ) invalid gene at index %d: %s", iSkip,
						i, GetGene(i).c_str());
				break;
			}
		if (fRet) {
			m_Data.Initialize(GetGenes(), GetExperiments());
			for (i = 0; i < m_Data.GetRows(); ++i)
				m_Data.Set(i, &vecdData[i * m_Data.GetColumns()]);
		}
	}
	//delete[] acBuf;

	return fRet;
}

bool CPCLImpl::OpenHelper() {
	unsigned char* pb;
	pb = m_abData;
	size_t i, j;
	m_vecstrFeatures.resize(*(uint32_t*) pb);
	pb += sizeof(uint32_t);
	for (i = 0; i < m_vecstrFeatures.size(); ++i) {
		string& strFeature = m_vecstrFeatures[i];
		strFeature.resize(*(uint32_t*) pb);
		pb += sizeof(uint32_t);
		for (j = 0; j < strFeature.size(); ++j)
			strFeature[j] = *pb++;
	}
	m_vecstrExperiments.resize(*(uint32_t*) pb);
	i = 0;
	pb += sizeof(uint32_t);
	for (i = 0; i < m_vecstrExperiments.size(); ++i) {
		string& strExperiment = m_vecstrExperiments[i];
		strExperiment.resize(*(uint32_t*) pb);
		pb += sizeof(uint32_t);
		for (j = 0; j < strExperiment.size(); ++j) {
			strExperiment[j] = *pb++;
		}
	}
	m_vecstrGenes.resize(*(uint32_t*) pb);
	pb += sizeof(uint32_t);
	for (i = 0; i < m_vecstrGenes.size(); ++i) {
		string& strGene = m_vecstrGenes[i];
		strGene.resize(*(uint32_t*) pb);
		pb += sizeof(uint32_t);
		for (j = 0; j < strGene.size(); ++j)
			strGene[j] = *pb++;
		m_mapstriGenes[strGene] = i;
	}
	return OpenMemmap(pb);

}

bool CPCLImpl::OpenMemmap(const unsigned char* pb) {
	size_t i;

	m_aadData = new float*[m_vecstrGenes.size()];
	m_aadData[0] = (float*) pb;
	for (i = 1; i < m_vecstrGenes.size(); ++i)
		m_aadData[i] = m_aadData[i - 1] + m_vecstrExperiments.size();
	m_Data.Initialize(m_vecstrGenes.size(), m_vecstrExperiments.size(), (float**) m_aadData);
	return true;
}

bool CPCLImpl::OpenExperiments(std::istream& istmInput, size_t iFeatures,
			       string& acStr, bool rTable) {
	const char* pc;
	string strToken;
	size_t iToken;
	
	Reset();
	if (!m_fHeader) {
		m_vecstrFeatures.resize(1 + iFeatures);
		return true;
	}
	getline(istmInput, acStr);
	if (!acStr.size()) {
		return false;
	}
	
	if(rTable){	  
	  m_vecstrFeatures.push_back("");
	}
	
	for (iToken = 0, pc = acStr.c_str(); (strToken = OpenToken(pc, &pc)).length()
			|| *pc; ++iToken)
	  if (iToken < iFeatures)
	    m_vecstrFeatures.push_back(strToken);
	  else if (iToken <= iFeatures && !rTable)
	    m_vecstrFeatures.push_back(strToken);
	  else
	    m_vecstrExperiments.push_back(strToken);
	
	return true;
}

bool CPCLImpl::OpenGene(std::istream& istmInput, std::vector<float>& vecdData,
		string& acStr) {
	const char* pc;
	char* pcEnd;
	string strToken;
	size_t iToken, iData, iBase, i;
	float d;
	map<string, size_t> mapstriValues;
	map<string, size_t>::iterator iterValue;
	const char* acLine = acStr.c_str();

	iBase = vecdData.size();
	getline(istmInput, acStr);
	if ((strToken = OpenToken(acLine)).empty())
		return false;
	if (strToken == "EWEIGHT")
		return true;
	if (!m_vecstrExperiments.empty())
		vecdData.resize(vecdData.size() + m_vecstrExperiments.size());
	for (iData = iToken = 0, pc = acLine; (strToken = OpenToken(pc, &pc)).length()
			|| *pc; ++iToken) {
		if (!iToken) {
			m_mapstriGenes[strToken] = m_vecstrGenes.size();
			m_vecstrGenes.push_back(strToken);
		} else if (iToken < m_vecstrFeatures.size())
			m_vecvecstrFeatures[iToken - 1].push_back(strToken);
		else if (!m_vecstrExperiments.empty() && (iData
				>= m_vecstrExperiments.size()))
			return false;
		else {
			d = CMeta::GetNaN();
			strToken = CMeta::Trim(strToken.c_str());
			if (strToken.length()) {
				d = (float) strtod(strToken.c_str(), &pcEnd);
				if (pcEnd != (strToken.c_str() + strToken.length())) {
					iterValue = mapstriValues.find(strToken);
					if (iterValue == mapstriValues.end()) {
						i = mapstriValues.size();
						mapstriValues[strToken] = i;
						d = i;
					} else
						d = iterValue->second;
				}
			}
			if (m_vecstrExperiments.empty())
				vecdData.push_back(d);
			else if ((i = (iBase + iData++)) >= vecdData.size())
				return false;
			else
				vecdData[i] = d;
		}
	}

	if (m_vecstrExperiments.empty())
		m_vecstrExperiments.resize(vecdData.size());
	else
		while (iData < m_vecstrExperiments.size())
			vecdData[iBase + iData++] = CMeta::GetNaN();

	return !!iToken;
}

/*!
 * \brief
 * Save the PCL's header row to the given text stream.
 *
 * \param ostm
 * Stream into which PCL header is saved.
 *
 * \param fCDT
 * If true, generate an initial CDT GENE column header.
 *
 * If called with fCDT false, saves standard PCL header and EWEIGHT rows to the given text stream using the
 * CPCL's current feature and experimental headers.  If fCDT is true, the output will instead be in CDT
 * format, with an extra initial column header for gene index identifiers.
 *
 * \see
 * Save
 */
void CPCL::SaveHeader(std::ostream& ostm, bool fCDT) const {
	size_t i;

	if (!m_fHeader)
		return;

/*	if (fCDT)
		ostm << c_szGID << '\t';
	ostm << m_vecstrFeatures[0];*/
	for (i = 1; i < m_vecstrFeatures.size(); ++i)
		ostm << '\t' << m_vecstrFeatures[i];
	for (i = 0; i < m_vecstrExperiments.size(); ++i)
		ostm << '\t' << m_vecstrExperiments[i];
	ostm << endl;

//	ostm << c_szEWEIGHT;
//	for (i = fCDT ? 0 : 1; i < m_vecstrFeatures.size(); ++i)
//		ostm << '\t';
//	for (i = 0; i < m_vecstrExperiments.size(); ++i)
//		ostm << '\t' << 1;
//	ostm << endl;
}

/*!
 * \brief
 * Save a single gene row to the given text stream.
 *
 * \param ostm
 * Stream into which gene row is saved.
 *
 * \param iGene
 * Gene index to be saved.
 *
 * \param iOriginal
 * If not equal to -1, generate an initial CDT GENE ID using the given original gene index.
 *
 * If called with iGene set to -1, saves a single row of a standard PCL to the given text stream using the
 * CPCL's current gene ID, features, and values for the row.  If a gene index is provided, the output will
 * instead be in CDT format, with an extra initial column of IDs of the form "GENE###", where the number
 * indicates the gene's original index in the pre-clustered file.
 *
 * \see
 * Save
 */
void CPCL::SaveGene(std::ostream& ostm, size_t iGene, size_t iOriginal) const {
	size_t i;
	float d;

	if (iOriginal != -1)
		ostm << c_szGENE << iOriginal << '\t';
	ostm << GetGene(iGene);
	for (i = 0; i < m_vecvecstrFeatures.size(); ++i)
		ostm << '\t' << m_vecvecstrFeatures[i][iGene];
	for (i = 0; i < m_vecstrExperiments.size(); ++i) {
		ostm << '\t';
		if (!CMeta::IsNaN(d = Get(iGene, i)))
			ostm << Get(iGene, i);
	}
	ostm << endl;
}

/*!
 * \brief
 * Save a PCL to the given text stream.
 *
 * \param ostm
 * Stream into which PCL file is saved.
 *
 * \param pveciGenes
 * If non-null, generate an initial CDT GENE column using the given original gene indices.
 *
 * If called without a vector of gene indices, saves a standard PCL to the given text stream using the CPCL's
 * current gene ID list, features, experiments, and data.  If a vector of gene indices is provided, the
 * output will instead be in CDT format, with an extra initial column of IDs of the form "GENE###", where the
 * number indicates the gene's original index in the pre-clustered file.
 *
 * \remarks
 * If pveciGenes is non-null, it must be of the same length as the PCL's gene list.
 *
 * \see
 * Open | CClustHierarchical
 */
//Save only certain genes
bool CPCL::Save(const char* szFile, const std::vector<size_t>* pveciGenes) const {
	ofstream ofsm;
	if (!strcmp(szFile + strlen(szFile) - 4, ".bin")) {
		ofsm.open(szFile, ios::binary);
		SaveBinary(ofsm);
	} else {
		ofsm.open(szFile);
		Save(ofsm, pveciGenes);
	}

}
void CPCL::Save(std::ostream& ostm, const std::vector<size_t>* pveciGenes) const {
	size_t i;

	SaveHeader(ostm, !!pveciGenes);

	for (i = 0; i < GetGenes(); ++i) {
		if (m_setiGenes.find(i) != m_setiGenes.end())
			continue;
		SaveGene(ostm, i, pveciGenes ? (*pveciGenes)[i] : -1);
	}
}

/*!
 * \brief
 * Save a PCL to the given binary stream.
 *
 * \param ostm
 * Stream into which PCL file is saved.
 *
 * \remarks
 * Most PCLs are saved as text files; binary files are more compact but very simple, encoding the size and
 * gene/feature/experiment values as simple character strings and saving data using the underlying
 * CFullMatrix::SaveBinary method.
 *
 * \see
 * SaveBinary
 */
void CPCL::SaveBinary(std::ostream& ostm) const {
	uint32_t iTmp;
	size_t i;

	iTmp = GetFeatures();
	ostm.write((const char*) &iTmp, sizeof(iTmp));
	for (i = 0; i < GetFeatures(); ++i)
		SaveString(ostm, GetFeature(i));
	iTmp = GetExperiments();
	ostm.write((const char*) &iTmp, sizeof(iTmp));
	for (i = 0; i < GetExperiments(); ++i)
		SaveString(ostm, GetExperiment(i));
	size_t iCount = 0;
	for (i = 0; i < GetGenes(); i++)
		if (!IsMasked(i))
			iCount++;
	iTmp = iCount;
	ostm.write((const char*) &iTmp, sizeof(iTmp));
	for (i = 0; i < GetGenes(); ++i)
		if (!IsMasked(i))
		SaveString(ostm, GetGene(i));

	for (i = 0; i < GetGenes(); ++i)
		if (!IsMasked(i))
		ostm.write((const char*) Get(i), GetExperiments() * sizeof(*Get(i)));
}

/*!
 * \brief
 * Load a PCL from the given binary stream.
 *
 * \param istm
 * Stream from which PCL file is loaded.
 *
 * \see
 * OpenBinary
 */
bool CPCL::OpenBinary(std::istream& istm) {
	uint32_t iTmp;
	size_t i;

	if (!istm.good())
		return false;
	

	Reset();
	istm.read((char*) &iTmp, sizeof(iTmp));
	m_vecstrFeatures.resize(iTmp);
	for (i = 0; i < m_vecstrFeatures.size(); ++i)
		OpenString(istm, m_vecstrFeatures[i]);

	istm.read((char*) &iTmp, sizeof(iTmp));
	m_vecstrExperiments.resize(iTmp);
	for (i = 0; i < m_vecstrExperiments.size(); ++i)
		OpenString(istm, m_vecstrExperiments[i]);

	istm.read((char*) &iTmp, sizeof(iTmp));
	m_vecstrGenes.resize(iTmp);
	for (i = 0; i < m_vecstrGenes.size(); ++i) {
		OpenString(istm, m_vecstrGenes[i]);
		m_mapstriGenes[m_vecstrGenes[i]] = i;
	}

	m_Data.Initialize(GetGenes(), GetExperiments());
	for (i = 0; i < m_Data.GetRows(); ++i)
		istm.read((char*) m_Data.Get(i), GetExperiments() * sizeof(*m_Data.Get(
				i)));

	return true;
}

/*!
 * \brief
 * Create a new PCL using the given genes, experiments, and features.
 *
 * \param vecstrGenes
 * Gene IDs to be used in the new PCL.
 *
 * \param vecstrExperiments
 * Experiment labels to be used in the new PCL.
 *
 * \param vecstrFeatures
 * Feature labels to be used in the new PCL (possibly empty).
 *
 * This creates an empty PCL (all entries missing) with the requested number of gene rows and experiment
 * columns; the given gene IDs and experiment labels are inserted into the appropriate header rows/columns.
 * The PCL will contain the given number of feature columns between the gene IDs and experimental values,
 * the values for which will be initialized to empty strings.
 */
void CPCL::Open(const std::vector<std::string>& vecstrGenes, const std::vector<
		std::string>& vecstrExperiments,
		const std::vector<std::string>& vecstrFeatures) {
	size_t i, j;

	Reset();
	if (vecstrFeatures.empty())
		m_vecstrFeatures.push_back("GID");
	else {
		m_vecstrFeatures.resize(vecstrFeatures.size());
		copy(vecstrFeatures.begin(), vecstrFeatures.end(),
				m_vecstrFeatures.begin());
		m_vecvecstrFeatures.resize(m_vecstrFeatures.size() - 1);
		for (i = 0; i < m_vecvecstrFeatures.size(); ++i)
			m_vecvecstrFeatures[i].resize(vecstrGenes.size());
	}

	m_vecstrGenes.resize(vecstrGenes.size());
	for (i = 0; i < GetGenes(); ++i)
		SetGene(i, vecstrGenes[i]);
	m_vecstrExperiments.resize(vecstrExperiments.size());
	for (i = 0; i < GetExperiments(); ++i)
		SetExperiment(i, vecstrExperiments[i]);

	m_Data.Initialize(GetGenes(), GetExperiments());
	for (i = 0; i < m_Data.GetRows(); ++i)
		for (j = 0; j < m_Data.GetColumns(); ++j)
			m_Data.Set(i, j, CMeta::GetNaN());
}

/*!
 * \brief
 * Create a new PCL using the given genes, experiments, and gene order.
 *
 * \param veciGenes
 * Order in which genes will be placed in the new PCL.
 *
 * \param vecstrGenes
 * Gene IDs to be used in the new PCL.
 *
 * \param vecstrExperiments
 * Experiment labels to be used in the new PCL.
 *
 * This creates an empty PCL (all entries missing) with the requested number of gene rows and experiment
 * columns; the given gene IDs and experiment labels are inserted into the appropriate header rows/columns.
 * However, the gene order will be the order of indices given in veciGenes.  For example, if veciGenes is
 * [1, 0, 2] and vecstrGenes contains [A, B, C], the order of genes within the new PCL will be [B, A, C].
 * Experiment order is unaffected and will be as given in vecstrExperiments.
 *
 * \remarks
 * veciGenes and vecstrGenes must be of the same length.
 *
 * \see
 * SortGenes
 */
void CPCL::Open(const std::vector<size_t>& veciGenes, const std::vector<
		std::string>& vecstrGenes,
		const std::vector<std::string>& vecstrExperiments) {
	size_t i, j;
	char ac[16];

	Reset();
	m_vecstrFeatures.resize(4);
	m_vecstrFeatures[0] = "GID";
	m_vecstrFeatures[1] = "YORF";
	m_vecstrFeatures[2] = "NAME";
	m_vecstrFeatures[3] = "GWEIGHT";
	m_vecvecstrFeatures.resize(m_vecstrFeatures.size() - 1);
	for (i = 0; i < m_vecvecstrFeatures.size(); ++i)
		m_vecvecstrFeatures[i].resize(veciGenes.size());
	for (i = 0; i < veciGenes.size(); ++i) {
		m_vecvecstrFeatures[0][i] = m_vecvecstrFeatures[1][i]
				= vecstrGenes[veciGenes[i]];
		m_vecvecstrFeatures[2][i] = "1";
	}

	m_vecstrGenes.resize(vecstrGenes.size());
	for (i = 0; i < m_vecstrGenes.size(); ++i) {
		m_vecstrGenes[i] = "GENE";
		sprintf_s(ac, "%d", veciGenes[i]);
		m_vecstrGenes[i] += ac;
		m_mapstriGenes[m_vecstrGenes[i]] = i;
	}
	m_vecstrExperiments.resize(vecstrExperiments.size());
	for (i = 0; i < m_vecstrExperiments.size(); ++i)
		m_vecstrExperiments[i] = vecstrExperiments[i];

	m_Data.Initialize(GetGenes(), GetExperiments());
	for (i = 0; i < m_Data.GetRows(); ++i)
		for (j = 0; j < m_Data.GetColumns(); ++j)
			m_Data.Set(i, j, CMeta::GetNaN());
}

/*!
 * \brief
 * Reorder the PCL's genes based on the order of the given indices.
 *
 * \param veciOrder
 * Index order in which genes should be placed.
 *
 * \returns
 * True if genes were reordered successfully.
 *
 * Reorders the PCL's gene rows based on the given indices.  For example, if the current row order is
 * [A, B, C] and the given indices are [1, 0, 2], the new row order will be [B, A, C].
 */
bool CPCL::SortGenes(const vector<size_t>& veciOrder) {
	size_t i;

	if (veciOrder.size() != m_Data.GetRows())
		return false;

	CMeta::Permute(m_Data.Get(), veciOrder);
	CMeta::Permute(m_vecstrGenes, veciOrder);
	for (i = 0; i < GetGenes(); ++i)
		m_mapstriGenes[GetGene(i)] = i;
	for (i = 0; i < m_vecvecstrFeatures.size(); ++i)
		CMeta::Permute(m_vecvecstrFeatures[i], veciOrder);

	return true;
}

/*!
 * \brief
 * Rank transform each row of the PCL in increasing order.
 *
 * Replaces all values in the PCL with their increasing ranks by row.  For example, a row containing
 * [-0.1, 0.1, 0.3] would be replaced with [0, 1, 2]; a row containing [0.5, 0.4, 0.3] would be replaced
 * with [2, 1, 0].
 *
 * \see
 * IMeasure::IsRank
 */
void CPCL::RankTransform() {
	size_t i, j, k;
	vector<size_t> veciRanks, veciCounts;

	veciRanks.resize(GetExperiments());
	veciCounts.resize(GetExperiments());
	for (i = 0; i < GetGenes(); ++i) {
		fill(veciRanks.begin(), veciRanks.end(), 0);
		for (j = 0; j < GetExperiments(); ++j) {
			if (CMeta::IsNaN(Get(i, j)))
				continue;
			for (k = 0; k < GetExperiments(); ++k) {
				if (CMeta::IsNaN(Get(i, k)))
					continue;
				if ((j != k) && (Get(i, k) < Get(i, j)))
					veciRanks[j]++;
			}
		}

		fill(veciCounts.begin(), veciCounts.end(), 0);
		for (j = 0; j < GetExperiments(); ++j)
			if (!CMeta::IsNaN(Get(i, j)))
				veciCounts[veciRanks[j]]++;

		for (j = 0; j < GetExperiments(); ++j)
			if (!CMeta::IsNaN(Get(i, j))) {
				k = veciRanks[j];
				// Closed form for sum(rank, rank + n) / n
				Set(i, j, k + ((veciCounts[k] + 1) / 2.0f));
			}
	}
}

/*!
 * \brief
 * Appends new, empty gene rows to the end of the PCL using the given gene IDs.
 *
 * \param vecstrGenes
 * Gene names to be appended to the PCL.
 *
 * \returns
 * True if the gene rows were appended successfully.
 *
 * \remarks
 * New rows will initially have empty strings for all features and missing values for all data.
 */
bool CPCL::AddGenes(const std::vector<std::string>& vecstrGenes) {
	size_t i, j, iStart;

	iStart = m_Data.GetRows();
	if (!m_Data.AddRows(vecstrGenes.size()))
		return false;
	for (i = iStart; i < m_Data.GetRows(); ++i)
		for (j = 0; j < m_Data.GetColumns(); ++j)
			m_Data.Set(i, j, CMeta::GetNaN());

	m_vecstrGenes.resize(m_vecstrGenes.size() + vecstrGenes.size());
	for (i = 0; i < vecstrGenes.size(); ++i)
		SetGene(iStart + i, vecstrGenes[i]);
	for (i = 0; i < m_vecvecstrFeatures.size(); ++i) {
		m_vecvecstrFeatures[i].resize(m_vecvecstrFeatures[i].size()
				+ vecstrGenes.size());
		if (m_vecstrFeatures[i + 1] == c_szNAME)
			for (j = 0; j < vecstrGenes.size(); ++j)
				m_vecvecstrFeatures[i][iStart + j] = vecstrGenes[j];
		else if (m_vecstrFeatures[i + 1] == c_szGWEIGHT)
			for (j = 0; j < vecstrGenes.size(); ++j)
				m_vecvecstrFeatures[i][iStart + j] = c_szOne;
	}

	return true;
}

/*!
 * \brief
 * Normalizes the PCL's values in the requested manner.
 *
 * \param eNormalize
 * Algorithm by which the PCL's values should be normalized.
 *
 * \see
 * ENormalize
 */
void CPCL::Normalize(ENormalize eNormalize) {
	size_t i, j, iCount;
	double dAve, dStd;
	float d, dMin, dMax;

	switch (eNormalize) {
	case ENormalizeMean:
		dMin = FLT_MAX;
		dAve = 0;
		for (iCount = i = 0; i < GetGenes(); ++i)
			for (j = 0; j < GetExperiments(); ++j)
				if (!CMeta::IsNaN(d = Get(i, j))) {
					if (d < dMin)
						dMin = d;
					iCount++;
					dAve += d;
				}
		if (!iCount)
			break;
		dAve = (dAve / iCount) - dMin;
		for (i = 0; i < GetGenes(); ++i)
			for (j = 0; j < GetExperiments(); ++j)
				if (!CMeta::IsNaN(d = Get(i, j)))
					Set(i, j, (float) ((d - dMin) / dAve));
		break;

	case ENormalizeZScore:
		dAve = dStd = 0;
		for (iCount = i = 0; i < GetGenes(); ++i)
			for (j = 0; j < GetExperiments(); ++j)
				if (!CMeta::IsNaN(d = Get(i, j))) {
					iCount++;
					dAve += d;
					dStd += d * d;
				}
		dAve /= iCount;
		dStd = sqrt((dStd / iCount) - (dAve * dAve));
		for (i = 0; i < GetGenes(); ++i)
			for (j = 0; j < GetExperiments(); ++j)
				if (!CMeta::IsNaN(d = Get(i, j)))
					Set(i, j, (float) ((d - dAve) / dStd));
		break;

	case ENormalizeRow:
		for (i = 0; i < GetGenes(); ++i) {
			dAve = CStatistics::Average(Get(i), Get(i) + GetExperiments());
			dStd = sqrt(CStatistics::Variance(Get(i),
					Get(i) + GetExperiments(), dAve));
			if (dStd)
				for (j = 0; j < GetExperiments(); ++j)
					Set(i, j, (float) ((Get(i, j) - dAve) / dStd));
		}
		break;

	case ENormalizeColumn:
	case ENormalizeColumnCenter:
	case ENormalizeColumnFraction:
		for (i = 0; i < GetExperiments(); ++i) {
			dAve = dStd = 0;
			for (iCount = j = 0; j < GetGenes(); ++j)
				if (!CMeta::IsNaN(d = Get(j, i))) {
					iCount++;
					dAve += d;
					dStd += d * d;
				}
			if (iCount) {
				if (eNormalize != ENormalizeColumnFraction) {
					dAve /= iCount;
					dStd = (dStd / iCount) - (dAve * dAve);
					dStd = (dStd <= 0) ? 1 : sqrt(dStd);
					if (eNormalize == ENormalizeColumnCenter)
						dStd = 1;
				}
				for (j = 0; j < GetGenes(); ++j)
					if (!CMeta::IsNaN(d = Get(j, i))) {
						if (eNormalize == ENormalizeColumnFraction)
							d /= (float) dAve;
						else
							d = (float) ((d - dAve) / dStd);
						Set(j, i, d);
					}
			}
		}
		break;

	case ENormalizeMinMax:
		dMin = FLT_MAX;
		dMax = -dMin;
		for (i = 0; i < GetGenes(); ++i)
			for (j = 0; j < GetExperiments(); ++j)
				if (!CMeta::IsNaN(d = Get(i, j))) {
					if (d < dMin)
						dMin = d;
					if (d > dMax)
						dMax = d;
				}
		if (!(dMax -= dMin))
			dMax = 1;
		for (i = 0; i < GetGenes(); ++i)
			for (j = 0; j < GetExperiments(); ++j)
				if (!CMeta::IsNaN(d = Get(i, j)))
					Set(i, j, (d - dMin) / dMax);
		break;
	case EMeanSubtractColumn:
		for (i = 0; i < GetExperiments(); ++i) {
			dAve = 0;
			for (iCount = j = 0; j < GetGenes(); ++j)
				if (!CMeta::IsNaN(d = Get(j, i))) {
					iCount++;
					dAve += d;

				}

			if (iCount) {
				dAve /= iCount;
				for (j = 0; j < GetGenes(); ++j)
					if (!CMeta::IsNaN(d = Get(j, i)))
						Set(j, i, (float) (d - dAve));
			}
		}
		break;
	}
}

/*!
 * \brief
 * Impute missing values using the knnimpute algorithm; optionally mask genes with too many missing values.
 *
 * \param iNeighbors
 * Number of nearest neighbors to use for imputation.
 *
 * \param dMinimumPresent
 * Fraction of conditions that must be present; genes with fewer than this many values present will be masked
 * instead of imputed.
 *
 * \param DatSimilarity
 * CDat similarity matrix from which nearest neighbors are determined.
 *
 * Imputes missing values in the PCL using the knnimpute algorithm of Troyanskaya et al, Bioinformatics 2001.
 * Briefly, for each gene with missing values, the k nearest neighbors (most similar genes) are found, and
 * the missing value is replaced with a weighted average of the values in the neighbors.  Optionally,
 * Impute can also completely mask genes that are missing too many values; genes with less than the
 * requested fraction of data present will be completely masked rather than imputed.
 *
 * \remarks
 * iNeighbors must be greater than zero, dMinimumPresent must be between zero and one inclusive, and
 * DatSimilarity must be of exactly the same size as the current PCL.
 */
void CPCL::Impute(size_t iNeighbors, float dMinimumPresent,
		const CDat& DatSimilarity) {
	vector<float> vecdMissing;
	vector<size_t> veciMissing;
	vector<SNeighbors> vecsNeighbors;
	size_t i, j, k, iCol, iMissing, iOne, iTwo;
	float d, dSum;
	const float* ad;

	vecdMissing.resize(GetGenes());
	for (i = 0; i < vecdMissing.size(); ++i) {
		for (iMissing = j = 0; j < GetExperiments(); ++j)
			if (CMeta::IsNaN(Get(i, j)))
				iMissing++;
		if ((vecdMissing[i] = (float) (GetExperiments() - iMissing)
				/ GetExperiments()) >= dMinimumPresent)
			veciMissing.push_back(i);
	}
	if (iNeighbors) {
		vecsNeighbors.resize(veciMissing.size());
		for (i = 0; i < vecsNeighbors.size(); ++i)
			vecsNeighbors[i].Initialize(Get(veciMissing[i]), GetExperiments(),
					iNeighbors);
		for (i = 0; i < veciMissing.size(); ++i) {
			if (!(i % 100))
				g_CatSleipnir().info(
						"CPCL::Impute( %d, %g ) finding neighbors for gene %d/%d",
						iNeighbors, dMinimumPresent, i, veciMissing.size());
			ad = Get(iOne = veciMissing[i]);
			for (j = (i + 1); j < veciMissing.size(); ++j) {
				d = DatSimilarity.Get(iOne, iTwo = veciMissing[j]);
				vecsNeighbors[i].Add(iTwo, d, Get(iTwo));
				vecsNeighbors[j].Add(iOne, d, ad);
			}
		}
	}

	for (iOne = i = 0; i < GetGenes(); ++i) {
		if (vecdMissing[i] < dMinimumPresent) {
			MaskGene(i);
			continue;
		}
		if (vecsNeighbors.empty())
			continue;
		{
			const SNeighbors& sGene = vecsNeighbors[iOne++];

			for (j = 0; j < GetExperiments(); ++j) {
				if (!CMeta::IsNaN(Get(i, j)))
					continue;

				iCol = sGene.GetColumn(j);
				for (d = dSum = 0, iMissing = k = 0; k < iNeighbors; ++k) {
					const pair<size_t, float>& prNeighbor =
							sGene.m_MatNeighbors.Get(k, iCol);

					if (prNeighbor.second == -FLT_MAX)
						continue;
					iMissing++;
					dSum += prNeighbor.second;
					d += Get(prNeighbor.first, j) * prNeighbor.second;
				}
				if (dSum)
					d /= dSum;
				Set(i, j, iMissing ? d : CMeta::GetNaN());
			}
		}
	}
}

/*!
 * \brief
 * Impute missing values using the knnimpute algorithm; optionally mask genes with too many missing values.
 *
 * \param iNeighbors
 * Number of nearest neighbors to use for imputation.
 *
 * \param dMinimumPresent
 * Fraction of conditions that must be present; genes with fewer than this many values present will be masked
 * instead of imputed.
 *
 * \param pMeasure
 * Similarity measure with which nearest neighbors are determined.
 *
 * \param fPrecompute
 * If true, precompute and cache in memory all pairwise similarities; otherwise, calculate these on the fly
 * from values in the PCL.
 *
 * Imputes missing values in the PCL using the knnimpute algorithm of Troyanskaya et al, Bioinformatics 2001.
 * Briefly, for each gene with missing values, the k nearest neighbors (most similar genes) are found, and
 * the missing value is replaced with a weighted average of the values in the neighbors.  Optionally,
 * Impute can also completely mask genes that are missing too many values; genes with less than the
 * requested fraction of data present will be completely masked rather than imputed.
 *
 * \remarks
 * iNeighbors must be greater than zero, dMinimumPresent must be between zero and one inclusive.  If
 * precomputation is requested, imputation will be faster, but memory proportional to the number of genes
 * squared will be required.  If it is not, imputation will be slower, but no additional memory is required.
 * As a rule of thumb, precomputation is desirable for anything less than ~25,000 genes (which will consume
 * roughly 1GB of memory; if you have RAM to spare, feel free to precompute at will).
 */
void CPCL::Impute(size_t iNeighbors, float dMinimumPresent,
		const IMeasure* pMeasure, bool fPrecompute) {
	CDat Dat;
	size_t i, j;

	if (fPrecompute) {
		Dat.Open(GetGeneNames());
		for (i = 0; i < Dat.GetGenes(); ++i)
			for (j = (i + 1); j < Dat.GetGenes(); ++j)
				Dat.Set(i, j, (float) pMeasure->Measure(Get(i),
						GetExperiments(), Get(j), GetExperiments()));
	} else
		Dat.Open(*this, pMeasure, false);
	Impute(iNeighbors, dMinimumPresent, Dat);
}

/*!
 * \brief
 * Resolves probesets by averaging agreeing probes and discarding conflicting ones.
 *
 * \param iSample
 * Number of random samples to generate for probeset null distribution.
 *
 * \param iBins
 * Number of bins for probeset distribution histogram.
 *
 * \param dBinSize
 * Size of bins for probeset distribution histogram in z-score units.
 *
 * \remarks
 * Performance can vary wildly depending on characteristics of the input PCL's probesets; use with caution.
 */
void CPCL::MedianMultiples(size_t iSample, size_t iBins, float dBinSize) {
	size_t i, j, k, iBin, iOne, iCutoff;
	CMeasureEuclidean Euclidean;
	vector<float> vecdEuclidean, vecdMapped, vecdBG, vecdFG, vecdMean;
	float d, dAve, dStd, dSample;
	vector<vector<size_t> > vecveciGenes;
	vector<size_t> veciMean;
	vector<string> vecstrAdd;

	{
		map<string, size_t> mapstriClasses;
		map<string, size_t>::iterator iterClass;

		for (i = 0; i < m_vecstrGenes.size(); ++i) {
			const string& strGene = m_vecstrGenes[i];

			if ((iterClass = mapstriClasses.find(strGene))
					== mapstriClasses.end()) {
				mapstriClasses[strGene] = vecveciGenes.size();
				vecveciGenes.push_back(vector<size_t> ());
				vecveciGenes.back().push_back(i);
			} else
				vecveciGenes[iterClass->second].push_back(i);
		}
	}

	dSample = 2.0f * iSample / GetGenes() / (GetGenes() - 1);
	dAve = dStd = 0;
	for (k = i = 0; i < GetGenes(); ++i)
		for (j = (i + 1); j < GetGenes(); ++j) {
			if (((float) rand() / RAND_MAX) >= dSample)
				continue;
			vecdEuclidean.push_back(d = (float) Euclidean.Measure(Get(i),
					GetExperiments(), Get(j), GetExperiments(),
					IMeasure::EMapNone));
			k++;
			dAve += d;
			dStd += d * d;
		}
	dAve /= k;
	dStd = sqrt((dStd / (k - 1)) - (dAve * dAve));
	if (!(iBins % 2))
		iBins++;

	vecdBG.resize(iBins);
	for (i = 0; i < vecdEuclidean.size(); ++i)
		vecdBG[MedianMultiplesBin(vecdEuclidean[i], dAve, dStd, iBins, dBinSize)]++;
	for (i = 0; i < vecdBG.size(); ++i)
		//		vecdBG[i] /= vecdEuclidean.size( );
		vecdBG[i] = (vecdBG[i] + 1) / (vecdEuclidean.size() + vecdBG.size());
	//*
	MedianMultiplesSmooth(2, vecdBG);
	//*/

	MedianMultiplesMapped(vecveciGenes, vecdMapped);
	vecdFG.resize(iBins);
	for (i = 0; i < vecdMapped.size(); ++i)
		vecdFG[MedianMultiplesBin(vecdMapped[i], dAve, dStd, iBins, dBinSize)]++;
	/*
	 dMax = 0;
	 for( iMax = i = 0; i < vecdBG.size( ); ++i )
	 if( vecdBG[i] > dMax ) {
	 dMax = vecdBG[i];
	 iMax = i; }
	 for( i = 0; i <= iMax; ++i )
	 vecdFG[i]++;
	 for( i = 0; i < vecdFG.size( ); ++i )
	 vecdFG[i] /= vecdMapped.size( ) + iMax + 1;
	 //*/
	//*/
	for (i = 0; i < vecdFG.size(); ++i)
		//		vecdFG[i] /= vecdMapped.size( );
		vecdFG[i] = (vecdFG[i] + 1) / (vecdMapped.size() + vecdFG.size());
	MedianMultiplesSmooth(2, vecdFG);
	/*
	 for( i = 0; i < vecdBG.size( ); ++i )
	 cerr << vecdBG[i] << '\t';
	 cerr << endl;
	 for( i = 0; i < vecdFG.size( ); ++i )
	 cerr << vecdFG[i] << '\t';
	 cerr << endl;
	 //*/
	for (iCutoff = 0; iCutoff < vecdFG.size(); ++iCutoff)
		if (vecdFG[iCutoff] < vecdBG[iCutoff])
			break;

	vecstrAdd.resize(1);
	veciMean.resize(GetExperiments());
	vecdMean.resize(GetExperiments());
	for (i = 0; i < vecveciGenes.size(); ++i) {
		const vector<size_t>& veciGenes = vecveciGenes[i];
		vector<size_t> veciAgree;

		if (veciGenes.size() < 2)
			continue;
		veciAgree.resize(veciGenes.size());
		for (j = 0; j < veciGenes.size(); ++j) {
			iOne = veciGenes[j];
			for (k = (j + 1); k < veciGenes.size(); ++k) {
				iBin = MedianMultiplesBin((float) Euclidean.Measure(Get(iOne),
						GetExperiments(), Get(veciGenes[k]), GetExperiments(),
						IMeasure::EMapNone), dAve, dStd, iBins, dBinSize);
				//				if( iBin <= iCutoff ) {
				if ((1.01 * vecdFG[iBin]) >= vecdBG[iBin]) {
					veciAgree[j]++;
					veciAgree[k]++;
				}
			}
		}

		fill(vecdMean.begin(), vecdMean.end(), 0.0f);
		fill(veciMean.begin(), veciMean.end(), 0);
		for (iBin = j = 0; j < veciGenes.size(); ++j)
			if ((2 * veciAgree[j]) >= (veciAgree.size() - 1)) {
				iBin++;
				iOne = veciGenes[j];
				for (k = 0; k < GetExperiments(); ++k)
					if (!CMeta::IsNaN(d = Get(iOne, k))) {
						veciMean[k]++;
						vecdMean[k] += d;
					}
			}

		for (j = 0; j < veciGenes.size(); ++j)
			MaskGene(veciGenes[j]);
		if ((iBin * 2) > veciAgree.size()) {
			for (j = 0; j < vecdMean.size(); ++j)
				vecdMean[j] = veciMean[j] ? (vecdMean[j] / veciMean[j])
						: CMeta::GetNaN();
			vecstrAdd[0] = GetGene(veciGenes[0]);
			AddGenes(vecstrAdd);
			Set(GetGenes() - 1, &vecdMean[0]);
		}
	}
}

void CPCLImpl::MedianMultiplesMapped(
		const vector<vector<size_t> >& vecveciGenes, vector<float>& vecdMapped) {
	size_t i, j, k;
	CMeasureEuclidean Euclidean;
	const float* adOne;

	for (i = 0; i < vecveciGenes.size(); ++i) {
		const vector<size_t>& veciGenes = vecveciGenes[i];

		if (veciGenes.size() < 2)
			continue;
		for (j = 0; j < veciGenes.size(); ++j) {
			adOne = m_Data.Get(veciGenes[j]);
			for (k = (j + 1); k < veciGenes.size(); ++k)
				vecdMapped.push_back((float) Euclidean.Measure(adOne,
						m_Data.GetColumns(), m_Data.Get(veciGenes[k]),
						m_Data.GetColumns(), IMeasure::EMapNone));
		}
	}
}

/**
 * Fills a PCL object from a tab text file (gene1 gene2 val)
 * gene1 is row, gene2 is col
 * used for creating a directed gene-gene network
 */
  bool CPCL::populate(const char* szFile, float dDefault){
  ifstream istm;
  const char*	pc;
  char*		pcTail;
  char*		acBuf;
  string		strToken, strCache, strValue;
  size_t		iOne, iTwo, i;
  float		dScore;
  
  istm.open( szFile );
  
  acBuf = new char[ c_iBufferSize ];
  while( istm.peek( ) != EOF ) {
    istm.getline( acBuf, c_iBufferSize - 1 );
    strToken = OpenToken( acBuf, &pc );
    if( !strToken.length( ) )
      break;
    //if( strToken == c_acComment )
    //  continue;
    if( strToken != strCache ) {
      strCache = strToken;			
      iOne = GetGene( strToken );
    }
    
    strToken = OpenToken( pc, &pc );
    if( !strToken.length( ) ) {
      delete[] acBuf;
      return false; }
    iTwo = GetGene( strToken );
    
    strValue = OpenToken( pc );
    if( !strValue.length( ) ) {
      if( CMeta::IsNaN( dScore = dDefault ) ) {
	delete[] acBuf;
	return false; } }
    else if( !( dScore = (float)strtod( strValue.c_str( ), &pcTail ) ) &&
	     ( pcTail != ( strValue.c_str( ) + strValue.length( ) ) ) ) {
      delete[] acBuf;
      return false; }
    
    Set(iOne, iTwo, dScore);
    
  }
  delete[] acBuf;
  
  return true;
}
  
}

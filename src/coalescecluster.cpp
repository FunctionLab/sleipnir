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
#include "coalescecluster.h"
#include "coalescei.h"
#include "coalescemotifs.h"
#include "pcl.h"
#include "halfmatrix.h"
#include "statistics.h"
#include "clusthierarchical.h"
#include "pst.h"

namespace Sleipnir {

    const char    CCoalesceClusterImpl::c_szMotifs[] = "_motifs.txt";
    const char    CCoalesceClusterImpl::c_szGenes[] = "Genes";
    const char    CCoalesceClusterImpl::c_szConditions[] = "Conditions";

/*!
 * \brief
 * Randomly initializes a new cluster from the given PCL by selecting a pair of correlated genes and
 * a surrounding seed of additional genes.
 *
 * \param PCL
 * Expression data from which cluster is seeded.
 *
 * \param Pot
 * Cluster initialized to the inverse of the current cluster, i.e. all genes not in the new cluster.
 *
 * \param vecsDatasets
 * Vector of dataset block structure to be used by the new cluster for subsequent selection of
 * significant conditions.
 *
 * \param setpriiSeeds
 * Set of previously failed seed pairs to be excluded for cluster initialization.
 *
 * \param vecdSeed
 * If non-empty, expression vector to be used as a cluster seed.
 *
 * \param iPairs
 * Maximum number of gene pairs to be sampled for seed pair discovery.
 *
 * \param dPValue
 * P-value threshhold for significant correlation.
 *
 * \param dProbability
 * Prior probability of gene inclusion.
 *
 * \param iThreads
 * Maximum number of simultaneous threads for gene pair sampling.
 *
 * \returns
 * True if the cluster was successfully initialized; false if no appropriate seed could be found.
 *
 * \throws <exception class>
 * Description of criteria for throwing this exception.
 *
 * A new cluster is initialized by selecting the most highly correlated gene pair from a random sample
 * of the input PCL that is below the significance threshhold.  The expression centroid of this pair is
 * then calculated, and all other genes significantly correlated with this centroid are subsequently
 * added.  Genes not selected for inclusion are added to the inverse cluster instead.  All significance
 * tests are appropriately Bonferroni corrected for multiple hypothesis testing.
 *
 * \remarks
 * Initial seed pair is the most highly correlated and significant gene pair sampled from the input PCL;
 * other genes significantly correlated with the resulting centroid are added subsequently to initialize
 * the cluster.  The provided dataset blocks are used for all subsequent condition significance tests
 * and covariance calculations in which the cluster is involved.
 */
    bool CCoalesceCluster::Initialize(const CPCL &PCL, CCoalesceCluster &Pot,
                                      const std::vector <SCoalesceDataset> &vecsDatasets,
                                      std::set <std::pair<size_t, size_t>> &setpriiSeeds,
                                      const vector<float> &vecdSeed, size_t iPairs, float dPValue, float dProbability,
                                      size_t iThreads) {
        size_t i;
        float dFraction;
        CPCL PCLCopy;

        m_setiDatasets.clear();
        m_setiGenes.clear();
        m_setsMotifs.clear();
        m_vecsDatasets.resize(vecsDatasets.size());
        Pot.m_vecsDatasets.resize(vecsDatasets.size());
        for (i = 0; i < m_vecsDatasets.size(); ++i) {
            m_vecsDatasets[i].m_psDataset = Pot.m_vecsDatasets[i].m_psDataset = &vecsDatasets[i];
            m_setiDatasets.insert(i);
        }
        m_vecdPriors.resize(PCL.GetGenes());
        fill(m_vecdPriors.begin(), m_vecdPriors.end(), 1 - dProbability);
        dFraction = min(1.0f, 2.0f * iPairs / PCL.GetGenes() / (PCL.GetGenes() - 1));

        PCLCopy.Open(PCL);
        PCLCopy.Normalize(CPCL::ENormalizeColumn);
        return ((vecdSeed.size() || AddSeedPair(PCLCopy, Pot, setpriiSeeds, dFraction, dPValue, iThreads)) &&
                AddCorrelatedGenes(PCLCopy, Pot, vecdSeed, dPValue));
    }

    void CCoalesceClusterImpl::Add(size_t iGene, CCoalesceCluster &Pot) {

        m_setiGenes.insert(iGene);
        Pot.m_setiGenes.erase(iGene);
    }

    void *CCoalesceClusterImpl::ThreadSeedPair(void *pData) {
        SThreadSeedPair *psData;
        double dR, dP;
        size_t i, j, iN, iExperiments, iPair, iPairs, iCur;
        pair <size_t, size_t> priiSeed;

        psData = (SThreadSeedPair *) pData;
        iExperiments = psData->m_pPCL->GetExperiments();
        psData->m_dMaxCorr = -(psData->m_dMinP = DBL_MAX);
        iPairs = psData->m_pPCL->GetGenes() * (psData->m_pPCL->GetGenes() - 1) / 2;
        for (psData->m_iOne = psData->m_iTwo = 0, iPair = psData->m_iOffset; iPair < iPairs;
             iPair += psData->m_iStep) {
            if (((float) rand() / RAND_MAX) > psData->m_dFraction)
                continue;
            for (i = 0, iCur = iPair; iCur >= (psData->m_pPCL->GetGenes() - 1 - i);
                 iCur -= psData->m_pPCL->GetGenes() - 1 - i++);
            j = i + iCur + 1;
            if (((dR = CMeasurePearson::Pearson(psData->m_pPCL->Get(i), iExperiments,
                                                psData->m_pPCL->Get(j), iExperiments, IMeasure::EMapNone, NULL, NULL,
                                                &iN)) < 0) ||
                CMeta::IsNaN(dR))
                continue;
            if (((dP = CStatistics::PValuePearson(dR, iN)) < psData->m_dMinP) ||
                (!dP && (dR > psData->m_dMaxCorr))) {
                priiSeed.first = i;
                priiSeed.second = j;
                if (psData->m_psetpriiSeeds->find(priiSeed) == psData->m_psetpriiSeeds->end()) {
                    psData->m_dMaxCorr = dR;
                    psData->m_dMinP = dP;
                    psData->m_iOne = i;
                    psData->m_iTwo = j;
                }
            }
        }

        return NULL;
    }

    bool CCoalesceClusterImpl::AddSeedPair(const CPCL &PCL, CCoalesceCluster &Pot,
                                           std::set <std::pair<size_t, size_t>> &setpriiSeeds, float dFraction,
                                           float dPValue, size_t iThreads) {
        size_t i, iOne, iTwo;
        double dMaxCorr, dMinP;
        pair <size_t, size_t> priiSeed;
        vector <pthread_t> vecpthdThreads;
        vector <SThreadSeedPair> vecsThreads;

        if (PCL.GetGenes() < 2) {
            g_CatSleipnir().error("CCoalesceClusterImpl::AddSeedPair( %g ) found no genes", dPValue);
            return false;
        }
        vecpthdThreads.resize(iThreads);
        vecsThreads.resize(vecpthdThreads.size());
        for (i = 0; i < vecsThreads.size(); ++i) {
            vecsThreads[i].m_iOffset = i;
            vecsThreads[i].m_iStep = vecsThreads.size();
            vecsThreads[i].m_pPCL = &PCL;
            vecsThreads[i].m_dFraction = dFraction;
            vecsThreads[i].m_psetpriiSeeds = &setpriiSeeds;
            if (pthread_create(&vecpthdThreads[i], NULL, ThreadSeedPair, &vecsThreads[i])) {
                g_CatSleipnir().error("CCoalesceClusterImpl::AddSeedPair( %g, %g, %d ) could not seed pair",
                                      dFraction, dPValue, iThreads);
                return false;
            }
        }
        dMaxCorr = -(dMinP = DBL_MAX);
        for (iOne = iTwo = i = 0; i < vecpthdThreads.size(); ++i) {
            pthread_join(vecpthdThreads[i], NULL);
            if ((vecsThreads[i].m_dMinP < dMinP) ||
                (!vecsThreads[i].m_dMinP && (vecsThreads[i].m_dMaxCorr > dMaxCorr))) {
                dMinP = vecsThreads[i].m_dMinP;
                dMaxCorr = vecsThreads[i].m_dMaxCorr;
                iOne = vecsThreads[i].m_iOne;
                iTwo = vecsThreads[i].m_iTwo;
            }
        }
        if ((dMinP * PCL.GetGenes() * (PCL.GetGenes() - 1) * dFraction * dFraction) <= (dPValue * 2)) {
            g_CatSleipnir().info("CCoalesceClusterImpl::AddSeedPair( %g, %g ) seeding: %s, %s, %g (p=%g)",
                                 dFraction, dPValue, PCL.GetGene(iOne).c_str(), PCL.GetGene(iTwo).c_str(), dMaxCorr,
                                 dMinP);
            priiSeed.first = iOne;
            priiSeed.second = iTwo;
            setpriiSeeds.insert(priiSeed);
            Add(iOne, Pot);
            Add(iTwo, Pot);
            return true;
        }

        g_CatSleipnir().notice("CCoalesceClusterImpl::AddSeedPair( %g, %g ) inadequate seed pair: %s, %s, %g (p=%g)",
                               dFraction, dPValue, PCL.GetGene(iOne).c_str(), PCL.GetGene(iTwo).c_str(), dMaxCorr,
                               dMinP);
        return false;
    }

    bool CCoalesceClusterImpl::AddCorrelatedGenes(const CPCL &PCL, CCoalesceCluster &Pot, const vector<float> &vecdSeed,
                                                  float dPValue) {
        size_t iGene, iN;
        double dR;

        CalculateCentroid(PCL);
        if (vecdSeed.size()) {
            if (vecdSeed.size() == m_vecdCentroid.size())
                copy(vecdSeed.begin(), vecdSeed.end(), m_vecdCentroid.begin());
            else
                g_CatSleipnir().error("CCoalesceClusterImpl::AddCorrelatedGenes( %g ) invalid seed provided, ignoring",
                                      dPValue);
        }
        for (iGene = 0; iGene < PCL.GetGenes(); ++iGene)
            if (!IsGene(iGene) &&
                ((dR = CMeasurePearson::Pearson(&m_vecdCentroid.front(), PCL.GetExperiments(),
                                                PCL.Get(iGene), PCL.GetExperiments(), IMeasure::EMapNone, NULL, NULL,
                                                &iN)) > 0) &&
                ((CStatistics::PValuePearson(dR, iN) * PCL.GetGenes()) <= dPValue))
                Add(iGene, Pot);

        return true;
    }

/*!
 * \brief
 * Updates the given motif score histograms based on the genes in the current cluster using the
 * provided per-gene motif scores.
 *
 * \param GeneScores
 * Per-gene motif scores from which score histograms are calculated.
 *
 * \param HistogramsCluster
 * Motif score histograms based on genes in the current cluster.
 *
 * \param pHistogramsPot
 * If non-null, motif score histograms based on genes not in the current cluster.
 *
 * COALESCE finds significant motifs for each cluster by determining which motifs have statistically
 * different distributions in the cluster versus the genomic background.  To do this, one histogram of
 * per-gene frequencies is constructed per motif; these histograms are in turn based on the frequencies
 * with which each motif appears in each gene.  For example, suppose we have three genes \c G1 through \c G3
 * and three motifs \c M1 through \c M3.  Based each gene's sequence, we determine that the motifs have
 * the following frequencies:
 * \code
 * Motif G1 G2 G3
 * M1    1  1  2
 * M2    2  0  2
 * M3    1  2  0
 * \endcode
 * If only genes \c G1 and \c G2 are in the cluster, we build a frequency histogram as follows:
 * \code
 * Motif 0  1  2
 * M1    0  2  0
 * M2    1  0  1
 * M3    0  1  1
 * \endcode
 * That is, in the cluster, there are zero genes within which \c M1 occurs zero times, two in which it
 * occurs once, zero in which it occurs twice, and so forth.  Each row (i.e. each motif's total histogram)
 * must sum to the number of genes in the cluster, and the resulting inverse histograms for genes not
 * in the cluster (i.e. \c G3) is:
 * \code
 * Motif 0  1  2
 * M1    0  0  1
 * M2    0  0  1
 * M3    1  0  0
 * \endcode
 * In COALESCE, additional complexity arises since motif frequencies are continuous (and must thus be
 * discretized in the histograms) and are calculated on a per-subsequence-type basis (giving rise to
 * multiple histograms per cluster).
 *
 * \see
 * Snapshot
 */
    void CCoalesceCluster::CalculateHistograms(const CCoalesceGeneScores &GeneScores,
                                               CCoalesceGroupHistograms &HistogramsCluster,
                                               CCoalesceGroupHistograms *pHistogramsPot) const {
        set<size_t>::const_iterator iterGene;
        size_t i;

        if (!GeneScores.GetMotifs())
            return;
        for (iterGene = GetGenes().begin(); iterGene != GetGenes().end(); ++iterGene)
            if (!binary_search(m_veciPrevGenes.begin(), m_veciPrevGenes.end(), *iterGene)) {
                HistogramsCluster.Add(GeneScores, *iterGene, false);
                if (pHistogramsPot)
                    pHistogramsPot->Add(GeneScores, *iterGene, true);
            }
        for (i = 0; i < m_veciPrevGenes.size(); ++i)
            if (GetGenes().find(m_veciPrevGenes[i]) == GetGenes().end()) {
                HistogramsCluster.Add(GeneScores, m_veciPrevGenes[i], true);
                if (pHistogramsPot)
                    pHistogramsPot->Add(GeneScores, m_veciPrevGenes[i], false);
            }
    }

/*!
 * \brief
 * Subtract the average expression value for each condition in the cluster from each gene in the cluster.
 *
 * \param PCL
 * Expression matrix from which cluster averages are subtracted.
 *
 * \param Pot
 * Inverse of genes in the cluster; used to determine the difference of the cluster's average from
 * the existing per-condition average.
 *
 * \remarks
 * This effectively masks the average effect of each condition in the cluster from the contained genes'
 * expression values so that the cluster won't be re-found later.  Actually subtracts the difference
 * between the cluster average and the overall average from each condition, since the overall average
 * need not be zero.
 */
    void CCoalesceCluster::Subtract(CPCL &PCL, const CCoalesceCluster &Pot) const {
        set<size_t>::const_iterator iterGene, iterDataset;
        size_t i, iCondition;
        float d, dAve;

        for (iterGene = GetGenes().begin(); iterGene != GetGenes().end(); ++iterGene)
            for (iterDataset = m_setiDatasets.begin(); iterDataset != m_setiDatasets.end(); ++iterDataset)
                for (i = 0; i < GetConditions(*iterDataset); ++i) {
                    iCondition = GetCondition(*iterDataset, i);
                    if (!CMeta::IsNaN(d = m_vecdCentroid[iCondition])) {
                        if (CMeta::IsNaN(dAve = Pot.m_vecdCentroid[iCondition]))
                            dAve = 0;
                        else
                            dAve = ((GetGenes().size() * d) + (Pot.GetGenes().size() * dAve)) /
                                   (GetGenes().size() + Pot.GetGenes().size());
                        PCL.Get(*iterGene, iCondition) -= d - dAve;
                    }
                }
    }

/*!
 * \brief
 * Subtract the average frequency of each motif in the cluster from the score for that motif in
 * each gene in the cluster.
 *
 * \param GeneScores
 * Per-gene motif frequency scores from which the cluster averages are subtracted.
 *
 * \remarks
 * This effectively masks the average effect of each motif in the cluster from the contained genes'
 * sequences so that they won't be re-found later.
 */
    void CCoalesceCluster::Subtract(CCoalesceGeneScores &GeneScores) const {
        set<size_t>::const_iterator iterGene;
        set<SMotifMatch>::const_iterator iterMotif;

        for (iterGene = GetGenes().begin(); iterGene != GetGenes().end(); ++iterGene)
            for (iterMotif = m_setsMotifs.begin(); iterMotif != m_setsMotifs.end(); ++iterMotif)
                GeneScores.Subtract(*iterMotif, *iterGene);
    }

    void CCoalesceClusterImpl::CalculateCentroid(const CPCL &PCL) {
        set<size_t>::const_iterator iterGene;
        size_t i, j;
        float d;

        m_veciCounts.resize(PCL.GetExperiments());
        fill(m_veciCounts.begin(), m_veciCounts.end(), 0);
        m_vecdCentroid.resize(PCL.GetExperiments());
        fill(m_vecdCentroid.begin(), m_vecdCentroid.end(), 0.0f);
        m_vecdStdevs.resize(PCL.GetExperiments());
        fill(m_vecdStdevs.begin(), m_vecdStdevs.end(), 0.0f);
        for (iterGene = GetGenes().begin(); iterGene != GetGenes().end(); ++iterGene)
            for (i = 0; i < m_veciCounts.size(); ++i)
                if (!CMeta::IsNaN(d = PCL.Get(*iterGene, i))) {
                    m_veciCounts[i]++;
                    m_vecdCentroid[i] += d;
                    m_vecdStdevs[i] += d * d;
                }
        for (i = 0; i < m_veciCounts.size(); ++i) {
            if ((j = m_veciCounts[i])) {
                m_vecdCentroid[i] /= j;
                m_vecdStdevs[i] = (m_vecdStdevs[i] / ((j == 1) ? 1 : (j - 1))) -
                                  (m_vecdCentroid[i] * m_vecdCentroid[i]);
                m_vecdStdevs[i] = (m_vecdStdevs[i] < 0) ? 0 : sqrt(m_vecdStdevs[i]);
            } else
                m_vecdCentroid[i] = CMeta::GetNaN();
            g_CatSleipnir().debug(
                    "CCoalesceClusterImpl::CalculateCentroid( ) condition %d: count %d, mean %g, stdev %g",
                    i, m_veciCounts[i], m_vecdCentroid[i], m_vecdStdevs[i]);
        }

        for (i = 0; i < m_vecsDatasets.size(); ++i) {
            SDataset &sDataset = m_vecsDatasets[i];

            if (sDataset.GetConditions() == 1)
                continue;
            sDataset.m_vecdCentroid.resize(sDataset.GetConditions());
            for (j = 0; j < sDataset.GetConditions(); ++j)
                sDataset.m_vecdCentroid[j] = m_vecdCentroid[sDataset.GetCondition(j)];
        }
    }

    void *CCoalesceClusterImpl::ThreadSelectCondition(void *pData) {
        vector<float> vecdDatasetCluster, vecdDatasetPot;
        vector <size_t> veciDatasetCluster, veciDatasetPot;
        size_t iDataset, iCondition, iCluster, iPot, iGene;
        double dP, dZ;
        float d;
        SThreadSelectCondition *psData;
        float *adCluster;
        float *adPot;

        psData = (SThreadSelectCondition *) pData;
        const vector <size_t> &veciPot = *psData->m_pveciPot;
        const vector <size_t> &veciCluster = *psData->m_pveciCluster;
        const CPCL &PCL = *psData->m_pPCL;

        adCluster = new float[veciCluster.size()];
        adPot = new float[veciPot.size()];
        for (iDataset = psData->m_iOffset; iDataset < psData->m_pvecsDatasets->size();
             iDataset += psData->m_iStep) {
            SDataset &sDataset = (*psData->m_pvecsDatasets)[iDataset];

            if (sDataset.GetConditions() == 1) {
                iCondition = sDataset.GetCondition(0);
                for (iCluster = iPot = iGene = 0; iGene < veciCluster.size(); ++iGene)
                    if (!CMeta::IsNaN(d = PCL.Get(veciCluster[iGene], iCondition)))
                        adCluster[iCluster++] = d;
                for (iGene = 0; iGene < veciPot.size(); ++iGene)
                    if (!CMeta::IsNaN(d = PCL.Get(veciPot[iGene], iCondition)))
                        adPot[iPot++] = d;
                if (!(iCluster && iPot))
                    continue;
                dZ = CStatistics::ZScore(adCluster, adCluster + iCluster, adPot, adPot + iPot);
                sDataset.m_dP = CStatistics::ZTest(dZ, min(iCluster, iPot));
                sDataset.m_dZ = (float) min(dZ, (double) FLT_MAX);
            } else {
                vecdDatasetCluster.resize(sDataset.GetConditions());
                fill(vecdDatasetCluster.begin(), vecdDatasetCluster.end(), 0.0f);
                vecdDatasetPot.resize(sDataset.GetConditions());
                fill(vecdDatasetPot.begin(), vecdDatasetPot.end(), 0.0f);
                veciDatasetCluster.resize(sDataset.GetConditions());
                fill(veciDatasetCluster.begin(), veciDatasetCluster.end(), 0);
                veciDatasetPot.resize(sDataset.GetConditions());
                fill(veciDatasetPot.begin(), veciDatasetPot.end(), 0);
                for (iGene = 0; iGene < veciCluster.size(); ++iGene)
                    for (iCondition = 0; iCondition < sDataset.GetConditions(); ++iCondition)
                        if (!CMeta::IsNaN(d = PCL.Get(veciCluster[iGene],
                                                      sDataset.GetCondition(iCondition)))) {
                            veciDatasetCluster[iCondition]++;
                            vecdDatasetCluster[iCondition] += d;
                        }
                for (iGene = 0; iGene < veciPot.size(); ++iGene)
                    for (iCondition = 0; iCondition < sDataset.GetConditions(); ++iCondition)
                        if (!CMeta::IsNaN(d = PCL.Get(veciPot[iGene],
                                                      sDataset.GetCondition(iCondition)))) {
                            veciDatasetPot[iCondition]++;
                            vecdDatasetPot[iCondition] += d;
                        }
                for (iCluster = iPot = iCondition = 0; iCondition < vecdDatasetCluster.size(); ++iCondition) {
                    if (veciDatasetCluster[iCondition] > iCluster)
                        iCluster = veciDatasetCluster[iCondition];
                    if (veciDatasetPot[iCondition] > iPot)
                        iPot = veciDatasetPot[iCondition];
                    if (veciDatasetCluster[iCondition])
                        vecdDatasetCluster[iCondition] /= veciDatasetCluster[iCondition];
                    if (veciDatasetPot[iCondition])
                        vecdDatasetPot[iCondition] /= veciDatasetPot[iCondition];
                }
                dP = CStatistics::MultivariateNormalCDF(vecdDatasetCluster, vecdDatasetPot,
                                                        sDataset.m_psDataset->m_MatSigmaChol, min(iCluster, iPot));
                if (dP > 0.5)
                    dP = 1 - dP;

                sDataset.m_dP = 2 * dP;
                dZ = 0;
                for (iCondition = 0; iCondition < sDataset.GetConditions(); ++iCondition)
                    if (sDataset.m_psDataset->m_vecdStdevs[iCondition]) {
                        d = (vecdDatasetCluster[iCondition] - vecdDatasetPot[iCondition]) /
                            sDataset.m_psDataset->m_vecdStdevs[iCondition];
                        dZ += d * d;
                    }
                dZ = sqrt(dZ);
                sDataset.m_dZ = (float) min(dZ, (double) FLT_MAX);
            }
        }
        delete[] adCluster;
        delete[] adPot;

        return NULL;
    }

/*!
 * \brief
 * Performs feature selection to include significant expression conditions in a converging cluster.
 *
 * \param PCL
 * Expression dataset from which significant datasets are selected.
 *
 * \param Pot
 * Inverse of current cluster used for genomic background calculations.
 *
 * \param iThreads
 * Maximum number of simultaneous threads for condition significance calculations.
 *
 * \param dPValue
 * P-value threshhold for condition significance.
 *
 * \param dZScore
 * Z-score effect size threshhold for condition significance.
 *
 * \returns
 * True if zero or more conditions were selected successfully, false otherwise.
 *
 * Selects zero or more conditions in which the cluster's current gene set is differentially expressed.
 * That is, in each selected condition, the average expression of genes in the cluster must differ from
 * the genomic background with at least the given significance and effect size threshholds.  If dataset
 * blocks were given at cluster initialization time, all conditions in the block are added (or not)
 * simultaneously using a multivariate significance test.
 *
 * \see
 * SelectMotifs | SelectGenes
 */
    bool CCoalesceCluster::SelectConditions(const CPCL &PCL, const CCoalesceCluster &Pot, size_t iThreads,
                                            float dPValue, float dZScore) {
        vector <pthread_t> vecpthdThreads;
        vector <SThreadSelectCondition> vecsThreads;
        size_t i;
        vector <size_t> veciCluster, veciPot;

        m_setiDatasets.clear();
        veciCluster.resize(GetGenes().size());
        copy(GetGenes().begin(), GetGenes().end(), veciCluster.begin());
        veciPot.resize(Pot.GetGenes().size());
        copy(Pot.GetGenes().begin(), Pot.GetGenes().end(), veciPot.begin());
        vecpthdThreads.resize(iThreads);
        vecsThreads.resize(vecpthdThreads.size());
        for (i = 0; i < vecsThreads.size(); ++i) {
            vecsThreads[i].m_iOffset = i;
            vecsThreads[i].m_iStep = vecsThreads.size();
            vecsThreads[i].m_pveciCluster = &veciCluster;
            vecsThreads[i].m_pveciPot = &veciPot;
            vecsThreads[i].m_pvecsDatasets = &m_vecsDatasets;
            vecsThreads[i].m_pPCL = &PCL;
            if (pthread_create(&vecpthdThreads[i], NULL, ThreadSelectCondition, &vecsThreads[i])) {
                g_CatSleipnir().error("CCoalesceCluster::SelectConditions( %d, %g ) could not select conditions",
                                      iThreads, dZScore);
                return false;
            }
        }
        for (i = 0; i < vecpthdThreads.size(); ++i)
            pthread_join(vecpthdThreads[i], NULL);

        dPValue /= m_vecsDatasets.size();
        for (i = 0; i < m_vecsDatasets.size(); ++i) {
            const SDataset &sDataset = m_vecsDatasets[i];

            if ((sDataset.m_dP < dPValue) && (fabs(sDataset.m_dZ) > dZScore)) {
                g_CatSleipnir().debug("CCoalesceCluster::SelectConditions( %d, %g ) selected dataset %d at %g, z=%g",
                                      iThreads, dZScore, i, sDataset.m_dP, sDataset.m_dZ);
                m_setiDatasets.insert(i);
            } else
                g_CatSleipnir().debug("CCoalesceCluster::SelectConditions( %d, %g ) rejected dataset %d at %g, z=%g",
                                      iThreads, dZScore, i, sDataset.m_dP, sDataset.m_dZ);
        }

        return true;
    }

    void *CCoalesceClusterImpl::ThreadSelectMotif(void *pData) {
        SThreadSelectMotif *psData;
        uint32_t iMotif;

        psData = (SThreadSelectMotif *) pData;
        for (iMotif = psData->m_iOffset; iMotif < (psData->m_pveciMotifs ? psData->m_pveciMotifs->size() :
                                                   psData->m_pMotifs->GetMotifs()); iMotif += psData->m_iStep)
            if (!AddSignificant(*psData->m_pMotifs, psData->m_pveciMotifs ? (*psData->m_pveciMotifs)[iMotif] :
                                                    iMotif, *psData->m_pHistsCluster, *psData->m_pHistsPot,
                                psData->m_dPValue, psData->m_dZScore,
                                psData->m_vecsMotifs))
                break;

        return NULL;
    }

/*!
 * \brief
 * Performs feature selection to include significant sequence motifs in a converging cluster.
 *
 * \param HistsCluster
 * Precalculated histogram of motif frequencies using genes in the cluster.
 *
 * \param HistsPot
 * Precalculated histogram of motif frequencies using genes not in the cluster.
 *
 * \param dPValue
 * P-value threshhold for motif significance.
 *
 * \param dZScore
 * Z-score effect size threshhold for motif significance.
 *
 * \param iMaxMotifs
 * Maximum number of motifs associated with a cluster; if more motifs are present, no new selection is
 * performed.
 *
 * \param iThreads
 * Maximum number of simultaneous threads for motif significance calculations.
 *
 * \param pMotifs
 * If non-null, motif library with which motifs are managed.
 *
 * \returns
 * True if zero or more significant motifs were selected; false otherwise.
 *
 * Selects zero or more motifs differentially over- or under-enriched in the genes currently in a
 * converging cluster.  Motifs are selected on a per-sequence-subtype basis, so a motif may be enriched,
 * for example, in an upstream but not downstream flank.  All motifs managed by the given library are tested
 * for significance, including all k-mers, reverse complements, and any currently constructed PSTs.
 * Significance testing is performed by z-scoring the within versus without frequency histograms, and all
 * p-values are Bonferroni corrected.
 *
 * \remarks
 * No motif selection can be performed if pMotifs is null; this will not generate an error, but only
 * expression biclustering (i.e. conditions and genes) will be performed.
 *
 * \see
 * CalculateHistograms | SelectConditions | SelectGenes
 */
    bool CCoalesceCluster::SelectMotifs(const CCoalesceGroupHistograms &HistsCluster,
                                        const CCoalesceGroupHistograms &HistsPot, float dPValue, float dZScore,
                                        size_t iMaxMotifs,
                                        size_t iThreads, const CCoalesceMotifLibrary *pMotifs) {
        uint32_t i, j;
        vector <pthread_t> vecpthdThreads;
        vector <SThreadSelectMotif> vecsThreads;
        vector <uint32_t> veciMotifs;

        if (!pMotifs)
            return true;
        if (m_setsMotifs.size() > iMaxMotifs) {
            set <uint32_t> setiMotifs;
            set<SMotifMatch>::const_iterator iterMotif;

            for (iterMotif = m_setsMotifs.begin(); iterMotif != m_setsMotifs.end(); ++iterMotif)
                setiMotifs.insert(iterMotif->m_iMotif);
            veciMotifs.resize(setiMotifs.size());
            copy(setiMotifs.begin(), setiMotifs.end(), veciMotifs.begin());
        }
        m_setsMotifs.clear();
        vecpthdThreads.resize(iThreads);
        vecsThreads.resize(vecpthdThreads.size());
        for (i = 0; i < vecsThreads.size(); ++i) {
            vecsThreads[i].m_iOffset = i;
            vecsThreads[i].m_iStep = vecsThreads.size();
            vecsThreads[i].m_pMotifs = pMotifs;
            vecsThreads[i].m_pHistsCluster = &HistsCluster;
            vecsThreads[i].m_pHistsPot = &HistsPot;
            vecsThreads[i].m_dPValue = dPValue;
            vecsThreads[i].m_dZScore = dZScore;
            vecsThreads[i].m_pveciMotifs = veciMotifs.empty() ? NULL : &veciMotifs;
            if (pthread_create(&vecpthdThreads[i], NULL, ThreadSelectMotif, &vecsThreads[i])) {
                g_CatSleipnir().error("CCoalesceCluster::SelectMotifs( %g, %d ) could not select motifs",
                                      dZScore, iThreads);
                return false;
            }
        }
        for (i = 0; i < vecpthdThreads.size(); ++i) {
            pthread_join(vecpthdThreads[i], NULL);
            for (j = 0; j < vecsThreads[i].m_vecsMotifs.size(); ++j)
                m_setsMotifs.insert(vecsThreads[i].m_vecsMotifs[j]);
        }

        return true;
    }

    bool CCoalesceClusterImpl::AddSignificant(const CCoalesceMotifLibrary &Motifs, uint32_t iMotif,
                                              const CCoalesceGroupHistograms &HistsCluster,
                                              const CCoalesceGroupHistograms &HistsPot, float dPValue,
                                              float dZScore, vector <SMotifMatch> &vecsMotifs) {
        size_t iTypeCluster, iTypePot;
        double dP, dAveOne, dAverage, dZ;
        CCoalesceSequencerBase::ESubsequence eSubsequence;

        dPValue /= Motifs.GetMotifs();
        for (iTypeCluster = 0; iTypeCluster < HistsCluster.GetTypes(); ++iTypeCluster) {
            const string &strTypeCluster = HistsCluster.GetType(iTypeCluster);

            if ((iTypePot = HistsPot.GetType(strTypeCluster)) == -1)
                continue;
            for (eSubsequence = CCoalesceSequencerBase::ESubsequenceBegin;
                 (size_t) eSubsequence < HistsCluster.GetSubsequences(iTypeCluster);
                 eSubsequence = (CCoalesceSequencerBase::ESubsequence) ((size_t) eSubsequence + 1)) {
                const CCoalesceHistogramSet<> &HistSetCluster = HistsCluster.Get(iTypeCluster, eSubsequence);
                const CCoalesceHistogramSet<> &HistSetPot = HistsPot.Get(iTypePot, eSubsequence);

                if ((HistSetCluster.GetMembers() <= iMotif) ||
                    (HistSetPot.GetMembers() <= iMotif) ||
                    !(HistSetCluster.GetTotal() && HistSetPot.GetTotal()))
                    continue;
                dP = HistSetCluster.CohensD(iMotif, HistSetPot, dAveOne, dAverage, dZ);
                if ((dP < dPValue) && (fabs(dZ) > dZScore)) {
                    SMotifMatch sMotif(iMotif, strTypeCluster, eSubsequence, (float) dZ, (float) (dAveOne -
                                                                                                  dAverage));

                    if (g_CatSleipnir().isInfoEnabled()) {
                        ostringstream ossm;

                        ossm << "CCoalesceClusterImpl::AddSignificant( " << iMotif << ", " << dZScore <<
                             " ) adding at " << dP << ":" << endl << sMotif.Save(&Motifs) << endl <<
                             "Cluster	" << HistSetCluster.Save(iMotif) << endl <<
                             "Pot	" << HistSetPot.Save(iMotif);
                        g_CatSleipnir().info(ossm.str().c_str());
                    }
                    vecsMotifs.push_back(sMotif);
                } else if (g_CatSleipnir().isDebugEnabled()) {
                    ostringstream ossm;

                    ossm << "CCoalesceClusterImpl::AddSignificant( " << iMotif << ", " << dZScore <<
                         " ) failed at " << dP << ":" << endl << SMotifMatch(iMotif, strTypeCluster, eSubsequence,
                                                                             (float) dZ,
                                                                             (float) (dAveOne - dAverage)).Save(&Motifs)
                         << endl <<
                         "Cluster	" << HistSetCluster.Save(iMotif) << endl <<
                         "Pot	" << HistSetPot.Save(iMotif);
                    g_CatSleipnir().debug(ossm.str().c_str());
                }
            }
        }

        return true;
    }

    void *CCoalesceClusterImpl::ThreadCentroid(void *pData) {
        SThreadCentroid *psData;

        psData = (SThreadCentroid *) pData;
        psData->m_pCluster->CalculateCentroid(psData->m_PCL);

        return NULL;
    }

    void *CCoalesceClusterImpl::ThreadSignificantGene(void *pData) {
        SThreadSignificantGene *psData;
        size_t i;

        psData = (SThreadSignificantGene *) pData;
        for (i = psData->m_iOffset; i < psData->m_pvecfSignificant->size(); i += psData->m_iStep)
            (*psData->m_pvecfSignificant)[i] = psData->m_pCluster->IsSignificant(i, *psData->m_pPCL,
                                                                                 *psData->m_pvecdStdevs,
                                                                                 psData->m_pMotifs,
                                                                                 *psData->m_pGeneScores,
                                                                                 *psData->m_pHistsCluster,
                                                                                 *psData->m_pHistsPot, *psData->m_pPot,
                                                                                 *psData->m_pveciDatasets,
                                                                                 psData->m_dBeta,
                                                                                 psData->m_iMinimum,
                                                                                 psData->m_dProbability);

        return NULL;
    }

/*!
 * \brief
 * Add and remove (im)probable genes to a converging cluster using Bayesian integration based on the
 * currently selected expression conditions and sequence motifs.
 *
 * \param PCL
 * PCL from which genes are selected.
 *
 * \param GeneScores
 * Per-gene motif scores used to determine each gene's probability based on motif frequencies.
 *
 * \param HistsCluster
 * Precalculated histogram of motif frequencies using genes in the cluster.
 *
 * \param HistsPot
 * Precalculated histogram of motif frequencies using genes not in the cluster.
 *
 * \param iMinimum
 * Minimum number of genes for which data must be present for a given motif for it to contribute to gene
 * probabilities.
 *
 * \param iThreads
 * Maximum number of simultaneous threads for gene probability calculations.
 *
 * \param Pot
 * Inverse of genes in the cluster; used to determine genomic background probabilities.
 *
 * \param dProbability
 * Probability threshhold above which genes are included in the cluster.
 *
 * \param pMotifs
 * If non-null, motif library managing significant sequence motifs.
 *
 * \returns
 * True if zero or more genes were included in the cluster, false otherwise.
 *
 * Adds and removes genes to/from a converging cluster based on the currently selected significant
 * expression conditions and sequence motifs.  The probability of gene membership in the cluster given
 * some data, <tt>P(g in C|D)</tt>, is calculated using Bayes rule as <tt>P(D|g in C)P(g in C)/(P(D|g in C) +
 * P(D|g notin C))</tt>.  All data features are assumed to be independent, such that <tt>P(D|g in C)</tt> is
 * a product over each significant condition/motif of A) the probability of the gene's expression value being
 * drawn from the cluster's distribution of expression values or B) the probability of its frequency for some
 * motif being drawn from the cluster's distribution of frequencies for that motif.  Distributions are assumed
 * to be normal, and frequencies are Laplace smoothed to avoid zeros.  The prior <tt>P(g in C)</tt> is used
 * to stabilize cluster convergence, such that it is equal to 1 if a gene was already included in the cluster
 * in a previous iteration and equal to dProbability if the gene was not previously included.
 *
 * \remarks
 * A prior is also applied to the two main data types, expression conditions and sequence motifs.  If X and M
 * are the sets of all possible conditions and motifs, respectively, and CX and CM are the significant
 * conditions and motifs in the cluster, then expression probabilities have prior |CX|/|C| and motifs have prior
 * |CM|/|M|.
 *
 * \see
 * SelectConditions | SelectMotifs
 */
    bool CCoalesceCluster::SelectGenes(const CPCL &PCL, const CCoalesceGeneScores &GeneScores,
                                       const CCoalesceGroupHistograms &HistsCluster,
                                       const CCoalesceGroupHistograms &HistsPot, size_t iMinimum,
                                       size_t iThreads, CCoalesceCluster &Pot, float dProbability,
                                       const CCoalesceMotifLibrary *pMotifs) {
        size_t i, iCluster, iPot;
        vector <pthread_t> vecpthdThreads;
        vector<bool> vecfSignificant;
        vector <SThreadSignificantGene> vecsThreads;
        vector <size_t> veciDatasets;
        float dSSCluster, dSSPot, dAve;
        set <size_t> setiMotifs;
        set<SMotifMatch>::const_iterator iterMotif;
        vector<float> vecdStdevs;

        vecpthdThreads.resize(iThreads);
        {
            SThreadCentroid sUs(this, PCL), sThem(&Pot, PCL);

            if (iThreads > 1) {
                if (pthread_create(&vecpthdThreads[0], NULL, ThreadCentroid, &sUs) ||
                    pthread_create(&vecpthdThreads[1], NULL, ThreadCentroid, &sThem)) {
                    g_CatSleipnir().error("CCoalesceCluster::SelectGenes( %d, %g ) could not calculate centroids",
                                          iThreads, dProbability);
                    return false;
                }
                for (i = 0; i < 2; ++i)
                    pthread_join(vecpthdThreads[i], NULL);
            } else if (ThreadCentroid(&sUs) || ThreadCentroid(&sThem)) {
                g_CatSleipnir().error("CCoalesceCluster::SelectGenes( %d, %g ) could not calculate centroids",
                                      iThreads, dProbability);
                return false;
            }
        }
        iCluster = GetGenes().size();
        iPot = Pot.GetGenes().size();
        vecdStdevs.resize(m_vecdStdevs.size());
        for (i = 0; i < vecdStdevs.size(); ++i) {
            dSSCluster = m_vecdStdevs[i];
            dSSCluster = iCluster * ((dSSCluster * dSSCluster) + (m_vecdCentroid[i] *
                                                                  m_vecdCentroid[i]));
            dSSPot = Pot.m_vecdStdevs[i];
            dSSPot = iPot * ((dSSPot * dSSPot) + (Pot.m_vecdCentroid[i] * Pot.m_vecdCentroid[i]));
            dAve = ((m_vecdCentroid[i] * iCluster) + (Pot.m_vecdCentroid[i] * iPot)) / PCL.GetGenes();
            vecdStdevs[i] = sqrt(((dSSCluster + dSSPot) / (iCluster + iPot)) - (dAve * dAve));
        }

        veciDatasets.resize(m_setiDatasets.size());
        copy(m_setiDatasets.begin(), m_setiDatasets.end(), veciDatasets.begin());
        vecfSignificant.resize(PCL.GetGenes());
        vecsThreads.resize(vecpthdThreads.size());
        for (i = 0; i < vecpthdThreads.size(); ++i) {
            vecsThreads[i].m_iOffset = i;
            vecsThreads[i].m_iStep = vecpthdThreads.size();
            vecsThreads[i].m_pvecfSignificant = &vecfSignificant;
            vecsThreads[i].m_pPCL = &PCL;
            vecsThreads[i].m_pMotifs = pMotifs;
            vecsThreads[i].m_pGeneScores = &GeneScores;
            vecsThreads[i].m_pHistsCluster = &HistsCluster;
            vecsThreads[i].m_pHistsPot = &HistsPot;
            vecsThreads[i].m_pCluster = this;
            vecsThreads[i].m_pPot = &Pot;
            vecsThreads[i].m_pveciDatasets = &veciDatasets;
            vecsThreads[i].m_pvecdStdevs = &vecdStdevs;
            vecsThreads[i].m_dBeta = m_setsMotifs.size() ? ((float) m_setsMotifs.size() /
                                                            (m_setiDatasets.size() + m_setsMotifs.size())) : 0.5f;
            vecsThreads[i].m_iMinimum = iMinimum;
            vecsThreads[i].m_dProbability = dProbability;
            if (pthread_create(&vecpthdThreads[i], NULL, ThreadSignificantGene, &vecsThreads[i])) {
                g_CatSleipnir().error("CCoalesceCluster::SelectGenes( %d, %g ) could not calculate significance",
                                      iThreads, dProbability);
                return false;
            }
        }
        for (i = 0; i < vecpthdThreads.size(); ++i)
            pthread_join(vecpthdThreads[i], NULL);
        m_setiGenes.clear();
        Pot.m_setiGenes.clear();
        for (i = 0; i < vecfSignificant.size(); ++i)
            if (vecfSignificant[i]) {
                m_setiGenes.insert(i);
                m_vecdPriors[i] = dProbability;
            } else {
                Pot.m_setiGenes.insert(i);
                m_vecdPriors[i] = 1 - dProbability;
            }

        return true;
    }

    bool CCoalesceClusterImpl::IsSignificant(size_t iGene, const CPCL &PCL, const vector<float> &vecdStdevs,
                                             const CCoalesceMotifLibrary *pMotifs,
                                             const CCoalesceGeneScores &GeneScores,
                                             const CCoalesceGroupHistograms &HistsCluster,
                                             const CCoalesceGroupHistograms &HistsPot,
                                             const CCoalesceCluster &Pot, const vector <size_t> &veciDatasets,
                                             float dBeta, size_t iMinimum,
                                             float dProbability) const {
        float d, dP, dLogPMotifsGivenIn, dLogPMotifsGivenOut, dLogPExpressionGivenIn, dLogPExpressionGivenOut;
        bool fClustered;

        fClustered = IsGene(iGene);
// We want P(g in C|S,E) =
// P(S,E|g in C)P(g in C)/P(S,E) =
// P(S|g in C)P(E|g in C)P(g in C)/(P(S,E|g in C)P(g in C) + P(S,E|g notin C)P(g notin C)) =
// P(S|g in C)P(E|g in C)P(g in C)/(P(S|g in C)P(E|g in C)P(g in C) + P(S|g notin C)P(E|g notin C)P(g notin C))
        if (!(CalculateProbabilityExpression(iGene, PCL, vecdStdevs, Pot, veciDatasets, fClustered,
                                             dLogPExpressionGivenIn, dLogPExpressionGivenOut) &&
              CalculateProbabilityMotifs(GeneScores, iGene,
                                         HistsCluster, HistsPot, fClustered, iMinimum, dLogPMotifsGivenIn,
                                         dLogPMotifsGivenOut)))
            return false;
        if ((d = m_vecdPriors[iGene]) >= 1)
            dP = 1;
        else if (d <= 0)
            dP = 0;
        else
            dP = 1 / (1 + exp((dBeta * (dLogPExpressionGivenOut - dLogPExpressionGivenIn) + (1 - dBeta) *
                                                                                            (dLogPMotifsGivenOut -
                                                                                             dLogPMotifsGivenIn)) /
                              (0.5f + 2 * pow(0.5f - dBeta, 2)) +
                              log((1 - d) / d)));

        if (g_CatSleipnir().isDebugEnabled()) {
            set<SMotifMatch>::const_iterator iterMotif;
            size_t iType;

            g_CatSleipnir().debug(
                    "CCoalesceClusterImpl::IsSignificant( %s ) is %g, prior %g beta %g, exp. p=%g vs. %g, seq. p=%g vs %g",
                    PCL.GetGene(iGene).c_str(), dP, d, dBeta, dLogPExpressionGivenIn, dLogPExpressionGivenOut,
                    dLogPMotifsGivenIn, dLogPMotifsGivenOut);
            for (iterMotif = m_setsMotifs.begin(); iterMotif != m_setsMotifs.end(); ++iterMotif)
                if ((iType = GeneScores.GetType(iterMotif->m_strType)) != -1)
                    g_CatSleipnir().debug("%g	%s", GeneScores.Get(iType, iterMotif->m_eSubsequence, iGene,
                                                                     iterMotif->m_iMotif),
                                          iterMotif->Save(pMotifs).c_str());
        }

        return (dP >= dProbability);
    }

    bool CCoalesceClusterImpl::CalculateProbabilityExpression(size_t iGene, const CPCL &PCL,
                                                              const vector<float> &vecdStdevs,
                                                              const CCoalesceCluster &Pot,
                                                              const vector <size_t> &veciDatasets,
                                                              bool fClustered, float &dLogPIn, float &dLogPOut) const {
        static const double c_dEpsilonZero = 1e-10;
        float dGene, dCluster, dPot;
        double dPCluster, dPPot, dStdev;
        size_t iDataset, iPot, iCondition;
        long double dPIn, dPOut;

        iPot = PCL.GetGenes() - GetGenes().size();
        dLogPIn = dLogPOut = 0;
        dPIn = dPOut = 1;
        for (iDataset = 0; iDataset < veciDatasets.size(); ++iDataset) {
            const SDataset &sDataset = m_vecsDatasets[veciDatasets[iDataset]];

            if (sDataset.GetConditions() == 1) {
                iCondition = sDataset.GetCondition(0);
                if (CMeta::IsNaN(dGene = PCL.Get(iGene, iCondition)) ||
                    CMeta::IsNaN(dCluster = m_vecdCentroid[iCondition]) ||
                    CMeta::IsNaN(dPot = Pot.m_vecdCentroid[iCondition]))
                    continue;

                if (fClustered)
                    dCluster = ((dCluster * GetGenes().size()) - dGene) /
                               (GetGenes().size() - 1);
                else
                    dPot = ((dPot * iPot) - dGene) / (iPot - 1);
                dStdev = max(c_dEpsilonZero, (double) vecdStdevs[iCondition]);
                dPCluster = max(c_dEpsilonZero, CStatistics::NormalPDF(dGene, dCluster, dStdev));
                dPPot = max(c_dEpsilonZero, CStatistics::NormalPDF(dGene, dPot, dStdev));
            } else {
                vector<float> vecdGene;

                vecdGene.resize(sDataset.GetConditions());
                for (iCondition = 0; iCondition < sDataset.GetConditions(); ++iCondition)
                    vecdGene[iCondition] = PCL.Get(iGene, sDataset.GetCondition(iCondition));
                dPCluster = min(1.0, max(c_dEpsilonZero, CStatistics::MultivariateNormalPDF(vecdGene,
                                                                                            sDataset.m_vecdCentroid,
                                                                                            sDataset.m_psDataset->m_dSigmaDetSqrt,
                                                                                            sDataset.m_psDataset->m_MatSigmaInv)));
                dPPot = min(1.0, max(c_dEpsilonZero, CStatistics::MultivariateNormalPDF(vecdGene,
                                                                                        Pot.m_vecsDatasets[veciDatasets[iDataset]].m_vecdCentroid,
                                                                                        sDataset.m_psDataset->m_dSigmaDetSqrt,
                                                                                        sDataset.m_psDataset->m_MatSigmaInv)));
            }
            dPIn *= dPCluster;
            dPOut *= dPPot;
            if ((dPIn < DBL_MIN) || (dPOut < DBL_MIN)) {
                dLogPIn += (float) log(dPIn);
                dLogPOut += (float) log(dPOut);
                dPIn = dPOut = 1;
            }
        }
        dLogPIn += (float) log(dPIn);
        dLogPOut += (float) log(dPOut);

        return true;
    }

    bool CCoalesceClusterImpl::CalculateProbabilityMotifs(const CCoalesceGeneScores &GeneScores, size_t iGene,
                                                          const CCoalesceGroupHistograms &HistsCluster,
                                                          const CCoalesceGroupHistograms &HistsPot,
                                                          bool fClustered, size_t iMinimum, float &dLogPIn,
                                                          float &dLogPOut) const {
        set<SMotifMatch>::const_iterator iterMotif;
        size_t iType, iCluster, iPot;
        unsigned short sCluster, sPot;
        double dPCluster, dPPot;
        float dGene;
        const float *adValues;
        long double dPIn, dPOut;

        dLogPIn = dLogPOut = 0;
        if (m_setsMotifs.empty())
            return true;

// TODO: this is the slowest part of the entire algorithm
        dPIn = dPOut = 1;
        for (iterMotif = m_setsMotifs.begin(); iterMotif != m_setsMotifs.end(); ++iterMotif) {
            if (((iType = GeneScores.GetType(iterMotif->m_strType)) == -1) ||
                !(adValues = GeneScores.Get(iType, iterMotif->m_eSubsequence, iGene)) ||
                (HistsCluster.GetType(iterMotif->m_strType) == -1) ||
                (HistsPot.GetType(iterMotif->m_strType) == -1))
                continue;
            dGene = adValues[iterMotif->m_iMotif];
            {
                const CCoalesceHistogramSet<> &HistSetCluster = HistsCluster.Get(iterMotif->m_strType,
                                                                                 iterMotif->m_eSubsequence);
                const CCoalesceHistogramSet<> &HistSetPot = HistsPot.Get(iterMotif->m_strType,
                                                                         iterMotif->m_eSubsequence);

                sCluster = HistSetCluster.Get(iterMotif->m_iMotif, dGene);
                iCluster = HistSetCluster.GetTotal();
                sPot = HistSetPot.Get(iterMotif->m_iMotif, dGene);
                iPot = HistSetPot.GetTotal();
                if ((iCluster < iMinimum) || (iPot < iMinimum))
                    continue;

                if (fClustered) {
                    if (sCluster) {
                        sCluster--;
                        iCluster--;
                    } else
                        g_CatSleipnir().warn(
                                "CCoalesceClusterImpl::CalculateProbabilityMotifs( %d, %d, %g, %g ) no motifs of %d in cluster: %g, %s\n%s",
                                iGene, fClustered, dLogPIn, dLogPOut, iCluster, dGene,
                                iterMotif->Save(NULL).c_str(), HistSetCluster.Save(iterMotif->m_iMotif).c_str());
                } else if (sPot) {
                    sPot--;
                    iPot--;
                } else
                    g_CatSleipnir().warn(
                            "CCoalesceClusterImpl::CalculateProbabilityMotifs( %d, %d, %g, %g ) no motifs of %d in pot: %g, %s\n%s",
                            iGene, fClustered, dLogPIn, dLogPOut, iPot, dGene, iterMotif->Save(NULL).c_str(),
                            HistSetPot.Save(iterMotif->m_iMotif).c_str());
                dPCluster = (double) (sCluster + 1) / (iCluster + HistSetCluster.GetEdges());
                dPPot = (double) (sPot + 1) / (iPot + HistSetPot.GetEdges());
            }
            dPIn *= dPCluster;
            dPOut *= dPPot;
            if ((dPIn < DBL_MIN) || (dPOut < DBL_MIN)) {
                dLogPIn += (float) log(dPIn);
                dLogPOut += (float) log(dPOut);
                dPIn = dPOut = 1;
            }
        }
        dLogPIn += (float) log(dPIn);
        dLogPOut += (float) log(dPOut);

        return true;
    }

/*!
 * \brief
 * Saves cluster in a pair of PCL and motif text files in the given directory.
 *
 * \param strDirectory
 * Directory in which cluster files are generated.
 *
 * \param iID
 * Integer ID assigned to the cluster.
 *
 * \param PCL
 * PCL from which gene and condition labels are read.
 *
 * \param pMotifs
 * If non-null, motif library with which cluster motifs are formatted.
 *
 * \returns
 * True if the cluster was saved successfully, false otherwise.
 *
 * Saves the cluster as a pair of text files (PCL and motifs) in the given directory.  Filenames are randomly
 * generated in the format \c c<id>_<rand>.pcl and \c c<id>_<rand>_motifs.txt.  The PCL file will contain only
 * genes included in the bicluster along with all conditions; conditions in the bicluster will be sorted to
 * the left and tagged with an initial \c *.  Motif files are of the format:
 * \code
 * <type1>	<subtype1>	<score1>	<tf1>:<match1>|<tf2>:<match2>|...|<tfy>:<matchx>
 * <pst1>
 * <pwm1>
 * <type2>	<subtype2>	<score2>	<tf1>:<match1>|<tf2>:<match2>|...|<tfy>:<matchx>
 * <pst2>
 * <pwm2>
 * ...
 * <typez>	<subtypez>	<scorez>	<tf1>:<match1>|<tf2>:<match2>|...|<tfy>:<matchx>
 * <pstz>
 * <pwmz>
 * \endcode
 *
 * \see
 * Open
 */
    bool CCoalesceCluster::Save(const std::string &strDirectory, size_t iID, const CPCL &PCL,
                                const CCoalesceMotifLibrary *pMotifs) const {
        ofstream ofsm;
        char *szTmp;
        CPCL PCLCopy;
        size_t iGeneFrom, iGeneTo, iExpFrom, iExpTo, iLength;
        string strBase;
        set<SMotifMatch>::const_iterator iterMotif;
        set <size_t> setiConditions;

        GetConditions(setiConditions);
        szTmp = new char[iLength = (strDirectory.length() + 16)];
        sprintf_s(szTmp, iLength - 1, "%s/c%04d_XXXXXX", strDirectory.empty() ? "." : strDirectory.c_str(),
                  iID);
        _mktemp_s(szTmp, iLength - 1);
        szTmp[iLength - 1] = 0;
        strBase = szTmp;
        delete[] szTmp;

        ofsm.open((strBase + c_szMotifs).c_str());
        if (!ofsm.is_open())
            return false;
        for (iterMotif = m_setsMotifs.begin(); iterMotif != m_setsMotifs.end(); ++iterMotif)
            ofsm << iterMotif->Save(pMotifs) << endl;
        ofsm.close();

        PCLCopy.Open(PCL);
// Reorder header
        for (iExpTo = iExpFrom = 0; iExpFrom < PCL.GetExperiments(); ++iExpFrom)
            if (setiConditions.find(iExpFrom) != setiConditions.end())
                PCLCopy.SetExperiment(iExpTo++, c_cStar + PCL.GetExperiment(iExpFrom));
        for (iExpFrom = 0; iExpFrom < PCL.GetExperiments(); ++iExpFrom)
            if (setiConditions.find(iExpFrom) == setiConditions.end())
                PCLCopy.SetExperiment(iExpTo++, PCL.GetExperiment(iExpFrom));
// Reorder genes
        for (iGeneTo = iGeneFrom = 0; iGeneFrom < PCL.GetGenes(); ++iGeneFrom)
            if (IsGene(iGeneFrom) && !SaveCopy(PCL, setiConditions, iGeneFrom, PCLCopy, iGeneTo++, false))
                return false;
        for (iGeneFrom = 0; iGeneFrom < PCL.GetGenes(); ++iGeneFrom)
            if (!IsGene(iGeneFrom))
                PCLCopy.MaskGene(iGeneTo++);

        ofsm.clear();
        ofsm.open((strBase + CPCL::GetExtension()).c_str());
        if (!ofsm.is_open())
            return false;
        PCLCopy.Save(ofsm);
        ofsm.close();

        return true;
    }

    bool CCoalesceClusterImpl::SaveCopy(const CPCL &PCLFrom, const std::set <size_t> &setiConditions,
                                        size_t iGeneFrom, CPCL &PCLTo, size_t iGeneTo, bool fClustered) const {
        size_t i, iExpTo, iExpFrom;

        PCLTo.SetGene(iGeneTo, PCLFrom.GetGene(iGeneFrom));
        for (i = 1; i < PCLFrom.GetFeatures(); ++i)
            PCLTo.SetFeature(iGeneTo, i, (string)((fClustered && (i == 1)) ? "*" : "") +
                                         PCLFrom.GetFeature(iGeneFrom, i));
        for (iExpTo = iExpFrom = 0; iExpFrom < PCLFrom.GetExperiments(); ++iExpFrom)
            if (setiConditions.find(iExpFrom) != setiConditions.end())
                PCLTo.Set(iGeneTo, iExpTo++, PCLFrom.Get(iGeneFrom, iExpFrom));
        for (iExpFrom = 0; iExpFrom < PCLFrom.GetExperiments(); ++iExpFrom)
            if (setiConditions.find(iExpFrom) == setiConditions.end())
                PCLTo.Set(iGeneTo, iExpTo++, PCLFrom.Get(iGeneFrom, iExpFrom));

        return true;
    }

/*!
 * \brief
 * Saves a textual representation of the cluster (including genes, conditions, and motifs) to the given
 * output stream.
 *
 * \param ostm
 * Output stream to which cluster description is saved.
 *
 * \param iID
 * Integer ID assigned to the cluster.
 *
 * \param PCL
 * PCL from which gene and condition labels are read.
 *
 * \param pMotifs
 * If non-null, motif library with which cluster motifs are formatted.
 *
 * \param dCutoffPWMs
 * Minimum information threshhold (in bits) for a PWM to be saved.
 *
 * \param dPenaltyGap
 * Alignment score penalty for gaps.
 *
 * \param dPenaltyMismatch
 * Alignment score penalty for mismatches.
 *
 * \param fNoRCs
 * If true, resolve the given motif into a single strand without reverse complements before saving PWM.
 *
 * Saves cluster in the following tab-delimited textual format:
 * \code
 * Cluster	<id>
 * Genes	<gene1>	<gene2>	...	<genew>
 * Conditions	<cond1>	<cond2>	...	<condx>
 * Motifs
 * <type1>	<subtype1>	<score1>	<tf1>:<match1>|<tf2>:<match2>|...|<tfy>:<matchy>
 * <pst1>
 * <pwm1>
 * <type2>	<subtype2>	<score2>	<tf1>:<match1>|<tf2>:<match2>|...|<tfy>:<matchy>
 * <pst2>
 * <pwm2>
 * ...
 * <typez>	<subtypez>	<scorez>	<tf1>:<match1>|<tf2>:<match2>|...|<tfy>:<matchy>
 * <pstz>
 * <pwmz>
 * \endcode
 *
 * \see
 * Open
 */
    void CCoalesceCluster::Save(std::ostream &ostm, size_t iID, const CPCL &PCL,
                                const CCoalesceMotifLibrary *pMotifs, float dCutoffPWMs, float dPenaltyGap,
                                float dPenaltyMismatch,
                                bool fNoRCs) const {
        set<size_t>::const_iterator iterID;
        set<SMotifMatch>::const_iterator iterMotif;
        size_t i;
        string strMotif;

        ostm << "Cluster\t" << iID << endl;
        ostm << c_szGenes;
        for (iterID = GetGenes().begin(); iterID != GetGenes().end(); ++iterID)
            ostm << '\t' << PCL.GetGene(*iterID);
        ostm << endl << c_szConditions;
        for (iterID = m_setiDatasets.begin(); iterID != m_setiDatasets.end(); ++iterID)
            for (i = 0; i < GetConditions(*iterID); ++i)
                ostm << '\t' << PCL.GetExperiment(GetCondition(*iterID, i));
        ostm << endl << "Motifs" << endl;
        for (iterMotif = GetMotifs().begin(); iterMotif != GetMotifs().end(); ++iterMotif)
            if (!(strMotif = iterMotif->Save(pMotifs, true, dCutoffPWMs, dPenaltyGap, dPenaltyMismatch,
                                             fNoRCs)).empty())
                ostm << strMotif << endl;
        ostm.flush();
    }

/*!
 * \brief
 * Opens a cluster based on the given PCL file and, if present, accompanying motif file.
 *
 * \param strPCL
 * Description of parameter strPCL.
 *
 * \param iSkip
 * Number of feature columns to skip between the gene IDs and first experimental column.
 *
 * \param PCL
 * PCL with which gene and condition labels are associated.
 *
 * \param pMotifs
 * If non-null, motif library with which cluster motifs are created.
 *
 * \returns
 * Number of successfully read cluster conditions, or -1 on failure.
 *
 * \remarks
 * Input PCL should be in the format generated by Save, i.e. only genes included in the bicluster, all
 * conditions, included conditions sorted to the left and tagged with an initial \c *.  For a given PCL
 * file \c &lt;filename>.pcl, motifs are read from \c &lt;filename>_motifs.txt.
 *
 * \see
 * Save
 */
    size_t CCoalesceCluster::Open(const std::string &strPCL, size_t iSkip, const CPCL &PCL,
                                  CCoalesceMotifLibrary *pMotifs) {
        CPCL PCLCluster;
        size_t i, j, iGene;
        string strMotifs;
        ifstream ifsm;

        Clear();
        if (!PCLCluster.Open(strPCL.c_str(), iSkip))
            return -1;

        for (i = 0; i < PCLCluster.GetExperiments(); ++i) {
            if (PCLCluster.GetExperiment(i)[0] != c_cStar)
                break;
            for (j = 0; j < PCL.GetExperiments(); ++j)
                if (PCL.GetExperiment(j) == (PCLCluster.GetExperiment(i).c_str() + 1)) {
                    m_setiDatasets.insert(j);
                    break;
                }
        }
        for (i = 0; i < PCLCluster.GetGenes(); ++i) {
            if ((iGene = PCL.GetGene(PCLCluster.GetGene(i))) == -1) {
                g_CatSleipnir().error("CCoalesceCluster::Open( %s, %i ) unrecognized gene: %s", strPCL.c_str(),
                                      iSkip, PCLCluster.GetGene(i).c_str());
                return -1;
            }
            m_setiGenes.insert(iGene);
        }

        if (!pMotifs || ((i = strPCL.rfind(CPCL::GetExtension())) !=
                         (strPCL.length() - strlen(CPCL::GetExtension()))))
            return PCLCluster.GetExperiments();
        strMotifs = strPCL.substr(0, i) + c_szMotifs;
        ifsm.open(strMotifs.c_str());
        if (!ifsm.is_open())
            return PCLCluster.GetExperiments();
        while (!ifsm.eof() && (ifsm.peek() != -1)) {
            SMotifMatch sMotif;

            if (!sMotif.Open(ifsm, *pMotifs)) {
                g_CatSleipnir().warn("CCoalesceCluster::Open( %s, %d ) could not open: %s", strPCL.c_str(),
                                     iSkip, strMotifs.c_str());
                return -1;
            }
            m_setsMotifs.insert(sMotif);
        }

        return PCLCluster.GetExperiments();
    }

/*!
 * \brief
 * Opens a cluster based on the textual representation (including genes, conditions, and motifs) in the given
 * input stream.
 *
 * \param istm
 * Input stream from which cluster is read.
 *
 * \param PCL
 * PCL with which gene and condition labels are associated.
 *
 * \param pMotifs
 * If non-null, motif library with which cluster motifs are created.
 *
 * \returns
 * True if cluster was opened successfully, false otherwise.
 *
 * \remarks
 * Opening can fail if the given stream contains gene or cluster IDs not present in the given PCL.
 * Motifs will be ignored if pMotifs is null.
 *
 * \see
 * Save
 */
    size_t CCoalesceCluster::Open(std::istream &istm, const CPCL &PCL, CCoalesceMotifLibrary *pMotifs) {
        string strBuffer;
        vector <string> vecstrLine;
        size_t i, j;
        vector <SMotifMatch> vecsMotifs;

        Clear();
        strBuffer.resize(CFile::GetBufferSize());
        istm.getline(&strBuffer[0], strBuffer.size());

        istm.getline(&strBuffer[0], strBuffer.size());
        CMeta::Tokenize(strBuffer.c_str(), vecstrLine);
        if (vecstrLine.empty() || (vecstrLine[0] != c_szGenes)) {
            g_CatSleipnir().error("CCoalesceCluster::Open( ) invalid line: %s", strBuffer.c_str());
            return -1;
        }
        for (i = 1; i < vecstrLine.size(); ++i) {
            if ((j = PCL.GetGene(vecstrLine[i])) == -1) {
                g_CatSleipnir().warn("CCoalesceCluster::Open( ) unrecognized gene: %s",
                                     vecstrLine[i].c_str());
                continue;
            }
            m_setiGenes.insert(j);
        }

        istm.getline(&strBuffer[0], strBuffer.size());
        vecstrLine.clear();
        CMeta::Tokenize(strBuffer.c_str(), vecstrLine);
        if (vecstrLine.empty() || (vecstrLine[0] != c_szConditions)) {
            g_CatSleipnir().error("CCoalesceCluster::Open( ) invalid line: %s", strBuffer.c_str());
            return -1;
        }
        for (i = 1; i < vecstrLine.size(); ++i) {
            if ((j = PCL.GetExperiment(vecstrLine[i])) == -1) {
                g_CatSleipnir().error("CCoalesceCluster::Open( ) unrecognized condition: %s",
                                      vecstrLine[i].c_str());
                return -1;
            }
            m_setiDatasets.insert(j);
        }

        istm.getline(&strBuffer[0], CFile::GetBufferSize());
        if (!CCoalesceMotifLibrary::Open(istm, vecsMotifs, pMotifs))
            return -1;
        for (i = 0; i < vecsMotifs.size(); ++i)
            m_setsMotifs.insert(vecsMotifs[i]);

        return m_setiDatasets.size();
    }

/*!
 * \brief
 * Creates a new cluster by merging the one or more clusters in the given hierarchy.
 *
 * \param Hierarchy
 * Hierarchy of cluster indices to be merged into the newly opened cluster.
 *
 * \param vecClusters
 * Vector of preexisting clusters, a superset of those in the given hierarchy.
 *
 * \param vecstrClusters
 * Vector of preexisting cluster IDs.
 *
 * \param dFraction
 * Minimum fraction of input cluster in which a gene must appear to be included in the output cluster.
 *
 * \param dCutoff
 * Edit distance threshhold beyond which merged motif alignments will be discarded.
 *
 * \param iCutoff
 * Maximum number of input motifs which will be merged exactly; larger numbers of motifs will be merged
 * heuristically.
 *
 * \param pMotifs
 * If non-null, motif library by which input and output motifs are managed.
 *
 * \returns
 * True of the new cluster was created successfully, false otherwise.
 *
 * \remarks
 * Currently only genes are filtered by dFraction; all conditions and motifs in any input cluster will be
 * included in the output cluster.
 *
 * \see
 * GetSimilarity | CClustHierarchical
 */
    bool CCoalesceCluster::Open(const CHierarchy &Hierarchy, const std::vector <CCoalesceCluster> &vecClusters,
                                const std::vector <std::string> &vecstrClusters, float dFraction, float dCutoff,
                                size_t iCutoff,
                                CCoalesceMotifLibrary *pMotifs) {
        map <size_t, size_t> mapiiGenes, mapiiDatasets;
        size_t i, iClusters;
        map<size_t, size_t>::const_iterator iterItem;
        vector < map < string, set < SMotifMatch > > > vecmapstrsetsMotifs;

        vecmapstrsetsMotifs.resize(CCoalesceSequencerBase::ESubsequenceEnd);
        g_CatSleipnir().notice("CCoalesceCluster::Open( %g ) merging clusters:", dFraction);
        if (!(iClusters = CCoalesceClusterImpl::Open(Hierarchy, vecClusters, vecstrClusters, mapiiGenes,
                                                     mapiiDatasets, vecmapstrsetsMotifs))) {
            g_CatSleipnir().error("CCoalesceCluster::Open( %g ) no clusters found", dFraction);
            return false;
        }

        Clear();
        for (iterItem = mapiiGenes.begin(); iterItem != mapiiGenes.end(); ++iterItem)
            if (((float) iterItem->second / iClusters) >= dFraction)
                m_setiGenes.insert(iterItem->first);
        for (iterItem = mapiiDatasets.begin(); iterItem != mapiiDatasets.end(); ++iterItem)
//		if( ( (float)iterItem->second / iClusters ) >= dFraction )
            m_setiDatasets.insert(iterItem->first);
        if (!pMotifs)
            return true;

        for (i = 0; i < vecmapstrsetsMotifs.size(); ++i) {
            const map <string, set<SMotifMatch>> &mapstrsetsMotifs = vecmapstrsetsMotifs[i];
            map < string, set < SMotifMatch > > ::const_iterator
            iterMotifs;

            for (iterMotifs = mapstrsetsMotifs.begin(); iterMotifs != mapstrsetsMotifs.end(); ++iterMotifs) {
                const set <SMotifMatch> &setsMotifs = iterMotifs->second;

                if (!((setsMotifs.size() < iCutoff) ?
                      CCoalesceClusterImpl::OpenMotifs(setsMotifs, *pMotifs, dCutoff) :
                      CCoalesceClusterImpl::OpenMotifsHeuristic(setsMotifs, *pMotifs, dCutoff, iCutoff)))
                    return false;
            }
        }

        return true;
    }

    bool CCoalesceClusterImpl::OpenMotifsHeuristic(const std::set <SMotifMatch> &setsMotifs,
                                                   CCoalesceMotifLibrary &Motifs, float dCutoff, size_t iCutoff) {
        vector <SMotifMatch> vecsMotifs;
        bool fDone;
        size_t i, iMotifs;
        set <SMotifMatch> setsMerged;

        iMotifs = setsMotifs.size();
        g_CatSleipnir().notice("CCoalesceClusterImpl::OpenMotifsHeuristic( %g ) resolving %d motifs", dCutoff,
                               iMotifs);

        vecsMotifs.resize(iMotifs);
        copy(setsMotifs.begin(), setsMotifs.end(), vecsMotifs.begin());
        do {
            fDone = true;
            sort(vecsMotifs.begin(), vecsMotifs.end());
            for (i = 0; (i + 1) < vecsMotifs.size(); ++i) {
                SMotifMatch &sOne = vecsMotifs[i];
                SMotifMatch &sTwo = vecsMotifs[i + 1];

                if ((sOne.m_iMotif == -1) || (sTwo.m_iMotif == -1))
                    break;
                if (Motifs.Align(sOne.m_iMotif, sTwo.m_iMotif, dCutoff) > dCutoff)
                    continue;
                if (sTwo.Open(sOne, sTwo, Motifs) == -1)
                    return false;
                if (Motifs.GetPST(sTwo.m_iMotif)->Integrate() > iCutoff)
                    Motifs.Simplify(sTwo.m_iMotif);
                fDone = false;
                iMotifs--;
                sOne.m_iMotif = -1;
            }
        } while (!fDone);

        for (i = 0; i < vecsMotifs.size(); ++i)
            if (vecsMotifs[i].m_iMotif != -1)
                setsMerged.insert(vecsMotifs[i]);

        return OpenMotifs(setsMerged, Motifs, dCutoff);
    }

    bool CCoalesceClusterImpl::OpenMotifs(const std::set <SMotifMatch> &setsMotifs, CCoalesceMotifLibrary &Motifs,
                                          float dCutoff) {
        vector <SMotifMatch> vecsMotifs;
        CDistanceMatrix MatSimilarity;
        CHierarchy *pHierMotifs;
        size_t i, j;
        bool fRet;

        g_CatSleipnir().notice("CCoalesceClusterImpl::OpenMotifs( %g ) resolving %d motifs", dCutoff,
                               setsMotifs.size());

        vecsMotifs.resize(setsMotifs.size());
        copy(setsMotifs.begin(), setsMotifs.end(), vecsMotifs.begin());
        MatSimilarity.Initialize(vecsMotifs.size());
        for (i = 0; i < vecsMotifs.size(); ++i)
            for (j = (i + 1); j < vecsMotifs.size(); ++j)
                MatSimilarity.Set(i, j, -Motifs.Align(vecsMotifs[i].m_iMotif,
                                                      vecsMotifs[j].m_iMotif, dCutoff));
        if (!(pHierMotifs = CClustHierarchical::Cluster(MatSimilarity)))
            return false;
        fRet = CCoalesceClusterImpl::OpenMotifs(Motifs, *pHierMotifs, vecsMotifs, dCutoff, m_setsMotifs);
        pHierMotifs->Destroy();
        return fRet;
    }

    size_t CCoalesceClusterImpl::Open(const CHierarchy &Hier, const std::vector <CCoalesceCluster> &vecClusters,
                                      const std::vector <std::string> &vecstrClusters,
                                      std::map <size_t, size_t> &mapiiGenes,
                                      std::map <size_t, size_t> &mapiiDatasets,
                                      TVecMapStrSetSMotifs &vecmapstrsetsMotifs) {
        set<size_t>::const_iterator iterFrom;
        map<size_t, size_t>::iterator iterTo;
        set<SMotifMatch>::const_iterator iterMotif;

        if (!Hier.IsGene())
            return (Open(Hier.Get(false), vecClusters, vecstrClusters, mapiiGenes, mapiiDatasets,
                         vecmapstrsetsMotifs) + Open(Hier.Get(true), vecClusters, vecstrClusters, mapiiGenes,
                                                     mapiiDatasets, vecmapstrsetsMotifs));

        const CCoalesceCluster &Cluster = vecClusters[Hier.GetID()];

        g_CatSleipnir().notice("CCoalesceClusterImpl::Open( ) cluster %s",
                               vecstrClusters[Hier.GetID()].c_str());
        for (iterFrom = Cluster.GetGenes().begin(); iterFrom != Cluster.GetGenes().end(); ++iterFrom)
            if ((iterTo = mapiiGenes.find(*iterFrom)) == mapiiGenes.end())
                mapiiGenes[*iterFrom] = 1;
            else
                iterTo->second++;
        for (iterFrom = Cluster.m_setiDatasets.begin(); iterFrom != Cluster.m_setiDatasets.end(); ++iterFrom)
            if ((iterTo = mapiiDatasets.find(*iterFrom)) == mapiiDatasets.end())
                mapiiDatasets[*iterFrom] = 1;
            else
                iterTo->second++;
        for (iterMotif = Cluster.m_setsMotifs.begin(); iterMotif != Cluster.m_setsMotifs.end(); ++iterMotif)
            vecmapstrsetsMotifs[iterMotif->m_eSubsequence][iterMotif->m_strType].insert(*iterMotif);

        return 1;
    }

    bool CCoalesceClusterImpl::OpenMotifs(CCoalesceMotifLibrary &Motifs, const CHierarchy &Hier,
                                          const std::vector <SMotifMatch> &vecsMotifs, float dCutoff,
                                          std::set <SMotifMatch> &setsMotifs) {

        if (Hier.IsGene() || (-Hier.GetSimilarity() < dCutoff)) {
            SMotifMatch sMotif;
            size_t iCount;

            sMotif.m_dZ = 0;
            if (sMotif.Open(Hier, vecsMotifs, Motifs, iCount = 0) == -1)
                return false;
            sMotif.m_dZ /= iCount;
            setsMotifs.insert(sMotif);
            return true;
        }

        return (OpenMotifs(Motifs, Hier.Get(false), vecsMotifs, dCutoff, setsMotifs) &&
                OpenMotifs(Motifs, Hier.Get(true), vecsMotifs, dCutoff, setsMotifs));
    }

/*!
 * \brief
 * Calculates a similarity score between the given and current clusters.
 *
 * \param Cluster
 * Cluster to be compared to the current cluster.
 *
 * \param iGenes
 * Total number of available genes.
 *
 * \param iDatasets
 * Total number of available datasets.
 *
 * \returns
 * Similarity score for the two clusters (minimum zero, increasing indicates greater similarity).
 *
 * \remarks
 * Used by COALESCE to postprocess modules.  Generally returns fraction of shared genes proportional
 * to the smaller of the two clusters; Jaccard index can also work well.
 */
    float CCoalesceCluster::GetSimilarity(const CCoalesceCluster &Cluster, size_t iGenes,
                                          size_t iDatasets) const {
        size_t iOverlapGenes;
        set<size_t>::const_iterator iterItem;
        float dRet;

        if (GetGenes().empty() || Cluster.GetGenes().empty())
            return 0;
        for (iOverlapGenes = 0, iterItem = GetGenes().begin(); iterItem != GetGenes().end(); ++iterItem)
            if (Cluster.IsGene(*iterItem))
                iOverlapGenes++;
// simple gene overlap works best with an ~0.75 cutoff
        dRet = (float) iOverlapGenes / min(GetGenes().size(), Cluster.GetGenes().size());
// jaccard index works best with an 0.25-0.33 cutoff
//	dRet = (float)iOverlapGenes / ( GetGenes( ).size( ) + Cluster.GetGenes( ).size( ) - iOverlapGenes );
        g_CatSleipnir().debug("CCoalesceCluster::GetSimilarity( %d, %d ): %d, %d, %d = %g", iGenes,
                              iDatasets, iOverlapGenes, GetGenes().size(), Cluster.GetGenes().size(), dRet);
        return dRet;
    }

/*!
 * \brief
 * Records the state of the cluster between convergence iterations.
 *
 * \param GeneScores
 * Per-gene sequence scores used to snapshot motif histograms.
 *
 * \param Histograms
 * Motif histograms to be updated using the current gene scores.
 *
 * \remarks
 * Should be called at the start of each COALESCE iteration, after new histograms are calculated but
 * before anything else.
 *
 * \see
 * IsConverged
 */
    void CCoalesceCluster::Snapshot(const CCoalesceGeneScores &GeneScores,
                                    CCoalesceGroupHistograms &Histograms) {

        Histograms.SetTotal(GeneScores, GetGenes());
        m_setiHistory.insert(GetHash());
        CCoalesceClusterImpl::Snapshot(m_setiDatasets, m_veciPrevDatasets);
        CCoalesceClusterImpl::Snapshot(m_setsMotifs, m_vecsPrevMotifs);
        CCoalesceClusterImpl::Snapshot(GetGenes(), m_veciPrevGenes);
    }

/*!
 * \brief
 * Labels all significant motifs in the cluster with any known TF motifs that match below the given
 * threshhold.
 *
 * \param Motifs
 * Motif library to be used for cluster and known motifs.
 *
 * \param eMatchType
 * Type of match to be performed: correlation, rmse, etc.
 *
 * \param dPenaltyGap
 * Alignment score penalty for gaps.
 *
 * \param dPenaltyMismatch
 * Alignment score penalty for mismatches.
 *
 * \param dPValue
 * P-value (or other score) threshhold below which known TFs must match.
 *
 * \returns
 * True if the labeling process succeeded (possibly with no matching labels), false otherwise.
 *
 * Labels any significant motifs in the cluster with zero or more matching known TF binding motifs.
 * These known labels are scored and will be read and written with subsequence Open/Save calls that
 * include this cluster's motifs.
 *
 * \remarks
 * For each significant motif in the cluster, scans any known TF motifs in Motifs for matches based
 * on one of several PWM matching algorithms (usually correlation p-value).
 *
 * \see
 * CCoalesceMotifLibrary::GetKnown
 */
    bool CCoalesceCluster::LabelMotifs(const CCoalesceMotifLibrary &Motifs, SMotifMatch::EType eMatchType,
                                       float dPenaltyGap, float dPenaltyMismatch, float dPValue) {
        set<SMotifMatch>::iterator iterMotif;

        if (!Motifs.GetKnowns())
            return true;
        for (iterMotif = m_setsMotifs.begin(); iterMotif != m_setsMotifs.end(); ++iterMotif) {
            SMotifMatch &sMotif = (SMotifMatch &) *iterMotif;
// For some obscure reason, Linux won't let me use this as a non-const iterator...
            if (!sMotif.Label(Motifs, eMatchType, dPenaltyGap, dPenaltyMismatch, dPValue))
                return false;
        }

        return true;
    }

}

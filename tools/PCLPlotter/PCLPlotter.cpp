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

static const char c_acNucleotides[] = "ACGT";

int main(int iArgs, char **aszArgs) {
    gengetopt_args_info sArgs;
    CPCL PCL, PCLBackground;
    vector<float> vecdSumsIn, vecdSumSqsIn, vecdSumsOut, vecdSumSqsOut;
    vector <size_t> veciCountsIn, veciCountsOut;
    size_t i, j, k;
    float d;
    set <size_t> setiGenes;

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }
    CMeta Meta(sArgs.verbosity_arg);

    if (!PCL.Open(sArgs.input_arg, sArgs.skip_arg)) {
        cerr << "Could not open: " << (sArgs.input_arg ? sArgs.input_arg : "stdin") << endl;
        return 1;
    }
    if (sArgs.background_arg && !PCLBackground.Open(sArgs.background_arg, sArgs.skip_arg)) {
        cerr << "Could not open: " << sArgs.background_arg << endl;
        return 1;
    }

    if (sArgs.fasta_arg) {
        CFASTA FASTA;
        map <string, size_t> mapstriTypes;
        map<string, size_t>::const_iterator iterType;
        size_t iRow, iColumn, iLength, iCount;
        set <size_t> setiForeground;
        CFullMatrix <size_t> MatSums, MatSumSqs, MatCounts;
        CFullMatrix <CHMM> MatHMMs;
        float dAve, dStd;

        if (!FASTA.Open(sArgs.fasta_arg)) {
            cerr << "Could not open: " << sArgs.fasta_arg << endl;
            return 1;
        }

        if (sArgs.motifs_arg) {
            ifstream ifsm;
            vector <SMotifMatch> vecsMotifs;
            CCoalesceMotifLibrary Motifs(sArgs.k_arg);
            SCoalesceModifiers sModifiers;
            SCoalesceModifierCache sCache(sModifiers);
            CPCL PCLScores;
            vector <string> vecstrMotifs;

            ifsm.open(sArgs.motifs_arg);
            if (!CCoalesceMotifLibrary::Open(ifsm, vecsMotifs, &Motifs)) {
                cerr << "Could not open: " << sArgs.motifs_arg << endl;
                return 1;
            }
            ifsm.close();

            vecstrMotifs.resize(vecsMotifs.size());
            for (i = 0; i < vecstrMotifs.size(); ++i)
                vecstrMotifs[i] = Motifs.GetMotif(vecsMotifs[i].m_iMotif);
            PCLScores.Open(PCL.GetGeneNames(), vecstrMotifs, vector<string>());
            PCLScores.Clear();

            for (iRow = 0; iRow < PCL.GetGenes(); ++iRow) {
                vector <SFASTASequence> vecsSequences;

                if (!(iRow % 100))
                    cerr << "Processing " << iRow << '/' << PCL.GetGenes() << endl;
                if ((i = FASTA.GetGene(PCL.GetGene(iRow))) == -1) {
                    for (i = 0; i < PCLScores.GetExperiments(); ++i)
                        PCLScores.Set(iRow, i, CMeta::GetNaN());
                    continue;
                }
                if (!FASTA.Get(i, vecsSequences))
                    return 1;
                for (iColumn = 0; iColumn < vecsMotifs.size(); ++iColumn) {
                    const SMotifMatch &sMotif = vecsMotifs[iColumn];

                    for (iLength = i = 0; i < vecsSequences.size(); ++i) {
                        const SFASTASequence &sSequence = vecsSequences[i];

                        if (sMotif.m_strType != sSequence.m_strType)
                            continue;
                        for (j = 0; j < sSequence.m_vecstrSequences.size(); ++j) {
                            const string &strSequence = sSequence.m_vecstrSequences[j];

                            if ((sMotif.m_eSubsequence != CCoalesceSequencerBase::ESubsequenceTotal) &&
                                (sMotif.m_eSubsequence != ((sSequence.m_fIntronFirst == !(j % 2)) ?
                                                           CCoalesceSequencerBase::ESubsequenceIntrons :
                                                           CCoalesceSequencerBase::ESubsequenceExons)))
                                continue;
                            iLength += strSequence.size();
                            PCLScores.Get(iRow, iColumn) += Motifs.GetMatch(strSequence,
                                                                            sMotif.m_iMotif, 0, sCache);
                        }
                    }
                    if (iLength)
                        PCLScores.Get(iRow, iColumn) /= iLength;
                }
            }

            PCLScores.Save(cout);
            return 0;
        }

        for (i = 0; i < FASTA.GetGenes(); ++i)
            if (PCL.GetGene(FASTA.GetGene(i)) != -1) {
                setiGenes.insert(i);
                setiForeground.insert(i);
            } else if (PCLBackground.GetGenes()) {
                if (PCLBackground.GetGene(FASTA.GetGene(i)) != -1)
                    setiGenes.insert(i);
            } else
                setiGenes.insert(i);
        MatSums.Initialize(FASTA.GetTypes().size(), 2);
        MatSums.Clear();
        MatSumSqs.Initialize(FASTA.GetTypes().size(), 2);
        MatSumSqs.Clear();
        MatCounts.Initialize(FASTA.GetTypes().size(), 2);
        MatCounts.Clear();
        MatHMMs.Initialize(FASTA.GetTypes().size(), 2);
        for (i = 0; i < MatHMMs.GetRows(); ++i)
            for (j = 0; j < MatHMMs.GetColumns(); ++j)
                MatHMMs.Get(i, j).Open(sArgs.degree_arg, c_acNucleotides);
        for (i = 0; i < FASTA.GetGenes(); ++i) {
            vector <SFASTASequence> vecsSequences;

            if (setiGenes.find(i) == setiGenes.end())
                continue;
            if (!FASTA.Get(i, vecsSequences))
                return 1;
            iColumn = (setiForeground.find(i) == setiForeground.end()) ? 1 : 0;
            for (j = 0; j < vecsSequences.size(); ++j) {
                const SFASTASequence &sSequence = vecsSequences[j];

                if ((iterType = mapstriTypes.find(sSequence.m_strType)) == mapstriTypes.end()) {
                    iRow = mapstriTypes.size();
                    mapstriTypes[sSequence.m_strType] = iRow;
                } else
                    iRow = iterType->second;
                MatCounts.Get(iRow, iColumn)++;
                for (iLength = k = 0; k < sSequence.m_vecstrSequences.size(); ++k) {
                    const string &strSequence = sSequence.m_vecstrSequences[k];

                    iLength += strSequence.length();
                    MatHMMs.Get(iRow, iColumn).Add(strSequence);
                }
                MatSums.Get(iRow, iColumn) += iLength;
                MatSumSqs.Get(iRow, iColumn) += iLength * iLength;
            }
        }

        for (iterType = mapstriTypes.begin(); iterType != mapstriTypes.end(); ++iterType) {
            cout << iterType->first;
            for (i = 0; i < MatCounts.GetColumns(); ++i) {
                iCount = max((size_t) 1, MatCounts.Get(iterType->second, i));
                dAve = (float) MatSums.Get(iterType->second, i) / iCount;
                dStd = sqrt(((float) MatSumSqs.Get(iterType->second, i) / iCount) - (dAve * dAve));
                cout << '\t' << dAve << '\t' << dStd;
            }
            cout << endl;

            for (i = 0; i < MatHMMs.GetColumns(); ++i)
                MatHMMs.Get(iterType->second, i).Save(cout);
        }
        return 0;
    }

    if (sArgs.genes_arg) {
        CGenome Genome;
        CGenes Genes(Genome);

        if (!Genes.Open(sArgs.genes_arg)) {
            cerr << "Could not open: " << sArgs.genes_arg << endl;
            return 1;
        }
        for (i = 0; i < Genes.GetGenes(); ++i)
            if ((j = PCL.GetGene(Genes.GetGene(i).GetName())) != -1)
                setiGenes.insert(j);
    }

    vecdSumsIn.resize(PCL.GetExperiments());
    vecdSumsOut.resize(PCL.GetExperiments());
    vecdSumSqsIn.resize(PCL.GetExperiments());
    vecdSumSqsOut.resize(PCL.GetExperiments());
    veciCountsIn.resize(PCL.GetExperiments());
    veciCountsOut.resize(PCL.GetExperiments());
    for (i = 0; i < PCL.GetGenes(); ++i) {
        bool fIn = sArgs.background_arg ||
                   (sArgs.genes_arg && (setiGenes.find(i) != setiGenes.end())) ||
                   (PCL.GetFeatures() && (PCL.GetFeature(i, 1)[0] == '*'));
        vector<float> &vecdSums = fIn ? vecdSumsIn : vecdSumsOut;
        vector<float> &vecdSumSqs = fIn ? vecdSumSqsIn : vecdSumSqsOut;
        vector <size_t> &veciCounts = fIn ? veciCountsIn : veciCountsOut;

        for (j = 0; j < PCL.GetExperiments(); ++j)
            if (!CMeta::IsNaN(d = PCL.Get(i, j))) {
                veciCounts[j]++;
                vecdSums[j] += d;
                vecdSumSqs[j] += d * d;
            }
    }

    if (sArgs.background_arg) {
        vector <size_t> veciExperiments;
        map <string, size_t> mapstriExperiments;
        map<string, size_t>::const_iterator iterExperiment;

        if (PCLBackground.GetExperiments() != PCL.GetExperiments()) {
            cerr << "Inconsistent experiments: " << PCL.GetExperiments() << " vs. " <<
                 PCLBackground.GetExperiments() << endl;
            return 1;
        }

        setiGenes.clear();
        for (i = 0; i < PCL.GetGenes(); ++i)
            if ((j = PCLBackground.GetGene(PCL.GetGene(i))) != -1)
                setiGenes.insert(j);
        for (i = 0; i < PCL.GetExperiments(); ++i)
            mapstriExperiments[PCLBackground.GetExperiment(i)] = i;
        veciExperiments.resize(PCL.GetExperiments());
        for (i = 0; i < veciExperiments.size(); ++i) {
            const string &strExperiment = PCL.GetExperiment(i);

            if ((iterExperiment = mapstriExperiments.find((strExperiment[0] == '*') ?
                                                          strExperiment.substr(1) : strExperiment)) ==
                mapstriExperiments.end()) {
                cerr << "Could not match experiment " << i << ": " << strExperiment << endl;
                return 1;
            }
            veciExperiments[i] = iterExperiment->second;
        }

        for (i = 0; i < PCLBackground.GetGenes(); ++i) {
            if (setiGenes.find(i) != setiGenes.end())
                continue;
            for (j = 0; j < veciExperiments.size(); ++j)
                if (!CMeta::IsNaN(d = PCLBackground.Get(i, veciExperiments[j]))) {
                    veciCountsOut[j]++;
                    vecdSumsOut[j] += d;
                    vecdSumSqsOut[j] += d * d;
                }
        }
    }

    for (i = 0; i < PCL.GetExperiments(); ++i) {
        if ((j = veciCountsIn[i])) {
            d = (vecdSumsIn[i] /= j);
            vecdSumSqsIn[i] = sqrt((vecdSumSqsIn[i] / j) - (d * d));
        }
        if ((j = veciCountsOut[i])) {
            d = (vecdSumsOut[i] /= j);
            vecdSumSqsOut[i] = sqrt((vecdSumSqsOut[i] / j) - (d * d));
        }
    }

    for (i = 0; i < PCL.GetExperiments(); ++i)
        cout << (i ? "\t" : "") << PCL.GetExperiment(i);
    cout << endl;
    for (i = 0; i < PCL.GetExperiments(); ++i)
        cout << (i ? "\t" : "") << vecdSumsIn[i];
    cout << endl;
    for (i = 0; i < PCL.GetExperiments(); ++i)
        cout << (i ? "\t" : "") << vecdSumSqsIn[i];
    cout << endl;
    for (i = 0; i < PCL.GetExperiments(); ++i)
        cout << (i ? "\t" : "") << vecdSumsOut[i];
    cout << endl;
    for (i = 0; i < PCL.GetExperiments(); ++i)
        cout << (i ? "\t" : "") << vecdSumSqsOut[i];
    cout << endl;


    return 0;
}

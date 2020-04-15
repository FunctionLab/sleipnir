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
#include <cmath>

#include "statistics.h"
#include "datapair.h"

int main(int iArgs, char **aszArgs) {
    gengetopt_args_info sArgs;
    CGenome Genome;
    CGenes Genes(Genome);
    ifstream ifsm;
    CDat Dat;
    size_t i, j, k;
    bool fModified;

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }

    CMeta Meta(sArgs.verbosity_arg, sArgs.random_arg);

    if (sArgs.genes_arg) {
        ifsm.open(sArgs.genes_arg);
        if (!Genes.Open(ifsm)) {
            cerr << "Could not open: " << sArgs.genes_arg << endl;
            return 1;
        }
        ifsm.close();
    }

    fModified = sArgs.normalize_flag || sArgs.subsample_arg;
    if (sArgs.input_arg) {
        if (!Dat.Open(sArgs.input_arg, sArgs.memmap_flag && !fModified,
                      sArgs.skip_arg, !!sArgs.zscore_flag, !!sArgs.duplicates_flag)) {
            cerr << "Could not open: " << sArgs.input_arg << endl;
            return 1;
        }
    } else if (!Dat.Open(cin, CDat::EFormatText, (float) HUGE_VAL, !!sArgs.duplicates_flag)) {
        cerr << "Could not open input" << endl;
        return 1;
    }

    if (sArgs.genelist_flag) {
        for (i = 0; i < Dat.GetGenes(); ++i)
            cout << Dat.GetGene(i) << endl;
        return 0;
    }

    if (sArgs.remap_arg) {
        static const size_t c_iBuffer = 1024;
        char acBuffer[c_iBuffer];
        vector <string> vecstrTokens;
        map <string, string> mapNames;
        map<string, string>::const_iterator iterName;

        ifsm.clear();
        ifsm.open(sArgs.remap_arg);
        if (!ifsm.is_open()) {
            cerr << "Could not open: " << sArgs.remap_arg << endl;
            return 1;
        }
        while (!ifsm.eof()) {
            ifsm.getline(acBuffer, c_iBuffer - 1);
            acBuffer[c_iBuffer - 1] = 0;
            vecstrTokens.clear();
            CMeta::Tokenize(acBuffer, vecstrTokens);
            if (vecstrTokens.empty())
                continue;
            if (vecstrTokens.size() < 2) {
                cerr << "Illegal remap line (" << vecstrTokens.size() << "): " << acBuffer << endl;
                return 1;
            }
            if (vecstrTokens[0] == vecstrTokens[1])
                continue;
            if ((iterName = mapNames.find(vecstrTokens[0])) == mapNames.end())
                mapNames[vecstrTokens[0]] = vecstrTokens[1];
            else if (iterName->second != vecstrTokens[1]) {
                cerr << "Ambiguous mapping: " << vecstrTokens[0] << " to " << iterName->second <<
                     ", " << vecstrTokens[1] << endl;
                return 1;
            }
        }

        for (i = 0; i < Dat.GetGenes(); ++i)
            if ((iterName = mapNames.find(Dat.GetGene(i))) != mapNames.end()) {
                if (Dat.GetGene(iterName->second) != -1) {
                    cerr << "Duplicate mapping: " << Dat.GetGene(i) << " to " <<
                         iterName->second << endl;
                    return 1;
                }
                Dat.SetGene(i, iterName->second);
            }
    }

    // should I set all non-missing values to defualt value?
    if (sArgs.dval_given) {
        for (i = 0; i < Dat.GetGenes(); ++i)
            for (j = (i + 1); j < Dat.GetGenes(); ++j) {
                if (CMeta::IsNaN(Dat.Get(i, j)))
                    continue;
                Dat.Set(i, j, sArgs.dval_arg);
            }
    }
    // should I add random noise from standard Normal?
    if (sArgs.noise_flag) {
        float d;
        for (i = 0; i < Dat.GetGenes(); ++i)
            for (j = (i + 1); j < Dat.GetGenes(); ++j)
                if (!CMeta::IsNaN(d = Dat.Get(i, j))) {
                    Dat.Set(i, j, d + CStatistics::SampleNormalStandard());
                }
    }
    // Rescale prior
    if (sArgs.prior_given && sArgs.newprior_given) {
        float r = sArgs.newprior_arg / sArgs.prior_arg;
        float rdiff = (1 - sArgs.newprior_arg) / (1 - sArgs.prior_arg);
        float d;
        for (i = 0; i < Dat.GetGenes(); ++i) {
            for (j = (i + 1); j < Dat.GetGenes(); ++j) {
                if (!CMeta::IsNaN(d = Dat.Get(i, j))) {
                    d = (d * r) / (d * r + (1 - d) * rdiff);
                    Dat.Set(i, j, d);
                }
            }
        }
    }
    // Convert posterior to log-likelihood ratio
    if (sArgs.prior_given && sArgs.logratio_flag) {
        float logprior = log(sArgs.prior_arg / (1 - sArgs.prior_arg));
        float d;
        for (i = 0; i < Dat.GetGenes(); ++i) {
            for (j = (i + 1); j < Dat.GetGenes(); ++j) {
                if (!CMeta::IsNaN(d = Dat.Get(i, j))) {
                    d = log(d / (1 - d)) - logprior;
                    Dat.Set(i, j, d);
                }
            }
        }

    }

    if (sArgs.randomize_flag)
        Dat.Randomize();
    if (sArgs.rank_flag)
        Dat.Rank();
    if (sArgs.normalize_flag || sArgs.zscore_flag)
        Dat.Normalize(sArgs.zscore_flag ? CDat::ENormalizeZScore : CDat::ENormalizeMinMax);

    if (sArgs.normalizeNPone_flag)
        Dat.Normalize(CDat::ENormalizeMinMaxNPone);

    if (sArgs.zero_flag || sArgs.dmissing_given) {
        for (i = 0; i < Dat.GetGenes(); ++i) {
            for (j = (i + 1); j < Dat.GetGenes(); ++j) {
                if (CMeta::IsNaN(Dat.Get(i, j))) {
                    if (sArgs.zero_flag) {
                        Dat.Set(i, j, 0);
                    } else {
                        Dat.Set(i, j, sArgs.dmissing_arg);
                    }
                }
            }
        }
    }

    if (sArgs.NegExp_flag) {
        float d;
        for (i = 0; i < Dat.GetGenes(); ++i) {
            for (j = (i + 1); j < Dat.GetGenes(); ++j) {
                if (CMeta::IsNaN(d = Dat.Get(i, j)))
                    continue;
                Dat.Set(i, j, exp(-d));
            }
        }
    }

    if (sArgs.abs_flag) {
        float d;
        for (i = 0; i < Dat.GetGenes(); ++i) {
            for (j = (i + 1); j < Dat.GetGenes(); ++j) {
                if (CMeta::IsNaN(d = Dat.Get(i, j)))
                    continue;
                Dat.Set(i, j, abs(d));
            }
        }
    }

    if (sArgs.flip_flag)
        Dat.Invert();

    if (sArgs.normalizeDeg_flag) {
        size_t iCutoff;
        float d;
        vector <size_t> veciCounts;
        vector<float> vecdTotals;

        veciCounts.resize(Dat.GetGenes());
        fill(veciCounts.begin(), veciCounts.end(), 0);
        vecdTotals.resize(Dat.GetGenes());
        fill(vecdTotals.begin(), vecdTotals.end(), 0.0f);

        for (iCutoff = i = 0; i < Dat.GetGenes(); ++i) {
            for (j = (i + 1); j < Dat.GetGenes(); ++j) {
                if (!CMeta::IsNaN(d = Dat.Get(i, j))) {
                    if (!sArgs.cutoff_given || (d >= sArgs.cutoff_arg)) {
                        // d = abs (d);
                        iCutoff++;
                        veciCounts[i]++;
                        veciCounts[j]++;
                        vecdTotals[i] += d;
                        vecdTotals[j] += d;
                        d *= d;
                    }
                }
            }
        }

        for (i = 0; i < Dat.GetGenes(); ++i) {
            for (j = (i + 1); j < Dat.GetGenes(); ++j) {
                if (!CMeta::IsNaN(d = Dat.Get(i, j))) {
                    if (!sArgs.cutoff_given || (d >= sArgs.cutoff_arg)) {
                        d = d * veciCounts[i] / sqrt(vecdTotals[i] * vecdTotals[j]);
                        Dat.Set(i, j, d);
                    }
                }
            }
        }
    }

    if (sArgs.normalizeLoc_flag) {
        size_t iCutoff;
        float d, zi, zj;
        vector <size_t> veciCounts;
        vector<float> vecdTotals, vecdAvgs, vecdSquares;

        veciCounts.resize(Dat.GetGenes());
        fill(veciCounts.begin(), veciCounts.end(), 0);
        vecdTotals.resize(Dat.GetGenes());
        fill(vecdTotals.begin(), vecdTotals.end(), 0.0f);
        vecdAvgs.resize(Dat.GetGenes());
        fill(vecdAvgs.begin(), vecdAvgs.end(), 0.0f);
        vecdSquares.resize(Dat.GetGenes());
        fill(vecdSquares.begin(), vecdSquares.end(), 0.0f);

        for (iCutoff = i = 0; i < Dat.GetGenes(); ++i) {
            for (j = (i + 1); j < Dat.GetGenes(); ++j) {
                if (!CMeta::IsNaN(d = Dat.Get(i, j))) {
                    if (!sArgs.cutoff_given || (d >= sArgs.cutoff_arg)) {
                        // d = abs (d);
                        iCutoff++;
                        veciCounts[i]++;
                        veciCounts[j]++;
                        vecdTotals[i] += d;
                        vecdTotals[j] += d;

                        d *= d;
                        vecdSquares[i] += d;
                        vecdSquares[j] += d;
                    }
                }
            }
        }

        for (i = 0; i < vecdSquares.size(); ++i) {
            d = vecdTotals[i] / veciCounts[i];
            vecdAvgs[i] = d;
            vecdSquares[i] = sqrt((vecdSquares[i] / veciCounts[i]) - (d * d));
        }

        for (i = 0; i < Dat.GetGenes(); ++i) {
            for (j = (i + 1); j < Dat.GetGenes(); ++j) {
                if (!CMeta::IsNaN(d = Dat.Get(i, j))) {
                    if (!sArgs.cutoff_given || (d >= sArgs.cutoff_arg)) {
                        if (vecdSquares[i] == 0) { zi = 0; }
                        else { zi = (d - vecdAvgs[i]) / vecdSquares[i]; }

                        if (vecdSquares[j] == 0) { zj = 0; }
                        else { zj = (d - vecdAvgs[j]) / vecdSquares[j]; }

                        d = (zi + zj) / sqrt(2);
                        Dat.Set(i, j, d);
                    }
                }
            }
        }
    }

    if (Genes.GetGenes())
        Dat.FilterGenes(Genes, CDat::EFilterInclude);
    if (sArgs.genex_arg)
        Dat.FilterGenes(sArgs.genex_arg, CDat::EFilterExclude);
    if (sArgs.genee_arg)
        Dat.FilterGenes(sArgs.genee_arg, CDat::EFilterEdge);
    if (sArgs.gexedges_arg)
        Dat.FilterGenes(sArgs.gexedges_arg, CDat::EFilterExEdge);


    if (sArgs.ccoeff_flag) {
        float sxy, sxxyy, sgxy;
        for (i = 0; i < Dat.GetGenes(); ++i) {
            sxy = 0.0;
            sxxyy = 0.0;
            sgxy = 0.0;
            for (j = 0; j < Dat.GetGenes(); ++j) {
                if (i == j)
                    continue;

                sxy = sxy + Dat.Get(i, j);
                sxxyy = sxxyy + (Dat.Get(i, j) * Dat.Get(i, j));

                for (k = j + 1; k < Dat.GetGenes(); ++k) {
                    if (i == k || j == k)
                        continue;
                    sgxy = sgxy + Dat.Get(i, j) * Dat.Get(i, k) * Dat.Get(j, k);
                }
            }
            cout << Dat.GetGene(i) << '\t' << sgxy / (sxy * sxy - sxxyy) << endl;
        }

        return 0;
    }

    if (sArgs.mar_flag) {
        float sxy, sxxyy;

        for (i = 0; i < Dat.GetGenes(); ++i) {
            sxy = 0;
            sxxyy = 0;

            for (j = 0; j < Dat.GetGenes(); ++j) {
                if (i == j)
                    continue;

                sxxyy = sxxyy + (Dat.Get(i, j) * Dat.Get(i, j));
                sxy = sxy + Dat.Get(i, j);
            }
            cout << Dat.GetGene(i) << '\t' << sxxyy / sxy << endl;
        }
        return 0;
    }

    if (sArgs.hubbiness_flag) {
        float sume;
        for (i = 0; i < Dat.GetGenes(); ++i) {
            sume = 0.0;
            for (j = 0; j < Dat.GetGenes(); ++j) {
                if (i == j)
                    continue;
                sume = sume + Dat.Get(i, j);
            }

            cout << Dat.GetGene(i) << '\t' << sume / (Dat.GetGenes() - 1) << endl;
        }

        return 0;
    }

    if (sArgs.paircount_flag) {
        size_t iTotal, iCutoff;
        float d, dAve, dStd;
        vector <size_t> veciCounts;
        vector<float> vecdTotals, vecdSquares;

        dAve = dStd = 0;
        veciCounts.resize(Dat.GetGenes());
        fill(veciCounts.begin(), veciCounts.end(), 0);
        vecdTotals.resize(Dat.GetGenes());
        fill(vecdTotals.begin(), vecdTotals.end(), 0.0f);
        vecdSquares.resize(Dat.GetGenes());
        fill(vecdSquares.begin(), vecdSquares.end(), 0.0f);
        for (iTotal = iCutoff = i = 0; i < Dat.GetGenes(); ++i)
            for (j = (i + 1); j < Dat.GetGenes(); ++j)
                if (!CMeta::IsNaN(d = Dat.Get(i, j))) {
                    if (!sArgs.cutoff_given || (d >= sArgs.cutoff_arg)) {
                        // d = abs (d);
                        dAve += d;
                        iCutoff++;
                        veciCounts[i]++;
                        veciCounts[j]++;
                        vecdTotals[i] += d;
                        vecdTotals[j] += d;
                        d *= d;
                        dStd += d;
                        vecdSquares[i] += d;
                        vecdSquares[j] += d;
                    }
                    iTotal++;
                }
        dAve /= iCutoff;
        dStd = sqrt((dStd / iCutoff) - (dAve * dAve));
        for (i = 0; i < vecdSquares.size(); ++i) {
            d = vecdTotals[i] / veciCounts[i];
            vecdSquares[i] = sqrt((vecdSquares[i] / veciCounts[i]) - (d * d));
        }

        cout << iTotal << endl;
        if (sArgs.cutoff_given)
            cout << iCutoff << endl;
        cout << dAve << '\t' << dStd << endl;
        for (i = 0; i < Dat.GetGenes(); ++i)
            cout << Dat.GetGene(i) << '\t' << vecdTotals[i] << '\t' << veciCounts[i] << '\t' <<
                 vecdSquares[i] << endl;
        return 0;
    }

    if (sArgs.cutoff_given)
        for (i = 0; i < Dat.GetGenes(); ++i)
            for (j = (i + 1); j < Dat.GetGenes(); ++j)
                if (Dat.Get(i, j) < sArgs.cutoff_arg)
                    Dat.Set(i, j, CMeta::GetNaN());

    if (sArgs.subsample_arg < 1)
        for (i = 0; i < Dat.GetGenes(); ++i)
            for (j = (i + 1); j < Dat.GetGenes(); ++j)
                if (!CMeta::IsNaN(Dat.Get(i, j)) &&
                    (((float) rand() / RAND_MAX) > sArgs.subsample_arg))
                    Dat.Set(i, j, CMeta::GetNaN());

    if (sArgs.edges_arg) {
        CDat DatLk1;
        vector <size_t> veciGenesOne;
        size_t iOne, iTwo;

        if (!DatLk1.Open(sArgs.edges_arg)) {
            cerr << "Could not open: " << sArgs.edges_arg << endl;
            return 1;
        }
        veciGenesOne.resize(Dat.GetGenes());
        for (i = 0; i < veciGenesOne.size(); ++i)
            veciGenesOne[i] = DatLk1.GetGene(Dat.GetGene(i));
        for (i = 0; i < Dat.GetGenes(); ++i) {
            if ((iOne = veciGenesOne[i]) == -1) {
                for (j = (i + 1); j < Dat.GetGenes(); ++j)
                    Dat.Set(i, j, CMeta::GetNaN());
                continue;
            }
            for (j = (i + 1); j < Dat.GetGenes(); ++j)
                if (((iTwo = veciGenesOne[j]) == -1) ||
                    CMeta::IsNaN(DatLk1.Get(iOne, iTwo)))
                    Dat.Set(i, j, CMeta::GetNaN());
        }
    }

    if (sArgs.exedges_arg) {
        CDat DatLk1;
        vector <size_t> veciGenesOne;
        size_t iOne, iTwo;

        if (!DatLk1.Open(sArgs.exedges_arg)) {
            cerr << "Could not open: " << sArgs.exedges_arg << endl;
            return 1;
        }

        veciGenesOne.resize(Dat.GetGenes());

        for (i = 0; i < veciGenesOne.size(); ++i)
            veciGenesOne[i] = DatLk1.GetGene(Dat.GetGene(i));

        for (i = 0; i < Dat.GetGenes(); ++i) {
            if ((iOne = veciGenesOne[i]) == -1)
                continue;
            for (j = (i + 1); j < Dat.GetGenes(); ++j) {
                if (((iTwo = veciGenesOne[j]) == -1) ||
                    CMeta::IsNaN(DatLk1.Get(iOne, iTwo)))
                    continue;
                Dat.Set(i, j, CMeta::GetNaN());
            }
        }
    }

    if (sArgs.summary_flag) {
        double sum, sq_sum, cTotal, mean, variance, d;

        sum = 0.0;
        sq_sum = 0.0;
        cTotal = 0.0;

        for (i = 0; i < Dat.GetGenes(); ++i)
            for (j = (i + 1); j < Dat.GetGenes(); ++j)
                if (!CMeta::IsNaN(d = Dat.Get(i, j))) {
                    sum += d;
                    sq_sum += (d * d);
                    cTotal += 1;
                }

        mean = sum / cTotal;
        variance = ((sq_sum) - cTotal * (mean * mean)) / (cTotal - 1);
        cout << mean << endl;
        cout << sqrt(variance) << endl;
        return 0;
    }

    if (sArgs.lookups1_arg) {
        CGenes GenesLk1(Genome);
        vector <size_t> veciGenesOne;
        size_t iOne, iTwo;
        float d;

        ifsm.clear();
        ifsm.open(sArgs.lookups1_arg);
        if (!GenesLk1.Open(ifsm)) {
            cerr << "Could not open: " << sArgs.lookups1_arg << endl;
            return 1;
        }
        ifsm.close();
        veciGenesOne.resize(GenesLk1.GetGenes());
        for (i = 0; i < veciGenesOne.size(); ++i)
            veciGenesOne[i] = Dat.GetGene(GenesLk1.GetGene(i).GetName());
        if (sArgs.lookups2_arg) {
            CGenes GenesLk2(Genome);
            vector <size_t> veciGenesTwo;

            ifsm.clear();
            ifsm.open(sArgs.lookups2_arg);
            if (!GenesLk2.Open(ifsm)) {
                cerr << "Could not open: " << sArgs.lookups2_arg << endl;
                return 1;
            }
            ifsm.close();
            veciGenesTwo.resize(GenesLk2.GetGenes());
            for (i = 0; i < veciGenesTwo.size(); ++i)
                veciGenesTwo[i] = Dat.GetGene(GenesLk2.GetGene(i).GetName());
            for (i = 0; i < veciGenesOne.size(); ++i) {
                if ((iOne = veciGenesOne[i]) == -1)
                    continue;
                for (j = 0; j < veciGenesTwo.size(); ++j)
                    if (((iTwo = veciGenesTwo[j]) != -1) && !CMeta::IsNaN(d = Dat.Get(iOne, iTwo)))
                        cout << Dat.GetGene(iOne) << '\t' << Dat.GetGene(iTwo) << '\t' << d << endl;
            }
        } else if (sArgs.lookup1_arg) {
            if ((iOne = Dat.GetGene(sArgs.lookup1_arg)) != -1)
                for (i = 0; i < veciGenesOne.size(); ++i)
                    if (((iTwo = veciGenesOne[i]) != -1) &&
                        !CMeta::IsNaN(d = Dat.Get(iOne, iTwo)))
                        cout << Dat.GetGene(iOne) << '\t' << Dat.GetGene(iTwo) << '\t' << d << endl;
        } else
            for (i = 0; i < veciGenesOne.size(); ++i)
                if ((iOne = veciGenesOne[i]) != -1)
                    for (j = (i + 1); j < veciGenesOne.size(); ++j)
                        if (((iTwo = veciGenesOne[j]) != -1) &&
                            !CMeta::IsNaN(d = Dat.Get(iOne, iTwo)))
                            cout << Dat.GetGene(iOne) << '\t' << Dat.GetGene(iTwo) << '\t' << d << endl;
        return 0;
    } else if (sArgs.lookup1_arg) {
        if ((i = Dat.GetGene(sArgs.lookup1_arg)) == -1) {
            cerr << "Unknown gene: " << sArgs.lookup1_arg << endl;
            return 1;
        }
        if (sArgs.lookup2_arg) {
            if ((j = Dat.GetGene(sArgs.lookup2_arg)) == -1) {
                cerr << "Unknown gene: " << sArgs.lookup2_arg << endl;
                return 1;
            }
            cout << Dat.GetGene(i) << '\t' << Dat.GetGene(j) << '\t' << Dat.Get(i, j) << endl;
        } else
            for (j = 0; j < Dat.GetGenes(); ++j) {
                float d;

                if ((i != j) && !CMeta::IsNaN(d = Dat.Get(i, j)))
                    cout << Dat.GetGene(i) << '\t' << Dat.GetGene(j) << '\t' << d << endl;
            }
    } else if (sArgs.table_flag) {
        for (i = 1; i < Dat.GetGenes(); ++i)
            cout << '\t' << Dat.GetGene(i);
        cout << endl;
        for (i = 0; (i + 1) < Dat.GetGenes(); ++i) {
            cout << Dat.GetGene(i);
            for (j = 0; j < i; ++j)
                cout << '\t';
            for (++j; j < Dat.GetGenes(); ++j) {
                float d;

                cout << '\t';
                if (!CMeta::IsNaN(d = Dat.Get(i, j)))
                    cout << d;
            }
            cout << endl;
        }
    } else if (sArgs.output_arg) {
        CDataPair datOut;
        datOut.Open(Dat);
        if (sArgs.quant_arg) {
            datOut.OpenQuants(sArgs.quant_arg);
            datOut.Quantize();
        }
        datOut.Save(sArgs.output_arg);
    } else {
        Dat.Save(cout, CDat::EFormatText);
        cout.flush();
    }

    return 0;
}

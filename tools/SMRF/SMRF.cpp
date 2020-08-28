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
#include <random>

static inline float enrange(float d, double dEpsilon) {
    float dTmp;

    if (d < dEpsilon)
        return (float) dEpsilon;
    if (d > (dTmp = (1 - (float) dEpsilon)))
        return dTmp;

    return d;
}

int main(int iArgs, char **aszArgs) {
    gengetopt_args_info sArgs;
    CPCL PCL;
    CDat Dat;
    size_t iIter, iCur, iPCL, iDat, i;
    float d, dCur, dLogIn, dLogOut;
    vector <size_t> veciPCL2Dat, veciShuffle;

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }
    CMeta Meta(sArgs.verbosity_arg, sArgs.random_arg);

    if (!Dat.Open(sArgs.graph_arg, !!sArgs.memmap_flag)) {
        cerr << "Could not open: " << sArgs.graph_arg << endl;
        return 1;
    }
    if (!PCL.Open(sArgs.pvalues_arg, sArgs.skip_arg)) {
        cerr << "Could not open: " << sArgs.pvalues_arg << endl;
        return 1;
    }
    veciPCL2Dat.resize(PCL.GetGenes());
    for (i = 0; i < veciPCL2Dat.size(); ++i)
        veciPCL2Dat[i] = Dat.GetGene(PCL.GetGene(i));

    veciShuffle.resize(PCL.GetGenes());
    for (i = 0; i < veciShuffle.size(); ++i)
        veciShuffle[i] = i;
    for (iIter = 0; iIter < (size_t) sArgs.iterations_arg; ++iIter) {
        cerr << "Iteration: " << iIter << '/' << sArgs.iterations_arg << endl;
        random_device rng;
        mt19937 urng(rng());
        shuffle(veciShuffle.begin(), veciShuffle.end(), urng);
        for (iCur = 0; iCur < veciShuffle.size(); ++iCur) {
            if (!(iCur % 1000))
                cerr << "Gene: " << iCur << '/' << veciShuffle.size() << endl;
            iPCL = veciShuffle[iCur];
            if ((iDat = veciPCL2Dat[iPCL]) == -1)
                continue;
            dLogIn = dLogOut = 0;
            for (i = 0; i < PCL.GetGenes(); ++i) {
                if (((sArgs.neighbors_arg < 1) && (((float) rand() / RAND_MAX) > sArgs.neighbors_arg)) ||
                    (i == iPCL) || (veciPCL2Dat[i] == -1) ||
                    CMeta::IsNaN(d = Dat.Get(iDat, veciPCL2Dat[i])))
                    continue;
                d = enrange(d, sArgs.epsilon_arg);
                dCur = enrange(PCL.Get(i, 0), sArgs.epsilon_arg);
                dLogIn += log((d * (1 - dCur)) + ((1 - d) * dCur));
                dLogOut += log((d * dCur) + ((1 - d) * (1 - dCur)));
            }
            d = enrange(PCL.Get(iPCL, 0), sArgs.epsilon_arg);
// P(gout|data) = P(data|gout)P(gout)/P(data) = P(n1|gout)P(n2|gout)...P(gout)/(P(data|gout)P(gout)+P(data|~gout)P(~gout))
// P(n1|gout) = P(n1)P(rel) + P(~n1)P(~rel)
            PCL.Set(iPCL, 0, 1 / (1 + ((1 - d) * exp((dLogIn - dLogOut) / PCL.GetGenes()) / d)));
        }
    }

    PCL.Save(cout);

    return 0;
}

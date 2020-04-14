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

static const char c_acXDSL[] = ".xdsl";
static const char c_acDSL[] = ".dsl";

int main(int iArgs, char **aszArgs) {
    gengetopt_args_info sArgs;
    CBayesNetMinimal BNDefault;
    vector <CBayesNetMinimal> vecBNs;
    ifstream ifsm;
    uint32_t iSize;
    size_t i;
    CBayesNetSmile BNSmile;
    CPCL PCLDatasets(false);
    string strDir, strFile;
    vector <string> vecstrGenes;

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }
    CMeta Meta(sArgs.verbosity_arg);

    if (!PCLDatasets.Open(sArgs.datasets_arg, 1)) {
        cerr << "Could not open: " << (sArgs.datasets_arg ? sArgs.datasets_arg : "stdin") << endl;
        return 1;
    }
    vecstrGenes.resize(PCLDatasets.GetGenes());
    for (i = 0; i < vecstrGenes.size(); ++i)
        vecstrGenes[i] = PCLDatasets.GetFeature(i, 1);

    ifsm.open(sArgs.input_arg, ios_base::binary);
    if (!BNDefault.Open(ifsm)) {
        cerr << "Could not read: " << sArgs.input_arg << endl;
        return 1;
    }
    ifsm.read((char *) &iSize, sizeof(iSize));
    vecBNs.resize(iSize);
    for (i = 0; i < vecBNs.size(); ++i)
        if (!vecBNs[i].Open(ifsm)) {
            cerr << "Could not read: " << sArgs.input_arg << " (" << i << ")" << endl;
            return 1;
        }
    ifsm.close();

    strDir = (string) sArgs.output_arg + '/';
    if (!BNSmile.Open(BNDefault, vecstrGenes))
        return 1;
    strFile = strDir + BNDefault.GetID() + (sArgs.xdsl_flag ? c_acXDSL : c_acDSL);
    cerr << "Saving: " << strFile << endl;
    BNSmile.Save(strFile.c_str());
    for (i = 0; i < vecBNs.size(); ++i) {
        if (!BNSmile.Open(vecBNs[i], vecstrGenes))
            return 1;
        strFile = strDir + vecBNs[i].GetID() + (sArgs.xdsl_flag ? c_acXDSL : c_acDSL);
        cerr << "Saving: " << strFile << endl;
        BNSmile.Save(strFile.c_str());
    }

    return 0;
}

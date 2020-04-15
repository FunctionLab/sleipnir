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

int main(int iArgs, char **aszArgs) {
    gengetopt_args_info sArgs;
    CPCL PCL;
    CDat Dat;
    int iRet;
    ofstream ofsm;

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }
    CMeta Meta(sArgs.verbosity_arg);

    if (iRet = CPCL::Distance(sArgs.input_arg, sArgs.skip_arg, sArgs.neighbors_arg ? sArgs.distance_arg : NULL,
                              false, false, !!sArgs.autocorrelate_flag, sArgs.genes_arg, CMeta::GetNaN(),
                              sArgs.limit_arg, PCL,
                              Dat)) {
        cmdline_parser_print_help();
        return iRet;
    }

    PCL.Impute(sArgs.neighbors_arg, (float) sArgs.missing_arg, Dat);

    if (sArgs.output_arg) {
        ofsm.open(sArgs.output_arg);
        PCL.Save(ofsm);
        ofsm.close();
    } else {
        PCL.Save(cout);
        cout.flush();
    }

    return 0;
}

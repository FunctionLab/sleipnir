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

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }
    CMeta Meta(sArgs.verbosity_arg);

    IMeasure::EMap eM = IMeasure::EMapCenter;

    if (sArgs.centering_flag == 0) {
        eM = IMeasure::EMapNone;
        fprintf(stderr, "INFO: centering is off\n");
    }


    if ((iRet = CPCL::Distance(sArgs.input_arg, sArgs.skip_arg, sArgs.weights_arg, sArgs.distance_arg,
                               sArgs.normalize_flag != 0,
                               sArgs.zscore_flag != 0,
                               sArgs.autocorrelate_flag != 0,
                               sArgs.genes_arg,
                               sArgs.cutoff_given ? (float) sArgs.cutoff_arg : CMeta::GetNaN(),
                               sArgs.limit_arg, PCL, Dat,
                               eM,
                               sArgs.freqweight_flag != 0, sArgs.alpha_arg, sArgs.threads_arg))) {
        cmdline_parser_print_help();
        return iRet;
    }
    if (sArgs.flip_flag)
        Dat.Invert();

    if (sArgs.quant_arg) {
        CDataPair datOut;
        datOut.Open(Dat);
        datOut.OpenQuants(sArgs.quant_arg);
        datOut.Quantize();
        datOut.Save(sArgs.output_arg);
    } else {
        Dat.Save(sArgs.output_arg);
    }

    return 0;
}

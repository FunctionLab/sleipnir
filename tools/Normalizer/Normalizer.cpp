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

const struct {
    const char *m_szName;
    CPCL::ENormalize m_eType;
} c_asTypesPCL[] = {
        {"columnz",   CPCL::ENormalizeColumn},
        {"rowz",      CPCL::ENormalizeRow},
        {"globalz",   CPCL::ENormalizeZScore},
        {"0to1",      CPCL::ENormalizeMinMax},
        {"colcenter", CPCL::ENormalizeColumnCenter},
        {"colfrac",   CPCL::ENormalizeColumnFraction},
        {NULL,        CPCL::ENormalizeNone}
};

const struct {
    const char *m_szName;
    CDat::ENormalize m_eType;
} c_asTypesDAT[] = {
        {"globalz", CDat::ENormalizeZScore},
        {"0to1",    CDat::ENormalizeMinMax},
        {"sigmoid", CDat::ENormalizeSigmoid},
        {"normcdf", CDat::ENormalizeNormCDF},
        {"pcc",     CDat::ENormalizePCC},
        {NULL,      CDat::ENormalizeNone}
};

int main(int iArgs, char **aszArgs) {
    CDat Dat;
    CPCL PCL;
    size_t i, j;
    ofstream ofsm;
    gengetopt_args_info sArgs;
    float d;
    istream *pistm;
    ifstream ifsm;
    ostream *postm;
    CPCL::ENormalize eTypePCL;
    CDat::ENormalize eTypeDAT;

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }
    CMeta Meta(sArgs.verbosity_arg);

    if (!strcmp("pcl", sArgs.itype_arg)) {
        if (sArgs.input_arg) {
            ifsm.open(sArgs.input_arg);
            pistm = &ifsm;
        } else
            pistm = &cin;
        if (!PCL.Open(*pistm, sArgs.skip_arg)) {
            cerr << "Could not open input: " << (sArgs.input_arg ? sArgs.input_arg : "standard input") << endl;
            return 1;
        }
        if (sArgs.input_arg)
            ifsm.close();

        if (!strcmp("medmult", sArgs.otype_arg))
            PCL.MedianMultiples();
        else {
            eTypePCL = CPCL::ENormalizeNone;
            for (i = 0; c_asTypesPCL[i].m_szName; ++i)
                if (!strcmp(c_asTypesPCL[i].m_szName, sArgs.otype_arg)) {
                    eTypePCL = c_asTypesPCL[i].m_eType;
                    break;
                }
            PCL.Normalize(eTypePCL);
        }

        if (sArgs.output_arg) {
            ofsm.open(sArgs.output_arg);
            postm = &ofsm;
        } else
            postm = &cout;
        PCL.Save(*postm);
        if (sArgs.output_arg)
            ofsm.close();
    } else if (!sArgs.input_arg) {
        cmdline_parser_print_help();
        return 1;
    } else {
        if (!Dat.Open(sArgs.input_arg)) {
            cerr << "Could not open input: " << sArgs.input_arg << endl;
            return 1;
        }

        eTypeDAT = CDat::ENormalizeNone;
        for (i = 0; c_asTypesDAT[i].m_szName; ++i)
            if (!strcmp(c_asTypesDAT[i].m_szName, sArgs.otype_arg)) {
                eTypeDAT = c_asTypesDAT[i].m_eType;
                break;
            }

        Dat.Normalize(eTypeDAT);
        if (sArgs.flip_flag)
            for (i = 0; i < Dat.GetGenes(); ++i)
                for (j = (i + 1); j < Dat.GetGenes(); ++j)
                    if (!CMeta::IsNaN(d = Dat.Get(i, j)))
                        Dat.Set(i, j, 1 - d);

        Dat.Save(sArgs.output_arg ? sArgs.output_arg : sArgs.input_arg);
    }

    return 0;
}

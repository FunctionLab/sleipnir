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
#include "statistics.h"
static const char c_szRBF[] = "rbf";
static const char c_szPolynomial[] = "poly";
static const char c_acDab[] = ".dab";
static const char c_acPcl[] = ".pcl";

int main(int iArgs, char** aszArgs) {
	gengetopt_args_info sArgs;
	CGenome Genome;
	CGenes Genes(Genome);
	ifstream ifsm;
	CPCL PCL;
	size_t i, j;
	bool fModified;
	char* file_ext = NULL;

	if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
		cmdline_parser_print_help();
		return 1;
	}
	CMeta Meta(sArgs.verbosity_arg);

	if (sArgs.input_arg) {
	  if (!PCL.Open(sArgs.input_arg, sArgs.skip_arg, !!sArgs.mmap_flag, sArgs.rPCL_flag)) {
			cerr << "Could not open: " << sArgs.input_arg << endl;
			return 1;
		}
	}
	// if coming from stdin, always assume it's non binary form

	else if (!PCL.Open(cin, sArgs.skip_arg)) {
		cerr << "Could not open input" << endl;
		return 1;
	}
	if (sArgs.Genes_flag) {
		for (size_t i = 0; i < PCL.GetGenes(); i++) {
			cout << i << '\t' << PCL.GetGene(i) << endl;
		}
	}
	else if (sArgs.transpose_flag){
	  for ( j = 0; j < PCL.GetGenes(); j++)
	    cout<<'\t'<<PCL.GetGene(j);
	  cout<<endl;
	  for ( i = 0; i < PCL.GetExperiments(); i++) {
	    cout<<PCL.GetExperiment(i);
	    for ( j = 0; j < PCL.GetGenes(); j++) {
	      cout << '\t' << PCL.Get(j,i);
	    }
	    cout<<endl;
	  }
	}
 else {
		//vector<size_t> vec_iGenes;
		ifstream ifsm;

		if (sArgs.genex_given) {
			ifsm.open(sArgs.genex_arg);
			if (!Genes.Open(ifsm)) {
				cerr << "Could not open: " << sArgs.genex_arg << endl;
				return 1;
			}

			for (i = 0; i < PCL.GetGenes(); i++)
				if (Genes.GetGene(PCL.GetGene(i)) != -1)
					PCL.MaskGene(i, true);
		}

		else if (sArgs.genes_given) {
			ifsm.open(sArgs.genes_arg);
			if (!Genes.Open(ifsm)) {
				cerr << "Could not open: " << sArgs.genes_arg << endl;
				return 1;
			}
			for (i = 0; i < PCL.GetGenes(); i++)
				if (Genes.GetGene(PCL.GetGene(i)) == -1)
					PCL.MaskGene(i, true);

		}
		ifsm.close();

		/*	if( Genes.GetGenes( ) )
		 PCL.FilterGenes( Genes, CPCL::EFilterInclude );
		 if( sArgs.genex_arg )
		 PCL.FilterGenes( sArgs.genex_arg, CPCL::EFilterExclude );
		 */
		//Normalize Zscore
		if (sArgs.normalize_flag)
			PCL.Normalize(CPCL::ENormalizeMinMax);
		if (sArgs.zrow_flag)
					PCL.Normalize(CPCL::ENormalizeRow);
		if (sArgs.zcol_flag)
					PCL.Normalize(CPCL::ENormalizeColumn);
		if (sArgs.scol_flag)
					PCL.Normalize(CPCL::EMeanSubtractColumn);
		else if (sArgs.normalize_flag)
			PCL.Normalize(CPCL::ENormalizeMinMax);
		if (sArgs.output_arg) {

				PCL.Save(sArgs.output_arg);
		} else {

				PCL.Save(cout, NULL);
			cout.flush();
		}

		return 0;

	}
}

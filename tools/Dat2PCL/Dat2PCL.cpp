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

int main(int iArgs, char** aszArgs) {
	gengetopt_args_info sArgs;
	CGenome Genome;
	CGenes Genes(Genome);       
	CPCL PCL;
	CDat Dat;
	size_t i, j;
	bool fModified;
	vector<string> features;
	
	if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
		cmdline_parser_print_help();
		return 1;
	}
	CMeta Meta(sArgs.verbosity_arg);
	
	if( !Dat.Open( sArgs.input_arg, false, 0, false, true)){
	  cerr << "Could not open: " << sArgs.input_arg << endl;
	  return 1; }
	
	PCL.Open( Dat.GetGeneNames( ), Dat.GetGeneNames( ), features);
	
	cerr << "Gene count: " << Dat.GetGenes() << endl;
	
	PCL.populate( sArgs.input_arg );
		
	if (sArgs.output_arg) {	  
	  PCL.Save(sArgs.output_arg);
	} else {
	  PCL.Save(cout, NULL);
	  cout.flush();
	}	
	return 0;
}

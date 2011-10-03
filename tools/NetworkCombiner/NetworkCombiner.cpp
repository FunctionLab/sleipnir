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

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	int					iRet;
	size_t				i, j, k, l;
	float d;
	DIR* dp;
	struct dirent* ep;
	CDat					DatOut, DatCur;
	vector<size_t>				veciGenesCur;	
	// store all input file names
	vector<string> input_files;
	
	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg );
	
	// now collect the data files from directory if given
	if(sArgs.directory_arg){
	  dp = opendir (sArgs.directory_arg);
	  if (dp != NULL){
	    while (ep = readdir (dp)){
	      // skip . .. files and temp files with ~
	      if (ep->d_name[0] == '.' || ep->d_name[strlen(ep->d_name)-1] == '~') 
		continue;
	      
	      // currently opens all files. Add filter here if want pick file extensions
	      input_files.push_back((string)sArgs.directory_arg + "/" + ep->d_name);	      
	    }
	    (void) closedir (dp);	    
	  }
	  else{
	    cerr << "Couldn't open the directory: " << sArgs.directory_arg << '\n';
	    return 1;
	  }
	}
	else{
	  input_files.resize( sArgs.inputs_num );
	  copy( sArgs.inputs, sArgs.inputs + sArgs.inputs_num, input_files.begin( ) );
	}
	
	// now iterate dat/dab networks
	for( i = 0; i < input_files.size( ); ++i ) {
	  // open first network, we will just add edge weights to this CDat
	  if( i == 0 ){
	    if( !DatCur.Open( input_files[ i ].c_str() ) ) {
	      cerr << "Couldn't open: " << input_files[ i ] << endl;
	      return 1; }
	    
	    DatOut.Open( DatCur );	    
	    cerr << "opened: " << input_files[ i ] << endl;
	    continue;	 
	  }
	  
	  if( !DatCur.Open( input_files[ i ].c_str() ) ) {
	    cerr << "Couldn't open: " << input_files[ i ] << endl;
	    return 1; }
	  cerr << "opened: " << input_files[ i ] << endl;
	  
	  if( sArgs.map_flag ){
	    // Get gene index match	  
	    veciGenesCur.clear();
	    veciGenesCur.resize(DatOut.GetGenes());
	    for( l = 0; l < DatOut.GetGenes(); l++){
	      veciGenesCur[ l ] = DatCur.GetGene( DatOut.GetGene(l) );
	      if( veciGenesCur[ l ] == -1 ){
		cerr << "ERROR: missing gene" << input_files[ l ] << endl;
		return 1;	      
	      }
	    }
	  }
	  
	  // now add edges to Dat
	  for( j = 0; j < DatOut.GetGenes( ); ++j )
	    for( k = ( j + 1 ); k < DatOut.GetGenes( ); ++k ) {
	      
	      if( sArgs.map_flag ){
		// we are assuming a fully connected network
		if( CMeta::IsNaN( d = DatCur.Get( veciGenesCur[ j ], veciGenesCur[ k ] ) ) ){
		  cerr << "ERROR: missing values" << input_files[ i ] << endl;
		  return 1;
		}
	      }
	      else{
		if( CMeta::IsNaN( d = DatCur.Get(  j, k ) ) ){
		  cerr << "ERROR: missing values" << input_files[ i ] << endl;
		  return 1;
		}
	      }
	      
	      DatOut.Set( j, k, DatOut.Get( j, k ) + d );
	    }
	}
	
	// now convert sum to mean
	for( j = 0; j < DatOut.GetGenes( ); ++j )
	  for( k = ( j + 1 ); k < DatOut.GetGenes( ); ++k )
	    DatOut.Set( j, k, DatOut.Get( j, k ) / input_files.size( ) );
	
	DatOut.Save( sArgs.output_arg );
	return iRet; 
}

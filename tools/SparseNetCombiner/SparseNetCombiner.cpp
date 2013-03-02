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

static const char   c_acDab[]   = ".dab";
static const char   c_acQDab[]   = ".qdab";
static const char   c_acDat[]   = ".dat";
static const char   c_acSVM[]   = ".svm";
static const char   c_acOUT[]   = ".out";

enum EMethod {
	EMethodBegin	= 0,
	EMethodMean		= EMethodBegin,
	EMethodMax		= EMethodMean + 1,
	EMethodEnd		= EMethodMax + 1
};

static const char*	c_aszMethods[]	= {
	"mean", "max", NULL
};

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	int					iRet;
	size_t				i, j, k, l;
	size_t    iRuns, inputNums;
	float d, dout;
	DIR* dp;
	struct dirent* ep;
	CDat					DatOut, DatTrack, DatCur;
	vector<size_t>				veciGenesCur;	
	// store all input file names
	vector<string> input_files;
	vector<string>		vecstrInputs, vecstrDatasets; 
	EMethod		eMethod;
	ifstream	ifsm;
	vector<float>	vecWeights;
	map<string, size_t> mapstriDatasets;
	
	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg );
	
	for( eMethod = EMethodBegin; eMethod < EMethodEnd; eMethod = (EMethod)( eMethod + 1 ) )
	  if( !strcmp( c_aszMethods[eMethod], sArgs.method_arg ) ){
	    cerr << "combine method: " << c_aszMethods[eMethod] << endl;
	    break;
	  }
	
	// now collect the data files from directory if given
	if(sArgs.directory_given){
	  dp = opendir (sArgs.directory_arg);
	  i = 0;
	  if (dp != NULL){
	    while (ep = readdir (dp)){
	      if(  strstr(ep->d_name, c_acDab) != NULL ||
		   strstr(ep->d_name, c_acQDab) != NULL ||
		   strstr(ep->d_name, c_acDat) != NULL
		   ){
		if(strstr(ep->d_name, c_acSVM) != NULL)
		  continue;
		if(strstr(ep->d_name, c_acOUT) != NULL)
		  continue;
		vecstrInputs.push_back((string)sArgs.directory_arg + "/" + ep->d_name);		
		vecstrDatasets.push_back(vecstrInputs[i]);		
		mapstriDatasets[ ep->d_name ] = i;
		++i;
	      }
	    }
	    (void) closedir (dp);
	    
	    inputNums = vecstrDatasets.size();
	    
	    // sort by ASCI 
	    //std::sort( vecstrInputs.begin(), vecstrInputs.end() );
	    
	    cerr << "Number of datasets: " << inputNums << '\n';
	  }
	  else{
	    cerr << "Couldn't open the directory: " << sArgs.directory_arg << '\n';
	    return 1;
	  }
	}
	
	// read dataset weight file
	if( sArgs.weight_given ) {
	  static const size_t	c_iBuffer	= 1024;
	  char					szBuffer[ c_iBuffer ];
	  string				strFile;
	  vector<string>			vecstrLine;
	  map<string, size_t>::const_iterator	iterDataset;
	  size_t  numDataset;
	  vecWeights.resize(vecstrInputs.size( ) );
	  
	  // set default weight
	  for(i=0; i < vecWeights.size(); i++){
	    vecWeights[i] = 0;
	  }
	  
	  ifsm.clear( );
	  ifsm.open( sArgs.weight_arg );
	  if( !ifsm.is_open( ) ) {
            cerr << "Could not open: " << sArgs.weight_arg << endl;
            return 1;
	  }
	  while( !ifsm.eof( ) ) {
            ifsm.getline( szBuffer, c_iBuffer - 1 );
            szBuffer[ c_iBuffer - 1 ] = 0;
            if( !szBuffer[ 0 ] )
	      continue;
            vecstrLine.clear( );
            CMeta::Tokenize( szBuffer, vecstrLine );
            if( vecstrLine.size( ) != 2 ) {
	      cerr << "Illegal weight line: " << szBuffer << endl;
	      return 1;
            }
	    
	    if( ( iterDataset = mapstriDatasets.find( vecstrLine[ 0 ] ) ) == mapstriDatasets.end( ) )
	      cerr << "Dataset in weights but not database: " << vecstrLine[ 0 ] << endl;
            else	      
	      vecWeights[ iterDataset->second ] = (float)atof( vecstrLine[ 1 ].c_str( ) );
	  }
	  ifsm.close( );
	  
	  // now only keep datasets with non-zero weighting
	  vecstrDatasets.clear();
	  numDataset = 0;
	  for(i = 0; i < vecstrInputs.size(); ++i){
	    if(vecWeights[i] == 0)
	      continue;
	    vecstrDatasets.push_back(vecstrInputs[i]);
	    vecWeights[numDataset] = vecWeights[i];
	    numDataset++;
	  }
	  vecstrDatasets.resize(numDataset);
	  vecWeights.resize(numDataset);
	  
	  cerr << "Total number of datasets combining with non-zero weights: " << numDataset << endl;
	}
	
	// now iterate dat/dab networks
	for( i = 0; i < vecstrDatasets.size( ); ++i ) {
	  
	  // open first network, we will just add edge weights to this CDat
	  if( i == 0 ){
	    if( !DatCur.Open( vecstrDatasets[ i ].c_str() ) ) {
	      cerr << "Couldn't open: " << vecstrDatasets[ i ] << endl;
	      return 1; }
	    
	    DatOut.Open( DatCur );	    	    
	    DatTrack.Open( DatCur );
	    
	    // this Dat is used to track various values (count, max)
	    for( j = 0; j < DatTrack.GetGenes( ); ++j )
	      for( k = ( j + 1 ); k < DatTrack.GetGenes( ); ++k ){
		if( CMeta::IsNaN( d = DatCur.Get( j, k)))
		  continue;
		
		if( sArgs.weight_given ){
		  d *= vecWeights[i];
		  // store weighted value (numerator)
		  DatOut.Set( j, k, d);
		}
		
		switch( eMethod ) {
		case EMethodMax:
		  DatOut.Set( j, k, d);
		  break;
		case EMethodMean:
		  // keep track of denominator
		  DatTrack.Set( j, k, 1.0);
		}
	      }
	    
	    cerr << "opened: " << vecstrDatasets[ i ] << endl;
	    continue;	 
	  }
	  
	  if( !DatCur.Open( vecstrDatasets[ i ].c_str() ) ) {
	    cerr << "Couldn't open: " << vecstrDatasets[ i ] << endl;
	    return 1; }
	  cerr << "opened: " << vecstrDatasets[ i ] << endl;
	  
	  if( sArgs.map_flag ){
	    // Get gene index match	  
	    cerr << "inside map flag" << endl;
	    veciGenesCur.clear();
	    veciGenesCur.resize(DatOut.GetGenes());
	    for( l = 0; l < DatOut.GetGenes(); l++){
	      veciGenesCur[ l ] = DatCur.GetGene( DatOut.GetGene(l) );
	      if( veciGenesCur[ l ] == -1 ){
		cerr << "ERROR: missing gene: " << DatOut.GetGene(l) << ", " << vecstrDatasets[ i ] << endl;
		return 1;
	      }
	    }
	  }
	  
	  // now add edges to Dat
	  for( j = 0; j < DatOut.GetGenes( ); ++j )
	    for( k = ( j + 1 ); k < DatOut.GetGenes( ); ++k ) {
	      
	      if( sArgs.map_flag ){
		if( CMeta::IsNaN( d = DatCur.Get( veciGenesCur[ j ], veciGenesCur[ k ] ) ) ){
		  continue;
		}
	      }
	      else{
		if( CMeta::IsNaN( d = DatCur.Get(  j, k ) ) ){
		  continue;
		}
	      }
	      
	      if( sArgs.weight_given )
		  d *= vecWeights[i];
	      
	      switch( eMethod ) {
	      case EMethodMax:
		if( CMeta::IsNaN( (dout = DatOut.Get( j, k ))) || d > dout )
		  DatOut.Set( j, k, d);
		break;
		
	      case EMethodMean:
		DatTrack.Set( j, k, DatTrack.Get( j, k ) + 1 );
		
		if( CMeta::IsNaN( (dout = DatOut.Get( j, k ))) )
		  DatOut.Set( j, k, d );
		else
		  DatOut.Set( j, k, DatOut.Get( j, k ) + d ); 	       
		break;
	      }
	    }
	}
	
	switch( eMethod ) {
	case EMethodMean:
	  // now convert sum to mean
	  for( j = 0; j < DatOut.GetGenes( ); ++j )
	    for( k = ( j + 1 ); k < DatOut.GetGenes( ); ++k )
	      DatOut.Set( j, k, DatOut.Get( j, k ) / DatTrack.Get( j, k ) );
	}
	
	DatOut.Save( sArgs.output_arg );
	return iRet; 
}
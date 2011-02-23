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


static const char   c_acDab[]   = ".dab";
static const char   c_acQDab[]   = ".qdab";

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info			sArgs;
	CDataPair					DatOne, DatTwo;
	size_t						iDatOne, iDatTwo;	
	size_t						i, j, iGeneOne, iGeneTwo, iCountOne, iCountTwo;
	size_t						iCountJoint, iJoint, iValueOne, iValueTwo;
	size_t    iRuns, inputNums;
	vector<size_t>				veciGenesOne, veciGenesTwo, veciOne, veciTwo, veciDefaults, veciSizes;
	vector<vector<size_t> >		vecveciJoint;
	vector<string>		vecstrInputs, vecstrDatasets; 
	map<string, size_t>			mapZeros;
	float						dOne, dJoint, dMI, dSubsample, dValueOne, dValueTwo;
	size_t  non_zero_bins;
	DIR* dp;
	struct dirent* ep;
	const char* fileOne;
	const char* fileTwo;

	set<string>               setstrGenes;
	set<string>::const_iterator iterGene;
	
	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg, sArgs.random_arg );
				
	ostream* posm;
	ofstream ofsm;
	ifstream ifsm;
	bool opened_DatOne;

	// open output file
	if (sArgs.output_given){
	  ofsm.open(sArgs.output_arg);
	  posm=&ofsm;
	}
	else{
	  posm=&cout;
	}
		
	if (sArgs.datasets_arg) {
		char acLine[ 1024 ];
		vector<string> vecstrTok;
		string datstr;

		ifsm.open(sArgs.datasets_arg);
		if ( !ifsm.is_open() ) {
			cerr << "Could not open: " << sArgs.datasets_arg << endl;
			return 1;
		}		
		while ( !ifsm.eof() ) {
			ifsm.getline( acLine, ARRAYSIZE(acLine) - 1 );
			acLine[ ARRAYSIZE(acLine) -1 ] = 0;
			if ( !acLine[0] )
				continue;
			vecstrTok.clear();
			//CMeta::Tokenize( acLine, vecstrTok );
			//if ( vecstrTok.size() < 2 ) {
			//	cerr << "Illegal line: " << acLine << endl;
			//	return 1;
			//}
			
			datstr = string(acLine);
			vecstrDatasets.push_back(datstr);		 
		}
	}


	// now collect the data files from directory if given
	if(sArgs.directory_arg){
	  dp = opendir (sArgs.directory_arg);
	  if (dp != NULL){
	    while (ep = readdir (dp)){
	      if(  strstr(ep->d_name, c_acDab) != NULL ||
		   strstr(ep->d_name, c_acQDab) != NULL
		   ){
		vecstrInputs.push_back((string)sArgs.directory_arg + "/" + ep->d_name);
	      }
	    }
	    (void) closedir (dp);
	    
	    inputNums = vecstrInputs.size();
	    cerr << "Number of datasets: " << inputNums << '\n';
	  }
	  else{
	    cerr << "Couldn't open the directory: " << sArgs.directory_arg << '\n';
	    return 1;
	  }
	}
	else{
	  inputNums = sArgs.inputs_num;
	}
	
	
	// read in zeros file if given
	if( sArgs.zeros_arg ) {
	  ifstream		ifsm;
	  vector<string>	vecstrZeros;
	  char			acLine[ 1024 ];
	  
	  ifsm.open( sArgs.zeros_arg );
	  if( !ifsm.is_open( ) ) {
	    cerr << "Could not open: " << sArgs.zeros_arg << endl;
	    return 1; }
	  while( !ifsm.eof( ) ) {
	    ifsm.getline( acLine, ARRAYSIZE(acLine) - 1 );
	    acLine[ ARRAYSIZE(acLine) - 1 ] = 0;
	    vecstrZeros.clear( );
	    CMeta::Tokenize( acLine, vecstrZeros );
	    if( vecstrZeros.empty( ) )
	      continue;
	    mapZeros[ vecstrZeros[ 0 ] ] = atoi( vecstrZeros[ 1 ].c_str( ) ); 
	  } 
	}
	
	// now set default values if given in zeros file use the given default value
	// if not, use random value. A value of -1 for a given dataset index in  veciDefaults 
	// represents to use random values
	veciDefaults.resize( inputNums );
	fill( veciDefaults.begin(), veciDefaults.end(), -1);
	if( sArgs.zeros_arg )
	  for( i = 0; i < veciDefaults.size( ); ++i ) {
	    map<string, size_t>::const_iterator	iterZero;
	    if( sArgs.inputs_num == 0 &&
		( iterZero = mapZeros.find( CMeta::Deextension( CMeta::Basename(vecstrInputs[ i ].c_str( ) ) ) ) ) != mapZeros.end( ) 
		){
	      veciDefaults[ i ] = iterZero->second;
	    }
	    else if( sArgs.inputs_num > 0 &&
		     ( iterZero = mapZeros.find( CMeta::Deextension( CMeta::Basename(sArgs.inputs[ i ]) ) ) ) != mapZeros.end( ) 
		     ){
	      veciDefaults[ i ] = iterZero->second;
	    }
	  }
	
	
	
	
	// now iterate through experiments to calculate MI
	for( iDatOne = iRuns = 0; iDatOne < ( sArgs.datasets_arg ? vecstrDatasets.size() : inputNums ); ++iDatOne){	  	  
	  opened_DatOne = false;
	  for( iDatTwo = ( sArgs.datasets_arg ? 0 : iDatOne );  iDatTwo < inputNums; ++iDatTwo, ++iRuns){

	    // ok skip runs not in range
	    if(sArgs.start_arg != -1 && iRuns < sArgs.start_arg ){
	      continue;
	    }
	    if(sArgs.end_arg != -1 && iRuns >= sArgs.end_arg){
	      continue;
	    }
	    
	    if (iRuns % 50 == 0)
	      cerr << "Processing " << iRuns << " among # of datasets: " << inputNums << '\n';
	    
	    // have I opned the first Data file?
	    if(!opened_DatOne){	      
	      // now open first Dat
	      if(sArgs.datasets_arg) {
		fileOne = vecstrDatasets[iDatOne].c_str();
	      }
	      else if(sArgs.directory_arg){
		fileOne = vecstrInputs[iDatOne].c_str();		
	      }
	      else{
		fileOne = sArgs.inputs[ iDatOne ];		
	      }
	      
	      if( !DatOne.Open(fileOne, false) ) {
		cerr << "Could not open: " << fileOne << '\n';
		return 1; }
	      
	      DatOne.Quantize();
	      opened_DatOne = true;
	      
	      veciOne.resize( DatOne.GetValues() );
	      
	      // laplace smoothingand also prevent devision by zeros
	      fill( veciOne.begin( ), veciOne.end( ), 1 );
	      iCountOne = DatOne.GetValues();
	      
	      // for entropy calculation
	      if(iDatTwo == iDatOne || !sArgs.union_flag){		
		// now get the count distribution of experiment i
		for(i = 0; i < DatOne.GetGenes(); ++i){
		  for( j = ( i + 1 ); j < DatOne.GetGenes(); ++j ) {
		    if( !CMeta::IsNaN( iValueOne = (size_t)DatOne.Get(i, j) ) ){
		      veciOne[ iValueOne ]++;
		      iCountOne++;
		    }
		  }
		}
		
	      }
	      
	      // When equal Only calculate self Entropy!
	      if( iDatTwo == iDatOne ){
		
		// compute the Entropy
		for( non_zero_bins = 0, dMI = 0,i = 0; i < veciOne.size( ); ++i )
		  if( dOne = (float)veciOne[ i ] / iCountOne ){
		    dMI += dOne * log( 1 / dOne );
		    ++non_zero_bins;
		  }
		
		//unbiased estimate correction
		dMI += ( non_zero_bins - 1 ) / ( 2.0f * ( iCountOne + iCountOne ) ); 
		dMI /= log( 2.0f );	    
		
		*posm << fileOne << '\t' << fileOne << '\t' << dMI << '\n';
		
		// ok now lets skip to next dataset do some MI calculation
		continue;
	      }
	    }
	    else if(sArgs.union_flag){
	      // reset since now you  have different number of genes	      
	      // laplace smoothingand also prevent devision by zeros
	      fill( veciOne.begin( ), veciOne.end( ), 1 );
	      
	      iCountOne = DatOne.GetValues();

	    }
	    
	    // now open Second Dat
	    if(sArgs.directory_arg){
	      fileTwo = vecstrInputs[iDatTwo].c_str( );
	    }
	    else{
	      fileTwo = sArgs.inputs[ iDatTwo ];
	    }
	    
	    if( !DatTwo.Open( fileTwo, false) ) {
	      cerr << "Could not open: " << fileTwo << '\n';
	      return 1; }
	    
	    
	    DatTwo.Quantize();	    
	    veciTwo.resize( DatTwo.GetValues() );
	    
	    // laplace smoothingand also prevent devision by zeros
	    fill( veciTwo.begin( ), veciTwo.end( ), 1 );
	    iCountTwo = DatTwo.GetValues();

	    // create & initialize joint bins
	    vecveciJoint.resize( veciOne.size( ) );
	    for( i = 0; i < vecveciJoint.size( ); ++i ) {
	      vecveciJoint[ i ].resize( veciTwo.size( ) );
	      fill( vecveciJoint[ i ].begin( ), vecveciJoint[ i ].end( ), 0 ); 
	    }
	    
	    // complie gene list from two dabs
	    setstrGenes.clear();
	    for(i = 0; i < DatOne.GetGenes(); i++)
	      setstrGenes.insert(DatOne.GetGene(i));	    
	    for(i = 0; i < DatTwo.GetGenes(); i++)
	      setstrGenes.insert(DatTwo.GetGene(i));	  
	    
	    veciGenesOne.clear();
	    veciGenesTwo.clear();
	    veciGenesOne.resize(setstrGenes.size());
	    veciGenesTwo.resize(setstrGenes.size());
	    
	    // store if a gene belongs to a dataset
	    for( iterGene = setstrGenes.begin( ),i = 0; iterGene != setstrGenes.end( ); ++iterGene,++i ){
	      veciGenesOne[ i ] = DatOne.GetGene(*iterGene);
	      veciGenesTwo[ i ] = DatTwo.GetGene(*iterGene);
	    }
	    
	    // when using only shared genes
	    if(!sArgs.union_flag){
	    
	      for( i = iCountJoint = 0; i < setstrGenes.size(); ++i ) {
		for( j = i + 1; j < setstrGenes.size(); ++j ) {
		  
		  // does this gene belong to DatTwo?
		  if(veciGenesTwo[ i ] != -1 && 
		     veciGenesTwo[ j ] != -1 &&
		     !CMeta::IsNaN( iValueTwo = (size_t)DatTwo.Get(veciGenesTwo[ i ], veciGenesTwo[ j ]) )){		  
		    
		    iCountTwo++;
		    veciTwo[ iValueTwo ]++;		  
		    
		    if( iValueTwo >= veciTwo.size( ) ){
		      cerr << iValueTwo << "\tj-1:" << j << endl;
		      return 1;
		    }
		  }
		  else{
		    // ok data set two doesn't have this gene pair
		    continue;
		  }
		  
		  // does this gene belong to DatOne?
		  if ( veciGenesOne[ i ] != -1 && 
		       veciGenesOne[ j ] != -1 && 
		       !CMeta::IsNaN( iValueOne = (size_t)DatOne.Get(veciGenesOne[ i ], veciGenesOne[ j ]))){
		    
		    iCountJoint++;		  
		    vecveciJoint[ iValueOne ][ iValueTwo ]++; 		  
		  }		
		}
	      }
	    }
	    else{ // when using union of genes of two datasets
	      
	      for( i =  iCountOne = iCountJoint = 0; i < setstrGenes.size(); ++i ) {
		for( j = i + 1; j < setstrGenes.size(); ++j ) {
		  
		  // does this gene belong to DatOne? If not, randomize it.
		  if ( veciGenesOne[ i ] == -1 || 
		       veciGenesOne[ j ] == -1 || 
		       CMeta::IsNaN( iValueOne = (size_t)DatOne.Get(veciGenesOne[ i ], veciGenesOne[ j ]))){
		    
		    if( veciDefaults[ iDatOne ] == -1 )
		      iValueOne = rand( ) % DatOne.GetValues( );
		    else
		      iValueOne = veciDefaults[ iDatOne ];
		  }
		  
		  // does this gene belong to DatTwo? If not, randomize it.
		  if(veciGenesTwo[ i ] == -1 || 
		     veciGenesTwo[ j ] == -1 ||
		     CMeta::IsNaN( iValueTwo = (size_t)DatTwo.Get(veciGenesTwo[ i ], veciGenesTwo[ j ]) )){		  
		    
		    if( veciDefaults[ iDatTwo ] == -1 )
		      iValueTwo = rand( ) % DatTwo.GetValues( );
		    else
		      iValueTwo = veciDefaults[ iDatTwo ];
		  }
		  
		  iCountJoint++;
		  iCountTwo++;
		  iCountOne++;
		  veciOne[ iValueOne ]++;
		  veciTwo[ iValueTwo ]++;		  
		  vecveciJoint[ iValueOne ][ iValueTwo ]++;		  
		}	      
	      }	      
	    }
	    
	    if(iCountJoint > 0){
	      // now calculate the MI between the two experiments
	      for( non_zero_bins = 0, dMI = 0,i = 0; i < veciOne.size( ); ++i ) {
		dOne = (float)veciOne[ i ] / iCountOne;
		for( j = 0; j < veciTwo.size( ); ++j )
		  if( iJoint = vecveciJoint[ i ][ j ] ) {
		    ++non_zero_bins;
		    dJoint = (float)iJoint / iCountJoint;
		    dMI += dJoint * log( dJoint * iCountTwo / dOne / veciTwo[ j ] );
		  } 
	      }
	      
	      //unbiased estimate correction
	      dMI += ( non_zero_bins - 1 ) / ( 2.0f * ( iCountOne + iCountTwo ) );
	      dMI = ( dMI < 0 ) ? 0 : ( dMI / log( 2.0f ) ); 
	    }
	    else{
	      dMI = 0;
	    }
	    	    
	    *posm << fileOne << '\t' << fileTwo << '\t' << dMI << '\n';	    	    
	  }  
	}
	
	return 0; 
}

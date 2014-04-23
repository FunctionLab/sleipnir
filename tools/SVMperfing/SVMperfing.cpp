#include <fstream>

#include <vector>
#include <queue>

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

using namespace SVMLight;

struct ParamStruct {
	vector<float> vecK, vecTradeoff;
	vector<size_t> vecLoss;
	vector<char*> vecNames;
};

ParamStruct ReadParamsFromFile(ifstream& ifsm, string outFile) {
	static const size_t c_iBuffer = 1024;
	char acBuffer[c_iBuffer];
	char* nameBuffer;
	vector<string> vecstrTokens;
	size_t extPlace;
	string Ext, FileName;
	if ((extPlace = outFile.find_first_of(".")) != string::npos) {
		FileName = outFile.substr(0, extPlace);
		Ext = outFile.substr(extPlace, outFile.size());
	} else {
		FileName = outFile;
		Ext = "";
	}
	ParamStruct PStruct;
	size_t index = 0;
	while (!ifsm.eof()) {
		ifsm.getline(acBuffer, c_iBuffer - 1);
		acBuffer[c_iBuffer - 1] = 0;
		vecstrTokens.clear();
		CMeta::Tokenize(acBuffer, vecstrTokens);
		if (vecstrTokens.empty())
			continue;
		if (vecstrTokens.size() != 3) {
			cerr << "Illegal params line (" << vecstrTokens.size() << "): "
					<< acBuffer << endl;
			continue;
		}
		if (acBuffer[0] == '#') {
			cerr << "skipping " << acBuffer << endl;
		} else {
			PStruct.vecLoss.push_back(atoi(vecstrTokens[0].c_str()));
			PStruct.vecTradeoff.push_back(atof(vecstrTokens[1].c_str()));
			PStruct.vecK.push_back(atof(vecstrTokens[2].c_str()));
			PStruct.vecNames.push_back(new char[c_iBuffer]);
			if (PStruct.vecLoss[index] == 4 || PStruct.vecLoss[index] == 5)
				sprintf(PStruct.vecNames[index], "%s_l%d_c%4.6f_k%4.3f%s",
						FileName.c_str(), PStruct.vecLoss[index],
						PStruct.vecTradeoff[index], PStruct.vecK[index],
						Ext.c_str());
			else
				sprintf(PStruct.vecNames[index], "%s_l%d_c%4.6f%s",
						FileName.c_str(), PStruct.vecLoss[index],
						PStruct.vecTradeoff[index], Ext.c_str());
			index++;
		}

	}
	return PStruct;
}

bool ReadModelFile(ifstream& ifsm, vector<float>& SModel) {
	static const size_t c_iBuffer = 1024;
	char acBuffer[c_iBuffer];
	char* nameBuffer;
	vector<string> vecstrTokens;
	size_t extPlace;
	string Ext, FileName;
	size_t index = 0;
	
	while (!ifsm.eof()) {
		ifsm.getline(acBuffer, c_iBuffer - 1);
		acBuffer[c_iBuffer - 1] = 0;
		vecstrTokens.clear();
		CMeta::Tokenize(acBuffer, vecstrTokens);
		if (vecstrTokens.empty())
			continue;
		if (vecstrTokens.size() > 1) {
			cerr << "Illegal model line (" << vecstrTokens.size() << "): "
					<< acBuffer << endl;
			continue;
		}
		if (acBuffer[0] == '#') {
			cerr << "skipping " << acBuffer << endl;
		} else {
		  SModel.push_back(atof(vecstrTokens[0].c_str()));
		}
		
		
	}
	ifsm.close();
	return true;
}

bool ReadGenesHoldoutFoldFile(ifstream& ifsm, map<string, size_t>& mapGene2Fold) {
	static const size_t c_iBuffer = 1024;
	char acBuffer[c_iBuffer];
	char* nameBuffer;
	vector<string> vecstrTokens;
	
	while (!ifsm.eof()) {
		ifsm.getline(acBuffer, c_iBuffer - 1);
		acBuffer[c_iBuffer - 1] = 0;
		vecstrTokens.clear();
		CMeta::Tokenize(acBuffer, vecstrTokens);
		if (vecstrTokens.empty())
			continue;
		if (vecstrTokens.size() != 2) {
			cerr << "Illegal line (" << vecstrTokens.size() << "): "
					<< acBuffer << endl;
			continue;
		}
		
		if (acBuffer[0] == '#') {
		  cerr << "skipping " << acBuffer << endl;
		} else {		  		  
		  mapGene2Fold[ vecstrTokens[0] ] = atoi( vecstrTokens[1].c_str() );
		}
	}
	return true;
}

// Read in the 
bool ReadProbParamFile(char* prob_file, float& A, float& B) {
	static const size_t c_iBuffer = 1024;
	char acBuffer[c_iBuffer];
	char* nameBuffer;
	vector<string> vecstrTokens;
	size_t i, extPlace;
	string Ext, FileName;
	size_t index = 0;
	ifstream ifsm;
	
	ifsm.open( prob_file );
	i = 0;
	while (!ifsm.eof()) {
		ifsm.getline(acBuffer, c_iBuffer - 1);
		acBuffer[c_iBuffer - 1] = 0;
		vecstrTokens.clear();
		CMeta::Tokenize(acBuffer, vecstrTokens);
		if (vecstrTokens.empty())
			continue;
		if (vecstrTokens.size() > 1) {
			cerr << "Illegal model line (" << vecstrTokens.size() << "): "
					<< acBuffer << endl;
			continue;
		}
		if (acBuffer[0] == '#') {
		  cerr << "skipping " << acBuffer << endl;
		} else {
		  if( i == 0 )		  
		    A = atof(vecstrTokens[0].c_str());
		  else if( i == 1 )
		    B = atof(vecstrTokens[0].c_str());
		  else{
		    cerr << "" << endl;
		    return false;
		  }		  
		  i++;
		}		
	}
	cerr << "Reading Prob file, A: " << A << ", B: " << B << endl;
	return true;
}

// Platt's binary SVM Probablistic Output
// Assume dec_values and labels have same dimensions and genes
static void sigmoid_train(CDat& dec_values, 
			  CDat& labels, 
			  float& A, float& B){
	double prior1=0, prior0 = 0;
	size_t i, j, idx, k;
	float d, lab;
	
	int max_iter=100;	// Maximal number of iterations
	double min_step=1e-10;	// Minimal step taken in line search
	double sigma=1e-12;	// For numerically strict PD of Hessian
	double eps=1e-5;
	vector<double> t;
	double fApB,p,q,h11,h22,h21,g1,g2,det,dA,dB,gd,stepsize;
	double newA,newB,newf,d1,d2;
	int iter; 
	
	// Negatives are values less than 0
	for(i = 0; i < dec_values.GetGenes(); i++)
	  for(j = (i+1); j < dec_values.GetGenes(); j++)
	    if (!CMeta::IsNaN(d = dec_values.Get(i, j)) && !CMeta::IsNaN(lab = labels.Get(i, j))  ){
	      if(lab > 0)
		prior1 += 1;
	      else if(lab < 0)
		prior0 += 1;	      
	    }
	
	// initialize size
	t.resize(prior0+prior1);
	
	// Initial Point and Initial Fun Value
	A=0.0; B=log((prior0+1.0)/(prior1+1.0));
	double hiTarget=(prior1+1.0)/(prior1+2.0);
	double loTarget=1/(prior0+2.0);			
	double fval = 0.0;
		
	for(idx = i = 0; idx < dec_values.GetGenes(); idx++)
	  for(j = (idx+1); j < dec_values.GetGenes(); j++)
	    if (!CMeta::IsNaN(d = dec_values.Get(idx, j)) && !CMeta::IsNaN(lab = labels.Get(idx, j))  ){
	      if (lab > 0 ) t[i]=hiTarget;
	      else t[i]=loTarget;
	      	      
	      fApB = d*A+B;
	      if (fApB>=0)
		fval += t[i]*fApB + log(1+exp(-fApB));
	      else
		fval += (t[i] - 1)*fApB +log(1+exp(fApB));	    
	      ++i;
	    }
	
	for (iter=0;iter<max_iter;iter++){
	  // Update Gradient and Hessian (use H' = H + sigma I)
	  h11=sigma; // numerically ensures strict PD
	  h22=sigma;
	  h21=0.0;g1=0.0;g2=0.0;
	  
	  for(i = idx = 0; idx < dec_values.GetGenes(); idx++)
	    for(j = (idx+1); j < dec_values.GetGenes(); j++)
	      if (!CMeta::IsNaN(d = dec_values.Get(idx, j)) && !CMeta::IsNaN(lab = labels.Get(idx, j))  ){			    	    
		fApB = d*A+B;		
		
		if (fApB >= 0){
		  p=exp(-fApB)/(1.0+exp(-fApB));
		  q=1.0/(1.0+exp(-fApB));
		}
		else{
		  p=1.0/(1.0+exp(fApB));
		  q=exp(fApB)/(1.0+exp(fApB));
		}
		d2=p*q;
		h11+=d*d*d2;
		h22+=d2;
		h21+=d*d2;
		d1=t[i]-p;
		g1+=d*d1;
		g2+=d1;
		
		++i;
	      }
	  
	  // Stopping Criteria
	  if (fabs(g1)<eps && fabs(g2)<eps)
	    break;
	  
	  // Finding Newton direction: -inv(H') * g
	  det=h11*h22-h21*h21;
	  dA=-(h22*g1 - h21 * g2) / det;
	  dB=-(-h21*g1+ h11 * g2) / det;
	  gd=g1*dA+g2*dB;
	  
	  stepsize = 1;		// Line Search
	  while (stepsize >= min_step){
	    newA = A + stepsize * dA;
	    newB = B + stepsize * dB;
	    
	    // New function value
	    newf = 0.0;
	    
	    for(i = idx = 0; idx < dec_values.GetGenes(); idx++)
	      for(j = (idx+1); j < dec_values.GetGenes(); j++)
		if (!CMeta::IsNaN(d = dec_values.Get(idx, j)) && !CMeta::IsNaN(lab = labels.Get(idx, j))  ){			    	    
		  fApB = d*newA+newB;
		  
		  if (fApB >= 0)
		    newf += t[i]*fApB + log(1+exp(-fApB));
		  else
		    newf += (t[i] - 1)*fApB +log(1+exp(fApB));
		  
		  ++i;
		}
	    
	    // Check sufficient decrease
	    if (newf<fval+0.0001*stepsize*gd){
	      A=newA;B=newB;fval=newf;
	      break;
	    }
	    else
	      stepsize = stepsize / 2.0;
	  }
	  
	  if (stepsize < min_step){
	    cerr << "Line search fails in two-class probability estimates: " << stepsize << ',' << min_step << endl;
	    break;
	  }
	}
	
	if (iter>=max_iter)
	  cerr << "Reaching maximal iterations in two-class probability estimates" << endl;	
}

static void sigmoid_predict(CDat& dec_values, float A, float B){
  size_t i, j;
  float d, fApB;
  
  for(i = 0; i < dec_values.GetGenes(); i++)
    for(j = (i+1); j < dec_values.GetGenes(); j++)
      if (!CMeta::IsNaN(d = dec_values.Get(i, j))){			    	    	
	fApB = d*A+B;
	// 1-p used later; avoid catastrophic cancellation
	if (fApB >= 0)
	  dec_values.Set(i,j, exp(-fApB)/(1.0+exp(-fApB)));
	else
	  dec_values.Set(i,j, 1.0/(1+exp(fApB)));
      }
}


int main(int iArgs, char** aszArgs) {
	gengetopt_args_info sArgs;	
	SVMLight::CSVMPERF SVM;
	
	size_t i, j, iGene, jGene, numpos, numneg, iSVM;
	ifstream ifsm;
	float d, sample_rate, dval;
	int   iRet;
	map<string, size_t>	mapstriZeros, mapstriDatasets;
	vector<string> vecstrDatasets;
	vector<bool> mapTgene;
	vector<bool> mapCgene;
	vector<size_t> mapTgene2fold;
	vector<int> tgeneCount;
	
	map<string, size_t> mapGene2Fold;

	DIR* dp;
	struct dirent* ep;	
	CGenome Genome;
        CGenes Genes(Genome);
	
	CGenome GenomeTwo;
        CGenes Context(GenomeTwo);
	
	CGenome GenomeThree;
        CGenes Allgenes(GenomeThree);
	
	CGenome GenomeFour;
        CGenes labelgenes(GenomeFour);
	
	if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
		cmdline_parser_print_help();
		return 1;
	}
	
	CMeta Meta( sArgs.verbosity_arg, sArgs.random_arg );
	
	SVM.SetVerbosity(sArgs.verbosity_arg);
	SVM.SetLossFunction(sArgs.error_function_arg);
	if (sArgs.k_value_arg > 1) {
		cerr << "k_value is >1. Setting default 0.5" << endl;
		SVM.SetPrecisionFraction(0.5);
	} else if (sArgs.k_value_arg <= 0) {
		cerr << "k_value is <=0. Setting default 0.5" << endl;
		SVM.SetPrecisionFraction(0.5);
	} else {
		SVM.SetPrecisionFraction(sArgs.k_value_arg);
	}

	
	if (sArgs.cross_validation_arg < 1){
	  cerr << "cross_valid is <1. Must be set at least 1" << endl;
	  return 1;
	}
	else if(sArgs.cross_validation_arg < 2){
	  cerr << "cross_valid is set to 1. No cross validation holdouts will be run." << endl;
	}
	
	SVM.SetTradeoff(sArgs.tradeoff_arg);
	if (sArgs.slack_flag)
		SVM.UseSlackRescaling();
	else
		SVM.UseMarginRescaling();
	
	
	if (!SVM.parms_check()) {
		cerr << "Sanity check failed, see above errors" << endl;
		return 1;
	}
	
	// read in the list of datasets
	if(sArgs.directory_arg ) {
	  dp = opendir (sArgs.directory_arg);
	  if (dp != NULL){
	    while (ep = readdir (dp)){
	      // skip . .. files and temp files with ~
	      if (ep->d_name[0] == '.' || ep->d_name[strlen(ep->d_name)-1] == '~') 
		continue;
	      
	      // currently opens all files. Add filter here if want pick file extensions
	      vecstrDatasets.push_back((string)sArgs.directory_arg + "/" + ep->d_name);	      
	    }
	    (void) closedir (dp);	    
	    	    
	    cerr << "Input Dir contrains # datasets: " << vecstrDatasets.size() << '\n';
	    // sort datasets in alphabetical order
	    std::sort(vecstrDatasets.begin(), vecstrDatasets.end());
	  }
	  else{
	    cerr << "Couldn't open the directory: " << sArgs.directory_arg << '\n';
	    return 1;
	  }	  
	}
	
	// read in the gene list
	if( sArgs.genes_arg ) {
	  ifsm.open( sArgs.genes_arg );
	  if( !labelgenes.Open( ifsm ) ) {
	    cerr << "Could not open: " << sArgs.genes_arg << endl;
	    return 1; }
	  ifsm.close( ); }
	
	// read target gene list
	if(sArgs.tgene_given ) {
	  ifstream ifsm;
	  ifsm.open(sArgs.tgene_arg);
	  
	  if (!Genes.Open(ifsm)) {
	    cerr << "Could not open: " << sArgs.tgene_arg << endl;
	    return 1;
	  }
	  ifsm.close();
	}

	// read context gene list
	if(sArgs.context_given ) {
	  ifstream ifsm;
	  ifsm.open(sArgs.context_arg);
	  
	  if (!Context.Open(ifsm)) {
	    cerr << "Could not open: " << sArgs.context_arg << endl;
	    return 1;
	  }
	  ifsm.close();
	}

	// read all gene list
	// IF given this flag predict for all gene pairs
	if(sArgs.allgenes_given ) {
	  ifstream ifsm;
	  ifsm.open(sArgs.allgenes_arg);
	  
	  if (!Allgenes.Open(ifsm)) {
	    cerr << "Could not open: " << sArgs.allgenes_arg << endl;
	    return 1;
	  }
	  ifsm.close();
	}

	// read the gene holdout fold
	if( sArgs.GenesHoldoutFold_given ){
	  ifstream ifsm;
	  ifsm.open( sArgs.GenesHoldoutFold_arg );
	  
	  if (!ReadGenesHoldoutFoldFile(ifsm, mapGene2Fold) ) {
	    cerr << "Could not open: " << sArgs.GenesHoldoutFold_arg << endl;
	    return 1;
	  }
	  
	  ifsm.close();	  	  
	}
	
	///######################
	// Chris added
	vector<SVMLight::SVMLabelPair*> vecLabels;
	CDat Labels;
	CDat Results;
	

	///######################
	/// Read in labels and apply various modifications based on user flags/args
	/// 
	if ( sArgs.labels_given ) {	  
	  if (!Labels.Open(sArgs.labels_arg, sArgs.mmap_flag)) {
	    cerr << "Could not open input labels Dat" << endl;
	    return 1;
	  }

	  if( labelgenes.GetGenes( ) )
	    Labels.FilterGenes( labelgenes, CDat::EFilterInclude );
	  	  
	  // random sample labels
	  if( sArgs.subsample_given ){
	    cerr << "Sub-sample labels to rate:" << sArgs.subsample_arg << endl;
	    for( i = 0; i < Labels.GetGenes( ); ++i )
	      for( j = ( i + 1 ); j < Labels.GetGenes( ); ++j )
		if( !CMeta::IsNaN( Labels.Get( i, j ) ) &&
		    ( ( (float)rand( ) / RAND_MAX ) > sArgs.subsample_arg ) )
		  Labels.Set( i, j, CMeta::GetNaN( ) );	    	    
	  }	  
	  
	  // set all NaN values to negatives
	  if( sArgs.nan2neg_given ){
	    cerr << "Set NaN labels dat as negatives" << endl;
	    
	    for(i = 0; i < Labels.GetGenes(); i++)
	      for(j = (i+1); j < Labels.GetGenes(); j++)
		if (CMeta::IsNaN(d = Labels.Get(i, j)))  
		  Labels.Set(i, j, -1);
	  }
	  
	  if( sArgs.tgene_given ){
	    mapTgene.resize(Labels.GetGenes());
	    
	    for(i = 0; i < Labels.GetGenes(); i++){
	      if(Genes.GetGene(Labels.GetGene(i)) == -1)
		mapTgene[i] = false;
	      else
		mapTgene[i] = true;
	    }
	    
	    // keep track of positive gene counts
	    tgeneCount.resize(Labels.GetGenes());
	    
	    // if given a target gene file
	    // Only keep eges that have only one gene in this targe gene list
	    if( sArgs.onetgene_flag ){
	      cerr << "Filtering to only include edges with one gene in gene file: " << sArgs.tgene_arg << endl;
	      
	      for(i = 0; i < Labels.GetGenes(); i++)
		for(j = (i+1); j < Labels.GetGenes(); j++)
		  if (!CMeta::IsNaN(d = Labels.Get(i, j))){
		    if(mapTgene[i] && mapTgene[j])
		      Labels.Set(i, j, CMeta::GetNaN());
		    else if(!mapTgene[i] && !mapTgene[j])
		      Labels.Set(i, j, CMeta::GetNaN());
		  }
	    }
	  } // if edgeholdout flag not given, we are doing gene holdout by default. 
	  // Since target gene list was not given we are using all genes in labels as target genes to cross holdout
	  else if( !sArgs.edgeholdout_flag ){
	    mapTgene.resize(Labels.GetGenes());
	    
	    // all genes are target genes
	    for(i = 0; i < Labels.GetGenes(); i++){
	      mapTgene[i] = true;
	    }
	    
	    // keep track of positive gene counts
	    tgeneCount.resize(Labels.GetGenes());
	  }
	  
	  //if given a context map the context genes
	  if( sArgs.context_given ){
	    mapCgene.resize(Labels.GetGenes());
	    
	    for(i = 0; i < Labels.GetGenes(); i++){
	      if(Context.GetGene(Labels.GetGene(i)) == -1)
		mapCgene[i] = false;
	      else
		mapCgene[i] = true;
	    }
	  }

	  // Set target prior
	  if(sArgs.prior_given){
	    numpos = 0;
	    numneg = 0;
	    for(i = 0; i < Labels.GetGenes(); i++)
	      for(j = (i+1); j < Labels.GetGenes(); j++)
		if (!CMeta::IsNaN(d = Labels.Get(i, j))){
		  if(d > 0){
		    ++numpos;}
		  else if(d < 0){
		    ++numneg;
		  }
		}

	    if( ((float)numpos / (numpos + numneg)) < sArgs.prior_arg){
	      
	      cerr << "Convert prior from orig: " << ((float)numpos / (numpos + numneg)) << " to target: " << sArgs.prior_arg << endl;
	      
	      sample_rate = ((float)numpos / (numpos + numneg)) / sArgs.prior_arg;
	      
	      // remove neg labels to reach prior
	      for(i = 0; i < Labels.GetGenes(); i++)
		for(j = (i+1); j < Labels.GetGenes(); j++)
		  if (!CMeta::IsNaN(d = Labels.Get(i, j)) && d < 0){
		    if((float)rand() / RAND_MAX  > sample_rate)
		      Labels.Set(i, j, CMeta::GetNaN());
		  }
	    }
	  }
	  
	  // output sample labels for eval/debug purpose
	  if(sArgs.OutLabels_given){
	    cerr << "save sampled labels as (1,0)s to: " << sArgs.OutLabels_arg << endl;
	    Labels.Normalize( CDat::ENormalizeMinMax );
	    Labels.Save(sArgs.OutLabels_arg);
	    return 0;
	  }
	  
	  // Exclude labels without context genes
	  if(sArgs.context_given && !sArgs.allContextPred_flag){
	    if( sArgs.touchContext_flag ){
	      Labels.FilterGenes( sArgs.context_arg, CDat::EFilterEdge );
	    }else{
	      Labels.FilterGenes( Context, CDat::EFilterInclude );
	    }
	  }
	  
	  // If not given a SVM model/models we are in learning mode, thus construct each SVMLabel object for label
	  if( !sArgs.model_given && !sArgs.modelPrefix_given ){
	    numpos = 0;
	    for(i = 0; i < Labels.GetGenes(); i++)
	      for(j = (i+1); j < Labels.GetGenes(); j++)
		if (!CMeta::IsNaN(d = Labels.Get(i, j))){
		  if (d != 0)  
		    vecLabels.push_back(new SVMLight::SVMLabelPair(d, i, j));
		  if(d > 0)
		    ++numpos;
		}
	    
	    // check to see if you have enough positives to learn from
	    if(sArgs.mintrain_given && sArgs.mintrain_arg > numpos){
	      cerr << "Not enough positive examples from: " << sArgs.labels_arg << " numpos: " << numpos << endl;
	      return 1;
	    }
	  }
	  
	  // save sampled labels 
	  if(sArgs.SampledLabels_given) {
	    cerr << "Save sampled labels file: " << sArgs.SampledLabels_arg << endl;
	    Labels.Save(sArgs.SampledLabels_arg);	    
	  }	  
	} /// Done with reading labels
	
	SVMLight::SAMPLE* pTrainSample;
	SVMLight::SAMPLE* pAllSample;
	vector<SVMLight::SVMLabelPair*> pTrainVector[sArgs.cross_validation_arg];
	vector<SVMLight::SVMLabelPair*> pTestVector[sArgs.cross_validation_arg];
	
	// ###################################
	//
	// Now conduct Learning if given a labels and no SVM models	
	//
	if (sArgs.output_given && sArgs.labels_given && !sArgs.modelPrefix_given && !sArgs.model_given) {
	  //do learning and classifying with cross validation
	  if( sArgs.cross_validation_arg > 1){
	    cerr << "setting cross validation holds" << endl;
	    
	    mapTgene2fold.resize(mapTgene.size());
	    
	    // assign target genes to there cross validation fold
	    if(sArgs.tgene_given || !sArgs.edgeholdout_flag){
	      for(i = 0; i < mapTgene.size(); i++){
		if(!mapTgene[i]){
		  mapTgene2fold[i] = -1; 
		  continue;
		}
		
		if( sArgs.GenesHoldoutFold_given ){
		  // Does not check if this gene has been assigned a random fold
		  mapTgene2fold[i] = mapGene2Fold[ Labels.GetGene(i) ];
		}else{
		  mapTgene2fold[i] = rand() % sArgs.cross_validation_arg;
		}
		
	      }
	      
	      // cross-fold by target gene
	      for (i = 0; i < sArgs.cross_validation_arg; i++) {
		cerr << "cross validation holds setup:" << i << endl;
		
		for (j = 0; j < vecLabels.size(); j++) {
		  //if( j % 1000 == 0)
		  //cerr << "cross validation push labels:" << j << endl;
		  
		  if(mapTgene[vecLabels[j]->iidx] || mapTgene[vecLabels[j]->jidx]){
		    
		    if(mapTgene2fold[vecLabels[j]->iidx] == i || mapTgene2fold[vecLabels[j]->jidx] == i){
		      
		      // only add if both genes are in context
		      if( sArgs.context_given  && 
			  !sArgs.allContextPred_flag  &&
			  !sArgs.touchContext_flag &&
			  ( !mapCgene[vecLabels[j]->iidx] || !mapCgene[vecLabels[j]->jidx])){

			if( !( sArgs.onlyPos_flag && vecLabels[j]->Target < 0 ) )
			  continue;
		      }
		      
		      if( sArgs.context_given  && 
			  !sArgs.allContextPred_flag  &&
			  sArgs.touchContext_flag &&
			  ( !mapCgene[vecLabels[j]->iidx] && !mapCgene[vecLabels[j]->jidx])){
			
			if( !( sArgs.onlyPos_flag && vecLabels[j]->Target < 0 ) )
			  continue;			
		      }

		      pTestVector[i].push_back(vecLabels[j]);
		    }else{		      		      

		      // only add if both genes are in context
		      if( sArgs.context_given  && 
			  !sArgs.touchContext_flag &&
			  ( !mapCgene[vecLabels[j]->iidx] || !mapCgene[vecLabels[j]->jidx])){

			if( !( sArgs.onlyPos_flag && vecLabels[j]->Target < 0 ) )
			  continue;

		      }
		      
		      // only add if both genes are in context
		      if( sArgs.context_given  && 
			  sArgs.touchContext_flag &&
			  ( !mapCgene[vecLabels[j]->iidx] && !mapCgene[vecLabels[j]->jidx])){
			
			if( !( sArgs.onlyPos_flag && vecLabels[j]->Target < 0 ) )
			  continue;
		      }
		      
		      pTrainVector[i].push_back(vecLabels[j]); 
		    }	
		    /*else if(mapTgene2fold[vecLabels[j]->iidx] != i && mapTgene2fold[vecLabels[j]->jidx] != i){		      		      
		      pTrainVector[i].push_back(vecLabels[j]); 
		      }*/    
		  }else{
		    cerr << "Error: edge exist without a target gene" << endl; 
		    return 1;
		  }
		}		
		
		cerr << "test,"<< i <<": " << pTestVector[i].size() << endl;
		int numpos = 0;
		int numpos_test = 0;
		for(j=0; j < pTrainVector[i].size(); j++)
		  if(pTrainVector[i][j]->Target > 0)
		    ++numpos;
		for(j=0; j < pTestVector[i].size(); j++)
		  if(pTestVector[i][j]->Target > 0)
		    ++numpos_test;
		
		if( numpos < 1 || (sArgs.mintrain_given && sArgs.mintrain_arg > numpos) ){						
		  cerr << "Not enough positive examples from fold: " << i  << " file: " << sArgs.labels_arg << " numpos: " << numpos << endl;
		  return 1;
		}
		
		cerr << "train,"<< i <<","<<numpos<<": " << pTrainVector[i].size() << endl;
		cerr << "test,"<< i <<","<<numpos_test<<": " <<  pTestVector[i].size() << endl;
	      }
	    }
	    else{ //randomly set eges into cross-fold
	      cerr << "Edge holdout" << endl;
	      /*
		if( sArgs.context_given ){
		cerr << "context not implemented yet for random edge holdout" << endl;
		return 1;
		}
	      */
	      
	      for (i = 0; i < sArgs.cross_validation_arg; i++) {
		pTestVector[i].reserve((size_t) vecLabels.size()
				       / sArgs.cross_validation_arg + sArgs.cross_validation_arg);
		pTrainVector[i].reserve((size_t) vecLabels.size()
					/ (sArgs.cross_validation_arg)
					* (sArgs.cross_validation_arg - 1)
					+ sArgs.cross_validation_arg);
		for (j = 0; j < vecLabels.size(); j++) {
		  if (j % sArgs.cross_validation_arg == i) {
		    pTestVector[i].push_back(vecLabels[j]);
		  }else{
		    if( sArgs.context_given ){
		      if( mapCgene[vecLabels[j]->iidx] && mapCgene[vecLabels[j]->jidx] )
			pTrainVector[i].push_back((vecLabels[j]));
		    }else{
		      pTrainVector[i].push_back((vecLabels[j]));
		    }
		  }
		}
	      	
		// print out number of examples
		int numpos = 0;
		int numpos_test = 0;
		for(j=0; j < pTrainVector[i].size(); j++)
		  if(pTrainVector[i][j]->Target > 0)
		    ++numpos;
		for(j=0; j < pTestVector[i].size(); j++)
		  if(pTestVector[i][j]->Target > 0)
		    ++numpos_test;		
		cerr << "train: "<< i <<" , "<<numpos<<": " << pTrainVector[i].size() << endl;
		cerr << "test: "<< i <<" , "<<numpos_test<<": " <<  pTestVector[i].size() << endl;
		
	      }      	      
	    }
	  }
	  else{ // if you have less than 2 fold cross, no cross validation is done, all train genes are used and predicted
	    
	    // no holdout so train is the same as test gene set
	    pTestVector[0].reserve((size_t) vecLabels.size() + sArgs.cross_validation_arg);
	    pTrainVector[0].reserve((size_t) vecLabels.size() + sArgs.cross_validation_arg);
	    
	    for (j = 0; j < vecLabels.size(); j++) {
	      pTestVector[0].push_back(vecLabels[j]);		      
	      pTrainVector[0].push_back(vecLabels[j]);		    
	    }
	  }
	  
	  // initalize the results
	  Results.Open(Labels.GetGeneNames(), true);
	  
	  /////
	  // Create feature vectors for all Label pairs using input datasets
	  //
	  cerr << "CreateDocs!"<< endl;
	  if(sArgs.normalizeZero_flag){
	    SVMLight::CSVMPERF::CreateDoc(vecstrDatasets,
					  vecLabels,
					  Labels.GetGeneNames(),
					  Sleipnir::CDat::ENormalizeMinMax);
	  }else if(sArgs.normalizeNPone_flag){
	    SVMLight::CSVMPERF::CreateDoc(vecstrDatasets,
					  vecLabels,
					  Labels.GetGeneNames(),
					  Sleipnir::CDat::ENormalizeMinMaxNPone);
	  }else if(sArgs.zscore_flag){
	    SVMLight::CSVMPERF::CreateDoc(vecstrDatasets,
					  vecLabels,
					  Labels.GetGeneNames(),
					  Sleipnir::CDat::ENormalizeZScore);
	  }else{
	    SVMLight::CSVMPERF::CreateDoc(vecstrDatasets,
					  vecLabels,
					  Labels.GetGeneNames());
	  }
	  
	  ////
	  // Start learning for each cross validation fold
	  for (i = 0; i < sArgs.cross_validation_arg; i++) {
	    std::stringstream sstm;
	    
	    // build up the output SVM model file name
	    if(sArgs.context_given){
	      std::string path(sArgs.context_arg);
	      size_t pos = path.find_last_of("/");
	      std::string cname;
	      if(pos != std::string::npos)
		cname.assign(path.begin() + pos + 1, path.end());
	      else
		cname = path;
	      
	      sstm << sArgs.output_arg << "." << sArgs.tradeoff_arg  << "." << cname << "." << i << ".svm";		      
	    }else
	      sstm << sArgs.output_arg << "." << sArgs.tradeoff_arg  << "." << i << ".svm";
	    
	    cerr << "Cross validation fold: " << i << endl;	    	    
	    pTrainSample = SVMLight::CSVMPERF::CreateSample(pTrainVector[i]);
	    
	    // Skip learning if SVM model file already exist
	    if( sArgs.skipSVM_flag && 
		access((char*)(sstm.str().c_str()), R_OK ) != -1
		){
	      //SVM.ReadModel((char*)(sstm.str().c_str()));
	      SVM.ReadModelLinear((char*)(sstm.str().c_str()));
	      cerr << "Using existing trained SVM model: " << sstm.str() << endl;
	    }else{
	      SVM.Learn(*pTrainSample);
	      cerr << "SVM model Learned" << endl;
	    }
	    
	    SVM.Classify(Results,
			 pTestVector[i]);
	    
	    cerr << "SVM model classified holdout" << endl;
	    
	    if( sArgs.savemodel_flag ){
	      SVM.WriteModel((char*)(sstm.str().c_str()));
	    }
	    
	    // DEBUG
	    SVMLight::CSVMPERF::FreeSample_leave_Doc(*pTrainSample);
	    free(pTrainSample);
	  }
	  
	  if( sArgs.prob_flag ){
	    cerr << "Converting prediction values to estimated probablity" << endl;
	    float A, B;
	    
	    // TODO add function to read in prob parameter file if already existing
	    sigmoid_train(Results, Labels, A, B);
	    sigmoid_predict(Results, A, B);
	  }
	  else if( sArgs.probCross_flag ){
	    float A, B;
	    size_t k, ctrain, itrain;
	    vector<SVMLight::SVMLabelPair*> probTrainVector;
	    
	    for (i = 0; i < sArgs.cross_validation_arg; i++) {		      
	      cerr << "Convert to probability for cross fold: " << i << endl;		      
	      
	      // construct prob file name
	      std::stringstream pstm;		
	      ofstream ofsm;		      
	      if(sArgs.context_given){
		std::string path(sArgs.context_arg);
		size_t pos = path.find_last_of("/");
		std::string cname;
		if(pos != std::string::npos)
		  cname.assign(path.begin() + pos + 1, path.end());
		else
		  cname = path;
		
		pstm << sArgs.output_arg << "." << sArgs.tradeoff_arg  << "." << cname  << "." << i << ".svm.prob";		      
	      }else
		pstm << sArgs.output_arg << "." << sArgs.tradeoff_arg  << "." << i << ".svm.prob";
	      
	      
	      if( sArgs.skipSVM_flag && 
		  access((char*)(pstm.str().c_str()), R_OK ) != -1
		  ){
		
		// read in parameter file
		if(!ReadProbParamFile((char*)(pstm.str().c_str()), A, B)){
		  cerr << "Failed to read Probablity parameter file: " << pstm.str() << endl;
		  return 1;
		}
		cerr << pstm.str() << ": read in values, A: " << A << ", B: " << B << endl;		
		
	      }else{		
		ctrain = 0;
		for (j = 0; j < sArgs.cross_validation_arg; j++) {		      
		  if(i == j)
		    continue;
		  ctrain += pTrainVector[j].size();
		}
		
		probTrainVector.resize(ctrain);
		itrain = 0;
		for (j = 0; j < sArgs.cross_validation_arg; j++) {		      
		  if(i == j)
		    continue;
		  for(k = 0; k < pTrainVector[j].size(); k++){
		    probTrainVector[itrain] = pTrainVector[j][k];
		    itrain += 1;
		  }
		}
		
		// train A,B sigmoid perameters
		SVM.sigmoid_train(Results, probTrainVector, A, B);		
	      }
	      
	      SVM.sigmoid_predict(Results, pTestVector[i], A, B);
	      
	      // open prob param file
	      if(sArgs.savemodel_flag){
		ofsm.open(pstm.str().c_str());
		ofsm << A << endl;
		ofsm << B << endl;
		ofsm.close();
	      }
	    }
	  }
	  
	  // only save cross-validated results when not predicting all genes
	  if( !sArgs.allgenes_given )
	    Results.Save(sArgs.output_arg);
	}
	
	///##############
	// If given all genes arg, this puts in prediction mode for all gene pairs from the gene list
	//
	if ( sArgs.allgenes_given && sArgs.output_given) { //read model and classify all	  
	  size_t iData;
	  vector<vector<float> > vecSModel;
	  vecSModel.resize(sArgs.cross_validation_arg);
	  ifstream mifsm;
	  vector<CDat* > vecResults;
	  vector<size_t> veciGene;
	  
	  cerr << "Predicting for all genes given." << endl;
	  
	  // open SVM models for prefix file
	  if(sArgs.modelPrefix_given){
	    for(i = 0; i < sArgs.cross_validation_arg; i++){
	      std::stringstream sstm;
	      
	      sstm << sArgs.modelPrefix_arg << "." << i << ".svm";	      
	      if( access((char*)(sstm.str().c_str()), R_OK) == -1 ){
		cerr << "ERROR: SVM model file cannot be opned: " << sstm.str() << endl;
		return 1;
	      }
	      
	      mifsm.open((char*)(sstm.str().c_str()));	      		      
	      ReadModelFile(mifsm, vecSModel[i]);
	    }
	  }else if( sArgs.model_given ){ // open SVM model from file
	    //vector<float> SModel;
	    vecSModel.resize(1);
	    
	    if( access(sArgs.model_arg, R_OK) == -1 ){
	      cerr << "ERROR: SVM model file cannot be opned: " << sArgs.model_arg << endl;
	      return 1;
	    }
	    
	    mifsm.open(sArgs.model_arg);
	    ReadModelFile(mifsm, vecSModel[0]);
	    
	    // DEBUG check if this is ok
	    sArgs.cross_validation_arg = 1;
	  }else{ 
	    // open SVM model file from which was just trained
	    for(i = 0; i < sArgs.cross_validation_arg; i++){
	      std::stringstream sstm;
	      
	      if(sArgs.context_given){
		std::string path(sArgs.context_arg);
		size_t pos = path.find_last_of("/");
		std::string cname;
		if(pos != std::string::npos)
		  cname.assign(path.begin() + pos + 1, path.end());
		else
		  cname = path;
		
		sstm << sArgs.output_arg << "." << sArgs.tradeoff_arg  << "." << cname << "." << i << ".svm";		      
	      }else
		sstm << sArgs.output_arg << "." << sArgs.tradeoff_arg  << "." << i << ".svm";
	      
	      
	      if( access((char*)(sstm.str().c_str()), R_OK) == -1 ){
		cerr << "ERROR: SVM model file cannot be opned: " << sstm.str() << endl;
		return 1;
	      }
	      
	      mifsm.open((char*)(sstm.str().c_str()));	      
	      ReadModelFile(mifsm, vecSModel[i]);
	    }
	  }
	  
	  // Initialize output for input all gene list
	  vecResults.resize(sArgs.cross_validation_arg);
	  for(i = 0; i < sArgs.cross_validation_arg; i++){
	    vecResults[ i ] = new CDat();
	    vecResults[ i ]->Open(Allgenes.GetGeneNames( ), true);
	  }
	  
	  	  
	  // Now iterate over all datasets to make predictions for all gene pairs
	  CDat wDat;		    
	  for(iData = 0; iData < vecstrDatasets.size(); iData++){	    
	    if(!wDat.Open(vecstrDatasets[iData].c_str(), sArgs.mmap_flag)) {
	      cerr << vecstrDatasets[iData].c_str() << endl;
	      cerr << "Could not open: " << vecstrDatasets[iData] << endl;
	      return 1;
	    }
	    
	    cerr << "Open: " << vecstrDatasets[iData] << endl;
	    
	    // normalize data file
	    if(sArgs.normalizeZero_flag){
	      cerr << "Normalize input [0,1] data" << endl;
	      wDat.Normalize( Sleipnir::CDat::ENormalizeMinMax );
	    }else if(sArgs.normalizeNPone_flag){
	      cerr << "Normalize input [-1,1] data" << endl;
	      wDat.Normalize( Sleipnir::CDat::ENormalizeMinMaxNPone );
	    }else if(sArgs.zscore_flag){
	      cerr << "Normalize input data to zscore" << endl;
	      wDat.Normalize( Sleipnir::CDat::ENormalizeZScore );
	    }
	    
	    // map result gene list to dataset gene list
	    veciGene.resize( vecResults[ 0 ]->GetGenes() );	    
	    for(i = 0; i < vecResults[ 0 ]->GetGenes(); i++){
	      veciGene[ i ] = wDat.GetGene( vecResults[ 0 ]->GetGene( i ) );
	    }
	    
	    // compute prediction component for this dataset/SVM model
	    for(i = 0; i < vecResults[ 0 ]->GetGenes(); i++){
	      iGene = veciGene[i];
	      
	      for(j = i+1; j < vecResults[ 0 ]->GetGenes(); j++){
		jGene = veciGene[j];
		
		// if no value, set feature to 0
		if( iGene == -1 || jGene == -1 || CMeta::IsNaN(d = wDat.Get(iGene, jGene) ) )
		  d = 0;
		
		// iterate each SVM model
		for(iSVM = 0; iSVM < sArgs.cross_validation_arg; iSVM++){
		  if( CMeta::IsNaN(dval = vecResults[ iSVM ]->Get(i, j)) )
		    vecResults[ iSVM ]->Set(i, j, (d * vecSModel[iSVM][iData])  );
		  else
		    vecResults[ iSVM ]->Set(i, j, (dval + (d * vecSModel[iSVM][iData])) );
		  
		}
		
	      }
	    }	    	    
	  }
	  
	  
	  // convert the SVM predictions for each model to Probablity if required
	  if( sArgs.prob_flag || sArgs.probCross_flag ){
	    // iterate over each SVM model and its output and convert to probablity
	    float A;
	    float B;
	    
	    cerr << "convert prediction dabs to probablity" << endl;
	    
	    for(iSVM = 0; iSVM < sArgs.cross_validation_arg; iSVM++){
	      // Read A,B perameters
	      std::stringstream pstm;
	      if(sArgs.modelPrefix_given){
		pstm << sArgs.modelPrefix_arg  << "." << iSVM << ".svm.prob";
	      }else if( sArgs.model_given ){ // open SVM model from file
		
		// DEBUG, this might need to looked at!!!
		pstm << sArgs.model_arg << ".prob";
		
	      }else{ // open SVM model file from which was just trained		  
		if(sArgs.context_given){
		  std::string path(sArgs.context_arg);
		  size_t pos = path.find_last_of("/");
		  std::string cname;
		  if(pos != std::string::npos)
		    cname.assign(path.begin() + pos + 1, path.end());
		  else
		    cname = path;
		  
		  pstm << sArgs.output_arg << "." << sArgs.tradeoff_arg  << "." << cname  << "." << iSVM << ".svm.prob";		      
		}else
		  pstm << sArgs.output_arg << "." << sArgs.tradeoff_arg  << "." << iSVM << ".svm.prob";
	      }
	      
	      // read in parameter file
	      if(!ReadProbParamFile((char*)(pstm.str().c_str()), A, B)){
		cerr << "Failed to read Probablity parameter file: " << pstm.str() << endl;
		return 1;
	      }
	      
	      cerr << pstm.str() << ", A: " << A << ",B: " << B << endl;
	  
	      // debug
	      //cerr << "known:" << vecResults[ iSVM ]->Get(x,y) << endl;
	      
	      // now convert the SVM model prediction values to probablity based on param A,B
	      sigmoid_predict( *(vecResults[ iSVM ]), A, B);
	      
	      //cerr << "known PROB:" << vecResults[ iSVM ]->Get(x,y) << endl;
	    }
	  }
	  
	  // filter results
	  if( sArgs.tgene_given ){
	    for(iSVM = 0; iSVM < sArgs.cross_validation_arg; iSVM++){
	      vecResults[ iSVM ]->FilterGenes( sArgs.tgene_arg, CDat::EFilterEdge );
	    }
	  }
	  
	  // Exclude pairs without context genes
	  if(sArgs.context_given && !sArgs.allContextPred_flag){
	    for(iSVM = 0; iSVM < sArgs.cross_validation_arg; iSVM++){
	      if( sArgs.touchContext_flag ){
		vecResults[ iSVM ]->FilterGenes( sArgs.context_arg, CDat::EFilterEdge );
	      }else{
		vecResults[ iSVM ]->FilterGenes( Context, CDat::EFilterInclude );
	      }
	    }
	  }
	  
	  // Take average of prediction value from each model
	  // The final results will be stored/overwritten into the first dab (i.e. vecResults[ 0 ])
	  for(i = 0; i < vecResults[ 0 ]->GetGenes(); i ++)
	    for(j = i+1; j < vecResults[ 0 ]->GetGenes(); j++){
	      if (CMeta::IsNaN(dval = vecResults[ 0 ]->Get(i, j)))
		continue;
	      
	      // start from the second
	      // Assume all SVM model prediction dabs have identical NaN locations
	      if(sArgs.aggregateMax_flag){ // take the maximum prediction value
		float maxval = dval;
		for(iSVM = 1; iSVM < sArgs.cross_validation_arg; iSVM++)
		  if( vecResults[ iSVM ]->Get(i, j) > dval )
		    dval = vecResults[ iSVM ]->Get(i, j);
		
		vecResults[ 0 ]->Set(i, j, dval);
	      }else{ // Average over prediction values
		for(iSVM = 1; iSVM < sArgs.cross_validation_arg; iSVM++)
		  dval += vecResults[ iSVM ]->Get(i, j);	      
		
		vecResults[ 0 ]->Set(i, j, (dval / sArgs.cross_validation_arg) );
	      }
	    }
	  
	  // Replace gene-pair prediction values with labels from the cross-validation result
	  // This basically replaces the prediction value with one prediction value with which this label was heldout
	  // Only do this if cross-validation was conducted or given a replacement dab
	  if( (!sArgs.NoCrossPredict_flag && sArgs.output_given && sArgs.labels_given && !sArgs.modelPrefix_given && !sArgs.model_given) || 
	      sArgs.CrossResult_given ){
	    
	    if( sArgs.CrossResult_given ){	      
	      if(!Results.Open(sArgs.CrossResult_arg, sArgs.mmap_flag)) {
		cerr << "Could not open: " << sArgs.CrossResult_arg << endl;
		return 1;
	      }
	    }
	    
	    cerr << "Label pairs set to cross-validation values" << endl;
	    
	    // map result gene list to dataset gene list
	    veciGene.resize( vecResults[ 0 ]->GetGenes() );	    
	    for(i = 0; i < vecResults[ 0 ]->GetGenes(); i++){
	      veciGene[ i ] = Results.GetGene( vecResults[ 0 ]->GetGene( i ) );
	    }
	    
	    // compute prediction component for this dataset/SVM model
	    for(i = 0; i < vecResults[ 0 ]->GetGenes(); i++){
	      if( (iGene = veciGene[i]) == -1 )
		continue;
	      
	      for(j = i+1; j < vecResults[ 0 ]->GetGenes(); j++){
		if( (jGene = veciGene[j]) == -1 )
		  continue;
		
		if( CMeta::IsNaN(d = Results.Get(iGene, jGene) ) )
		  continue;
		
		vecResults[ 0 ]->Set(i, j, d);
	      }
	    }
	  }
	  
	  // now save the averged prediction values
	  vecResults[ 0 ]->Save(sArgs.output_arg);
	}
	
}


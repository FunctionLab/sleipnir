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


// Platt's binary SVM Probablistic Output
// Assume dec_values and labels have same dimensions and genes
static void sigmoid_train(CDat& dec_values, 
			  CDat& labels, 
			  float& A, float& B){
	double prior1=0, prior0 = 0;
	size_t i, j, idx;
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
	
	size_t i, j, iGene, jGene, numpos, numneg;
	ifstream ifsm;
	float d, sample_rate;
	int   iRet;
	map<string, size_t>	mapstriZeros, mapstriDatasets;
	vector<string> vecstrDatasets;
	vector<bool> mapTgene;
	vector<bool> mapCgene;
	vector<size_t> mapTgene2fold;
	vector<int> tgeneCount;
	
	DIR* dp;
	struct dirent* ep;	
	CGenome Genome;
        CGenes Genes(Genome);
	
	CGenome GenomeTwo;
        CGenes Context(GenomeTwo);
	
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

	//  cout << "there are " << vecLabels.size() << " labels processed" << endl;
	/*
	if (sArgs.input_given) {
		if (!PCL.Open(sArgs.input_arg, sArgs.skip_arg, sArgs.mmap_flag)) {
		cerr << "Could not open input PCL" << endl;
		return 1;
		}
	}
	*/
	
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
	  }
	  else{
	    cerr << "Couldn't open the directory: " << sArgs.directory_arg << '\n';
	    return 1;
	  }	  
	}
	
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
		
	///######################
	// Chris added
	vector<SVMLight::SVMLabelPair*> vecLabels;
	/*********
	 ** PCL gold standard later
	 
	CPCL Labels;
	if (sArgs.labels_given) {	  
	  if (!Labels.Open(sArgs.labels_arg, sArgs.skip_arg, sArgs.mmap_flag)) {
	    cerr << "Could not open input PCL" << endl;
	    return 1;
	  }
	  
	  for(size_t i = 0; i < PCL.GetGenes(); i++)
	    for(size_t j = 0; j < PCL.GetExperiments(); j++)
	      if (!CMeta::IsNaN(d = PCL.Get(i, j)))  
		vecLabels.push_back(SVMLight::SVMLabelPairPair(d, i, j));
	  
	}
	****/
	CDat Labels;
	CDat Results;
	if (sArgs.labels_given) {	  
	  if (!Labels.Open(sArgs.labels_arg, sArgs.mmap_flag)) {
	    cerr << "Could not open input labels Dat" << endl;
	    return 1;
	  }
	  
	  // set all NaN values to negatives
	  if( sArgs.nan2neg_given ){
	    cerr << "Set NaN labels dat as negatives" << endl;
	    
	    for(i = 0; i < Labels.GetGenes(); i++)
	      for(j = (i+1); j < Labels.GetGenes(); j++)
		if (CMeta::IsNaN(d = Labels.Get(i, j)))  
		  Labels.Set(i, j, -1);
	  }
	  
	  // Only keep edges that have one gene in this given list
	  if( sArgs.geneq_arg ){
	    Labels.FilterGenes( sArgs.geneq_arg, CDat::EFilterExEdge );	    
	    Labels.FilterGenes( sArgs.geneq_arg, CDat::EFilterEdge );
	  }
	  
	  // if given a target gene file
	  // Only keep eges that have only one gene in this targe gene list
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
	    
	    for(i = 0; i < Labels.GetGenes(); i++)
	      for(j = (i+1); j < Labels.GetGenes(); j++)
		if (!CMeta::IsNaN(d = Labels.Get(i, j))){
		  if(mapTgene[i] && mapTgene[j])
		    Labels.Set(i, j, CMeta::GetNaN());
		  else if(!mapTgene[i] && !mapTgene[j])
		    Labels.Set(i, j, CMeta::GetNaN());
		}
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
	
	SVMLight::SAMPLE* pTrainSample;
	SVMLight::SAMPLE* pAllSample;
	vector<SVMLight::SVMLabelPair*> pTrainVector[sArgs.cross_validation_arg];
	vector<SVMLight::SVMLabelPair*> pTestVector[sArgs.cross_validation_arg];
	//vector<SVMLight::Result> AllResults;
	//vector<SVMLight::Result> tmpAllResults;
	
	if (sArgs.model_given && sArgs.labels_given) { //learn once and write to file
	  //DEBUG implement
	  //pTrainSample = CSVMPERF::CreateSample(vecstrDatasets, vecLabels, Labels.GetGeneNames());
	  //SVM.Learn(*pTrainSample);
	  //SVM.WriteModel(sArgs.model_arg);
	} else if (sArgs.model_given && sArgs.output_given) { //read model and classify all
	  
	  // DEBUG need to figure out which/where to select all genes to predict
	  // Also, need to paralize the full predictions
	  
	  SVM.ReadModel(sArgs.model_arg);
	  //AllResults = SVM.Classify(PCL, vecAllLabels);
	  
	  // DEBUG output predictions
	  // sArgs.output_arg
	  
	} else if (sArgs.output_given && sArgs.labels_given) {
		//do learning and classifying with cross validation
	        if( sArgs.cross_validation_arg > 1){
		  cerr << "setting cross validation holds" << endl;
		  
		  mapTgene2fold.resize(mapTgene.size());
		  
		  // assign target genes to there cross validation fold
		  if(sArgs.tgene_given){
		    for(i = 0; i < mapTgene.size(); i++){
		      if(!mapTgene[i]){
			mapTgene2fold[i] = -1; 
			continue;
		      }
		      //cerr << "what's up?" << endl;
		      mapTgene2fold[i] = rand() % sArgs.cross_validation_arg;
		    }
		    
		    // cross-fold by target gene
		    for (i = 0; i < sArgs.cross_validation_arg; i++) {
		      cerr << "cross validation holds setup:" << i << endl;
		      
		      // keep track of positive gene counts
		      if(sArgs.balance_flag){
			cerr << "Set up balance: " << i << endl;
			for(j = 0; j < Labels.GetGenes(); j++)
			  tgeneCount[j] = 0;
			
			for(j = 0; j < vecLabels.size(); j++)
			  if(vecLabels[j]->Target > 0){
			    ++(tgeneCount[vecLabels[j]->iidx]);
			    ++(tgeneCount[vecLabels[j]->jidx]);
			  }
			
			if(sArgs.bfactor_given)
			  for(j = 0; j < vecLabels.size(); j++)
			    if(tgeneCount[vecLabels[j]->jidx] < 500)
			      tgeneCount[vecLabels[j]->jidx] = sArgs.bfactor_arg*tgeneCount[vecLabels[j]->jidx];
		      }
		      
		      for (j = 0; j < vecLabels.size(); j++) {
			//if( j % 1000 == 0)
			//cerr << "cross validation push labels:" << j << endl;
			
			// assume only one gene is a target gene in a edge
			if(mapTgene[vecLabels[j]->iidx]){
			  if(vecLabels[j]->Target < 0){
			    --(tgeneCount[vecLabels[j]->iidx]);
			  }
			  
			  if(mapTgene2fold[vecLabels[j]->iidx] == i)			    
			    pTestVector[i].push_back(vecLabels[j]);
			  else{
			    //cerr << tgeneCount[vecLabels[j]->iidx] << endl;
			    
			    if( sArgs.balance_flag && vecLabels[j]->Target < 0 && tgeneCount[vecLabels[j]->iidx] < 0){
			      continue;
			    }
			    
			    // only add if both genes are in context
			    if( sArgs.context_given  && ( !mapCgene[vecLabels[j]->iidx] || !mapCgene[vecLabels[j]->jidx]))
			      continue;
			    
			    pTrainVector[i].push_back(vecLabels[j]); 
			  }
			}
			else if(mapTgene[vecLabels[j]->jidx]){
			  if(vecLabels[j]->Target < 0)
			    --(tgeneCount[vecLabels[j]->jidx]);
			  
			  if(mapTgene2fold[vecLabels[j]->jidx] == i)
			    pTestVector[i].push_back(vecLabels[j]);
			  else{
			    //cerr << tgeneCount[vecLabels[j]->jidx] << endl;
			    
			    if( sArgs.balance_flag && vecLabels[j]->Target < 0 && tgeneCount[vecLabels[j]->jidx] < 0){
			      continue;
			    }
			    
			    // only add if both genes are in context
			    if( sArgs.context_given && ( !mapCgene[vecLabels[j]->iidx] || !mapCgene[vecLabels[j]->jidx]))
			      continue;
			    
			    pTrainVector[i].push_back(vecLabels[j]); 
			  }
			}
			else{
			  cerr << "Error: edge exist without a target gene" << endl; return 1;
			}
		      }
		      
		      cerr << "test,"<< i <<": " << pTestVector[i].size() << endl;
		      int numpos = 0;
		      for(j=0; j < pTrainVector[i].size(); j++)
			if(pTrainVector[i][j]->Target > 0)
			  ++numpos;
		      
		      if( numpos < 1 || (sArgs.mintrain_given && sArgs.mintrain_arg > numpos) ){						
			cerr << "Not enough positive examples from fold: " << i  << " file: " << sArgs.labels_arg << " numpos: " << numpos << endl;
			return 1;
		      }
		      
		      cerr << "train,"<< i <<","<<numpos<<": " << pTrainVector[i].size() << endl;
		      
		    }
		  }
		  else{ //randomly set eges into cross-fold
		    if( sArgs.context_given ){
		      cerr << "context not implemented yet for random edge holdout" << endl;
		      return 1;
		    }

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
			} else {
			  pTrainVector[i].push_back((vecLabels[j]));
			}
		      }
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
		
		//if (sArgs.all_flag) {
		  // DEBUG, probably don't allow all flag
		  // orig for classify all genes
		//}
		if (sArgs.params_given) { //reading paramters from file
			ifsm.close();
			ifsm.clear();
			ifsm.open(sArgs.params_arg);
			if (!ifsm.is_open()) {
				cerr << "Could not open: " << sArgs.params_arg << endl;
				return 1;
			}
			ParamStruct PStruct;
			string outFile(sArgs.output_arg);
			PStruct = ReadParamsFromFile(ifsm, outFile);

			size_t iParams;
			ofstream ofsm;
			SVMLight::SAMPLE * ppTrainSample[sArgs.cross_validation_arg];
			
			//build all the samples since they are being reused
			/*
			for (i = 0; i < sArgs.cross_validation_arg; i++)
			  ppTrainSample[i] = SVMLight::CSVMPERF::CreateSample(vecstrDatasets,
									      pTrainVector[i]);
			
			for (iParams = 0; iParams < PStruct.vecTradeoff.size(); iParams++) {
			  // DEBUG need to implement			  
			}
			*/
		} else { //run once
		  
		  // create sample for all labels
		  // DEBUG need to implement
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
		  }else{
		    SVMLight::CSVMPERF::CreateDoc(vecstrDatasets,
						  vecLabels,
						  Labels.GetGeneNames());
		  }
		  
		  for (i = 0; i < sArgs.cross_validation_arg; i++) {
		    cerr << "Cross validation: " << i << endl;
		    
		    // DEBUG need to implement
		    pTrainSample = SVMLight::CSVMPERF::CreateSample(pTrainVector[i]);
		    cerr << "Cross Validation Trial " << i << endl;		    
		    
		    SVM.Learn(*pTrainSample);
		    cerr << "Learned" << endl;
		    SVM.Classify(Results,
				 pTestVector[i]);
		    //cerr << "Classified " << tmpAllResults.size() << " examples"
		    //<< endl;
		    //AllResults.insert(AllResults.end(), tmpAllResults.begin(),
		    //				      tmpAllResults.end());
		    //tmpAllResults.resize(0);
		    ///if (sArgs.all_flag) {
		      // DEBUG, probably don't allow all flag
		      // orig for classify all genes				  
		    //}
		    
		    // MAIRA, Also, why did you not start at i==0???
		    
		    if(sArgs.savemodel_flag){
		      // i cross validation
		      // sArgs.output_arg		  
		      std::stringstream sstm;
		      
		      if(sArgs.context_given){
			std::string path(sArgs.context_arg);
			size_t pos = path.find_last_of("/");
			std::string cname;
			if(pos != std::string::npos)
			  cname.assign(path.begin() + pos + 1, path.end());
			else
			  cname = path;
			
			sstm << sArgs.output_arg << "." << sArgs.tradeoff_arg  << "." << i << "." << cname << ".svm";		      
		      }else
			sstm << sArgs.output_arg << "." << sArgs.tradeoff_arg  << "." << i << ".svm";
		      SVM.WriteModel((char*)(sstm.str().c_str()));
		    }
		    
		    // DEBUG
		    SVMLight::CSVMPERF::FreeSample_leave_Doc(*pTrainSample);
		    free(pTrainSample);
		  }
		  
		  //if (sArgs.all_flag) { //add the unlabeled results
		  // DEBUG, probably don't allow all flag
		  // orig for classify all genes
		  //}
		  
		  if(sArgs.prob_flag){
		    cerr << "Converting prediction values to estimated probablity" << endl;
		    float A, B;
		    sigmoid_train(Results, Labels, A, B);
		    sigmoid_predict(Results, A, B);
		  }
		  
		  Results.Save(sArgs.output_arg);
		  return 0;
		}
	} else {
		cerr << "More options are needed" << endl;
	}

}


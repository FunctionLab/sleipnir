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
#ifndef SEEKQUERY_H
#define SEEKQUERY_H

#include "stdafx.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_sort_float.h>
#include <gsl/gsl_sort_vector_float.h>
#include <gsl/gsl_permute_vector_float.h>

namespace Sleipnir {

enum PartitionMode{
	LEAVE_ONE_IN = 0,
	LEAVE_ONE_OUT = LEAVE_ONE_IN + 1,
	CUSTOM_PARTITION = LEAVE_ONE_OUT + 1
};

class CSeekQuery{
public:
	CSeekQuery(){
		crossValGenes = NULL;
		queryGenes.clear();
		iNumFold = 0;
		iFoldSize = 0;
	}
	~CSeekQuery(){
		if(crossValGenes!=NULL){
			delete[] crossValGenes;
		}
		queryGenes.clear();
		iNumFold = 0;
		iFoldSize = 0;
	}

	bool InitializeQuery(vector<char> query){
		size_t i;
		for(i=0; i<query.size(); i++){
			if(query[i]==1){
				queryGenes.push_back(query[i]);
			}
		}
		queryGenes.resize(queryGenes.size());
		return true;
	}

	vector<int>& GetQuery(){
		return queryGenes;
	}

	vector<int>& GetCVQuery(size_t i){
		return crossValGenes[i];
	}

	bool CreateCVPartitions(gsl_rng *rnd, enum PartitionMode p, size_t iFold=-1){
		//must have run initializequery beforehand
		if(p!=LEAVE_ONE_IN && p!=LEAVE_ONE_OUT && p!=CUSTOM_PARTITION){
			cerr << "Error, unknown partition mode" << endl;
			return false;
		}
		qSize = queryGenes.size();
		size_t fold_size = 0;
		if(iFold==-1){
			if(p==LEAVE_ONE_IN){
				iFold = qSize;
				fold_size = 1;
			}else if(p==LEAVE_ONE_OUT){
				iFold = qSize;
				fold_size = qSize-1;
			}else{
				cerr << "Error, must specify number of folds if CustomPartition mode" << endl;
				return false;
			}
		}else{
			if(p==LEAVE_ONE_IN || p==LEAVE_ONE_OUT){
				cerr << "Error, specified number of folds, so this must NOT be LEAVE_ONE_OUT or LEAVE_ONE_IN" << endl;
				return false;
			}
			fold_size = qSize / iFold;
			if(qSize % iFold > 0){
				fold_size++;
			}

		}
		iNumFold = iFold;
		iFoldSize = fold_size;
		crossValGenes = new vector<int>[iNumFold];

		size_t i, j, k;
		int *q_b = (int*)malloc(qSize*sizeof(int));
		for(i=0; i<qSize; i++){
			q_b[i] = queryGenes[i];
		}

		gsl_ran_shuffle(rnd, q_b, qSize, sizeof(int));

		if(p==LEAVE_ONE_IN || p==CUSTOM_PARTITION){
			k = 0;
			for(i=0; i<iFold; i++){
				for(j=0; j<iFoldSize; j++){
					if(k==qSize) continue;
					crossValGenes[i].push_back(q_b[k]);
					k++;
				}
				crossValGenes[i].resize(crossValGenes[i].size());
			}
		}else if(p==LEAVE_ONE_OUT){
			int current_index = -1;
			for(i=0; i<iFold; i++){
				for(j=0; j<iFoldSize; j++){
					current_index = (i+j) % qSize;
					crossValGenes[i].push_back(q_b[current_index]);
				}
				crossValGenes[i].resize(crossValGenes[i].size());
			}
		}

		free(q_b);
		return true;
	}

private:
	vector<int> queryGenes;
	vector<int> *crossValGenes;
	size_t iNumFold;
	size_t iFoldSize;
	size_t qSize;

};

}
#endif

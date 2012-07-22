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
#include "seekquery.h"
#include "seekreader.h"


namespace Sleipnir {

CSeekQuery::CSeekQuery(){
	crossValGenes = NULL;
	queryGenes.clear();
	queryGenePresence.clear();
	iNumFold = 0;
	iFoldSize = 0;
}
CSeekQuery::~CSeekQuery(){
	if(crossValGenes!=NULL){
		delete[] crossValGenes;
	}
	queryGenePresence.clear();
	queryGenes.clear();
	iNumFold = 0;
	iFoldSize = 0;
}

bool CSeekQuery::InitializeQuery(const vector<char> &query){
	ushort i;
	queryGenePresence.resize(query.size());

	for(i=0; i<query.size(); i++){
		if(query[i]==1) queryGenes.push_back(i);
		queryGenes[i] = query[i];
	}
	queryGenes.resize(queryGenes.size());
	return true;
}

bool CSeekQuery::InitializeQuery(const vector<ushort> &query,
	const ushort &iGenes){
	ushort i;
	queryGenePresence.resize(iGenes);
	fill(queryGenePresence.begin(), queryGenePresence.end(), (char) 0);
	for(i=0; i<query.size(); i++){
		queryGenes.push_back(query[i]);
		queryGenePresence[query[i]] = 1;
	}
	queryGenes.resize(queryGenes.size());
	return true;
}

ushort CSeekQuery::GetNumFold() const{
	return iNumFold;
}

vector<char>& CSeekQuery::GetQueryPresence(){
	return queryGenePresence;
}

vector<ushort>& CSeekQuery::GetQuery(){
	return queryGenes;
}

vector<ushort>& CSeekQuery::GetCVQuery(ushort &i){
	return crossValGenes[i];
}

bool CSeekQuery::CreateCVPartitions(const gsl_rng *rnd,
		const enum PartitionMode &p, const ushort iFold){
	//must have run initializequery beforehand
	if(p!=LEAVE_ONE_IN && p!=LEAVE_ONE_OUT && p!=CUSTOM_PARTITION){
		cerr << "Error, unknown partition mode" << endl;
		return false;
	}
	qSize = queryGenes.size();
	ushort fold_size = 0;
	ushort iFoldx = iFold;
	if(CSeekTools::IsNaN(iFold)){
		if(p==LEAVE_ONE_IN){
			iFoldx = qSize;
			fold_size = 1;
		}else if(p==LEAVE_ONE_OUT){
			iFoldx = qSize;
			fold_size = qSize-1;
		}else{
			cerr << "Error, must specify number of folds if \
					CustomPartition mode" << endl;
			return false;
		}
	}else{
		if(p==LEAVE_ONE_IN || p==LEAVE_ONE_OUT){
			cerr << "Error, specified number of folds, so this must NOT be \
					LEAVE_ONE_OUT or LEAVE_ONE_IN" << endl;
			return false;
		}
		fold_size = qSize / iFoldx;
		if(qSize % iFoldx > 0){
			fold_size++;
		}
	}
	iNumFold = iFoldx;
	iFoldSize = fold_size;
	crossValGenes = new vector<ushort>[iNumFold];
	//printf("Fold size %d %d\n", iNumFold, iFoldSize);

	ushort i, j, k;
	ushort *q_b = (ushort*)malloc(qSize*sizeof(ushort));
	for(i=0; i<qSize; i++){
		q_b[i] = queryGenes[i];
		//printf("%d ", q_b[i]);
	}
	//printf("\n");
	//getchar();
	gsl_ran_shuffle(rnd, q_b, qSize, sizeof(ushort));

	if(p==LEAVE_ONE_IN || p==CUSTOM_PARTITION){
		k = 0;
		for(i=0; i<iNumFold; i++){
			for(j=0; j<iFoldSize; j++){
				if(k==qSize) continue;
				crossValGenes[i].push_back(q_b[k]);
				k++;
			}
			crossValGenes[i].resize(crossValGenes[i].size());
		}
	}else if(p==LEAVE_ONE_OUT){
		ushort current_index = -1;
		for(i=0; i<iNumFold; i++){
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
}

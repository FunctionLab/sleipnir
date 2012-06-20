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
#include "seekmap.h"
#include "stdafx.h"
#include "datapair.h"
#include "seekwriter.h"
#include "seekreader.h"

namespace Sleipnir {

bool CSeekWriter::GetGeneAverage(CDataPair &Dat, vector<string> &vecstrGenes,
		vector<float> &vecResult){
	/* assume datapair is already opened */
	size_t i, j;
	vector<size_t> veciGenes;
	veciGenes.clear();
	veciGenes.resize(vecstrGenes.size());
	for( i = 0; i < vecstrGenes.size( ); ++i )
		veciGenes[ i ] = Dat.GetGene( vecstrGenes[i] );

	CSeekTools::InitVector(vecResult, vecstrGenes.size(), CMeta::GetNaN());
	for(i=0; i<vecstrGenes.size(); i++){
		size_t s = veciGenes[i];
		if(s==-1) continue;
		float *v = Dat.GetFullRow(s);
		float sum = 0;
		int num = 0;
		for(j=0; j<vecstrGenes.size(); j++){
			size_t t = veciGenes[j];
			if(t==-1) continue;
			if(CMeta::IsNaN(v[t])) continue;
			sum+=v[t];
			num++;
		}
		vecResult[i] = sum / (float) num;
		free(v);
	}
	return true;
}

bool CSeekWriter::GetGenePresence(CDataPair &Dat, vector<string> &vecstrGenes,
		vector<char> &vecResult){
	/* assume datapair is already opened */
	size_t i, j;
	vector<size_t> veciGenes;
	veciGenes.clear();
	veciGenes.resize(vecstrGenes.size());
	for( i = 0; i < vecstrGenes.size( ); ++i )
		veciGenes[ i ] = Dat.GetGene( vecstrGenes[i] );

	CSeekTools::InitVector(vecResult, vecstrGenes.size(), (char) 0);

	for(i=0; i<vecstrGenes.size(); i++){
		if(veciGenes[i]==-1) continue;
		vecResult[i]=1;
	}
	return true;
}

}

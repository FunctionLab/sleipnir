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
#include "seekwriter.h"
#include "seekreader.h"

namespace Sleipnir {

bool CSeekWriter::GetGeneAverage(CDataPair &Dat,
	const vector<string> &vecstrGenes,
	vector<float> &vecResult, bool logit){

	/* assume datapair is already opened */
	ushort i, j;
	vector<ushort> veciGenes;
	veciGenes.clear();
	veciGenes.resize(vecstrGenes.size());
	for( i = 0; i < vecstrGenes.size( ); ++i )
		veciGenes[ i ] = Dat.GetGene( vecstrGenes[i] );

	CSeekTools::InitVector(vecResult, vecstrGenes.size(), CMeta::GetNaN());
	for(i=0; i<vecstrGenes.size(); i++){
		ushort s = veciGenes[i];
		if(CSeekTools::IsNaN(s)) continue;
		float *v = Dat.GetFullRow(s);
		float sum = 0;
		ushort num = 0;
		for(j=0; j<vecstrGenes.size(); j++){
			ushort t = veciGenes[j];
			if(CSeekTools::IsNaN(t)) continue;
			if(CMeta::IsNaN(v[t])) continue;
			if(logit){
				sum+=log(v[t]) - log((float) (1.0-v[t]));
			}else{
				sum+=v[t];
			}
			num++;
		}
		vecResult[i] = sum / (float) num;
		free(v);
	}
	return true;
}

bool CSeekWriter::GetGenePresence(CDataPair &Dat,
	const vector<string> &vecstrGenes,
	vector<char> &vecResult){
	/* assume datapair is already opened */
	ushort i, j;
	vector<ushort> veciGenes;
	veciGenes.clear();
	veciGenes.resize(vecstrGenes.size());
	for( i = 0; i < vecstrGenes.size( ); ++i )
		veciGenes[ i ] = Dat.GetGene( vecstrGenes[i] );

	CSeekTools::InitVector(vecResult, vecstrGenes.size(), (char) 0);

	for(i=0; i<vecstrGenes.size(); i++){
		if(CSeekTools::IsNaN(veciGenes[i])) continue;
		vecResult[i]=1;
	}
	return true;
}

bool CSeekWriter::GetDatasetSinfo(CDataPair &Dat,
	float &mean, float &stdev){
	ushort i, j;
	mean = CMeta::GetNaN();
	stdev = CMeta::GetNaN();

	ushort iGenes = Dat.GetGenes();

	unsigned int num = 0;
	float sum = 0;

	for(i=0; i<iGenes; i++){
		ushort s = i;
		float *v = Dat.GetFullRow(s);
		for(j=0; j<iGenes; j++){
			if(CMeta::IsNaN(v[j])) continue;
			sum+=v[j];
			num++;
		}
	}

	if(num==0) return true;

	mean = sum / (float) num;
	float diff = 0;
	for(i=0; i<iGenes; i++){
		ushort s = i;
		float *v = Dat.GetFullRow(s);
		for(j=0; j<iGenes; j++){
			if(CMeta::IsNaN(v[j])) continue;
			diff += (v[j] - mean) * (v[j] - mean);
		}
	}
	diff /= (float) num;
	stdev = sqrt(diff);

	return true;
}

}

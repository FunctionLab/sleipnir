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

#include "seekbasic.h"

namespace Sleipnir {

enum PartitionMode{
	LEAVE_ONE_IN = 0,
	LEAVE_ONE_OUT = LEAVE_ONE_IN + 1,
	CUSTOM_PARTITION = LEAVE_ONE_OUT + 1
};

class CSeekQuery{
public:
	CSeekQuery();
	~CSeekQuery();

	bool InitializeQuery(const vector<char>&);
	bool InitializeQuery(const vector<ushort>&, const ushort &);

	ushort GetNumFold() const;
	const vector<ushort>& GetQuery() const;
	const vector<char>& GetQueryPresence() const;
	const vector<ushort>& GetCVQuery(ushort&) const;
	bool CreateCVPartitions(const gsl_rng*, \
		const enum PartitionMode &, const ushort=-1);

private:
	vector<ushort> queryGenes;
	vector<char> queryGenePresence;

	vector<ushort> *crossValGenes;
	ushort iNumFold;
	ushort iFoldSize;
	ushort qSize;

};

}
#endif

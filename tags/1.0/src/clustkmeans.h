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
#ifndef CLUSTKMEANS_H
#define CLUSTKMEANS_H

#include "clustkmeansi.h"

namespace Sleipnir {

class IMeasure;

/*!
 * \brief
 * Utility class containing static k-means clustering methods.
 */
class CClustKMeans : protected CClustKMeansImpl {
public:
	static bool Cluster( const CDataMatrix& MatData, const IMeasure* pMeasure, size_t iK,
		std::vector<uint16_t>& vecsClusters, const CDataMatrix* pMatWeights = NULL );
};

}

#endif // CLUSTQTC_H

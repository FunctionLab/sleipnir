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
#ifndef STATISTICSI_H
#define STATISTICSI_H

#include "fullmatrix.h"
#include "mathb.h"
#include "meta.h"

namespace Sleipnir {

    class CStatisticsImpl : protected CMath {
    protected:
        static double Chi2CDF(double, double, double, double);

        static double IncompleteGamma(double, double);

        static double Normal01CDF(double);

        static double GammaLog(double);

        static double IncompleteBeta(double, double, double);

        static double IncompleteBetaCF(double, double, double);

        static double ModifiedBesselI(size_t, double);

        static bool MatrixLUSubstitute(CDataMatrix &, const std::vector <size_t> &, std::vector<float> &);

        static double EpsilonDouble() {
            double dRet;

            for (dRet = 1; 1 < (1 + dRet); dRet /= 2);

            return (2 * dRet);
        }

        static bool MatrixMultiply(const std::vector<float> &vecdLeft, const CDataMatrix &MatRight,
                                   std::vector<float> &vecdOut) {
            size_t i, j;

            if (vecdLeft.size() != MatRight.GetRows())
                return false;

            vecdOut.resize(MatRight.GetColumns());
            for (i = 0; i < vecdOut.size(); ++i) {
                vecdOut[i] = 0;
                for (j = 0; j < vecdLeft.size(); ++j)
                    vecdOut[i] += vecdLeft[j] * MatRight.Get(j, i);
            }

            return true;
        }

        static double MatrixMultiply(const std::vector<float> &vecdLeft, const std::vector<float> &vecdRight) {
            size_t i;
            double dRet;

            if (vecdLeft.size() != vecdRight.size())
                return CMeta::GetNaN();

            for (dRet = 0, i = 0; i < vecdLeft.size(); ++i)
                dRet += vecdLeft[i] * vecdRight[i];

            return dRet;
        }

        template<class tType>
        static bool SumsSkip(tType Value) {

            return false;
        }

        static bool SumsSkip(float dValue) {

            return CMeta::IsNaN(dValue);
        }

        static bool SumsSkip(double dValue) {

            return CMeta::IsNaN(dValue);
        }

        template<class tType>
        static void Sums(tType Begin, tType End, double *pdSum, double *pdSumSq, size_t *piN) {
            tType Cur;

            if (pdSum)
                *pdSum = 0;
            if (pdSumSq)
                *pdSumSq = 0;
            if (piN)
                *piN = 0;
            for (Cur = Begin; Cur != End; ++Cur) {
                if (SumsSkip(*Cur))
                    continue;
                if (piN)
                    *piN += 1;
                if (pdSum)
                    *pdSum += *Cur;
                if (pdSumSq)
                    *pdSumSq += *Cur * *Cur;
            }
        }
    };

}

#endif // STATISTICSI_H

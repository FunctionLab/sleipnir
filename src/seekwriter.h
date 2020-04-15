/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a part of SEEK (Search-based exploration of expression compendium)
* which is authored and maintained by: Qian Zhu (qzhu@princeton.edu)
*
* If you use this file, please cite the following publication:
* Qian Zhu, Aaron K Wong, Arjun Krishnan, Miriam R Aure, Alicja Tadych, 
* Ran Zhang, David C Corney, Casey S Greene, Lars A Bongo, 
* Vessela N Kristensen, Moses Charikar, Kai Li & Olga G Troyanskaya
* "Targeted exploration and analysis of large cross-platform human 
* transcriptomic compendia" Nat Methods (2015)
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library for development, or use any other Sleipnir executable
* tools, please also cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#ifndef SEEKWRITER_H
#define SEEKWRITER_H

#include "seekbasic.h"
#include "seekmap.h"
#include "seekevaluate.h"
#include "sparsematrix.h"
#include "datapair.h"
#include "seekreader.h"
#include "strassen.h"

namespace Sleipnir {

    class CSeekWriter {
    public:

        //should be either unsigned short or utype
        template<class tType>
        static bool ReadSeekSparseMatrixHeader(const char *fileName,
                                               CSeekIntIntMap &m) {
            FILE *f = fopen(fileName, "rb");
            if (f == NULL) {
                cerr << "File not found" << endl;
                return false;
            }
            size_t j;
            tType val, numPresent;
            int ret;

            //m need to be initialized to size vecstrGenes.size() first!
            ret = fread((char *) (&numPresent), 1, sizeof(numPresent), f);
            for (j = 0; j < numPresent; j++) {
                ret = fread((char *) (&val), 1, sizeof(val), f); //val stores the present gene
                m.Add((utype) val);
            }
            fclose(f);
            return true;
        }

        //compatibility
        template<class tType>
        static bool ReadSeekSparseMatrix(const char *fileName,
                                         CSparseFlatMatrix<float> &mat, CSeekIntIntMap &m, const int maxRank,
                                         const float rbp_p, const vector <string> &vecstrGenes) {

            FILE *f = fopen(fileName, "rb");
            if (f == NULL) {
                cerr << "File not found" << endl;
                return false;
            }

            size_t i, j;
            tType numGenes, numPresent, val;
            int ret;

            mat.Initialize(vecstrGenes.size());
            ret = fread((char *) (&numPresent), 1, sizeof(numPresent), f);
            for (j = 0; j < numPresent; j++) {
                ret = fread((char *) (&val), 1, sizeof(val), f); //val = gene ID
                m.Add((utype) val);
                mat.InitializeRow(val, maxRank * 2); //initial capacity
            }
            ret = fread((char *) (&numGenes), 1, sizeof(numGenes), f);

            vector<float> rbp_score;
            rbp_score.resize(maxRank);
            for (i = 0; i < maxRank; i++)
                rbp_score[i] = (1.0 - rbp_p) * pow(rbp_p, (float) i);

            for (i = 0; i < numGenes; i++) {
                tType id, id2;  //gene ID
                unsigned short numEntries, val; //rank
                ret = fread((char *) (&id), 1, sizeof(id), f);
                ret = fread((char *) (&numEntries), 1, sizeof(numEntries), f);
                for (j = 0; j < numEntries; j++) {
                    ret = fread((char *) (&id2), 1, sizeof(id2), f);
                    ret = fread((char *) (&val), 1, sizeof(val), f);
                    tType first = id;
                    tType second = id2;
                    //mat is a full matrix
                    mat.Add(first, second, rbp_score[val]);
                    mat.Add(second, first, rbp_score[val]);
                }
            }
            fclose(f);

            mat.Organize();
            size_t ii, jj;
            const vector <utype> &allRGenes = m.GetAllReverse();
            fprintf(stderr, "Begin calculating row sum\n");

            vector<float> vecSum;
            CSeekTools::InitVector(vecSum, vecstrGenes.size(), (float) 0);
            for (ii = 0; ii < m.GetNumSet(); ii++) {
                i = (size_t) allRGenes[ii];
                const vector <CPair<float>> &row = mat.GetRow(i);
                for (jj = 0; jj < row.size(); jj++)
                    vecSum[i] += row[jj].v;
            }

            vector<float> vecSqrtSum;
            CSeekTools::InitVector(vecSqrtSum, vecstrGenes.size(), (float) 0);

            for (ii = 0; ii < m.GetNumSet(); ii++) {
                i = (size_t) allRGenes[ii];
                if (vecSum[i] == 0) continue;
                vecSqrtSum[i] = sqrtf(vecSum[i]);
            }

            fprintf(stderr, "Begin normalization using row sum\n");
            float rv;
#pragma omp parallel for \
        shared(m, mat, allRGenes, vecSqrtSum) \
        private(ii, i, j, rv) schedule(dynamic)
            for (ii = 0; ii < m.GetNumSet(); ii++) {
                i = (size_t) allRGenes[ii];
                vector<CPair<float> >::iterator row_it;
                for (row_it = mat.RowBegin(i); row_it != mat.RowEnd(i); row_it++) {
                    j = (size_t) row_it->i;
                    rv = row_it->v;
                    if (vecSqrtSum[i] == 0 || vecSqrtSum[j] == 0) continue;
                    row_it->v = rv / vecSqrtSum[i] / vecSqrtSum[j];
                }
            }
            return true;
        }

        //compatibility
        template<class tType>
        static bool RemoveDominant(CSparseFlatMatrix<float> &mat,
                                   CSeekIntIntMap &m, const vector <string> &vecstrGenes) {

            size_t i, j;
            vector <vector<float>> tmp_mat, tmp2;
            tmp_mat.resize(vecstrGenes.size());
            tmp2.resize(vecstrGenes.size());
            for (i = 0; i < vecstrGenes.size(); i++) {
                tmp_mat[i].resize(vecstrGenes.size());
                tmp2[i].resize(vecstrGenes.size());
                for (j = 0; j < vecstrGenes.size(); j++) {
                    tmp_mat[i][j] = 0;
                    tmp2[i][j] = 0;
                }
            }

            size_t ii, jj;
            const vector <utype> &allRGenes = m.GetAllReverse();
            for (ii = 0; ii < m.GetNumSet(); ii++) {
                i = (size_t) allRGenes[ii];
                vector<CPair<float> >::iterator row_it;
                for (row_it = mat.RowBegin(i); row_it != mat.RowEnd(i); row_it++) {
                    j = (size_t) row_it->i;
                    tmp_mat[i][j] = row_it->v;
                }
            }

            int TOP = 100;
            fprintf(stderr, "Started removing dominant...\n");
            for (ii = 0; ii < m.GetNumSet(); ii++) {
                i = (size_t) allRGenes[ii];
                vector <CPair<float>> vp;
                vp.resize(m.GetNumSet());
                for (jj = 0; jj < m.GetNumSet(); jj++) {
                    j = (size_t) allRGenes[jj];
                    vp[jj].i = (utype) j;
                    vp[jj].v = tmp_mat[i][j];
                }
                nth_element(vp.begin(), vp.begin() + TOP, vp.end(), CDescendingValue<float>());
                sort(vp.begin(), vp.begin() + TOP, CDescendingValue<float>());

                //top 100
                size_t k, l;
                size_t max_ind = 0;
                float max_val = -1;
                for (k = 0; k < TOP; k++) {
                    size_t this_g = vp[k].i;
                    float v = 0;
                    for (l = 0; l < TOP; l++) {
                        size_t other_g = vp[l].i;
                        if (this_g == other_g) continue;
                        v += tmp_mat[this_g][other_g];
                    }
                    if (v > max_val) {
                        max_val = v;
                        max_ind = k;
                    }
                }
                for (jj = 0; jj < m.GetNumSet(); jj++) {
                    j = (size_t) allRGenes[jj];
                    tmp2[i][j] = tmp_mat[i][j] - tmp_mat[j][max_ind];
                    if (tmp2[i][j] < 0)
                        tmp2[i][j] = 0;
                }
            }

            fprintf(stderr, "Started re-normalizing matrix...\n");

            for (ii = 0; ii < m.GetNumSet(); ii++) {
                i = (size_t) allRGenes[ii];
                for (jj = ii + 1; jj < m.GetNumSet(); jj++) {
                    j = (size_t) allRGenes[jj];
                    float m = max(tmp2[i][j], tmp2[j][i]);
                    tmp2[i][j] = m;
                    tmp2[j][i] = m;
                }
            }

            vector<float> vecSum;
            CSeekTools::InitVector(vecSum, vecstrGenes.size(), (float) 0);
            for (ii = 0; ii < m.GetNumSet(); ii++) {
                i = (size_t) allRGenes[ii];
                for (jj = 0; jj < m.GetNumSet(); jj++) {
                    j = (size_t) allRGenes[jj];
                    vecSum[i] += tmp2[i][j];
                }
            }

            vector<float> vecSqrtSum;
            CSeekTools::InitVector(vecSqrtSum, vecstrGenes.size(), (float) 0);
            for (ii = 0; ii < m.GetNumSet(); ii++) {
                i = (size_t) allRGenes[ii];
                if (vecSum[i] == 0) continue;
                vecSqrtSum[i] = sqrtf(vecSum[i]);
            }

            for (ii = 0; ii < m.GetNumSet(); ii++) {
                i = (size_t) allRGenes[ii];
                vector<CPair<float> >::iterator row_it;
                for (row_it = mat.RowBegin(i); row_it != mat.RowEnd(i); row_it++) {
                    j = (size_t) row_it->i;
                    if (vecSqrtSum[i] == 0 || vecSqrtSum[j] == 0) continue;
                    row_it->v = tmp2[i][j] / vecSqrtSum[i] / vecSqrtSum[j];
                    //fprintf(stderr, "%d %d %.3e\n", i, j, row_it->v);
                }
            }
            return true;
        }

        //compatibility
        template<class tType>
        static bool WriteSparseMatrix(CDataPair &Dat,
                                      vector <map<tType, unsigned short>> &umat,
                                      const vector <string> &vecstrGenes, const char *fileName) {

            FILE *f = fopen(fileName, "wb");
            if (f == NULL) {
                fprintf(stderr, "File not found %s\n", fileName);
                return false;
            }

            size_t i, j;
            vector <tType> veciGenes;
            veciGenes.resize(vecstrGenes.size());
            for (i = 0; i < vecstrGenes.size(); i++)
                veciGenes[i] = (tType) Dat.GetGeneIndex(vecstrGenes[i]);

            CSeekIntIntMap mm(vecstrGenes.size());
            for (i = 0; i < vecstrGenes.size(); i++)
                if (veciGenes[i] != (tType) -1)
                    mm.Add((utype) i);

            tType numPresent = (tType) mm.GetNumSet();
            //1 tType
            fwrite((char *) (&numPresent), 1, sizeof(numPresent), f);
            const vector <utype> &allR = mm.GetAllReverse();
            //numPresent tType
            for (i = 0; i < numPresent; i++) {
                tType pr = (tType) allR[i];
                fwrite((char *) (&pr), 1, sizeof(pr), f);
            }

            tType numGenes = 0;
            for (i = 0; i < vecstrGenes.size(); i++) {
                if (umat[i].size() == 0) continue;
                numGenes++;
            }
            //1 tType
            fwrite((char *) (&numGenes), 1, sizeof(numGenes), f);

            for (i = 0; i < vecstrGenes.size(); i++) {
                unsigned short numEntries = umat[i].size(); //should be 1000
                if (numEntries == 0)
                    continue;
                //1 tType
                tType gi = (tType) i;
                fwrite((char *) (&gi), 1, sizeof(gi), f);
                //1 unsigned short
                fwrite((char *) (&numEntries), 1, sizeof(numEntries), f);
                typename map<tType, unsigned short>::iterator it;
                for (it = umat[i].begin(); it != umat[i].end(); it++) {
                    tType first = it->first;
                    unsigned short second = it->second;
                    //1 tType
                    fwrite((char *) (&first), 1, sizeof(first), f);
                    //1 unsigned short
                    fwrite((char *) (&second), 1, sizeof(second), f);
                }
            }
            fclose(f);
            return true;
        }

        //compatiblity
        template<class tType>
        static bool GetSparseRankMatrix(CDataPair &Dat,
                                        vector <map<tType, unsigned short>> &umat, int maxRank,
                                        const vector <string> &vecstrGenes) {

            size_t i, j;
            vector <tType> veciGenes;
            veciGenes.resize(vecstrGenes.size());
            for (i = 0; i < vecstrGenes.size(); ++i)
                veciGenes[i] = (tType) Dat.GetGeneIndex(vecstrGenes[i]);
            umat.resize(vecstrGenes.size());
            for (i = 0; i < vecstrGenes.size(); i++)
                umat[i] = map<tType, unsigned short>();

            fprintf(stderr, "Start reading DAB...\n");
            tType s, t;
            for (i = 0; i < vecstrGenes.size(); i++) {
                if ((s = veciGenes[i]) == (tType) -1) continue;
                if (i % 1000 == 0) fprintf(stderr, "Start reading gene %d...\n", i);

                float *v = Dat.GetFullRow(s);
                vector <AResultFloat> vv;
                vv.resize(vecstrGenes.size());

                for (j = 0; j < vecstrGenes.size(); j++) {
                    vv[j].i = (utype) j;
                    vv[j].f = -9999;
                    if ((t = veciGenes[j]) == (tType) -1) continue;
                    if (CMeta::IsNaN(v[t])) continue;
                    vv[j].f = v[t];
                }
                nth_element(vv.begin(), vv.begin() + maxRank, vv.end());
                sort(vv.begin(), vv.begin() + maxRank);

                for (j = 0; j < vecstrGenes.size(); j++) {
                    if (j < maxRank) {
                        tType first = (tType) i;
                        tType second = (tType) vv[j].i;
                        if ((tType) i >= (tType) vv[j].i) {
                            first = (tType) vv[j].i;
                            second = (tType) i;
                        }
                        typename map<tType, unsigned short>::iterator it;
                        if ((it = umat[first].find(second)) == umat[first].end())
                            umat[first][second] = (unsigned short) j;
                        else
                            umat[first][second] = std::min(it->second, (unsigned short) j);
                    }
                }
                delete[] v;
            }
            fprintf(stderr, "Finished reading DAB\n");
            return true;
        }

        //To be used after NormalizeDAB
        template<class tType>
        static bool ConvertToSparseMatrix(CDataPair &Dat,
                                          vector <map<tType, unsigned short>> &umat,
                                          const vector <string> &vecstrGenes, const float cutoff_val) {

            size_t i, j;
            vector <tType> veciGenes;
            veciGenes.resize(vecstrGenes.size());
            for (i = 0; i < vecstrGenes.size(); ++i)
                veciGenes[i] = (tType) Dat.GetGeneIndex(vecstrGenes[i]);
            umat.resize(vecstrGenes.size());
            for (i = 0; i < vecstrGenes.size(); i++)
                umat[i] = map<tType, unsigned short>();

            tType s, t;
            for (i = 0; i < vecstrGenes.size(); i++) {
                if ((s = veciGenes[i]) == (tType) -1) continue;
                if (i % 1000 == 0) fprintf(stderr, "Start reading gene %d...\n", i);

                for (j = i + 1; j < vecstrGenes.size(); j++) {
                    if ((t = veciGenes[j]) == (tType) -1) continue;
                    float r = Dat.Get(s, t);
                    if (CMeta::IsNaN(r)) continue;
                    if (r > cutoff_val)
                        umat[i][j] = (unsigned short) (r * 100.0);
                }
            }
            fprintf(stderr, "Finished reading DAB\n");
            return true;
        }

        //to be used for sparse matrix created from cutting-off z-scores
        template<class tType>
        static bool ReadSeekSparseMatrix(const char *fileName,
                                         CSparseFlatMatrix<float> &mat, CSeekIntIntMap &m,
                                         const vector <string> &vecstrGenes, const int initialCapacity,
                                         const float exponent) {

            if (exponent < 1.0) {
                fprintf(stderr, "Exponent must be >=1.0\n");
                return false;
            }

            FILE *f = fopen(fileName, "rb");
            if (f == NULL) {
                cerr << "File not found" << endl;
                return false;
            }

            size_t i, j;
            tType numGenes, numPresent, val;
            int ret;

            mat.Initialize(vecstrGenes.size());
            ret = fread((char *) (&numPresent), 1, sizeof(numPresent), f);
            for (j = 0; j < numPresent; j++) {
                ret = fread((char *) (&val), 1, sizeof(val), f); //val = gene ID
                m.Add((utype) val);
                mat.InitializeRow(val, initialCapacity); //initial capacity
            }
            ret = fread((char *) (&numGenes), 1, sizeof(numGenes), f);

            for (i = 0; i < numGenes; i++) {
                tType id, id2;  //gene ID
                unsigned short numEntries, val; //z-scores * 100.0
                ret = fread((char *) (&id), 1, sizeof(id), f);
                ret = fread((char *) (&numEntries), 1, sizeof(numEntries), f);
                for (j = 0; j < numEntries; j++) {
                    ret = fread((char *) (&id2), 1, sizeof(id2), f);
                    ret = fread((char *) (&val), 1, sizeof(val), f);
                    tType first = id;
                    tType second = id2;
                    float fval = (float) val / 100.0;
                    if (exponent > 1.0)
                        fval = pow(fval, exponent);
                    mat.Add(first, second, fval);
                    mat.Add(second, first, fval);
                }
            }
            fclose(f);

            mat.Organize();
            return true;
        }

        //===============================================================
        //not currently used
        static bool ReadSparseMatrix(const char *fileName,
                                     vector <map<utype, float>> &mat,
                                     CSeekIntIntMap &m, const int maxRank, const float rbp_p,
                                     const vector <string> &vecstrGenes);

        //not currently used
        static bool ProductNorm(const vector <map<utype, float>> &mat1,
                                const vector <map<utype, float>> &mat2, const CSeekIntIntMap &m1,
                                const CSeekIntIntMap &m2, vector <map<utype, float>> &re);

        //================================================================
        static bool SumSparseMatrix(CSparseFlatMatrix<float> &mat1,
                                    CSparseFlatHalfMatrix<float> &res, const CSeekIntIntMap &mi, const float w);

        static bool SumSparseMatrix(CSparseFlatMatrix<float> &mat1,
                                    CHalfMatrix<float> &res, const CSeekIntIntMap &mi, const float w);

        static bool NormalizeDAB(CDataPair &Dat, const vector <string> &vecstrGenes,
                //bool cutoff, float cutoff_val,
                                 bool expTransform, bool divideNorm, bool subtractNorm);

        static bool GetGeneAverage(CDataPair &Dat,
                                   const vector <string> &vecstrGenes,
                                   vector<float> &vecResult, bool logit = false, float top_percent = 1.0);

        static bool GetGenePresence(CDataPair &Dat,
                                    const vector <string> &vecstrGenes,
                                    vector<char> &vecResult);

        static bool GetDatasetSinfo(CDataPair &Dat, float &mean,
                                    float &stdev);

        static void TopologicalOverlap(CDataPair &Dat,
                                       const vector <string> &vecstrGenes);
    };

}
#endif

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
#ifndef FULLMATRIX_H
#define FULLMATRIX_H

#undef int64_t

#include <stdint.h>
#include "meta.h"
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

namespace Sleipnir {

/*!
 * \brief
 * An asymmetric two-dimensional matrix.
 * 
 * \param tType
 * Type of element contained by the matrix.
 * 
 * \remarks
 * Essentially nothing except a lack of bounds checking for efficiency's sake separates this from every
 * other 2D matrix implementation known to man.
 * 
 * \see
 * CCompactMatrix | CHalfMatrix
 */
    template<class tType>
    class CFullMatrix {
    public:
        CFullMatrix() : m_iR(0), m_iC(0), m_aaData(NULL) {}

        virtual ~CFullMatrix() {

            Reset();
        }

        /*!
         * \brief
         * Empties the matrix and deallocates all associated memory.
         */
        void Reset() {
            size_t i;

            if (m_aaData && m_fMemory) {
                for (i = 0; i < GetRows(); ++i)
                    delete[] m_aaData[i];
                delete[] m_aaData;
            }

            m_iR = m_iC = 0;
            m_aaData = NULL;
            m_fMemory = false;
        }

        /*!
         * \brief
         * Return the two-dimensional array backing the matrix.
         *
         * \returns
         * Two-dimensional array backing the matrix.
         */
        tType **const Get() const {

            return m_aaData;
        }

        /*!
         * \brief
         * Return a single row of the matrix.
         *
         * \param iY
         * Matrix row.
         *
         * \returns
         * Requested matrix row.
         *
         * \remarks
         * For efficiency, no bounds checking is performed.  The given row must be smaller than the size.  The
         * returned array will contain a number of elements equal to the number of columns.
         *
         * \see
         * Set
         */
        tType *Get(size_t iY) const {

            return m_aaData[iY];
        }

        /*!
         * \brief
         * Returns the value at the requested matrix position.
         *
         * \param iY
         * Matrix row.
         *
         * \param iX
         * Matrix column.
         *
         * \returns
         * Value at the requested matrix position.
         *
         * \remarks
         * For efficiency, no bounds checking is performed.  The given row and column must be smaller than the
         * size.
         *
         * \see
         * Set
         */
        tType &Get(size_t iY, size_t iX) const {

            return m_aaData[iY][iX];
        }

        /*!
         * \brief
         * Set the value at the requested matrix position.
         *
         * \param iY
         * Matrix row.
         *
         * \param iX
         * Matrix column.
         *
         * \param Value
         * Value to store.
         *
         * \remarks
         * For efficiency, no bounds checking is performed.  The given row and column must be smaller than the
         * size.
         *
         * \see
         * Get
         */
        void Set(size_t iY, size_t iX, const tType &Value) {

            m_aaData[iY][iX] = Value;
        }

        /*!
         * \brief
         * Set a single row of the matrix.
         *
         * \param iY
         * Matrix row.
         *
         * \param aValues
         * Data to be copied into the requested row.
         *
         * \remarks
         * For efficiency, no bounds checking is performed.  The given row must be smaller than the size,
         * and the given data array must be at least as long as the number of columns.
         *
         * \see
         * Get
         */
        void Set(size_t iY, const tType *aValues) {

            memcpy(m_aaData[iY], aValues, GetColumns() * sizeof(*aValues));
        }

        /*!
         * \brief
         * Create a new matrix using a reference to the given matrix.
         *
         * \param Mat
         * Matrix to duplicate in the created matrix.
         *
         * \remarks
         * A reference is held to the given matrix; it is not copied, so it should not be destroyed before
         * the newly initialized matrix.
         */
        void Initialize(const CFullMatrix &Mat) {

            Initialize(Mat.GetRows(), Mat.GetColumns(), Mat.m_aaData);
        }

        /*!
         * \brief
         * Create a new matrix of the requested size and, optionally, referencing the given data.
         *
         * \param iR
         * Matrix rows.
         *
         * \param iC
         * Matrix columns.
         *
         * \param aaData
         * If non-null, the memory that will back the newly created matrix.
         *
         * \remarks
         * A reference is held to the given memory; it is not copied, so it should not be destroyed before
         * the newly initialized matrix.
         */
        virtual void Initialize(size_t iR, size_t iC, tType **aaData = NULL) {
            size_t i;

            Reset();
            m_iR = iR;
            m_iC = iC;

            if (aaData != NULL) {
                m_aaData = aaData;
                m_fMemory = false; //reference a given memory
            } else {
                m_fMemory = true; //create a new memory
                m_aaData = new tType *[GetRows()];
                for (i = 0; i < GetRows(); ++i)
                    m_aaData[i] = new tType[GetColumns()];
            }
        }

        /*!
         * \brief
         * Return the number of rows in the matrix.
         *
         * \returns
         * Number of rows in the matrix.
         *
         * \see
         * GetColumns
         */
        size_t GetRows() const {

            return m_iR;
        }

        /*!
         * \brief
         * Return the number of columns in the matrix.
         *
         * \returns
         * Number of columns in the matrix.
         *
         * \see
         * GetRows
         */
        size_t GetColumns() const {

            return m_iC;
        }

        /*!
         * \brief
         * Loads a matrix from the given stream.
         *
         * \param istm
         * Stream from which matrix is loaded.
         *
         * \param fBinary
         * If true, matrix and stream are in binary format.
         *
         * \returns
         * True if matrix was loaded successfully.
         *
         * \remarks
         * A binary matrix is saved in row-major order and should have been produced using Save.  A text matrix
         * is a simple tab-delimited file from which matrix elements are read using >>.
         */
        bool Open(std::istream &istm, bool fBinary) {

            return (fBinary ? OpenBinary(istm) : OpenText(istm));
        }

        /*!
         * \brief
         * Loads a matrix from the given filename, or standard input if null.
         *
         * \param szFile
         * File from which binary matrix is opened; if null, use text from standard input.
         *
         * \returns
         * True if open was successful, false otherwise.
         *
         * \see
         * Save
         */
        bool Open(const char *szFile) {

            if (szFile) {
                std::ifstream ifsm;

                ifsm.open(szFile, ios_base::binary);
                return Open(ifsm, true);
            }

            return Open(std::cin, false);
        }

        /*!
         * \brief
         * Creates a copy of the given matrix's data.
         *
         * \param Mat
         * Matrix from which data is copied.
         *
         * \returns
         * True if matrix was copied successfully.
         *
         * \remarks
         * Matrix will be resized if dimensions do not match exactly.
         *
         * \see
         * Initialize
         */
        bool Copy(const CFullMatrix &Mat) {
            size_t i;

            if ((GetRows() != Mat.GetRows()) || (GetColumns() != Mat.GetColumns()))
                Initialize(Mat.GetRows(), Mat.GetColumns());

            for (i = 0; i < GetRows(); ++i)
                memcpy(Get(i), Mat.Get(i), GetColumns() * sizeof(*Get(i)));

            return true;
        }

        /*!
         * \brief
         * Saves a matrix to the given stream in either binary or tab-delimited text format.
         *
         * \param ostm
         * Stream to which matrix is saved.
         *
         * \param fBinary
         * If true, matrix is saved in binary format; otherwise, matrix is saved as tab-delimited text.
         *
         * \returns
         * True if matrix was saved successfully.
         *
         * \remarks
         * A binary matrix is saved in row-major order; a text matrix is a simple tab-delimited file from which
         * matrix elements are read using >>.
         *
         * \see
         * Open
         */
        bool Save(std::ostream &ostm, bool fBinary) const {

            return (fBinary ? SaveBinary(ostm) : SaveText(ostm));
        }

        /*!
         * \brief
         * Saves a matrix to the given file in binary form; if null, save text to standard output.
         *
         * \param szFile
         * File in which binary matrix is saved; if null, text is saved to standard output.
         *
         * \returns
         * True if matrix was saved successfully.
         *
         * \see
         * Open
         */
        bool Save(const char *szFile) const {

            if (szFile) {
                std::ofstream ofsm;

                ofsm.open(szFile, ios_base::binary);
                return Save(ofsm, true);
            }

            return Save(cout, false);
        }

        /*!
         * \brief
         * Saves a matrix to the given stream in either binary or tab-delimited text format.
         *
         * \param szFile
         * File to which matrix is saved.
         *
         * \param fBinary
         * If true, matrix is saved in binary format; otherwise, matrix is saved as tab-delimited text.
         *
         * \returns
         * True if matrix was saved successfully.
         *
         * \remarks
         * A binary matrix is saved in row-major order; a text matrix is a simple tab-delimited file from which
         * matrix elements are read using >>.
         *
         * \see
         * Open
         */
        bool Save(const char *szFile, bool fBinary) const {
            std::ofstream ofsm;

            ofsm.open(szFile, (fBinary ? ios_base::binary : ios_base::out) | ios_base::out);
            return Save(ofsm, fBinary);
        }

        /*!
         * \brief
         * Sets all entries of the matrix to 0 without changing its size.
         */
        void Clear() {
            size_t i;

            if (!m_aaData)
                return;

            for (i = 0; i < GetRows(); ++i)
                memset(m_aaData[i], 0, sizeof(*m_aaData[i]) * GetColumns());
        }

        /*!
         * \brief
         * Increases the number of rows in the matrix by the requested amount.
         *
         * \param iR
         * Number of rows to append to the matrix.
         *
         * \param fClear
         * If true, set all values in the new rows to zero.
         *
         * \returns
         * True if the rows were appended successfully.
         *
         * \remarks
         * On success, GetRows will be increased by the requested amount.  The news rows are not initialized
         * to any particular value.
         */
        bool AddRows(size_t iR, bool fClear = false) {
            tType **m_aaNew;
            size_t i;

            if (!m_fMemory)
                return false;

            m_aaNew = new tType *[GetRows() + iR];
            memcpy(m_aaNew, m_aaData, GetRows() * sizeof(*m_aaData));
            delete[] m_aaData;
            m_aaData = m_aaNew;
            for (i = 0; i < iR; ++i) {
                m_aaNew[GetRows() + i] = new tType[GetColumns()];
                if (fClear)
                    memset(m_aaNew[GetRows() + i], 0, GetColumns() * sizeof(**m_aaNew));
            }
            m_iR += iR;

            return true;
        }

        /*!
         * \brief
         * Performs matrix multiplication by applying the given vector on the right of the current matrix.
         *
         * \param vecRight
         * Vector by which current matrix is multiplied on the right.
         *
         * \returns
         * True on success, false if the given vector is of an inappropriate size.
         *
         * \remarks
         * vecRight must have the same number of entries as the current matrix has columns.
         *
         * \see
         * Multiply
         */
        bool Multiply(std::vector <tType> &vecRight) {
            size_t i, j;

            if (vecRight.size() != GetColumns())
                return false;

            for (i = 0; i < GetColumns(); ++i)
                for (j = 0; j < GetRows(); ++j)
                    Set(j, i, vecRight[i] * Get(j, i));

            return true;
        }

        /*!
         * \brief
         * Performs matrix multiplication by applying the given matrix on the right of the current matrix.
         *
         * \param MatRight
         * Matrix by which current matrix is multiplied on the right.
         *
         * \param fTranspose
         * If true, transpose MatRight prior to multiplication.
         *
         * \returns
         * True on success, false if MatRight is of inappropriate dimensions.
         *
         * \remarks
         * If fTranspose is false, MatRight must have the same number of rows as the current matrix has
         * columns; if it is true, MatRight and the current matrix must have the same number of columns.
         *
         * \see
         * Multiply
         */
        bool Multiply(const CFullMatrix &MatRight, bool fTranspose = false) {
            CFullMatrix MatTmp;
            size_t i, j, k, iRR, iRC;

            iRR = fTranspose ? MatRight.GetColumns() : MatRight.GetRows();
            if (GetColumns() != iRR)
                return false;

            MatTmp.Initialize(GetRows(), iRC = fTranspose ? MatRight.GetRows() : MatRight.GetColumns());
            MatTmp.Clear();
            for (i = 0; i < iRC; ++i)
                for (j = 0; j < GetRows(); ++j)
                    for (k = 0; k < iRR; ++k)
                        MatTmp.Get(j, i) += Get(j, k) * (fTranspose ? MatRight.Get(i, k) :
                                                         MatRight.Get(k, i));

            return Copy(MatTmp);
        }

        /*!
         * \brief
         * Performs singular value decomposition of the current matrix.
         *
         * \param MatU
         * Output U matrix of basis vectors.
         *
         * \param MatV
         * Output V matrix of analysis vectors.
         *
         * \param vecS
         * Output diagonal of Sigma matrix of singular values.
         *
         * \returns
         * True on success, false if the matrix cannot be decomposed.
         *
         * Outputs a singular value decomposition of the current matrix M such that M = U * Sigma * V'.
         * For a matrix of expression values, V behaves as eigengenes, U as eigenconditions, and S
         * quantifies the variability (or weight) attributable to each eigengene.
         *
         * \remarks
         * Matrix must not contain any missing values.  Implementation courtesy of AJ Sedgewick.
         */
        bool SVD(CFullMatrix<tType> &MatU, CFullMatrix<tType> &MatV, std::vector <tType> &vecS) const {
            size_t i, j, k, iNU, iCT, iRT;
            vector <tType> vecE, vecWork;
            CFullMatrix<tType> MatA;
            bool fU, fV;
            ESVD eCase;
            tType F;
            int iMarker, iP, iPP;

            iNU = min(GetRows(), GetColumns());
            vecS.resize(min(GetRows() + 1, GetColumns()));
            MatU.Initialize(GetRows(), iNU);
            MatU.Clear();
            MatV.Initialize(GetColumns(), GetColumns());

            vecE.resize(GetColumns());
            vecWork.resize(GetRows());
            MatA.Copy(*this);

            // Reduce A to bidiagonal form, storing the diagonal elements
            // in s and the super-diagonal elements in e.
            fU = fV = true;
            iCT = min(GetRows() - 1, GetColumns());
            iRT = max((size_t) 0, min(GetColumns() - 2, GetRows()));
            for (k = 0; k < max(iCT, iRT); ++k) {
                if (k < iCT) {
                    // Compute the transformation for the k-th column and
                    // place the k-th diagonal in s[k].
                    // Compute 2-norm of k-th column without under/overflow.
                    vecS[k] = 0;
                    for (i = k; i < GetRows(); ++i)
                        vecS[k] = Hypotenuse(vecS[k], MatA.Get(i, k));
                    if (vecS[k]) {
                        if (MatA.Get(k, k) < 0)
                            vecS[k] *= -1;
                        for (i = k; i < GetRows(); ++i)
                            MatA.Get(i, k) /= vecS[k];
                        MatA.Get(k, k) += 1;
                    }
                    vecS[k] *= -1;
                }

                for (j = (k + 1); j < GetColumns(); ++j) {
                    if ((k < iCT) && vecS[k]) {
                        tType T;

                        // Apply the transformation.
                        T = 0;
                        for (i = k; i < GetRows(); ++i)
                            T += MatA.Get(i, k) * MatA.Get(i, j);
                        T /= -MatA.Get(k, k);
                        for (i = k; i < GetRows(); ++i)
                            MatA.Get(i, j) += T * MatA.Get(i, k);
                    }

                    // Place the k-th row of A into e for the
                    // subsequent calculation of the row transformation.
                    vecE[j] = MatA.Get(k, j);
                }

                if (fU && (k < iCT))
                    // Place the transformation in U for subsequent back
                    // multiplication.
                    for (i = k; i < GetRows(); ++i)
                        MatU.Set(i, k, MatA.Get(i, k));
                if (k < iRT) {
                    // Compute the k-th row transformation and place the
                    // k-th super-diagonal in e[k].
                    // Compute 2-norm without under/overflow.
                    vecE[k] = 0;
                    for (i = (k + 1); i < GetColumns(); ++i)
                        vecE[k] = Hypotenuse(vecE[k], vecE[i]);
                    if (vecE[k]) {
                        if (vecE[k + 1] < 0)
                            vecE[k] *= -1;
                        for (i = (k + 1); i < GetColumns(); ++i)
                            vecE[i] /= vecE[k];
                        vecE[k + 1] += 1;
                    }
                    vecE[k] *= -1;

                    if (((k + 1) < GetRows()) && vecE[k]) {
                        // Apply the transformation.
                        for (i = (k + 1); i < GetRows(); ++i)
                            vecWork[i] = 0;
                        for (j = (k + 1); j < GetColumns(); ++j)
                            for (i = (k + 1); i < GetRows(); ++i)
                                vecWork[i] += vecE[j] * MatA.Get(i, j);
                        for (j = (k + 1); j < GetColumns(); ++j) {
                            tType T;

                            T = -vecE[j] / vecE[k + 1];
                            for (i = (k + 1); i < GetRows(); ++i)
                                MatA.Get(i, j) += T * vecWork[i];
                        }
                    }

                    if (fV)
                        // Place the transformation in V for subsequent
                        // back multiplication.
                        for (i = (k + 1); i < GetColumns(); ++i)
                            MatV.Set(i, k, vecE[i]);
                }
            }

            // Set up the final bidiagonal matrix or order p.
            iP = (int) min(GetColumns(), GetRows() + 1);
            if (iCT < GetColumns())
                vecS[iCT] = MatA.Get(iCT, iCT);
            if (GetRows() < (size_t) iP)
                vecS[iP - 1] = 0;
            if ((iRT + 1) < (size_t) iP)
                vecE[iRT] = MatA.Get(iRT, iP - 1);
            vecE[iP - 1] = 0;

            // If required, generate U.
            if (fU) {
                for (j = iCT; j < iNU; ++j) {
                    for (i = 0; i < GetRows(); ++i)
                        MatU.Set(i, j, 0);
                    MatU.Set(j, j, 1);
                }
                for (size_t iTmp = 0; iTmp < iCT; ++iTmp) {
                    k = iCT - 1 - iTmp;
                    if (vecS[k]) {
                        for (j = (k + 1); j < iNU; ++j) {
                            tType T;

                            T = 0;
                            for (i = k; i < GetRows(); ++i)
                                T += MatU.Get(i, k) * MatU.Get(i, j);
                            T /= -MatU.Get(k, k);
                            for (i = k; i < GetRows(); ++i)
                                MatU.Get(i, j) += T * MatU.Get(i, k);
                        }
                        for (i = k; i < GetRows(); ++i)
                            MatU.Get(i, k) *= -1;
                        MatU.Get(k, k) += 1;
                        for (i = 0; (i + 1) < k; ++i)
                            MatU.Set(i, k, 0);
                    } else {
                        for (i = 0; i < GetRows(); ++i)
                            MatU.Set(i, k, 0);
                        MatU.Set(k, k, 1);
                    }
                }
            }

            // If required, generate V.
            if (fV) {
                for (size_t iTmp = 0; iTmp < GetColumns(); ++iTmp) {
                    k = GetColumns() - 1 - iTmp;
                    if ((k < iRT) && vecE[k])
                        for (j = (k + 1); j < iNU; ++j) {
                            tType T;

                            T = 0;
                            for (i = (k + 1); i < GetColumns(); ++i)
                                T += MatV.Get(i, k) * MatV.Get(i, j);
                            T /= -MatV.Get(k + 1, k);
                            for (i = (k + 1); i < GetColumns(); ++i)
                                MatV.Get(i, j) += T * MatV.Get(i, k);
                        }
                    for (i = 0; i < GetColumns(); ++i)
                        MatV.Set(i, k, 0);
                    MatV.Set(k, k, 1);
                }
            }

            // Main iteration loop for the singular values.
            iPP = iP - 1;
            while (iP) {
                // Here is where a test for too many iterations would go.
                // This section of the program inspects for
                // negligible elements in the s and e arrays.  On
                // completion the variables kase and k are set as follows.
                for (iMarker = (iP - 2); iMarker >= 0; --iMarker)
                    if (fabs(vecE[iMarker]) <= (std::numeric_limits<tType>::epsilon() *
                                                (fabs(vecS[iMarker]) + fabs(vecS[iMarker + 1])))) {
                        vecE[iMarker] = 0;
                        break;
                    }
                if (iMarker == (iP - 2))
                    eCase = ESVDConverged;
                else {
                    int iKS;

                    for (iKS = (iP - 1); iKS > iMarker; --iKS) {
                        tType T;

                        T = ((iKS != iP) ? fabs(vecE[iKS]) : 0) +
                            ((iKS != (iMarker + 1)) ? fabs(vecE[iKS - 1]) : 0);
                        if (fabs(vecS[iKS]) <= (std::numeric_limits<tType>::epsilon() * T)) {
                            vecS[iKS] = 0;
                            break;
                        }
                    }
                    if (iKS == iMarker)
                        eCase = ESVDSmallE;
                    else if (iKS == (iP - 1))
                        eCase = ESVDSmallSE;
                    else {
                        eCase = ESVDSmallS;
                        iMarker = iKS;
                    }
                }
                iMarker++;

                switch (eCase) {
                    // Deflate negligible s(p).
                    case ESVDSmallSE:
                        F = vecE[iP - 2];
                        vecE[iP - 2] = 0;
                        for (int iTmp = (iP - 2); iTmp >= iMarker; --iTmp) {
                            tType T, CS, SN;

                            T = Hypotenuse(vecS[iTmp], F);
                            CS = vecS[iTmp] / T;
                            SN = F / T;
                            vecS[iTmp] = T;
                            if (iTmp != iMarker) {
                                F = -SN * vecE[iTmp - 1];
                                vecE[iTmp - 1] = CS * vecE[iTmp - 1];
                            }
                            if (fV)
                                for (i = 0; i < GetColumns(); ++i) {
                                    T = (CS * MatV.Get(i, iTmp)) + (SN * MatV.Get(i, iP - 1));
                                    MatV.Set(i, iP - 1, (-SN * MatV.Get(i, iTmp)) +
                                                        (CS * MatV.Get(i, iP - 1)));
                                    MatV.Set(i, iTmp, T);
                                }
                        }
                        break;

                        // Split at negligible s(k).
                    case ESVDSmallS:
                        F = vecE[iMarker - 1];
                        vecE[iMarker - 1] = 0;
                        for (j = iMarker; (int) j < iP; ++j) {
                            tType T, CS, SN;

                            T = Hypotenuse(vecS[j], F);
                            CS = vecS[j] / T;
                            SN = F / T;
                            vecS[j] = T;
                            F = -SN * vecE[j];
                            vecE[j] *= CS;
                            if (fU)
                                for (i = 0; i < GetRows(); ++i) {
                                    T = (CS * MatU.Get(i, j)) + (SN * MatU.Get(i, iMarker - 1));
                                    MatU.Set(i, iMarker - 1, (-SN * MatU.Get(i, j)) +
                                                             (CS * MatU.Get(i, iMarker - 1)));
                                    MatU.Set(i, j, T);
                                }
                        }
                        break;

                        // Perform one qr step.
                    case ESVDSmallE:
                        // Calculate the shift.
                        tType Scale, SP, SPM1, EPM1, SK, EK, B, C, Shift, G;

                        Scale = Max(fabs(vecS[iP - 1]), fabs(vecS[iP - 2]), fabs(vecE[iP - 2]),
                                    fabs(vecS[iMarker]), fabs(vecE[iMarker]));
                        SP = vecS[iP - 1] / Scale;
                        SPM1 = vecS[iP - 2] / Scale;
                        EPM1 = vecE[iP - 2] / Scale;
                        SK = vecS[iMarker] / Scale;
                        EK = vecE[iMarker] / Scale;
                        B = (((SPM1 + SP) * (SPM1 - SP)) + (EPM1 * EPM1)) / 2;
                        C = SP * SP * EPM1 * EPM1;
                        Shift = 0;
                        if (B || C) {
                            Shift = sqrt((B * B) + C);
                            if (B < 0)
                                Shift *= -1;
                            Shift = C / (B + Shift);
                        }
                        F = ((SK + SP) * (SK - SP)) + Shift;
                        G = SK * EK;

                        // Chase zeros.
                        for (j = iMarker; (int) (j + 1) < iP; ++j) {
                            tType T, CS, SN;

                            T = Hypotenuse(F, G);
                            CS = F / T;
                            SN = G / T;
                            if (j != iMarker)
                                vecE[j - 1] = T;
                            F = (CS * vecS[j]) + (SN * vecE[j]);
                            vecE[j] = (CS * vecE[j]) - (SN * vecS[j]);
                            G = SN * vecS[j + 1];
                            vecS[j + 1] = CS * vecS[j + 1];
                            if (fV)
                                for (i = 0; i < GetColumns(); ++i) {
                                    T = (CS * MatV.Get(i, j)) + (SN * MatV.Get(i, j + 1));
                                    MatV.Set(i, j + 1, (-SN * MatV.Get(i, j)) + (CS * MatV.Get(i, j + 1)));
                                    MatV.Set(i, j, T);
                                }
                            T = Hypotenuse(F, G);
                            CS = F / T;
                            SN = G / T;
                            vecS[j] = T;
                            F = (CS * vecE[j]) + (SN * vecS[j + 1]);
                            vecS[j + 1] = (-SN * vecE[j]) + (CS * vecS[j + 1]);
                            G = SN * vecE[j + 1];
                            vecE[j + 1] = CS * vecE[j + 1];
                            if (fU && ((j + 1) < GetRows()))
                                for (i = 0; i < GetRows(); ++i) {
                                    T = (CS * MatU.Get(i, j)) + (SN * MatU.Get(i, j + 1));
                                    MatU.Set(i, j + 1, (-SN * MatU.Get(i, j)) + (CS * MatU.Get(i, j + 1)));
                                    MatU.Set(i, j, T);
                                }
                        }
                        vecE[iP - 2] = F;
                        break;

                        // Convergence.
                    case ESVDConverged:
                        // Make the singular values positive.
                        if (vecS[iMarker] < 0) {
                            vecS[iMarker] = (vecS[iMarker] < 0) ? -vecS[iMarker] : 0;
                            if (fV)
                                for (i = 0; (int) i <= iPP; ++i)
                                    MatV.Get(i, iMarker) *= -1;
                        }

                        // Order the singular values.
                        while (iMarker < iPP) {
                            tType T;

                            if (vecS[iMarker] >= vecS[iMarker + 1])
                                break;
                            std::swap(vecS[iMarker], vecS[iMarker + 1]);
                            if (fV && ((iMarker + 1) < (int) GetColumns()))
                                for (i = 0; i < GetColumns(); ++i) {
                                    T = MatV.Get(i, iMarker + 1);
                                    MatV.Set(i, iMarker + 1, MatV.Get(i, iMarker));
                                    MatV.Set(i, k, T);
                                }
                            if (fU && ((iMarker + 1) < (int) GetRows()))
                                for (i = 0; i < GetRows(); ++i) {
                                    T = MatU.Get(i, iMarker + 1);
                                    MatU.Set(i, iMarker + 1, MatU.Get(i, iMarker));
                                    MatU.Set(i, iMarker, T);
                                }
                            iMarker++;
                        }
                        iP--;
                        break;
                }
            }

            return true;
        }

    private:
        static const size_t c_iBuffer = 131072;

        enum ESVD {
            ESVDSmallSE,    // s(p) and e[k-1] are negligible and k<p
            ESVDSmallS,        // s(k) is negligible and k<p
            ESVDSmallE,        // e[k-1] is negligible, k<p, s(k), ..., s(p) are not negligible (qr step)
            ESVDConverged    // e(p-1) is negligible (convergence).
        };

        static tType Max(const tType &A, const tType &B, const tType &C, const tType &D, const tType &E) {

            return max(A, max(B, max(C, max(D, E))));
        }

        static tType Hypotenuse(const tType &A, const tType &B) {

            return sqrt((A * A) + (B * B));
        }

        bool OpenBinary(std::istream &istm) {
            size_t i;
            uint32_t iSize;

            Reset();
            if (!istm.good() || istm.eof())
                return false;
            istm.read((char *) &iSize, sizeof(iSize));
            m_iR = iSize;
            istm.read((char *) &iSize, sizeof(iSize));
            m_iC = iSize;
            m_fMemory = true;
            m_aaData = new tType *[GetRows()];
            for (i = 0; i < GetRows(); ++i) {
                m_aaData[i] = new tType[GetColumns()];
                istm.read((char *) m_aaData[i], (std::streamsize)(sizeof(*m_aaData[i]) * GetColumns()));
            }

            return true;
        }

        bool OpenText(std::istream &istm) {
            char szLine[c_iBuffer];
            std::vector <tType> vecRow;
            std::vector < tType * > vecpRows;
            std::stringstream sstm;
            tType Cur;
            tType *aCur;
            size_t i;

            Reset();
            m_fMemory = true;
            while (istm.peek() != EOF) {
                istm.getline(szLine, c_iBuffer - 1);
                if (!szLine[0])
                    break;
                sstm.clear();
                sstm.str(szLine);
                vecRow.clear();
                while (sstm.peek() != EOF) {
                    sstm >> Cur;
                    vecRow.push_back(Cur);
                }
                if (!GetColumns())
                    m_iC = vecRow.size();
                else if (GetColumns() != vecRow.size())
                    return false;

                aCur = new tType[GetColumns()];
                for (i = 0; i < GetColumns(); ++i)
                    aCur[i] = vecRow[i];
                vecpRows.push_back(aCur);
            }

            m_aaData = new tType *[m_iR = vecpRows.size()];
            for (i = 0; i < GetRows(); ++i)
                m_aaData[i] = vecpRows[i];

            return true;
        }

        bool SaveBinary(std::ostream &ostm) const {
            size_t i;
            uint32_t iSize;

            iSize = (uint32_t) GetRows();
            ostm.write((const char *) &iSize, sizeof(iSize));
            iSize = (uint32_t) GetColumns();
            ostm.write((const char *) &iSize, sizeof(iSize));
            for (i = 0; i < GetRows(); ++i)
                ostm.write((const char *) m_aaData[i], (std::streamsize)(sizeof(*m_aaData[i]) *
                                                                         GetColumns()));

            return true;
        }

        bool SaveText(std::ostream &ostm) const {
            size_t i, j;

            if (!(GetRows() && GetColumns()))
                return false;

            for (i = 0; i < GetRows(); ++i) {
                ostm << m_aaData[i][0];
                for (j = 1; j < GetColumns(); ++j)
                    ostm << '\t' << m_aaData[i][j];
                ostm << std::endl;
            }

            return true;
        }

        size_t m_iR;
        size_t m_iC;
        bool m_fMemory;
        tType **m_aaData;
    };

/*!
 * \brief
 * A full matrix of four-byte floats, used for almost all Sleipnir continuously valued asymmetric matrices.
 */
    typedef CFullMatrix<float> CDataMatrix;

}

#endif // FULLMATRIX_H

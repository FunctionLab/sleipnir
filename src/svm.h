#ifndef SVM_H
#define SVM_H

#ifndef NO_SVM_LIGHT

#include "svmi.h"

namespace Sleipnir {

class CSVM : CSVMImpl {
public:
	/*!
	 * \brief
	 * Type of kernel used by the SVM.
	 */
	enum EKernel {
		/*!
		 * \brief
		 * Linear kernel.
		 */
		EKernelLinear		= 0,
		/*!
		 * \brief
		 * Polynomial kernel.
		 */
		EKernelPolynomial	= EKernelLinear + 1,
		/*!
		 * \brief
		 * Radial basis function kernel.
		 */
		EKernelRBF			= EKernelPolynomial + 1
	};

	bool OpenAlphas( std::istream& istm );
	bool Open( std::istream& istm );
	bool Save( std::ostream& ostm ) const;
	bool Learn( const CPCL& PCL, const CGenes& GenesPositive );
	bool Learn( const CPCL& PCL, const CGenes& GenesPositive, const CGenes& GenesNegative );
	bool Evaluate( const CPCL& PCL, std::vector<float>& vecdResults ) const;

	bool Learn( const char* szData ) {
		SData	sData;

		sData.m_eType = SData::EFile;
		sData.m_uData.m_szFile = szData;

		return CSVMImpl::Learn( sData ); }

	bool Learn( const CPCLSet& PCLs, const CDataPair& Answers ) {
		SData	sData;

		sData.m_eType = SData::EPCLs;
		sData.m_uData.m_pPCLs = &PCLs;
		sData.m_uAnswers.m_pAnswers = &Answers;

		return CSVMImpl::Learn( sData ); }

	bool Learn( const IDataset* pData, const CDataPair& Answers ) {
		SData	sData;

		sData.m_eType = SData::EData;
		sData.m_uData.m_pData = pData;
		sData.m_uAnswers.m_pAnswers = &Answers;

		return CSVMImpl::Learn( sData ); }

	bool Evaluate( const char* szFile, CDat& DatOut ) const {
		SData	sData;

		sData.m_eType = SData::EFile;
		sData.m_uData.m_szFile = szFile;

		return CSVMImpl::Evaluate( sData, NULL, DatOut ); }

	bool Evaluate( const CPCLSet& PCLs, CDat& DatOut ) const {
		SData	sData;

		sData.m_eType = SData::EPCLs;
		sData.m_uData.m_pPCLs = &PCLs;

		return CSVMImpl::Evaluate( sData, NULL, DatOut ); }

	bool Evaluate( const IDataset* pData, CDat& DatOut ) const {
		SData	sData;

		sData.m_eType = SData::EData;
		sData.m_uData.m_pData = pData;

		return CSVMImpl::Evaluate( sData, NULL, DatOut ); }

	bool Evaluate( const CPCLSet& PCLs, const CGenes& GenesIn, CDat& DatOut ) const {
		SData	sData;

		sData.m_eType = SData::EPCLs;
		sData.m_uData.m_pPCLs = &PCLs;

		return CSVMImpl::Evaluate( sData, &GenesIn, DatOut ); }

	bool Evaluate( const IDataset* pData, const CGenes& GenesIn, CDat& DatOut ) const {
		SData	sData;

		sData.m_eType = SData::EData;
		sData.m_uData.m_pData = pData;

		return CSVMImpl::Evaluate( sData, &GenesIn, DatOut ); }

	void SetIterations( size_t iIterations ) {

		m_sLearn.maxiter = iIterations; }

	void SetCache( size_t iMegabytes ) {

		m_sLearn.kernel_cache_size = iMegabytes; }

	void SetTradeoff( float dTradeoff ) {

		m_sLearn.svm_c = dTradeoff; }

	void SetGamma( float dGamma ) {

		m_sKernel.rbf_gamma = dGamma; }

	void SetDegree( size_t iDegree ) {

		m_sKernel.poly_degree = iDegree; }

	void SetKernel( EKernel eKernel ) {

		m_sKernel.kernel_type = eKernel; }

	void SetVerbosity( size_t iVerbosity ) {

		verbosity = iVerbosity; }
};

}

#endif // NO_SVM_LIGHT

#endif // SVM_H

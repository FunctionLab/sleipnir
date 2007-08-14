#include "stdafx.h"
#include "svm.h"
#include "pclset.h"
#include "dataset.h"
#include "meta.h"
#include "genome.h"

extern "C" {
KERNEL_CACHE* kernel_cache_init( long, long );
void kernel_cache_cleanup( KERNEL_CACHE* );
void svm_learn_classification( DOC**, double*, long, long, LEARN_PARM*, KERNEL_PARM*,
	KERNEL_CACHE*, MODEL*, double* );
void svm_learn_regression( DOC**, double*, long, long, LEARN_PARM*, KERNEL_PARM*,
	KERNEL_CACHE**, MODEL* );
void svm_learn_ranking( DOC**, double*, long, long, LEARN_PARM*, KERNEL_PARM*,
	KERNEL_CACHE**, MODEL* );
void svm_learn_optimization( DOC**, double*, long, long, LEARN_PARM*, KERNEL_PARM*,
	KERNEL_CACHE*, MODEL*, double* );
}

namespace libBioUtils {

SWORD	CSVMImpl::s_asWords[ CSVMImpl::c_iWords ];

CSVMImpl::SLearn::SLearn( ) {

	predfile[ 0 ] = 0;
	alphafile[ 0 ] = 0;
	biased_hyperplane = 1;
	sharedslack = 0;
	remove_inconsistent = 0;
	skip_final_opt_check = 0;
	svm_maxqpsize = 10;
	svm_newvarsinqp = 0;
	svm_iter_to_shrink = 100;
	maxiter = 100000;
	kernel_cache_size = 40;
	svm_c = 0;
	eps = 0.1;
	transduction_posratio = -1.0;
	svm_costratio = 1;
	svm_costratio_unlab = 1;
	svm_unlabbound = 1e-5;
	epsilon_crit = 0.001;
	epsilon_a = 1e-15;
	compute_loo = 0;
	rho = 1;
	xa_depth = 0;
	type = CLASSIFICATION; }

CSVMImpl::SKernel::SKernel( ) {

	kernel_type = 0;
	poly_degree = 3;
	rbf_gamma = 1;
	coef_lin = 1;
	coef_const = 1;
	custom[ 0 ] = 0; }

CSVMImpl::CSVMImpl( ) : m_apDocs(NULL), m_iDocs(0), m_adAlphas(NULL), m_iAlphas(0),
	m_pModel(NULL), m_adLabels(NULL) {

	verbosity = 2; }

CSVMImpl::~CSVMImpl( ) {

	Reset( true, true, true ); }

void CSVMImpl::Reset( bool fData, bool fModel, bool fAlphas ) {
	size_t	i;

	if( fModel && m_pModel ) {
		free_model( m_pModel, 0 );
		m_pModel = NULL; }
	if( fAlphas && m_adAlphas ) {
		free( m_adAlphas );
		m_adAlphas = NULL; }
	if( fData ) {
		if( m_apDocs ) {
			for( i = 0; i < m_iDocs; ++i )
				free_example( m_apDocs[ i ], 1 );
			delete[] m_apDocs;
			m_apDocs = NULL; }
		if( m_adLabels ) {
			delete[] m_adLabels;
			m_adLabels = NULL; } } }

size_t CSVMImpl::GetWords( const SData& sData ) const {
	size_t	i, iRet;

	switch( sData.m_eType ) {
		case SData::EPCLs:
			for( iRet = i = 0; i < sData.m_uData.m_pPCLs->GetPCLs( ); ++i )
				iRet += sData.m_uData.m_pPCLs->Get( i ).GetExperiments( );
			return ( iRet * 2 );

		case SData::EData:
			return sData.m_uData.m_pData->GetExperiments( );

		case SData::EFile:
			return m_iWords;

		case SData::EPCL:
			return sData.m_uData.m_pPCL->GetExperiments( ); }

	return -1; }

DOC* CSVMImpl::CreateDoc( const SData& sData, size_t iOne, size_t iTwo, size_t iDoc ) const {
	SWORD*	asWords;
	size_t	i, j, iWord, iWords;
	float	d;
	DOC*	pRet;

	iWords = GetWords( sData );
	asWords = ( iWords >= c_iWords ) ? new SWORD[ iWords + 1 ] : s_asWords;
	for( i = 0; i < iWords; ++i )
		asWords[ i ].wnum = i + 1;
	asWords[ i ].wnum = 0;
	if( sData.m_eType == SData::EPCLs ) {
		const CPCLSet&	PCLs	= *sData.m_uData.m_pPCLs;

		for( iWord = i = 0; i < PCLs.GetPCLs( ); ++i ) {
			for( j = 0; j < PCLs.Get( i ).GetExperiments( ); ++j ) {
				if( CMeta::IsNaN( d = PCLs.Get( i, iOne, j ) ) )
					d = 0;
				assert( ( iWord + j ) < iWords );
				asWords[ iWord + j ].weight = d;
				if( CMeta::IsNaN( d = PCLs.Get( i, iTwo, j ) ) )
					d = 0;
				assert( ( iWord + ( iWords / 2 ) + j ) < iWords );
				asWords[ iWord + ( iWords / 2 ) + j ].weight = d; }
			iWord += PCLs.Get( i ).GetExperiments( ); } }
	else {
		const IDataset*	pData	= sData.m_uData.m_pData;

		for( i = 0; i < pData->GetExperiments( ); ++i ) {
			if( CMeta::IsNaN( d = pData->GetContinuous( iOne, iTwo, i ) ) )
				d = 0;
			asWords[ i ].weight = d; } }

	pRet = create_example( iDoc, 0, 0, 1, create_svector( asWords, "", 1 ) );
	if( asWords != s_asWords )
		delete[] asWords;
	return pRet; }

DOC* CSVMImpl::CreateDoc( const SData& sData, size_t iGene ) const {
	SWORD*	asWords;
	size_t	i, iWords;
	DOC*	pRet;

	iWords = GetWords( sData );
	asWords = ( iWords >= c_iWords ) ? new SWORD[ iWords + 1 ] : s_asWords;
	for( i = 0; i < iWords; ++i )
		asWords[ i ].wnum = i + 1;
	asWords[ i ].wnum = 0;

	for( i = 0; i < iWords; ++i )
		asWords[ i ].weight = sData.m_uData.m_pPCL->Get( iGene, i );
	pRet = create_example( i, 0, 0, 1, create_svector( asWords, "", 1 ) );

	if( asWords != s_asWords )
		delete[] asWords;
	return pRet; }

bool CSVM::OpenAlphas( istream& istm ) {
	static const size_t	c_iBuf	= 1024;
	char			szBuf[ c_iBuf ];
	vector<float>	vecdAlphas;
	size_t			i;

	Reset( false, false, true );
	while( istm.peek( ) != EOF ) {
		istm.getline( szBuf, c_iBuf - 1 );
		vecdAlphas.push_back( (float)atof( szBuf ) ); }
	m_adAlphas = new double[ m_iAlphas = vecdAlphas.size( ) ];
	for( i = 0; i < m_iAlphas; ++i )
		m_adAlphas[ i ] = vecdAlphas[ i ];

	return true; }

bool CSVMImpl::Initialize( const SData& sData ) {
	size_t			i, j, iOne, iTwo, iDoc;
	vector<size_t>	veciGenes;
	float			d;

	Reset( true, false, false );
	if( sData.m_eType == SData::EFile ) {
		read_documents( (char*)sData.m_uData.m_szFile, &m_apDocs, &m_adLabels,
			(long*)&m_iWords, (long*)&m_iDocs );
		return true; }
	if( sData.m_eType == SData::EPCL ) {
		for( m_iDocs = i = 0; i < sData.m_uData.m_pPCL->GetGenes( ); ++i )
			if( !sData.m_uData.m_pPCL->IsMasked( i ) && ( !sData.m_pNegative ||
				sData.m_uAnswers.m_pGenes->IsGene( sData.m_uData.m_pPCL->GetGene( i ) ) ||
				sData.m_pNegative->IsGene( sData.m_uData.m_pPCL->GetGene( i ) ) ) )
				m_iDocs++;
		m_apDocs = new DOC*[ m_iDocs ];
		m_adLabels = new double[ m_iDocs ];
		for( i = j = 0; i < sData.m_uData.m_pPCL->GetGenes( ); ++i ) {
			const string&	strGene	= sData.m_uData.m_pPCL->GetGene( i );

			if( !sData.m_uData.m_pPCL->IsMasked( i ) ) {
				if( !sData.m_pNegative ) {
					m_apDocs[ j ] = CreateDoc( sData, i );
					m_adLabels[ j++ ] = sData.m_uAnswers.m_pGenes->IsGene( strGene ) ? 1 : -1; }
				else if( sData.m_uAnswers.m_pGenes->IsGene( strGene ) ) {
					m_apDocs[ j ] = CreateDoc( sData, i );
					m_adLabels[ j++ ] = 1; }
				else if( sData.m_pNegative->IsGene( strGene ) ) {
					m_apDocs[ j ] = CreateDoc( sData, i );
					m_adLabels[ j++ ] = -1; } } }
		return true; }

	veciGenes.resize( ( sData.m_eType == SData::EPCLs ) ?
		sData.m_uData.m_pPCLs->GetGenes( ) : sData.m_uData.m_pData->GetGenes( ) );
	for( i = 0; i < veciGenes.size( ); ++i )
		veciGenes[ i ] = sData.m_uAnswers.m_pAnswers->GetGene( ( sData.m_eType == SData::EPCLs ) ?
			sData.m_uData.m_pPCLs->GetGene( i ) : sData.m_uData.m_pData->GetGene( i ) );
	for( m_iDocs = i = 0; i < veciGenes.size( ); ++i )
		if( ( iOne = veciGenes[ i ] ) != -1 )
			for( j = ( i + 1 ); j < veciGenes.size( ); ++j )
				if( ( ( iTwo = veciGenes[ j ] ) != -1 ) &&
					!CMeta::IsNaN( sData.m_uAnswers.m_pAnswers->Get( iOne, iTwo ) ) )
					m_iDocs++;
	m_apDocs = new DOC*[ m_iDocs ];
	m_adLabels = new double[ m_iDocs ];

	for( iDoc = i = 0; i < veciGenes.size( ); ++i )
		if( ( iOne = veciGenes[ i ] ) != -1 )
			for( j = ( i + 1 ); j < veciGenes.size( ); ++j )
				if( ( ( iTwo = veciGenes[ j ] ) != -1 ) &&
					!CMeta::IsNaN( d = sData.m_uAnswers.m_pAnswers->Get( iOne, iTwo ) ) ) {
					m_adLabels[ iDoc ] = d ? 1 : -1;
					m_apDocs[ iDoc++ ] = CreateDoc( sData, i, j, iDoc ); }
	assert( iDoc == m_iDocs );

	return true; }

bool CSVM::Learn( const char* szData ) {
	SData	sData;

	sData.m_eType = SData::EFile;
	sData.m_uData.m_szFile = szData;

	return CSVMImpl::Learn( sData ); }

bool CSVM::Learn( const CPCLSet& PCLs, const CDataPair& Answers ) {
	SData	sData;

	sData.m_eType = SData::EPCLs;
	sData.m_uData.m_pPCLs = &PCLs;
	sData.m_uAnswers.m_pAnswers = &Answers;

	return CSVMImpl::Learn( sData ); }

bool CSVM::Learn( const IDataset* pData, const CDataPair& Answers ) {
	SData	sData;

	sData.m_eType = SData::EData;
	sData.m_uData.m_pData = pData;
	sData.m_uAnswers.m_pAnswers = &Answers;

	return CSVMImpl::Learn( sData ); }

bool CSVM::Learn( const CPCL& PCL, const CGenes& GenesPos, const CGenes& GenesNeg ) {
	SData	sData;

	sData.m_eType = SData::EPCL;
	sData.m_uData.m_pPCL = &PCL;
	sData.m_uAnswers.m_pGenes = &GenesPos;
	sData.m_pNegative = GenesNeg.GetGenes( ) ? &GenesNeg : NULL;

	return CSVMImpl::Learn( sData ); }

bool CSVMImpl::Learn( const SData& sData ) {
	KERNEL_CACHE*	pCache;
	size_t			iWords;

	Reset( false, true, false );
	m_pModel = (MODEL*)calloc( 1, sizeof(*m_pModel) );
	if( !Initialize( sData ) )
		return false;
	iWords = GetWords( sData );

	pCache = ( m_sKernel.kernel_type == LINEAR ) ? NULL :
		kernel_cache_init( m_iDocs, m_sLearn.kernel_cache_size );
	switch( m_sLearn.type ) {
		case CLASSIFICATION:
			svm_learn_classification( m_apDocs, m_adLabels, m_iDocs, iWords,
				(LEARN_PARM*)&m_sLearn, (KERNEL_PARM*)&m_sKernel, pCache, m_pModel,
				m_adAlphas );
			break;

		case REGRESSION:
			svm_learn_regression( m_apDocs, m_adLabels, m_iDocs, iWords,
				(LEARN_PARM*)&m_sLearn, (KERNEL_PARM*)&m_sKernel, &pCache, m_pModel );
			break;

		case RANKING:
			svm_learn_ranking( m_apDocs, m_adLabels, m_iDocs, iWords, (LEARN_PARM*)&m_sLearn,
				(KERNEL_PARM*)&m_sKernel, &pCache, m_pModel );
			break;

		case OPTIMIZATION:
			svm_learn_optimization( m_apDocs, m_adLabels, m_iDocs, iWords,
				(LEARN_PARM*)&m_sLearn, (KERNEL_PARM*)&m_sKernel, pCache, m_pModel,
				m_adAlphas );
			break; }

	if( pCache )
		kernel_cache_cleanup( pCache );

	return true; }

bool CSVM::Save( ostream& ostm ) const {
	size_t		i, j;
	SVECTOR*	pVec;

	if( !m_pModel )
		return false;

	ostm << "SVM-light Version " << VERSION << endl;
	ostm << m_pModel->kernel_parm.kernel_type << " # kernel type" << endl;
	ostm << m_pModel->kernel_parm.poly_degree << " # kernel parameter -d" << endl;
	ostm << m_pModel->kernel_parm.rbf_gamma << " # kernel parameter -g" << endl;
	ostm << m_pModel->kernel_parm.coef_lin << " # kernel parameter -s" << endl;
	ostm << m_pModel->kernel_parm.coef_const << " # kernel parameter -r" << endl;
	ostm << m_pModel->kernel_parm.custom << "# kernel parameter -u" << endl;
	ostm << m_pModel->totwords << " # highest feature index" << endl;
	ostm << m_pModel->totdoc << " # number of training documents" << endl;
 
	for( i = j = 1; i < (size_t)m_pModel->sv_num; ++i )
		for( pVec = m_pModel->supvec[ i ]->fvec; pVec; pVec = pVec->next )
			j++;
	ostm << j << " # number of support vectors plus 1" << endl;
	ostm << m_pModel->b <<
		" # threshold b, each following line is a SV (starting with alpha*y)" << endl;

	for( i = 1; i < (size_t)m_pModel->sv_num; ++i )
		for( pVec = m_pModel->supvec[ i ]->fvec; pVec; pVec = pVec->next ) {
			ostm << ( m_pModel->alpha[ i ] * pVec->factor ) << ' ';
			for( j = 0; pVec->words[ j ].wnum; ++j )
				ostm << pVec->words[ j ].wnum << ':' << pVec->words[ j ].weight << ' ';
			ostm << '#' << endl; }
//			ostm << '#' << pVec->userdefined << endl; }

	return true; }

bool CSVM::Evaluate( const char* szFile, CDat& DatOut ) const {
	SData	sData;

	sData.m_eType = SData::EFile;
	sData.m_uData.m_szFile = szFile;

	return CSVMImpl::Evaluate( sData, NULL, DatOut ); }

bool CSVM::Evaluate( const CPCLSet& PCLs, CDat& DatOut ) const {
	SData	sData;

	sData.m_eType = SData::EPCLs;
	sData.m_uData.m_pPCLs = &PCLs;

	return CSVMImpl::Evaluate( sData, NULL, DatOut ); }

bool CSVM::Evaluate( const IDataset* pData, CDat& DatOut ) const {
	SData	sData;

	sData.m_eType = SData::EData;
	sData.m_uData.m_pData = pData;

	return CSVMImpl::Evaluate( sData, NULL, DatOut ); }

bool CSVM::Evaluate( const CPCLSet& PCLs, const CGenes& GenesIn, CDat& DatOut ) const {
	SData	sData;

	sData.m_eType = SData::EPCLs;
	sData.m_uData.m_pPCLs = &PCLs;

	return CSVMImpl::Evaluate( sData, &GenesIn, DatOut ); }

bool CSVM::Evaluate( const IDataset* pData, const CGenes& GenesIn, CDat& DatOut ) const {
	SData	sData;

	sData.m_eType = SData::EData;
	sData.m_uData.m_pData = pData;

	return CSVMImpl::Evaluate( sData, &GenesIn, DatOut ); }

bool CSVMImpl::Evaluate( const SData& sData, const CGenes* pGenesIn, CDat& DatOut ) const {
	size_t			i, j, iGenes;
	DOC*			pDoc;
	bool			fGene;

	if( !m_pModel )
		return false;
	if( m_pModel->kernel_parm.kernel_type == 0 )
		add_weight_vector_to_linear_model( m_pModel );

	if( sData.m_eType == SData::EFile )
		return EvaluateFile( sData.m_uData.m_szFile, DatOut );

	fGene = true;
	iGenes = ( sData.m_eType == SData::EPCLs ) ? sData.m_uData.m_pPCLs->GetGenes( ) :
		sData.m_uData.m_pData->GetGenes( );
	for( i = 0; i < iGenes; ++i ) {
		const string&	strGeneOne	= ( sData.m_eType == SData::EPCLs ) ?
			sData.m_uData.m_pPCLs->GetGene( i ) : sData.m_uData.m_pData->GetGene( i );

		if( !( i % 10 ) )
			g_CatBioUtils.notice( "CSVMImpl::Evaluate( ) gene %d/%d", i, iGenes );
		if( pGenesIn )
			fGene = pGenesIn->IsGene( strGeneOne );
		for( j = ( i + 1 ); j < iGenes; ++j ) {
			const string&	strGeneTwo	= ( sData.m_eType == SData::EPCLs ) ?
				sData.m_uData.m_pPCLs->GetGene( j ) : sData.m_uData.m_pData->GetGene( j );

			if( !fGene && !pGenesIn->IsGene( strGeneTwo ) )
				continue;
			if( !( pDoc = CreateDoc( sData, i, j, 0 ) ) )
				return false;
			DatOut.Set( i, j, (float)( m_pModel->kernel_parm.kernel_type ?
				classify_example( m_pModel, pDoc ) :
				classify_example_linear( m_pModel, pDoc ) ) );
			free_example( pDoc, 1 ); } }

	return true; }

bool CSVM::Evaluate( const CPCL& PCL, vector<float>& vecdOut, bool fMasks ) const {
	size_t	i;
	DOC*	pDoc;
	SData	sData;

	if( !m_pModel )
		return false;
	if( m_pModel->kernel_parm.kernel_type == 0 )
		add_weight_vector_to_linear_model( m_pModel );

	sData.m_eType = SData::EPCL;
	sData.m_uData.m_pPCL = &PCL;
	for( i = 0; i < PCL.GetGenes( ); ++i ) {
		if( !( i % 1000 ) )
			g_CatBioUtils.notice( "CSVMImpl::Evaluate( ) gene %d/%d", i, PCL.GetGenes( ) );
		if( fMasks && !PCL.IsMasked( i ) )
			continue;

		if( !( pDoc = CreateDoc( sData, i ) ) )
			return false;
		vecdOut.push_back( (float)( m_pModel->kernel_parm.kernel_type ?
			classify_example( m_pModel, pDoc ) :
			classify_example_linear( m_pModel, pDoc ) ) );
		free_example( pDoc, 1 ); }

	return true; }

bool CSVMImpl::EvaluateFile( const char* szFile, CDat& DatOut ) const {
	static const size_t	c_iSize	= 512;
	char			szGene[ c_iSize ];
	SWORD*			asWords;
	char*			pc;
	ifstream		ifsm;
	vector<string>	vecstrGenes;
	uint32_t		i, j, k, iDocs, iWords, iGenes;
	float*			ad;
	DOC*			pDoc;

	ifsm.open( szFile, ios_base::binary );
	if( !ifsm.is_open( ) )
		return false;
	ifsm.read( (char*)&iWords, sizeof(iWords) );
	ifsm.read( (char*)&iDocs, sizeof(iDocs) );
	ifsm.seekg( iDocs * ( ( ( iWords + 1 ) * sizeof(float) ) +
		( 3 * sizeof(iDocs) ) ), ios_base::cur );
	ifsm.read( (char*)&iGenes, sizeof(iGenes) );
	vecstrGenes.resize( iGenes );
	for( i = 0; i < iGenes; ++i ) {
		for( pc = szGene; ; ++pc ) {
			ifsm.read( pc, 1 );
			if( !*pc )
				break; }
		vecstrGenes[ i ] = szGene; }
	DatOut.Open( vecstrGenes );

	asWords = ( iWords >= c_iWords ) ? new SWORD[ iWords + 1 ] : s_asWords;
	for( i = 0; i < iWords; ++i )
		asWords[ i ].wnum = i + 1;
	asWords[ i ].wnum = 0;

	ad = new float[ iWords + 1 ];
	ifsm.seekg( 2 * sizeof(iDocs), ios_base::beg );
	for( i = 0; i < iDocs; ++i ) {
		if( !( i % 1000 ) )
			g_CatBioUtils.notice( "CSVMImpl::EvaluateFile( %s ) pair %d/%d", szFile, i,
				iDocs );
		ifsm.read( (char*)ad, ( iWords + 1 ) * sizeof(*ad) );
		for( j = 0; j < iWords; ++j )
			asWords[ j ].weight = ad[ j + 1 ];
		pDoc = create_example( i, 0, 0, 1, create_svector( asWords, "", 1 ) );
		ifsm.read( (char*)&j, sizeof(j) );
		ifsm.read( (char*)&j, sizeof(j) );
		ifsm.read( (char*)&k, sizeof(k) );
		DatOut.Set( j, k, (float)( m_pModel->kernel_parm.kernel_type ?
			classify_example( m_pModel, pDoc ) :
			classify_example_linear( m_pModel, pDoc ) ) );
		free_example( pDoc, 1 ); }
	delete[] ad;

	if( asWords != s_asWords )
		delete[] asWords;

	return true; }

bool CSVM::Open( istream& istm ) {
	static const size_t	c_iBuf	= 2048;
	char			szBuf[ c_iBuf ];
	vector<string>	vecstrLine, vecstrToken;
	SWORD*			asWords;
	size_t			i, j;

	Reset( false, true, true );
	m_pModel = (MODEL*)calloc( 1, sizeof(*m_pModel) );

	istm.getline( szBuf, c_iBuf - 1 );
	istm >> m_pModel->kernel_parm.kernel_type;
	istm.getline( szBuf, c_iBuf - 1 );
	istm >> m_pModel->kernel_parm.poly_degree;
	istm.getline( szBuf, c_iBuf - 1 );
	istm >> m_pModel->kernel_parm.rbf_gamma;
	istm.getline( szBuf, c_iBuf - 1 );
	istm >> m_pModel->kernel_parm.coef_lin;
	istm.getline( szBuf, c_iBuf - 1 );
	istm >> m_pModel->kernel_parm.coef_const;
	istm.getline( szBuf, c_iBuf - 1 );
	istm.getline( szBuf, c_iBuf - 1 );
	CMeta::Tokenize( szBuf, vecstrLine, "#", true );
	if( vecstrLine.size( ) > 1 )
		strcpy_s( m_pModel->kernel_parm.custom, 49, vecstrLine[ 0 ].c_str( ) );
	istm >> m_pModel->totwords;
	istm.getline( szBuf, c_iBuf - 1 );
	istm >> m_pModel->totdoc;
	istm.getline( szBuf, c_iBuf - 1 );
	istm >> m_pModel->sv_num;
	istm.getline( szBuf, c_iBuf - 1 );
	istm >> m_pModel->b;
	istm.getline( szBuf, c_iBuf - 1 );

	m_pModel->supvec = (DOC**)malloc( m_pModel->sv_num * sizeof(*m_pModel->supvec) );
	m_pModel->alpha = (double*)malloc( m_pModel->sv_num * sizeof(*m_pModel->alpha) );
	m_pModel->index = NULL;
	m_pModel->lin_weights = NULL;

	asWords = new SWORD[ m_pModel->totwords + 1 ];
	asWords[ m_pModel->totwords ].wnum = 0;
	for( i = 1; i < (size_t)m_pModel->sv_num; ++i ) {
		istm.getline( szBuf, c_iBuf - 1 );
		vecstrLine.clear( );
		CMeta::Tokenize( szBuf, vecstrLine, CMeta::c_szWS, true );
		if( vecstrLine.size( ) != ( m_pModel->totwords + 2 ) ) {
			g_CatBioUtils.error( "CSVM::Open( ) wanted %d words but only found %d on line: %s",
				( m_pModel->totwords + 2 ), vecstrLine.size( ), szBuf );
			delete[] asWords;
			return false; }
		m_pModel->alpha[ i ] = atof( vecstrLine[ 0 ].c_str( ) );
		for( j = 1; ( j + 1 ) < vecstrLine.size( ); ++j ) {
			vecstrToken.clear( );
			CMeta::Tokenize( vecstrLine[ j ].c_str( ), vecstrToken, ":", true );
			if( vecstrToken.size( ) != 2 ) {
				g_CatBioUtils.error( "CSVM::Open( ) found illegal token \"%s\" on line: %s",
					vecstrLine[ j ].c_str( ), szBuf );
				delete[] asWords;
				return false; }
			asWords[ j - 1 ].wnum = atoi( vecstrToken[ 0 ].c_str( ) );
			asWords[ j - 1 ].weight = (float)atof( vecstrToken[ 1 ].c_str( ) ); }
		m_pModel->supvec[ i ] = create_example( -1, 0, 0, 0, create_svector( asWords, "",
			1 ) ); }

	delete[] asWords;
	return true; }

void CSVM::SetIterations( size_t iIter ) {

	m_sLearn.maxiter = iIter; }

void CSVM::SetCache( size_t iCache ) {

	m_sLearn.kernel_cache_size = iCache; }

void CSVM::SetTradeoff( float dTradeoff ) {

	m_sLearn.svm_c = dTradeoff; }

void CSVM::SetGamma( float dGamma ) {

	m_sKernel.rbf_gamma = dGamma; }

void CSVM::SetDegree( size_t iDegree ) {

	m_sKernel.poly_degree = iDegree; }

void CSVM::SetKernel( EKernel eKernel ) {

	m_sKernel.kernel_type = eKernel; }

void CSVM::SetVerbosity( size_t iVerbosity ) {

	verbosity = iVerbosity; }

}

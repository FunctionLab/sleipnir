#include "stdafx.h"
#include "svm.h"
#include "pclset.h"
#include "datapair.h"
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

CSVMImpl::SLearn::SLearn( ) {

	strcpy( predfile, "" );
	strcpy( alphafile, "" );
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
  strcpy( custom, "" ); }

size_t CSVMImpl::GetWords( const CPCLSet& Data ) {
	size_t	i, iRet;

	for( iRet = i = 0; i < Data.GetPCLs( ); ++i )
		iRet += Data.GetExperiments( i );

	return ( iRet * 2 ); }

DOC* CSVMImpl::CreateDoc( const CPCLSet& Data, size_t iOne, size_t iTwo, size_t iDoc ) {
	static const size_t	c_iWords	= 1024;
	static SWORD		l_asWords[ c_iWords ];
	SWORD*	asWords;
	size_t	i, j, iWord, iWords;
	float	d;
	DOC*	pRet;

	asWords = ( ( iWords = GetWords( Data ) ) >= c_iWords ) ? new SWORD[ iWords + 1 ] :
		l_asWords;

	for( i = 0; i < iWords; ++i )
		asWords[ i ].wnum = i + 1;
	asWords[ i ].wnum = 0;
	for( iWord = i = 0; i < Data.GetPCLs( ); ++i ) {
		for( j = 0; j < Data.GetExperiments( i ); ++j ) {
			if( CMeta::IsNaN( d = Data.Get( i, iOne, j ) ) )
				d = 0;
			assert( ( iWord + j ) < iWords );
			asWords[ iWord + j ].weight = d;
			if( CMeta::IsNaN( d = Data.Get( i, iTwo, j ) ) )
				d = 0;
			assert( ( iWord + ( iWords / 2 ) + j ) < iWords );
			asWords[ iWord + ( iWords / 2 ) + j ].weight = d; }
		iWord += Data.GetExperiments( i ); }

	pRet = create_example( iDoc, 0, 0, 1, create_svector( asWords, "", 1 ) );
	if( asWords != l_asWords )
		delete[] asWords;
	return pRet; }

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

bool CSVM::OpenAlphas( istream& istm ) {
	static const size_t	c_iBuf	= 1024;
	char	szBuf[ c_iBuf ];
	vector<float>	vecdAlphas;

	Reset( false, false, true );
	while( istm.peek( ) != EOF ) {
		istm.getline( szBuf, c_iBuf - 1 );
		vecdAlphas.push_back( (float)atof( szBuf ) ); }
	m_adAlphas = new double[ m_iAlphas = vecdAlphas.size( ) ];
	copy( vecdAlphas.begin( ), vecdAlphas.end( ), m_adAlphas );

	return true; }

bool CSVMImpl::Initialize( const CPCLSet& Data, const CDataPair& Answers ) {
	size_t			i, j, iOne, iTwo, iDoc;
	vector<size_t>	veciGenes;
	float			d;

	Reset( true, false, false );
	veciGenes.resize( Data.GetGenes( ) );
	for( i = 0; i < veciGenes.size( ); ++i )
		veciGenes[ i ] = Answers.GetGene( Data.GetGene( i ) );
	for( m_iDocs = i = 0; i < Data.GetGenes( ); ++i )
		if( ( iOne = veciGenes[ i ] ) != -1 )
			for( j = ( i + 1 ); j < Data.GetGenes( ); ++j )
				if( ( ( iTwo = veciGenes[ j ] ) != -1 ) &&
					!CMeta::IsNaN( Answers.Get( iOne, iTwo ) ) )
					m_iDocs++;
	m_apDocs = new DOC*[ m_iDocs ];
	m_adLabels = new double[ m_iDocs ];

	for( iDoc = i = 0; i < Data.GetGenes( ); ++i )
		if( ( iOne = veciGenes[ i ] ) != -1 )
			for( j = ( i + 1 ); j < Data.GetGenes( ); ++j )
				if( ( ( iTwo = veciGenes[ j ] ) != -1 ) &&
					!CMeta::IsNaN( d = Answers.Get( iOne, iTwo ) ) ) {
					m_adLabels[ iDoc ] = d ? 1 : -1;
					m_apDocs[ iDoc++ ] = CreateDoc( Data, i, j, iDoc ); }
	assert( iDoc == m_iDocs );

	return true; }

bool CSVM::Learn( const CPCLSet& Data, const CDataPair& Answers ) {
	KERNEL_CACHE*	pCache;

	Reset( false, true, false );
	m_pModel = (MODEL*)calloc( 1, sizeof(*m_pModel) );
	if( !Initialize( Data, Answers ) )
		return false;

	pCache = ( m_sKernel.kernel_type == LINEAR ) ? NULL :
		kernel_cache_init( m_iDocs, m_sLearn.kernel_cache_size );
	switch( m_sLearn.type ) {
		case CLASSIFICATION:
			svm_learn_classification( m_apDocs, m_adLabels, m_iDocs, GetWords( Data ),
				(LEARN_PARM*)&m_sLearn, (KERNEL_PARM*)&m_sKernel, pCache, m_pModel,
				m_adAlphas );
			break;

		case REGRESSION:
			svm_learn_regression( m_apDocs, m_adLabels, m_iDocs, GetWords( Data ),
				(LEARN_PARM*)&m_sLearn, (KERNEL_PARM*)&m_sKernel, &pCache, m_pModel );
			break;

		case RANKING:
			svm_learn_ranking( m_apDocs, m_adLabels, m_iDocs, GetWords( Data ), (LEARN_PARM*)&m_sLearn,
				(KERNEL_PARM*)&m_sKernel, &pCache, m_pModel );
			break;

		case OPTIMIZATION:
			svm_learn_optimization( m_apDocs, m_adLabels, m_iDocs, GetWords( Data ),
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

bool CSVM::Evaluate( const CPCLSet& Data, const CGenes& GenesIn, CDat& DatOut ) const {

	return CSVMImpl::Evaluate( Data, &GenesIn, DatOut ); }

bool CSVM::Evaluate( const CPCLSet& Data, CDat& DatOut ) const {

	return CSVMImpl::Evaluate( Data, NULL, DatOut ); }

bool CSVMImpl::Evaluate( const CPCLSet& Data, const CGenes* pGenesIn, CDat& DatOut ) const {
	size_t	i, j;
	DOC*	pDoc;
	bool	fGene;

	if( !m_pModel )
		return false;
	if( m_pModel->kernel_parm.kernel_type == 0 )
		add_weight_vector_to_linear_model( m_pModel );

	fGene = true;
	for( i = 0; i < Data.GetGenes( ); ++i ) {
		if( !( i % 10 ) )
			g_CatBioUtils.notice( "CSVMImpl::Evaluate( ) gene %d/%d", i, Data.GetGenes( ) );
		if( pGenesIn )
			fGene = pGenesIn->IsGene( Data.GetGene( i ) );
		for( j = ( i + 1 ); j < Data.GetGenes( ); ++j ) {
			if( !fGene && !pGenesIn->IsGene( Data.GetGene( j ) ) )
				continue;
			if( !( pDoc = CreateDoc( Data, i, j, 0 ) ) )
				return false;
			DatOut.Set( i, j, (float)( m_pModel->kernel_parm.kernel_type ?
				classify_example( m_pModel, pDoc ) :
				classify_example_linear( m_pModel, pDoc ) ) );
			free_example( pDoc, 1 ); } }

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
		strcpy( m_pModel->kernel_parm.custom, vecstrLine[ 0 ].c_str( ) );
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

void CSVM::SetKernel( EKernel eKernel ) {

	m_sKernel.kernel_type = eKernel; }

}

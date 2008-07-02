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
#ifndef BNSERVER_H
#define BNSERVER_H

class CDot;

struct SBNServerData {
	const IOntology**				m_apOntologies;
	const CBayesNetMinimal&			m_BNDefault;
	const CDatabase&				m_Database;
	const CGenome&					m_Genome;
	size_t							m_iLimit;
	const CDataMatrix&				m_MatBackgrounds;
	const CDataMatrix&				m_MatBetweenCC;
	const CDataMatrix&				m_MatBetweenDC;
	const CDataMatrix&				m_MatBetweenDD;
	const CDataMatrix&				m_MatParameters;
	const CDataMatrix&				m_MatWithinC;
	const CDataMatrix&				m_MatWithinD;
	string							m_strFiles;
	string							m_strGraphviz;
	const vector<CBayesNetMinimal>&	m_vecBNs;
	const vector<float>&			m_vecdPriors;
	const vector<size_t>&			m_veciContexts;
	const vector<size_t>&			m_veciDiseases;
	const vector<vector<size_t> >&	m_vecveciContexts;
	const vector<vector<size_t> >&	m_vecveciDiseases;

	SBNServerData( 
		const IOntology**				apOntologies,
		const CBayesNetMinimal&			BNDefault,
		const CDatabase&				Database,
		const CGenome&					Genome,
		size_t							iLimit,
		const CDataMatrix&				MatBackgrounds,
		const CDataMatrix&				MatBetweenCC,
		const CDataMatrix&				MatBetweenDC,
		const CDataMatrix&				MatBetweenDD,
		const CDataMatrix&				MatParameters,
		const CDataMatrix&				MatWithinC,
		const CDataMatrix&				MatWithinD,
		const string&					strFiles,
		const string&					strGraphviz,
		const vector<CBayesNetMinimal>&	vecBNs,
		const vector<float>&			vecdPriors,
		const vector<size_t>&			veciContexts,
		const vector<size_t>&			veciDiseases,
		const vector<vector<size_t> >&	vecveciContexts,
		const vector<vector<size_t> >&	vecveciDiseases
	) :
		m_apOntologies		(apOntologies),
		m_BNDefault			(BNDefault),
		m_Database			(Database),
		m_Genome			(Genome),
		m_iLimit			(iLimit),
		m_MatBackgrounds	(MatBackgrounds),
		m_MatBetweenCC		(MatBetweenCC),
		m_MatBetweenDC		(MatBetweenDC),
		m_MatBetweenDD		(MatBetweenDD),
		m_MatParameters		(MatParameters),
		m_MatWithinC		(MatWithinC),
		m_MatWithinD		(MatWithinD),
		m_strFiles			(strFiles),
		m_strGraphviz		(strGraphviz),
		m_vecBNs			(vecBNs),
		m_vecdPriors		(vecdPriors),
		m_veciContexts		(veciContexts),
		m_veciDiseases		(veciDiseases),
		m_vecveciContexts	(vecveciContexts),
		m_vecveciDiseases	(vecveciDiseases)
		{ }
};

class CBNServer : public IServerClient {
public:
	static bool Get( size_t, size_t, float*, const Sleipnir::CDatabase&,
		const std::vector<Sleipnir::CBayesNetMinimal>&, const Sleipnir::CBayesNetMinimal& );

	CBNServer( SOCKET, const string&, const SBNServerData& );
	~CBNServer( );

	IServerClient* NewInstance( SOCKET, uint32_t, uint16_t );
	void Destroy( );
	bool ProcessMessage( const std::vector<unsigned char>& );
	bool GenerateNetworkIcons( ) const;
	bool GenerateAssociations( const char*, size_t );

private:
	typedef size_t (CBNServer::*TPFNProcessor)( const std::vector<unsigned char>&, size_t );

	enum EGraphOutput {
		EGraphOutputFile,
		EGraphOutputSocket,
		EGraphOutputNamed
	};

	enum ESetType {
		ESetGenes	= 0,
		ESetContext	= ESetGenes + 1,
		ESetDisease	= ESetContext + 1
	};

	static const size_t			c_iDegree			= 1;
	static const size_t			c_iValues			= 4;
	static const size_t			c_iOverestimate		= 100;
	static const TPFNProcessor	c_apfnProcessors[];
	static const size_t			c_iProcessors;
	static const float			c_dCutoff;
	static const float			c_adColorMin[];
	static const float			c_adColorMax[];

	template<class tType>
	static bool SetPointer( tType* pPointer, tType Value ) {

		if( pPointer )
			*pPointer = Value;

		return !!pPointer; }

	template<class tType>
	static bool PushPointer( std::vector<tType>* pvecPointer, tType Value ) {

		if( pvecPointer )
			pvecPointer->push_back( Value );

		return !!pvecPointer; }

	template<class tType>
	static bool PushPointer( tType* pPointer, std::vector<tType>& vecValues, tType Value ) {

		if( pPointer )
			vecValues.push_back( Value );

		return !!pPointer; }

	template<class tType>
	static void SetArray( tType* ad, size_t iChunk, size_t iOffset, float dA, float dB, float dC, float dD ) {
		float	adParams[]	= {dA, dB, dC, dD};
		size_t	i;

		for( i = 0; i < ARRAYSIZE(adParams); ++i )
			ad[ ( i * iChunk ) + iOffset ] = adParams[ i ]; }

	static float GetMatrix( const CDataMatrix& Mat, size_t iR, size_t iC ) {

		return ( ( ( iR < Mat.GetRows( ) ) && ( iC < Mat.GetColumns( ) ) ) ? Mat.Get( iR, iC ) :
			CMeta::GetNaN( ) ); }

	static float GetBetween( const CDataMatrix& Mat, size_t iR, size_t iOne, size_t iTwo,
		size_t iChunk ) {
		float	d;

		return ( CMeta::IsNaN( d = GetMatrix( Mat, iR, ( iOne * iChunk ) + iTwo ) ) ? 0 : d ); }

	template<class tType>
	static bool Winsorize( std::vector<tType>& vecValues ) {

		return CStatistics::Winsorize( vecValues, ( vecValues.size( ) / 10 ) + 1 ); }

	// Utility
	bool Get( size_t, size_t, float* = NULL );
	bool Get( size_t, const std::vector<size_t>&, size_t, float* );
	bool GetGenes( const std::vector<size_t>&, size_t, float );
	bool GetWithin( const std::vector<size_t>&, size_t, float*, std::vector<float>* ) const;
	// Association processing
	bool GetAssociationsSet( unsigned char, const std::vector<size_t>&, size_t ) const;
	bool GetAssociationsDC( unsigned char, unsigned char, size_t, size_t, bool = false ) const;
	bool GetAssociation( size_t, const std::vector<unsigned char>&, const std::vector<size_t>&, size_t, bool,
		float*, float*, float*, std::vector<float>*, std::vector<float>*, float ) const;
	bool GetAssociation( const std::vector<size_t>&, const std::vector<size_t>&, size_t, float&,
		float&, float& ) const;
	// Graph processing
	bool GraphCreate( const std::vector<size_t>&, size_t, size_t, float, std::vector<bool>&,
		std::vector<size_t>&, Sleipnir::CDat& ) const;
	bool GraphWrite( const Sleipnir::CDat&, const std::vector<size_t>&, const std::vector<size_t>&,
		const std::vector<bool>&, size_t, EGraphOutput ) const;
	bool SelectNeighborsPixie( const std::vector<size_t>&, const std::vector<bool>&, size_t, size_t,
		const Sleipnir::CDataMatrix&, std::vector<size_t>& ) const;
	bool SelectNeighborsRatio( const std::vector<size_t>&, const std::vector<bool>&, size_t, size_t,
		const Sleipnir::CDataMatrix&, std::vector<size_t>& ) const;
	bool SendGenes( const std::vector<size_t>&, const std::vector<size_t>& ) const;
	// Message processors
	size_t ProcessInference( const std::vector<unsigned char>&, size_t );
	size_t ProcessData( const std::vector<unsigned char>&, size_t );
	size_t ProcessGraph( const std::vector<unsigned char>&, size_t );
	size_t ProcessContexts( const std::vector<unsigned char>&, size_t );
	size_t ProcessTermFinder( const std::vector<unsigned char>&, size_t );
	size_t ProcessDiseases( const std::vector<unsigned char>&, size_t );
	size_t ProcessGenes( const std::vector<unsigned char>&, size_t );
	size_t ProcessAssociation( const std::vector<unsigned char>&, size_t );
	size_t ProcessAssociations( const std::vector<unsigned char>&, size_t );

	size_t GetGenes( ) const {

		return m_sData.m_Database.GetGenes( ); }

	size_t GetContexts( ) const {

		return m_sData.m_vecveciContexts.size( ); }

	size_t GetDiseases( ) const {

		return m_sData.m_vecveciDiseases.size( ); }

	float GetBackground( size_t iContext, size_t iGene ) const {

		return ( ( ( iContext < m_sData.m_MatBackgrounds.GetRows( ) ) &&
			( iGene < m_sData.m_MatBackgrounds.GetColumns( ) ) ) ?
			m_sData.m_MatBackgrounds.Get( iContext, iGene ) : 1 ); }

	size_t InitializeDiseases( ) {
		size_t	iRet;

		iRet = c_iValues * GetDiseases( );
		if( !m_adDiseases )
			m_adDiseases = new float[ iRet ];

		return iRet; }

	size_t InitializeGenes( ) {
		size_t	iRet;

		iRet = c_iValues * GetGenes( );
		if( !m_adGenes )
			m_adGenes = new float[ iRet ];

		return iRet; }

	size_t InitializeContexts( ) {
		size_t	iRet;

		iRet = c_iValues * GetContexts( );
		if( !m_adContexts )
			m_adContexts = new float[ iRet ];

		return iRet; }

	const CBayesNetMinimal& GetBN( size_t iContext ) const {

		return ( ( iContext && GetContexts( ) ) ? m_sData.m_vecBNs[ ( iContext - 1 ) % GetContexts( ) ] :
			m_sData.m_BNDefault ); }

	const std::string& GetGene( size_t iGene ) const {

		return m_sData.m_Database.GetGene( iGene ); }

	size_t GetGene( const std::string& strGene ) const {

		return m_sData.m_Database.GetGene( strGene ); }

	float GetFraction( size_t iSize ) const {

		return ( ( iSize > m_sData.m_iLimit ) ? ( (float)m_sData.m_iLimit / iSize ) : 1 ); }

	bool IsFraction( float dFraction ) const {

		return ( ( dFraction < 1 ) && ( ( (float)rand( ) / RAND_MAX ) > dFraction ) ); }

	size_t GetContext( unsigned char bDiseases, size_t iContext, size_t iCurrent ) const {

		return ( ( bDiseases || ( iContext != -1 ) ) ? iContext : ( iCurrent + 1 ) ); }

	float GetPValue( float dBetween, float dBackground, float dWithin, size_t iContext, size_t iSmall,
		size_t iBig, size_t iCount = 0, bool fZ = false ) const {
//		static const size_t	c_iCutoff	= 75;
		const float*	adParams;
		float			dValue, dStd, dRet, dA, dB, dC;

		dValue = GetPrior( iContext ) * dBetween / dBackground / dWithin;
		if( CMeta::IsNaN( dValue ) )
			return 1;
		if( iSmall > iBig )
			swap( iBig, iSmall );
//		if( iBig > c_iCutoff ) {
//			iSmall = max( (size_t)( iSmall * ( (float)c_iCutoff / iBig ) ), (size_t)1 );
//			iBig = c_iCutoff; }
		adParams = GetParameters( iContext );
		dA = ( ( adParams[ 0 ] * iSmall ) + adParams[ 1 ] ) / ( iSmall + 1 );
		dB = adParams[ 2 ];
		dC = ( ( adParams[ 3 ] * iSmall ) + adParams[ 4 ] ) / ( iSmall + adParams[ 5 ] );
		dStd = ( ( dA * iBig ) + dB ) / ( iBig + dC );
		dRet = fZ ? ( ( dValue - 1 ) / dStd ) : ( 1 - (float)CStatistics::NormalCDF( dValue, 1, dStd ) );
/*
cerr << iContext << ":	" << iBig << '\t' << iSmall << endl;
cerr << dBetween << '\t' << dBackground << '\t' << dWithin << '\t' << GetPrior( iContext ) << endl;
for( size_t i = 0; i < 6; ++i )
cerr << ( i ? "\t" : "" ) << adParams[ i ];
cerr << endl << dA << '\t' << dB << '\t' << dC << '\t' << dStd << endl;
cerr << dValue << ":	" << dRet << endl;
//*/
		if( !fZ && iCount )
			dRet *= iCount;

		return dRet; }

	const float* GetParameters( size_t iContext ) const {

		return m_sData.m_MatParameters.Get( iContext ); }

	float GetPrior( size_t iContext ) const {

// Alternative: calculate priors in the current context rather than the global context
		iContext = 0;
		return m_sData.m_vecdPriors[ iContext ]; }

	const CDatabase& GetDatabase( ) const {

		return m_sData.m_Database; }

	const string& GetFiles( ) const {

		return m_sData.m_strFiles; }

	const vector<size_t>& GetContext( size_t iContext ) const {

		return m_sData.m_vecveciContexts[ iContext ]; }

	const vector<size_t>& GetDisease( size_t iDisease ) const {

		return m_sData.m_vecveciDiseases[ iDisease ]; }

	CGenome& GetGenome( ) const {

		return *(CGenome*)&m_sData.m_Genome; }

	const IOntology* GetOntology( size_t iOntology ) const {
		size_t	i;

		for( i = 0; ( i < iOntology ) && m_sData.m_apOntologies[ i ]; ++i );
		return m_sData.m_apOntologies[ i ]; }

	size_t GetDatachunk( ) const {

		return ( ( GetDatabase( ).GetDatasets( ) + 1 ) / 2 ); }

	const vector<size_t>& GetDiseaseGenes( ) const {

		return m_sData.m_veciDiseases; }

	const vector<vector<size_t> >& GetGeneSets( unsigned char bDiseases ) const {

		return ( bDiseases ? m_sData.m_vecveciDiseases : m_sData.m_vecveciContexts ); }

	float GetWithinContext( size_t iContext, size_t iSet ) const {

// Alternative: calculate within scores in the current context rather than the global context
		iContext = 0;
		return GetMatrix( m_sData.m_MatWithinC, iContext, iSet ); }

	float GetWithinDisease( size_t iContext, size_t iSet ) const {

// Alternative: calculate within scores in the current context rather than the global context
		iContext = 0;
		return GetMatrix( m_sData.m_MatWithinD, iContext, iSet ); }

	float GetWithin( unsigned char bDiseases, size_t iContext, size_t iSet ) const {

		return ( bDiseases ? GetWithinDisease( iContext, iSet ) : GetWithinContext( iContext, iSet ) ); }

	float GetBetween( size_t iContext, unsigned char bDiseaseOne, size_t iOne, unsigned char bDiseaseTwo,
		size_t iTwo ) const {

		return ( bDiseaseOne ?
			( bDiseaseTwo ? GetBetweenDD( iContext, iOne, iTwo ) : GetBetweenDC( iContext, iOne, iTwo ) ) :
			( bDiseaseTwo ? GetBetweenDC( iContext, iTwo, iOne ) : GetBetweenCC( iContext, iOne, iTwo ) ) ); }

	float GetBetweenCC( size_t iContext, size_t iOne, size_t iTwo ) const {

		return GetBetween( m_sData.m_MatBetweenCC, iContext, iOne, iTwo, GetContexts( ) ); }

	float GetBetweenDD( size_t iContext, size_t iOne, size_t iTwo ) const {

		return GetBetween( m_sData.m_MatBetweenDD, iContext, iOne, iTwo, GetDiseases( ) ); }

	float GetBetweenDC( size_t iContext, size_t iDisease, size_t iProcess ) const {

		return GetBetween( m_sData.m_MatBetweenDC, iContext, iDisease, iProcess, GetContexts( ) ); }

	float*					m_adContexts;
	float*					m_adDiseases;
	float*					m_adGenes;
	SOCKET					m_iSocket;
	const SBNServerData&	m_sData;
	string					m_strConnection;
};

#endif // BNSERVER_H

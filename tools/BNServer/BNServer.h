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

class CBNServer : public IServerClient {
public:
	static bool Get( size_t, size_t, float*, const Sleipnir::CDatabase&,
		const std::vector<Sleipnir::CBayesNetMinimal>&, const Sleipnir::CBayesNetMinimal& );

	CBNServer( const Sleipnir::CBayesNetMinimal&, const std::vector<Sleipnir::CBayesNetMinimal>&,
		const std::vector<std::vector<size_t> >&, SOCKET, const Sleipnir::CDatabase&, const string&,
		const char*, const char*, const Sleipnir::CDataMatrix&, const Sleipnir::CGenome&,
		const Sleipnir::IOntology**, const std::vector<std::vector<size_t> >&, const std::vector<size_t>&,
		const std::vector<size_t>&, size_t );
	~CBNServer( );

	IServerClient* NewInstance( SOCKET, uint32_t, uint16_t );
	void Destroy( );
	bool ProcessMessage( const std::vector<unsigned char>& );
	bool GenerateNetworkIcons( ) const;

private:
	typedef size_t (CBNServer::*TPFNProcessor)( const std::vector<unsigned char>&, size_t );

	enum EGraphOutput {
		EGraphOutputFile,
		EGraphOutputSocket,
		EGraphOutputNamed
	};

	static const size_t			c_iDegree			= 1;
	static const TPFNProcessor	c_apfnProcessors[];
	static const size_t			c_iProcessors;
	static const float			c_dCutoff;
	static const float			c_adColorMin[];
	static const float			c_adColorMax[];

	// Utility
	bool Get( size_t, size_t, float* = NULL );
	bool Get( size_t, const std::vector<size_t>&, size_t, float* );
	bool GetGenes( const std::vector<size_t>&, size_t );
	bool GetAssociation( const std::vector<unsigned char>&, float, const std::vector<size_t>&, size_t, bool,
		float*, float*, std::vector<float>*, std::vector<float>* ) const;
	bool GetAssociation( const std::vector<size_t>&, const std::vector<size_t>&, size_t, float&,
		float& ) const;
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

		return m_Database.GetGenes( ); }

	size_t GetContexts( ) const {

		return m_vecveciContexts.size( ); }

	size_t GetDiseases( ) const {

		return m_vecveciDiseases.size( ); }

	float GetBackground( size_t iContext, size_t iGene ) const {

		return ( ( m_MatBackgrounds.GetColumns( ) && m_MatBackgrounds.GetRows( ) ) ?
			m_MatBackgrounds.Get( iContext, iGene ) : 1 ); }

	size_t InitializeDiseases( ) {
		size_t	iRet;

		iRet = 2 * m_vecveciDiseases.size( );
		if( !m_adDiseases )
			m_adDiseases = new float[ iRet ];

		return iRet; }

	size_t InitializeGenes( ) {
		size_t	iRet;

		iRet = 2 * GetGenes( );
		if( !m_adGenes )
			m_adGenes = new float[ iRet ];

		return iRet; }

	size_t InitializeContexts( ) {
		size_t	iRet;

		iRet = 2 * m_vecBNs.size( );
		if( !m_adContexts )
			m_adContexts = new float[ iRet ];

		return iRet; }

	const CBayesNetMinimal& GetBN( size_t iContext ) const {

		return ( ( iContext && m_vecBNs.size( ) ) ? m_vecBNs[ ( iContext - 1 ) % m_vecBNs.size( ) ] :
			m_BNDefault ); }

	const std::string& GetGene( size_t iGene ) const {

		return m_Database.GetGene( iGene ); }

	size_t GetGene( const std::string& strGene ) const {

		return m_Database.GetGene( strGene ); }

	float GetFraction( size_t iSize ) const {

		return ( ( iSize > m_iLimit ) ? ( (float)m_iLimit / iSize ) : 1 ); }

	bool IsFraction( float dFraction ) const {

		return ( ( dFraction < 1 ) && ( ( (float)rand( ) / RAND_MAX ) > dFraction ) ); }

	size_t GetContext( unsigned char bDiseases, size_t iContext, size_t iCurrent ) const {

		return ( ( bDiseases || ( iContext != -1 ) ) ? iContext : ( iCurrent + 1 ) ); }

	const CBayesNetMinimal&			m_BNDefault;
	const vector<CBayesNetMinimal>&	m_vecBNs;
	SOCKET							m_iSocket;
	float*							m_adGenes;
	float*							m_adContexts;
	float*							m_adDiseases;
	const CDatabase&				m_Database;
	string							m_strConnection;
	string							m_strGraphviz;
	string							m_strFiles;
	const vector<vector<size_t> >&	m_vecveciContexts;
	const CDataMatrix&				m_MatBackgrounds;
	const CGenome&					m_Genome;
	const IOntology**				m_apOntologies;
	size_t							m_iOntologies;
	const vector<vector<size_t> >	m_vecveciDiseases;
	const vector<size_t>&			m_veciDiseases;
	const vector<size_t>&			m_veciContexts;
	size_t							m_iLimit;
};

#endif // BNSERVER_H

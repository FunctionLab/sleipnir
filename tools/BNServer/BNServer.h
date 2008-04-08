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
		const Sleipnir::CCompactFullMatrix&, SOCKET, const Sleipnir::CDatabase&, const string&,
		const char*, const char*, const Sleipnir::CDataMatrix&, const Sleipnir::CGenome&,
		const Sleipnir::IOntology**, const std::vector<std::vector<size_t> >& );
	~CBNServer( );

	IServerClient* NewInstance( SOCKET, uint32_t, uint16_t );
	void Destroy( );
	bool ProcessMessage( const std::vector<unsigned char>& );

private:
	typedef size_t (CBNServer::*TPFNProcessor)( const std::vector<unsigned char>&, size_t );

	static const size_t			c_iDegree			= 1;
	static const TPFNProcessor	c_apfnProcessors[];
	static const size_t			c_iProcessors;
	static const float			c_dCutoff;
	static const float			c_adColorMin[];
	static const float			c_adColorMax[];

	bool Get( size_t, size_t, float* = NULL );
	bool Get( size_t, const std::vector<size_t>&, size_t, float* );
	bool GetContext( size_t, const std::vector<unsigned char>&, size_t );
	bool GetDisease( size_t, size_t, const std::vector<unsigned char>&, size_t );
	bool GraphCreate( const std::vector<size_t>&, size_t, size_t, std::vector<bool>&, std::vector<size_t>&,
		Sleipnir::CDat& ) const;
	bool GraphWrite( const Sleipnir::CDat&, const std::vector<size_t>&, const std::vector<size_t>&,
		const std::vector<bool>&, size_t, bool ) const;
	bool SelectNeighborsPixie( const std::vector<size_t>&, const std::vector<bool>&, size_t, size_t,
		const Sleipnir::CDataMatrix&, std::vector<size_t>& ) const;
	bool SelectNeighborsRatio( const std::vector<size_t>&, const std::vector<bool>&, size_t, size_t,
		const Sleipnir::CDataMatrix&, std::vector<size_t>& ) const;
	bool SendGenes( const std::vector<size_t>&, const std::vector<size_t>& ) const;
	size_t ProcessInference( const std::vector<unsigned char>&, size_t );
	size_t ProcessData( const std::vector<unsigned char>&, size_t );
	size_t ProcessGraph( const std::vector<unsigned char>&, size_t );
	size_t ProcessContexts( const std::vector<unsigned char>&, size_t );
	size_t ProcessTermFinder( const std::vector<unsigned char>&, size_t );
	size_t ProcessDiseases( const std::vector<unsigned char>&, size_t );

	size_t GetGenes( ) const {

		return m_Database.GetGenes( ); }

	float GetBackground( size_t iContext, size_t iGene ) const {

		return ( ( m_MatBackgrounds.GetColumns( ) && m_MatBackgrounds.GetRows( ) ) ?
			m_MatBackgrounds.Get( iContext, iGene ) : 1 ); }

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
	const CCompactFullMatrix&		m_MatContexts;
	const CDataMatrix&				m_MatBackgrounds;
	const CGenome&					m_Genome;
	const IOntology**				m_apOntologies;
	size_t							m_iOntologies;
	const vector<vector<size_t> >	m_vecveciDiseases;
};

#endif // BNSERVER_H

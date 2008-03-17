#ifndef BNSERVER_H
#define BNSERVER_H

class CDot;

class CBNServer : public IServerClient {
public:
	static bool Get( size_t, size_t, float*, const CDatabase&, const vector<CBayesNetMinimal>&,
		const CBayesNetMinimal& );

	CBNServer( const CBayesNetMinimal&, const vector<CBayesNetMinimal>&, const CCompactFullMatrix&, SOCKET,
		const CDatabase&, const string&, const char*, const char*, const CDataMatrix&, const CGenome&,
		const IOntology** );
	~CBNServer( );

	IServerClient* NewInstance( SOCKET, uint32_t, uint16_t );
	void Destroy( );
	bool ProcessMessage( const vector<unsigned char>& );

private:
	typedef size_t (CBNServer::*TPFNProcessor)( const vector<unsigned char>&, size_t );

	static const size_t			c_iDegree			= 1;
	static const TPFNProcessor	c_apfnProcessors[];
	static const size_t			c_iProcessors;
	static const float			c_dCutoff;
	static const float			c_adColorMin[];
	static const float			c_adColorMax[];

	bool Get( size_t, size_t, float* = NULL );
	bool Get( size_t, const vector<size_t>&, size_t, float* );
	bool Get( const vector<unsigned char>&, size_t );
	bool GraphCreate( const vector<size_t>&, size_t, size_t, vector<bool>&, vector<size_t>&, CDat& ) const;
	bool GraphWrite( const CDat&, const vector<size_t>&, const vector<size_t>&, const vector<bool>&,
		size_t, bool ) const;
	bool SelectNeighborsPixie( const vector<size_t>&, const vector<bool>&, size_t, size_t, const CDataMatrix&,
		vector<size_t>& ) const;
	bool SelectNeighborsRatio( const vector<size_t>&, const vector<bool>&, size_t, size_t, const CDataMatrix&,
		vector<size_t>& ) const;
	bool SendGenes( const vector<size_t>&, const vector<size_t>& ) const;
	size_t ProcessInference( const vector<unsigned char>&, size_t );
	size_t ProcessData( const vector<unsigned char>&, size_t );
	size_t ProcessGraph( const vector<unsigned char>&, size_t );
	size_t ProcessContexts( const vector<unsigned char>&, size_t );
	size_t ProcessTermFinder( const vector<unsigned char>&, size_t );

	size_t GetGenes( ) const {

		return m_Database.GetGenes( ); }

	const CBayesNetMinimal&			m_BNDefault;
	const vector<CBayesNetMinimal>&	m_vecBNs;
	SOCKET							m_iSocket;
	float*							m_adGenes;
	float*							m_adContexts;
	const CDatabase&				m_Database;
	string							m_strConnection;
	string							m_strGraphviz;
	string							m_strFiles;
	const CCompactFullMatrix&		m_MatContexts;
	const CDataMatrix&				m_MatBackgrounds;
	const CGenome&					m_Genome;
	const IOntology**				m_apOntologies;
	size_t							m_iOntologies;
};

#endif // BNSERVER_H

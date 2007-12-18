#ifndef BNSERVER_H
#define BNSERVER_H

class CDot;

class CBNServer : public IServerClient {
public:
	CBNServer( const CBayesNetMinimal&, const vector<CBayesNetMinimal>&, const CCompactFullMatrix&, SOCKET,
		const CDatabase&, const string&, const char*, const char* );
	~CBNServer( );

	IServerClient* NewInstance( SOCKET, uint32_t, uint16_t, const CPropertyFile* );
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
	bool PixieCreate( const vector<size_t>&, size_t, size_t, vector<bool>&, CDat& ) const;
	bool PixieGraph( const CDat&, const vector<bool>&, bool ) const;
	bool PixieGraphWrite( const CDat&, const vector<bool>&, const CDot&, ostream& ) const;
	size_t ProcessInference( const vector<unsigned char>&, size_t );
	size_t ProcessData( const vector<unsigned char>&, size_t );
	size_t ProcessGraph( const vector<unsigned char>&, size_t );
	size_t ProcessContexts( const vector<unsigned char>&, size_t );

	const CBayesNetMinimal&			m_BNDefault;
	const vector<CBayesNetMinimal>&	m_vecBNs;
	size_t							m_iGenes;
	SOCKET							m_iSocket;
	float*							m_adGenes;
	float*							m_adContexts;
	const CDatabase&				m_Database;
	string							m_strConnection;
	string							m_strGraphviz;
	string							m_strFiles;
	const CCompactFullMatrix&		m_MatContexts;
};

#endif // BNSERVER_H

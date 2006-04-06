#ifndef BAYESNETFNI_H
#define BAYESNETFNI_H

namespace libBioUtils {

class CBayesNetFNNode {
protected:
	friend class CBayesNetFN;
	friend class CBayesNetFNImpl;

	static const char	c_szType[];

	static CBayesNetFNNode* Open( DSL_node* );

	const std::string& GetName( ) const;
	unsigned char GetParameters( ) const;
	void Reverse( );
	bool Save( DSL_node* ) const;
	bool Learn( const std::vector<size_t>& );

	virtual const char* GetType( ) const = 0;
	virtual void Randomize( ) = 0;
	virtual CBayesNetFNNode* New( DSL_node* ) const = 0;
	virtual bool Learn( const IDataset*, size_t, bool ) = 0;
	virtual bool Evaluate( float, std::vector<float>& ) const = 0;

	virtual bool IsContinuous( ) const {

		return true; }

	std::string			m_strName;
	const char*			m_szType;
	CFullMatrix<float>	m_Params;
};

class CBayesNetFNNodeDiscrete : protected CBayesNetFNNode {
protected:
	friend class CBayesNetFNNode;

	void Randomize( );
	bool Learn( const IDataset*, size_t, bool );
	bool Evaluate( float, std::vector<float>& ) const;

	CBayesNetFNNode* New( DSL_node* pNode ) const {

		return new CBayesNetFNNodeDiscrete( ); }

	const char* GetType( ) const {

		return "discrete"; }

	bool IsContinuous( ) const {

		return false; }
};

class CBayesNetFNNodeGaussian : protected CBayesNetFNNode {
protected:
	friend class CBayesNetFNNode;

	static const size_t	c_iMu		= 0;
	static const size_t	c_iSigma	= 1;

	void Randomize( );
	bool Learn( const IDataset*, size_t, bool );
	bool Evaluate( float, std::vector<float>& ) const;

	CBayesNetFNNode* New( DSL_node* pNode ) const {

		return new CBayesNetFNNodeGaussian( ); }

	const char* GetType( ) const {

		return "gaussian"; }
};

class CBayesNetFNNodeBeta : protected CBayesNetFNNode {
protected:
	friend class CBayesNetFNNode;

	static const size_t	c_iMin		= 0;
	static const size_t	c_iMax		= 1;
	static const size_t	c_iAlpha	= 2;
	static const size_t	c_iBeta		= 3;

	void Randomize( );
	bool Learn( const IDataset*, size_t, bool );
	bool Evaluate( float, std::vector<float>& ) const;

	CBayesNetFNNode* New( DSL_node* pNode ) const {

		return new CBayesNetFNNodeBeta( ); }

	const char* GetType( ) const {

		return "beta"; }
};

class CBayesNetFNNodeExponential : protected CBayesNetFNNode {
protected:
	friend class CBayesNetFNNode;

	static const size_t	c_iMin	= 0;
	static const size_t	c_iBeta	= 1;

	void Randomize( );
	bool Learn( const IDataset*, size_t, bool );
	bool Evaluate( float, std::vector<float>& ) const;

	CBayesNetFNNode* New( DSL_node* pNode ) const {

		return new CBayesNetFNNodeExponential( ); }

	const char* GetType( ) const {

		return "exponential"; }
};

class CBayesNetFNImpl : protected CBayesNetImpl {
protected:
	CBayesNetFNImpl( );
	~CBayesNetFNImpl( );

	void Reset( );
	bool Evaluate( const IDataset*, CDat*, std::vector<std::vector<float> >*, bool ) const;
	bool Evaluate( const IDataset*, size_t, size_t, bool, std::vector<float>& ) const;

	size_t				m_iNodes;
	CBayesNetFNNode**	m_apNodes;
	bool				m_fSmileNet;
	DSL_network			m_SmileNet;
};

}

#endif // BAYESNETFNI_H

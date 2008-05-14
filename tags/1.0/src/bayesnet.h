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
#ifndef BAYESNET_H
#define BAYESNET_H

#include "bayesnetint.h"
#include "bayesneti.h"

namespace Sleipnir {

class IDataset;

#ifndef NO_SMILE

/*!
 * \brief
 * Implements IBayesNet for networks using the SMILE library from the U. Pittsburgh Decision Systems Lab.
 * 
 * CBayesNetSmile loads and saves Bayes nets from DSL/XDSL files and performs Bayesian inference using
 * the SMILE library from the University of Pittsburgh Decision Systems Laboratory.  While SMILE is used
 * for internal representation of the Bayes net and for inference in many cases, Sleipnir implements
 * several optimizations.  Networks detected to have naive structures are learned and evaluated using
 * more efficient maximum likelihood methods, and Sleipnir implements its own EM and ELR parameter
 * learning algorithms.  Naive SMILE-based Bayes nets can be converted to extremely efficient
 * CBayesNetMinimal objects, and if the PNL library is present, they can also be converted to CBayesNetPNL
 * representations.
 * 
 * \remarks
 * Should minimally support any network type allowed by SMILE; only tested using discrete networks with
 * hierarchical structure.  Default values for individual nodes are stored in the "zero" property of the
 * SMILE network (can be visualized in the "User Properties" pane in GeNIe); for example, to provide a
 * default value of 2 for some node, give it a user property named "zero" with value "2".
 */
class CBayesNetSmile : public CBayesNetSmileImpl, public IBayesNet {
public:
	CBayesNetSmile( bool fGroup = true );

	bool Open( const std::vector<std::string>& vecstrFiles, size_t iValues );
#ifdef PNL_ENABLED
	bool Convert( CBayesNetPNL& BNPNL ) const;
#endif // PNL_ENABLED
	bool Open( const IDataset* pDataset, const std::vector<std::string>& vecstrNames,
		const std::vector<size_t>& veciDefaults );
	bool Open( const CBayesNetSmile& BNPrior, const std::vector<CBayesNetSmile*>& vecpBNs );
	bool Open( const CBayesNetMinimal& BNMinimal, const std::vector<std::string>& vecstrNames );
	float Evaluate( size_t iNode, unsigned char bValue ) const;
	unsigned char GetDefault( size_t iNode ) const;

	/*!
	 * \brief
	 * Provide a Bayes net of identical structure from which default parameter values can be obtained.
	 * 
	 * \param Defaults
	 * Bayes net with identical structure whose parameters will be used when insufficient data is available
	 * during parameter learning.
	 * 
	 * If a default network is provided, the IBayesNet::Learn method will use that network's probability
	 * distribution for any parameter column in which fewer than CBayesNetSmileImpl::c_iMinimum examples are
	 * present.  This prevents error introduced by methods such as Laplace smoothing when too few
	 * examples are present to estimate a reasonable maximum likelihood probability distribution.
	 */
	void SetDefault( const CBayesNetSmile& Defaults ) {

		m_pDefaults = &Defaults; }

	bool Learn( const IDataset* pDataset, size_t iIterations, bool fZero = false, bool fELR = false );
	bool Evaluate( const std::vector<unsigned char>& vecbDatum, std::vector<float>& vecdResults,
		bool fZero = false, size_t iNode = 0, bool fIgnoreMissing = false ) const;
	bool Evaluate( const CPCLPair& PCLData, CPCL& PCLResults, bool fZero = false,
		int iAlgorithm = DSL_ALG_BN_LAURITZEN ) const;
	void GetNodes( std::vector<std::string>& vecstrNodes ) const;
	void Randomize( );
	void Randomize( size_t iNode );
	void Reverse( size_t iNode );

	bool Open( const char* szFile ) {

		return ( m_fSmileNet = !m_SmileNet.ReadFile( szFile ) ); }

	bool Save( const char* szFile ) const {

		return ( m_fSmileNet ? !((CBayesNetSmile*)this)->m_SmileNet.WriteFile( szFile ) : false ); }

	bool GetCPT( size_t iNode, CDataMatrix& MatCPT ) const {

		return CBayesNetSmileImpl::GetCPT( m_SmileNet.GetNode( (int)iNode ), MatCPT ); }

	unsigned char GetValues( size_t iNode ) const {

		return m_SmileNet.GetNode( (int)iNode )->Definition( )->GetNumberOfOutcomes( ); }

	bool IsContinuous( size_t iNode ) const {
		UNUSED_PARAMETER(iNode);

		return IsContinuous( ); }

	bool IsContinuous( ) const {

		return CBayesNetSmileImpl::IsContinuous( ); }

	bool Evaluate( const IDataset* pDataset, std::vector<std::vector<float> >& vecvecdResults,
		bool fZero ) const {

		return CBayesNetSmileImpl::Evaluate( pDataset, NULL, &vecvecdResults, fZero ); }

	bool Evaluate( const IDataset* pDataset, CDat& DatResults, bool fZero ) const {

		return CBayesNetSmileImpl::Evaluate( pDataset, &DatResults, NULL, fZero ); }
};

/*!
 * \brief
 * Implements IBayesNet for networks using custom node types.
 * 
 * CBayesNetFN can be used to construct Bayes nets using arbitrary node types.  These are usually only
 * theoretically sound in a naive structure, but in such a case, any node type can be used for which parameters
 * can be estimated from data: discrete, Gaussian, Beta, Exponential, etc.  These networks are stored using
 * a SMILE network in a DSL/XDSL file, but the semantics of each node's parameters are dependent on the node
 * type.
 * 
 * \remarks
 * Node types are stored in the properties for each node (accessible through the "User Properties" pane in
 * GeNIe).  To determine node type, CBayesNetFN checks the "type" property of each node; if no such property
 * is available, the node is a standard discrete distribution (can also be given the value "discrete").
 * Other allowable values include "gaussian", "beta", "exponential", and "mog".
 * 
 * \see
 * CBayesNetFNNode
 */
class CBayesNetFN : CBayesNetFNImpl, public IBayesNet {
public:
	bool Open( const char* szFile );
	bool Save( const char* szFile ) const;
	bool Learn( const IDataset* pDataset, size_t iIterations, bool fZero = false, bool fELR = false );
	bool Evaluate( const std::vector<unsigned char>& vecbDatum, std::vector<float>& vecdResults,
		bool fZero = false, size_t iNode = 0, bool fIgnoreMissing = false ) const;
	void GetNodes( std::vector<std::string>& vecstrNodes ) const;
	unsigned char GetValues( size_t iNode ) const;
	bool IsContinuous( ) const;

	bool Evaluate( const IDataset* pDataset, std::vector<std::vector<float> >& vecvecdResults,
		bool fZero ) const {

		return CBayesNetFNImpl::Evaluate( pDataset, NULL, &vecvecdResults, fZero ); }

	bool Evaluate( const IDataset* pDataset, CDat& DatResults, bool fZero ) const {

		return CBayesNetFNImpl::Evaluate( pDataset, &DatResults, NULL, fZero ); }

	bool IsContinuous( size_t iNode ) const {

		return m_apNodes[ iNode ]->IsContinuous( ); }

	void Randomize( ) {
		size_t	i;

		for( i = 0; i < m_iNodes; ++i )
			Randomize( i ); }

	void Randomize( size_t iNode ) {

		m_apNodes[ iNode ]->Randomize( ); }

	void Reverse( size_t iNode ) {

		m_apNodes[ iNode ]->Reverse( ); }

	bool GetCPT( size_t iNode, CDataMatrix& MatCPT ) const {

		return CBayesNetSmileImpl::GetCPT( m_SmileNet.GetNode( (int)iNode ), MatCPT ); }

	bool Evaluate( const CPCLPair& PCLData, CPCL& PCLResults, bool fZero, int iAlgorithm ) const {
		UNUSED_PARAMETER(PCLData);
		UNUSED_PARAMETER(PCLResults);
		UNUSED_PARAMETER(fZero);
		UNUSED_PARAMETER(iAlgorithm);

		return false; }
};

#endif // NO_SMILE

#ifdef PNL_ENABLED

/*!
 * \brief
 * Implements IBayesNet for networks using the PNL library from Intel.
 * 
 * CBayesNetPNL loads and saves Bayes nets from XML files and performs Bayesian inference using the PNL
 * library from Intel.  PNL functionality is disabled by default in Sleipnir, since the PNL library is
 * very large, difficult to build, and much slower than SMILE for most tasks.  It has the benefit of
 * supporting more flexible node types, particularly continuous nodes.  Most methods not immediately
 * related to Bayesian learning and inference are not implemented in this class.
 * 
 * \remarks
 * Due to PNL's finickiness, this class has been minimally tested.  User beware!
 */
class CBayesNetPNL : public CBayesNetPNLImpl, public IBayesNet {
public:
	CBayesNetPNL( bool fGroup = true );

	bool Open( const char* szFile );
	bool Save( const char* szFile ) const;
	bool Learn( const IDataset* pDataset, size_t iIterations, bool fZero = false, bool fELR = false );

	void GetNodes( std::vector<std::string>& vecstrNodes ) const {
		UNUSED_PARAMETER(vecstrNodes); }

	bool IsContinuous( size_t iNode ) const {
		UNUSED_PARAMETER(iNode);

		return IsContinuous( ); }

	bool IsContinuous( ) const {

		return CBayesNetPNLImpl::IsContinuous( ); }

	bool Evaluate( const IDataset* pDataset, std::vector<std::vector<float> >& vecvecdResults,
		bool fZero ) const {

		return CBayesNetPNLImpl::Evaluate( pDataset, NULL, &vecvecdResults, fZero ); }

	bool Evaluate( const IDataset* pDataset, CDat& DatResults, bool fZero ) const {

		return CBayesNetPNLImpl::Evaluate( pDataset, &DatResults, NULL, fZero ); }

	void Randomize( ) { }

	void Randomize( size_t iNode ) {
		UNUSED_PARAMETER(iNode); }

	void Reverse( size_t iNode ) {
		UNUSED_PARAMETER(iNode); }

	virtual bool Evaluate( const std::vector<unsigned char>& vecbDatum, std::vector<float>& vecdResults,
		bool fZero = false, size_t iNode = 0, bool fIgnoreMissing = false ) const {

		return false; }

	unsigned char GetValues( size_t iNode ) const {
		UNUSED_PARAMETER(iNode);

		return 0; }

	bool GetCPT( size_t iNode, CDataMatrix& MatCPT ) const {
		UNUSED_PARAMETER(iNode);
		UNUSED_PARAMETER(MatCPT);

		return false; }

	bool Evaluate( const CPCLPair& PCLData, CPCL& PCLResults, bool fZero, int iAlgorithm ) const {
		UNUSED_PARAMETER(PCLData);
		UNUSED_PARAMETER(PCLResults);
		UNUSED_PARAMETER(fZero);
		UNUSED_PARAMETER(iAlgorithm);

		return false; }
};

#endif // PNL_ENABLED

/*!
 * \brief
 * Implements a heavily optimized discrete naive Bayesian classifier.
 * 
 * CBayesNetMinimal provides a custom implementation of a discrete naive Bayesian classifier heavily
 * optimized for rapid inference.  The intended use is to learn an appropriate network and parameters offline
 * using one of the more complex Bayes net implementations.  The resulting network can then be converted to a
 * minimal form and used for online (realtime) inference.  A minimal Bayes net always consists of one output
 * (class) node and zero or more data nodes, all discrete and taking one or more different values.
 */
class CBayesNetMinimal : public CBayesNetMinimalImpl {
public:
#ifndef NO_SMILE
	bool Open( const CBayesNetSmile& BNSmile );
#endif // NO_SMILE
	bool Open( std::istream& istm );
	bool OpenCounts( const char* szFileCounts, const std::map<std::string, size_t>& mapstriNodes,
		const std::vector<unsigned char>& vecbDefaults, const std::vector<float>& vecdAlphas,
		size_t iPseudocounts = -1, const CBayesNetMinimal* pBNDefault = NULL );
	void Save( std::ostream& ostm ) const;
	float Evaluate( const std::vector<unsigned char>& vecbDatum, size_t iOffset = 0 ) const;
	bool Evaluate( const std::vector<unsigned char>& vecbData, float* adResults, size_t iGenes,
		size_t iStart = 0 ) const;

	/*!
	 * \brief
	 * Return the conditional probability table matrix for the indicated node.
	 * 
	 * \param iNode
	 * Index of node whose CPT is returned (zero-based).
	 * 
	 * \returns
	 * CPT of requested node.
	 * 
	 * \remarks
	 * Classes are row-oriented, data values are column-oriented.  The root node is index zero, remaining
	 * data nodes begin at index one.
	 */
	const CDataMatrix& GetCPT( size_t iNode ) const {

		return ( iNode ? m_vecNodes[ iNode - 1 ].m_MatCPT : m_MatRoot ); }

	/*!
	 * \brief
	 * Return the total number of nodes in the Bayes net.
	 * 
	 * \returns
	 * Number of nodes in the Bayes net (including root and data nodes).
	 * 
	 * \remarks
	 * Includes root/class node an non-root/data nodes.
	 */
	size_t GetNodes( ) const {

		return ( m_vecNodes.size( ) + 1 ); }

	/*!
	 * \brief
	 * Sets the string identifier of the network.
	 * 
	 * \param strID
	 * String identifier for the network.
	 * 
	 * \remarks
	 * ID is not used internally and is purely for human convenience.
	 */
	void SetID( const std::string& strID ) {

		m_strID = strID; }

	/*!
	 * \brief
	 * Returns the string identifier of the network.
	 * 
	 * \returns
	 * String identifier for the network.
	 * 
	 * \remarks
	 * ID is not used internally and is purely for human convenience.
	 */
	const std::string& GetID( ) const {

		return m_strID; }

	/*!
	 * \brief
	 * Returns the default value (if no input is provided) for the requested node.
	 * 
	 * \param iNode
	 * Node for which default value is returned.
	 * 
	 * \returns
	 * Default value for the requested node (-1 if none).
	 * 
	 * \remarks
	 * Requested node index must be less than the number of nodes (beginning with the root node at index 0).
	 */
	const unsigned char GetDefault( size_t iNode ) const {

		return ( iNode ? m_vecNodes[ iNode - 1 ].m_bDefault : 0xFF ); }
};

}

#endif // BAYESNET_H

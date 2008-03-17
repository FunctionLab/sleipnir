#ifndef CLUSHIERARCHICAL_H
#define CLUSHIERARCHICAL_H

#include "clusthierarchicali.h"

namespace Sleipnir {

class CPCL;
class IMeasure;

/*!
 * \brief
 * Represents a simple node in a binary tree.
 * 
 * Generated by CClustHierarchical::Cluster, a CHierarchy is an extremely rudimentary representation of an
 * binary tree intended to be serialized to disk as the GTR file in a CDT/GTR pair.  Each node either zero
 * or two children, a unique integer identifier within the tree, and a similarity score indicating its height
 * within the tree.
 */
class CHierarchy : public CHierarchyImpl {
public:
	CHierarchy( size_t iID, float dSimilarity, const CHierarchy* pLeft, const CHierarchy* pRight );

	void GetGenes( std::vector<size_t>& veciGenes ) const;
	float SortChildren( const std::vector<float>& vecdScores );

	/*!
	 * \brief
	 * Save the hierarchy to the given stream in GTR format.
	 * 
	 * \param ostm
	 * Output stream into which the hierarchy is saved.
	 * 
	 * \param iGenes
	 * Total number of leaf nodes in the hierarchy.
	 * 
	 * \remarks
	 * iGenes can be calculated from the hierarchy; it is included as an input solely for convenience
	 * purposes, since the genes must be output in original order (not traversal order) to satisfy GTR
	 * file formatting requirements.
	 */
	void Save( std::ostream& ostm, size_t iGenes ) const {
		size_t	i;

		for( i = 0; ( i + 1 ) < iGenes; ++i )
			CHierarchyImpl::Save( ostm, i ); }

	/*!
	 * \brief
	 * Safety method to delete a hierarchy.
	 * 
	 * \remarks
	 * Included to avoid the necessity of directly deleting something allocated within a library method.
	 */
	void Destroy( ) {

		delete this; }

	/*!
	 * \brief
	 * Returns this node's height within the hierarchy.
	 * 
	 * \returns
	 * The current node's height within the hierarchy.
	 */
	float GetSimilarity( ) const {

		return m_dScore; }

	/*!
	 * \brief
	 * Returns true if the current node is a leaf node (i.e. represents a gene in the hierarchy).
	 * 
	 * \returns
	 * True if the current node is a leaf (has no children).
	 */
	bool IsGene( ) const {

		return CHierarchyImpl::IsGene( ); }

	/*!
	 * \brief
	 * Returns the current node's unique ID within the hierarchy.
	 * 
	 * \returns
	 * The current node's ID.
	 * 
	 * \remarks
	 * Leaf node IDs generally correspond to gene indices within the pre-clustered PCL; internal node IDs are
	 * arbitrary unique values.
	 */
	size_t GetID( ) const {

		return m_iID; }

	/*!
	 * \brief
	 * Returns the current node's left or right child.
	 * 
	 * \param fRight
	 * If true, return the right (second) child; otherwise, return the left (first).
	 * 
	 * \returns
	 * One of the current node's two children.
	 * 
	 * \remarks
	 * Do not call for leaf nodes.
	 * 
	 * \see
	 * IsLeaf
	 */
	const CHierarchy& Get( bool fRight ) const {

		return *( fRight ? m_pRight : m_pLeft ); }
};

/*!
 * \brief
 * Utility class containing static hierarchical clustering methods.
 */
class CClustHierarchical : CClustHierarchicalImpl {
public:
	static CHierarchy* Cluster( const CDistanceMatrix& MatSimilarities );
	static CHierarchy* Cluster( const CDistanceMatrix& MatSimilarities,
		const std::vector<bool>& vecfIncluded );
};

}

#endif // CLUSHIERARCHICAL_H

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
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, and Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#ifndef DOT_H
#define DOT_H

class CDot {
public:
	CDot( const CDat& Dat ) : m_Dat(Dat) {

		m_pProperties = new boost::dynamic_properties( boost::ignore_other_properties ); }

	~CDot( ) {

		delete m_pProperties; }

	bool Open( const char* );
	bool Save( ostream&, const std::vector<bool>&, size_t ) const;

private:
	typedef boost::property<boost::vertex_name_t, string, boost::property<boost::vertex_index1_t, string,
		boost::property<boost::vertex_index2_t, string, boost::property<boost::vertex_attribute_t, string> > > >
																				TAttributesVertex;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
		TAttributesVertex, boost::property<boost::edge_attribute_t, string> >
																				TGraph;
	typedef TGraph::edge_descriptor												TEdge;
	typedef TGraph::vertex_descriptor											TVertex;

	static const char	c_szCutoffBox[];
	static const char	c_szHeader00[];
	static const char	c_szHeader01[];
	static const char	c_szHeader02[];
	static const float	c_dEdgeOpacity;
	static const float	c_dScale;

	bool SaveEdge( ostream&, const TEdge& ) const;
	bool SaveVertex( ostream&, const TVertex&, const std::vector<bool>& ) const;

	const CDat&					m_Dat;
	TGraph						m_Graph;
	boost::dynamic_properties*	m_pProperties;
};

#endif // DOT_H

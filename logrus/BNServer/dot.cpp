#include "stdafx.h"
#include "dot.h"

using namespace boost;

const float	CDot::c_dEdgeOpacity	= 0.13f;
const char	CDot::c_szHeader[]		=
	"<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n"
	"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n"
	"<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n"
//	"	xmlns:a3=\"http://ns.adobe.com/AdobeSVGViewerExtensions/3.0/\"\n"
//	"	a3:scriptImplementation=\"Adobe\"\n"
//	"	onload=\"init( )\"\n"
//	"	viewBox=\"10 10 256 256"\"\n"
	"	>\n"
// If we used scripts for the viewbox, they'd go here.
	"	<g id=\"graph\" class=\"graph\" style=\"font-family:Times-Roman;font-size:12pt\">\n"
	"		<title>G</title>\n";

bool CDot::Open( const char* szDot ) {
	ifstream	ifsm;

	m_pProperties->property( "node_id", get( vertex_name, m_Graph ) );
	m_pProperties->property( "pos", get( vertex_attribute, m_Graph ) );
	m_pProperties->property( "width", get( vertex_index1, m_Graph ) );
	m_pProperties->property( "height", get( vertex_index2, m_Graph ) );
	m_pProperties->property( "pos", get( edge_attribute, m_Graph ) );
	ifsm.open( szDot );
	return ( ifsm.is_open( ) && read_graphviz( ifsm, m_Graph, *m_pProperties ) ); }

bool CDot::Save( ostream& ostm, const vector<bool>& vecfQuery ) const {
	graph_traits<TGraph>::edge_iterator		iterEdge, iterEdgeEnd;
	graph_traits<TGraph>::vertex_iterator	iterVertex, iterVertexEnd;

	ostm << c_szHeader;

	ostm << "		<g id=\"edges\" class=\"edgeList\">" << endl;
	for( tie( iterEdge, iterEdgeEnd ) = edges( m_Graph ); iterEdge != iterEdgeEnd; ++iterEdge )
		if( !SaveEdge( ostm, *iterEdge ) )
			return false;
	ostm << "		</g>" << endl;

	ostm << "		<g id=\"nodes\" class=\"nodeList\">" << endl;
	for( tie( iterVertex, iterVertexEnd ) = vertices( m_Graph ); iterVertex != iterVertexEnd; ++iterVertex )
		if( !SaveVertex( ostm, *iterVertex, vecfQuery ) )
			return false;
	ostm << "		</g>" << endl;

	ostm << "	</g>" << endl << "</svg>" << endl;

	return true; }

bool CDot::SaveEdge( ostream& ostm, const TEdge& Edge ) const {
	static const char	c_acTabs[]	= "			";
	static const size_t	c_iBuffer	= 16;
	char			acBuffer[ c_iBuffer ];
	size_t			i, iHead, iTail;
	float			dWeight;
	string			strPath;
	vector<string>	vecstrPath;

	iHead = source( Edge, m_Graph );
	iTail = target( Edge, m_Graph );
	if( iTail < iHead ) {
		i = iHead;
		iHead = iTail;
		iTail = i; }
	if( CMeta::IsNaN( dWeight = m_Dat.Get( iHead, iTail ) ) )
		return false;
	sprintf_s( acBuffer, "%d_%d", iHead, iTail );

	strPath = get( "pos", *m_pProperties, Edge );
	CMeta::Tokenize( strPath.c_str( ), vecstrPath, " " );
	if( vecstrPath.size( ) < 1 )
		return false;
	strPath = "M";
	strPath += vecstrPath[ 0 ];
	for( i = 1; i < vecstrPath.size( ); ++i )
		strPath += ( ( i == 1 ) ? "C" : " " ) + vecstrPath[ i ];

	ostm <<
		c_acTabs << "<g id=\"edge" << acBuffer << "\" class=\"edge\" stroke-opacity=\"" <<
			c_dEdgeOpacity << "\">" << endl <<
		c_acTabs << "	<title>" << iHead << "--" << iTail << "--" << dWeight << "</title>" << endl <<
		c_acTabs << "	<text display=\"none\" id=\"edge" << acBuffer << "_confidence\">" << dWeight <<
			"</text>" << endl <<
		c_acTabs << "	<g id=\"edge" << acBuffer << "_nodes\" display=\"none\">" << endl <<
		c_acTabs << "		<text>node" << iHead << "</text>" << endl <<
		c_acTabs << "		<text>node" << iTail << "</text>" << endl <<
		c_acTabs << "	</g>" << endl <<
		c_acTabs << "	<path style=\"fill:none;stroke:#" << CDat::GetColor( dWeight ) <<
			";stroke-width:3;\" d=\"" + strPath + "\" />" << endl <<
		c_acTabs << "</g>" << endl;

	return true; }

bool CDot::SaveVertex( ostream& ostm, const TVertex& Vertex, const vector<bool>& vecfQuery ) const {
	static const float	c_dScale	= 36;
	static const char	c_acTabs[]	= "			";
	static const size_t	c_iBuffer	= 16;
	char									acBuffer[ c_iBuffer ];
	graph_traits<TGraph>::out_edge_iterator	iterEdge, iterEdgeEnd;
	string									strID, strPos, strCX, strCY, strRX, strRY;
	string									strForeground, strBackground;
	size_t									i;
	float									dWidth, dHeight;
	vector<string>							vecstrPos;

	strID = get( "node_id", *m_pProperties, Vertex );
	dWidth = (float)atof( get( "width", *m_pProperties, Vertex ).c_str( ) );
	dHeight = (float)atof( get( "height", *m_pProperties, Vertex ).c_str( ) );
	strPos = get( "pos", *m_pProperties, Vertex );
	CMeta::Tokenize( strPos.c_str( ), vecstrPos, "," );
	if( vecstrPos.size( ) != 2 )
		return false;
	strCX = vecstrPos[ 0 ];
	strCY = vecstrPos[ 1 ];
	sprintf_s( acBuffer, "%d", (size_t)( c_dScale * dWidth ) );
	strRX = acBuffer;
	sprintf_s( acBuffer, "%d", (size_t)( c_dScale * dHeight ) );
	strRY = acBuffer;

	ostm <<
		c_acTabs << "<g id=\"node" + strID + "\" class=\"node\">" << endl <<
		c_acTabs << "	<title>" + m_Dat.GetGene( Vertex ) + "</title>" << endl <<
		c_acTabs << "	<g id=\"node" + strID + "_edges\" display=\"none\">" << endl;
	for( tie( iterEdge, iterEdgeEnd ) = out_edges( Vertex, m_Graph ); iterEdge != iterEdgeEnd; ++iterEdge ) {
		size_t	iHead, iTail;

		iHead = source( *iterEdge, m_Graph );
		iTail = target( *iterEdge, m_Graph );
		if( iTail < iHead ) {
			i = iHead;
			iHead = iTail;
			iTail = i; }
		ostm << c_acTabs << "		<text>edge" << iHead << '_' << iTail << "</text>" << endl; }
	ostm << c_acTabs << "	</g>" << endl;

	strForeground = "black";
	strBackground = vecfQuery[ Vertex ] ? "lightgrey" : "white";
	ostm <<
		c_acTabs << "	<ellipse cx=\"" << strCX << "\" cy=\"" << strCY << "\" rx=\"" << strRX <<
			"\" ry=\"" << strRY << "\" style=\"fill:" << strBackground << ";stroke:black\" />" << endl <<
		c_acTabs << "	<text id=\"node" << strID << "_title\" text-anchor=\"middle\" style=\"fill:" <<
			strForeground << ";stroke:" << strForeground << "\" x=\"" << strCX << "\" y=\"" <<
			( atof( strCY.c_str( ) ) + 5 ) << "\">" << m_Dat.GetGene( Vertex ) << "</text>" << endl <<
		c_acTabs << "</g>" << endl;

	return true; }

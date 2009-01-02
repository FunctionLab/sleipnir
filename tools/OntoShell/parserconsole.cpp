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
#include "stdafx.h"
#include "parserconsole.h"

const char*					CParserConsole::SArgs::c_aszFlags[]		=
	{ CParserConsole::c_szGenes, CParserConsole::c_szLong, CParserConsole::c_szSibs,
	CParserConsole::c_szZeroes, CParserConsole::c_szBonferroni,
	CParserConsole::c_szRecursive, CParserConsole::c_szBackground, NULL };
const char					CParserConsole::c_szDotDotDot[]			= "...";
const char					CParserConsole::c_szBackground[]		= "-k";
const char					CParserConsole::c_szBonferroni[]		= "-b";
const char					CParserConsole::c_szGenes[]				= "-g";
const char					CParserConsole::c_szLong[]				= "-l";
const char					CParserConsole::c_szSibs[]				= "-s";
const char					CParserConsole::c_szZeroes[]			= "-a";
const char					CParserConsole::c_szRecursive[]			= "-r";
const char					CParserConsole::c_szStar[]				= "*";
const char					CParserConsole::c_szHelpHelp[]			= "Commands:\n"
	"cat <gene>+                Displays information on individual genes.\n"
	"cd [path]                  Display or change current term.\n"
	"find <filename> [p] [bkg]  Runs term finder on the given gene list.\n"
	"help [command]             Provides help on command syntax.\n"
	"ls [path]                  List parents, children, and annotations.\n"
	"parentage <onto> <file>    For terms in onto, list parents in the given set.";
const CParserConsole::TPFnParser	CParserConsole::c_apfnParsers[]	=
	{ &CParserConsole::ParseCat, &CParserConsole::ParseCd, &CParserConsole::ParseFind,
	&CParserConsole::ParseHelp, &CParserConsole::ParseLs, &CParserConsole::ParseParentage, NULL };
const char*					CParserConsole::c_aszHelps[]			= {
	"cat [-l] [-r] [path]<gene>+\n\n"
	"Displays the name, synonyms, and annotations for the given gene(s).\n"
	"Annotations are listed per ontology as available, with ontology term glosses\n"
	"abbreviated unless the -l flag is given.  A * will list all genes annotated\n"
	"to the given location or, in combination with the -r flag, to it and its\n"
	"descendants.",
	"cd [path]\n\n"
	"With no argument, cd displays the current path - either an ontology term, an\n"
	"ontology name, or the root marker.  When given a path, cd changes the current\n"
	"term to that path's target.  As in DOS and Unix, paths can contain . and ..\n"
	"characters to indicate the current term or its parent.  In nodes with\n"
	"multiple parents, the parent ID must be specified explicitly.  / serves as\n"
	"the path separator and root marker.",
	"find <filename> [p=0.05] [bkg] [-l] [-b] [-g] [-a] [-s] [-k]\n\n"
	"Performs a hypergeometric test over each ontology using the given gene list.\n"
	"Only terms with probability less than p are displayed, with a default p-value\n"
	"of 0.05.  The total number of possible genes is assumed to be the entire\n"
	"genome.  If given, the second file is used as a background distribution.  The\n"
	"other optional flags are:\n"
	"-l  Long listings; deactivates term gloss abbreviation.\n"
	"-b  Bonferroni correction; deactivates Bonferroni correction.\n"
	"-g  Genes; display genes associated with each ontology term.\n"
	"-a  All listings; display additional information for ontology terms.\n"
	"-s  Siblings; deactivates child annotations during analysis.\n"
	"-k  Background; uses whole genome background in place of ontology background.",
	CParserConsole::c_szHelpHelp,
	"ls [-l] [-a] [-g] [-s] [-r] [path]\n\n"
	"With no arguments, the ls command displays the parents, children, and gene\n"
	"annotations of the current term.  Given a path, it displays the same\n"
	"information for that target instead.  The four optional flags are:\n"
	"-l  Long listings; deactivates term gloss abbreviation.\n"
	"-a  All listings; includes terms with zero gene annotations.\n"
	"-g  Genes; deactives gene listings.\n"
	"-s  Siblings; deactivates parent and child listings.\n"
	"-r  Recursive; descend into child nodes.",
	"parentage [-a] <ontology> <filename>\n\n"
	"Loads an ontology slim from the given filename.  Then, for each term in\n"
	"the indicated ontology, outputs the zero or more parents of that term that\n"
	"fall within the given set.  This \"bubbles up\" the ontology to the level\n"
	"given in the input slim file.  Optional flags are:\n"
	"-a  All listings; include terms with no parents in the slim.",
	NULL };

CParserConsole::SArgs::SArgs( ) : m_fGenes(m_afFlags[ 0 ]), m_fLong(m_afFlags[ 1 ]),
	m_fSibs(m_afFlags[ 2 ]), m_fZeroes(m_afFlags[ 3 ]), m_fBonferroni(m_afFlags[ 4 ]),
	m_fRecursive(m_afFlags[ 5 ]), m_fBackground(m_afFlags[ 6 ]) {

	m_fGenes = false;
	m_fLong = false;
	m_fSibs = true;
	m_fZeroes = false;
	m_fBonferroni = true;
	m_fRecursive = false;
	m_fBackground = false; }

bool CParserConsole::SArgs::Parse( const string& strArg ) {
	size_t	i;

	for( i = 0; SArgs::c_aszFlags[ i ]; ++i )
		if( strArg == SArgs::c_aszFlags[ i ] ) {
			m_afFlags[ i ] = !m_afFlags[ i ];
			return true; }

	return false; }

void CParserConsole::PrintLink( const IOntology* pOnto, size_t iNode, char cType,
	const SArgs& sArgs ) {
	size_t	iGenes, iCount;
	string	strID, strGloss;

	if( !( iCount = pOnto->GetGenes( iNode, true ) ) && !sArgs.m_fZeroes )
		return;

	cout << cType << ' ' << ( strID = pOnto->GetID( iNode ) );
	PrintSpaces( c_iWidthID - strID.size( ) );
	iGenes = pOnto->GetGenes( iNode );
	PrintNumber( iGenes, c_iWidthGenes );
	PrintNumber( iCount - iGenes, c_iWidthGenes );
	PrintGloss( pOnto->GetGloss( iNode ), c_iWidthGloss, sArgs.m_fLong );
	cout << endl; }

void CParserConsole::PrintNumber( size_t iNumber, size_t iWidth ) {
	size_t	iUsed;

	cout << (unsigned int)iNumber;
	iUsed = iNumber ? (size_t)log10( (float)iNumber ) : 0;
	PrintSpaces( iWidth - iUsed ); }

void CParserConsole::PrintSpaces( size_t iSpaces ) {
	size_t	i;

	for( i = 0; i < iSpaces; ++i )
		cout << ' '; }

void CParserConsole::PrintAnnotation( const IOntology* pOnto, size_t iNode,
	const SArgs& sArgs, const STermFound* psFound ) {
	char	szBuf[ 128 ];
	string	strID;
	size_t	iWidth;

	iWidth = c_iWidthGloss + c_iWidthGenes;
	cout << ( strID = pOnto->GetID( iNode ) );
	PrintSpaces( c_iWidthID - strID.size( ) );
	if( psFound ) {
		iWidth -= c_iWidthGenes;
		sprintf_s( szBuf, "%g", psFound->m_dP );
		cout << szBuf;
		PrintSpaces( c_iWidthP - strlen( szBuf ) );
		sprintf_s( szBuf, "%-4d %-4d %-4d %-4d ", psFound->m_iHitsTerm, psFound->m_iSizeTerm,
			psFound->m_iHitsTotal, psFound->m_iSizeTotal );
		iWidth -= strlen( szBuf );
		cout << szBuf; }
	PrintGloss( pOnto->GetGloss( iNode ), iWidth, sArgs.m_fLong );
	cout << endl; }

void CParserConsole::PrintGloss( string strGloss, size_t iWidth, bool fLong ) {

	if( ( strGloss.length( ) > iWidth ) && !fLong ) {
		strGloss.resize( iWidth );
		strGloss += c_szDotDotDot; }
	cout << strGloss; }

void CParserConsole::PrintGene( const CGene& Gene, const SArgs& sArgs ) {
	size_t				i, j;
	const IOntology*	pOnto;

	cout << Gene.GetName( );
	if( Gene.GetSynonyms( ) ) {
		cout << " (" << Gene.GetSynonym( 0 );
		for( i = 1; i < Gene.GetSynonyms( ); ++i )
			cout << ',' << Gene.GetSynonym( i );
		cout << ')'; }
	if( Gene.GetDubious( ) )
		cout << " Dubious";
	if( Gene.GetRNA( ) )
		cout << " RNA";
	cout << endl;
	if( Gene.GetGloss( ).length( ) )
		cout << Gene.GetGloss( ) << endl;
	for( i = 0; i < Gene.GetOntologies( ); ++i ) {
		pOnto = Gene.GetOntology( i );
		cout << pOnto->GetID( ) << ':' << endl;
		PrintAnnotation( pOnto, Gene.GetAnnotation( i, 0 ), sArgs );
		for( j = 1; j < Gene.GetAnnotations( i ); ++j )
			PrintAnnotation( pOnto, Gene.GetAnnotation( i, j ), sArgs ); } }

void CParserConsole::PrintGenes( const vector<const CGene*>& vecpGenes, size_t iWidth, const CGenes* pGenes ) {
	size_t			i, iCol, iCols, iSpaces;
	vector<string>	vecstrGenes;

	iSpaces = 1;
	i = FormatGenes( vecpGenes, vecstrGenes, pGenes );
	if( !iWidth )
		iWidth = i;
	iCols = ( iWidth >= c_iWidthScreen ) ? 1 : ( c_iWidthScreen / iWidth );
	for( iCol = i = 0; i < vecpGenes.size( ); ++i,iCol %= iCols ) {
		PrintSpaces( iSpaces );
		cout << vecstrGenes[ i ];
		if( ++iCol == iCols ) {
			iSpaces = 1;
			cout << endl; }
		else
			iSpaces = iWidth - vecstrGenes[ i ].length( ); }
	if( iCol )
		cout << endl; }

size_t CParserConsole::FormatGenes( const vector<const CGene*>& vecpGenes,
	vector<string>& vecstrGenes, const CGenes* pGenes ) {
	size_t	i, j, iRet;

	vecstrGenes.resize( vecpGenes.size( ) );
	for( iRet = i = 0; i < vecpGenes.size( ); ++i ) {
		vecstrGenes[ i ] = ( pGenes && pGenes->IsGene( vecpGenes[ i ]->GetName( ) ) ) ? "*" : "";
		vecstrGenes[ i ] += vecpGenes[ i ]->GetName( );
		if( vecpGenes[ i ]->GetSynonyms( ) ) {
			vecstrGenes[ i ] += "(" + vecpGenes[ i ]->GetSynonym( 0 );
			for( j = 1; j < vecpGenes[ i ]->GetSynonyms( ); ++j )
				vecstrGenes[ i ] += "," + vecpGenes[ i ]->GetSynonym( j );
			vecstrGenes[ i ] += ")"; }
		if( vecpGenes[ i ]->GetRNA( ) )
			vecstrGenes[ i ] += "'";
		if( vecpGenes[ i ]->GetDubious( ) )
			vecstrGenes[ i ] += "!";
		if( vecstrGenes[ i ].length( ) > iRet )
			iRet = vecstrGenes[ i ].length( ); }
	sort( vecstrGenes.begin( ), vecstrGenes.end( ) );

	return ++iRet; }

CParserConsole::CParserConsole( const IOntology** apOntologies, const CGenome& Genome ) :
	CParser( apOntologies, Genome ) {

	m_sLocation.m_pOnto = NULL;
	m_sLocation.m_iNode = -1; }

CParser::SLocation CParserConsole::GetLocation( const string& strLoc, bool fLast ) const {

	return CParser::GetLocation( m_vecpOntologies, strLoc, fLast, &m_sLocation ); }

bool CParserConsole::ProcessLine( const char* szLine ) {
	vector<string>	vecstrLine;
	size_t			i;
	string			strLine;
	const char*		pcPrev;
	const char*		pcNext;

	if( !szLine )
		return false;

	for( pcPrev = szLine; pcPrev && *pcPrev; pcPrev = pcNext ) {
		vecstrLine.clear( );
		if( pcNext = strchr( pcPrev, c_cSemicolon ) )
			strLine.assign( pcPrev, pcNext++ - pcPrev );
		else
			strLine.assign( pcPrev );
		CMeta::Tokenize( strLine.c_str( ), vecstrLine, CMeta::c_szWS, true );
		if( vecstrLine.empty( ) )
			continue;
		if( vecstrLine[ 0 ][ 0 ] == c_cShell ) {
			if( !ParseShell( strLine ) )
				return false;
			continue; }
		for( i = 0; c_aszParsers[ i ]; ++i )
			if( !strcmp( vecstrLine[ 0 ].c_str( ), c_aszParsers[ i ] ) )
				break;
		if( !c_aszParsers[ i ] ) {
			cout << "Unknown command: " << strLine << endl;
			return false; }
		if( !(this->*c_apfnParsers[ i ])( vecstrLine ) )
			return false; }

	return true; }

bool CParserConsole::ParseCat( const vector<string>& vecstrLine ) {
	size_t				i, j;
	vector<string>		vecstrGenes;
	SArgs				sArgs;
	string				strGene, strPath;
	SLocation			sLoc;
	vector<SLocation>	vecVisited;

	for( i = 1; i < vecstrLine.size( ); ++i )
		if( !sArgs.Parse( vecstrLine[ i ] ) )
			vecstrGenes.push_back( vecstrLine[ i ] );
	if( !vecstrGenes.size( ) ) {
		cout << "Cat, no genes given" << endl;
		return false; }

	for( i = 0; i < vecstrGenes.size( ); ++i ) {
		strPath.clear( );
		strGene = vecstrGenes[ i ];
		if( ( j = strGene.rfind( c_cSep ) ) != -1 ) {
			strPath = strGene.substr( 0, j );
			strGene = strGene.substr( j + 1 ); }
		if( strGene == c_szStar ) {
			if( !Recurse( GetLocation( strPath ), sArgs.m_fRecursive, sArgs.m_fZeroes,
				vecVisited ) ) {
				cout << "cat, illegal location: " << strPath << endl;
				return false; }
			PrintGenes( vecVisited, sArgs ); }
		else if( ( j = m_Genome.GetGene( strGene ) ) == -1 )
			cout << "cat, unknown gene: " << strGene << endl;
		else
			PrintGene( m_Genome.GetGene( j ), sArgs ); }

	return true; }

void CParserConsole::PrintGenes( const vector<SLocation>& vecVisited,
	const SArgs& sArgs ) const {
	TSetPGenes					setpGenes;
	TSetPGenes::const_iterator	iterGene;

	CParser::CollectGenes( vecVisited, setpGenes );
	for( iterGene = setpGenes.begin( ); iterGene != setpGenes.end( ); ++iterGene )
		PrintGene( **iterGene, sArgs ); }

bool CParserConsole::ParseCd( const vector<string>& vecstrLine ) {
	SLocation	sLoc;

	if( vecstrLine.size( ) < 2 ) {
		cout << m_sLocation.ToString( true ) << endl;
		return true; }

	sLoc = GetLocation( vecstrLine[ 1 ] );
	if( !sLoc.IsValid( ) ) {
		cout << "cd, illegal location: " << vecstrLine[ 1 ] << endl;
		return false; }
	m_sLocation = sLoc;

	return true; }

bool CParserConsole::ParseFind( const vector<string>& vecstrLine ) {
	CGenes					Genes( (CGenome&)m_Genome ), GenesBkg( (CGenome&)m_Genome );
	ifstream				ifsm;
	size_t					i, j, k, l, iWidth, iTotal;
	vector<STermFound>		vecsTerms;
	vector<size_t>			veciOnto;
	string					strFile, strP, strBkg;
	SArgs					sArgs;
	const IOntology*		pOnto;
	vector<const CGene*>	vecpGenes;
	vector<string>			vecstrGenes;
	float					dP	= 0.05f;

	if( vecstrLine.size( ) < 2 )
		return false;
	for( i = 1; i < vecstrLine.size( ); ++i ) {
		if( sArgs.Parse( vecstrLine[ i ] ) )
			continue;
		if( !strFile.length( ) )
			strFile = vecstrLine[ i ];
		else if( !strP.length( ) )
			strP = vecstrLine[ i ];
		else if( !strBkg.length( ) )
			strBkg = vecstrLine[ i ]; }
	ifsm.open( strFile.c_str( ) );
	if( !( ifsm.is_open( ) && Genes.Open( ifsm, false ) ) ) {
		cout << "find, can't open file: " << strFile << endl;
		return false; }
	ifsm.close( );
	if( strP.length( ) )
		dP = (float)atof( strP.c_str( ) );
	if( strBkg.length( ) ) {
		ifsm.clear( );
		ifsm.open( strBkg.c_str( ) );
		if( !( ifsm.is_open( ) && GenesBkg.Open( ifsm ) ) ) {
			cout << "find, can't open background: " << strBkg << endl;
			return false; }
		ifsm.close( ); }

	CParser::TermFinder( Genes, dP, GenesBkg, sArgs.m_fBonferroni, sArgs.m_fSibs, sArgs.m_fBackground,
		veciOnto, vecsTerms );

	for( i = j = 0; i < m_vecpOntologies.size( ); ++i ) {
		pOnto = m_vecpOntologies[ i ];
		if( j >= veciOnto[ i ] )
			continue;
		cout << pOnto->GetID( ) << ':' << endl;

		l = j;
		if( !sArgs.m_fGenes ) {
			vecpGenes.clear( );
			for( ; j < veciOnto[ i ]; ++j ) {
				iTotal = pOnto->GetGenes( vecsTerms[ j ].m_iID, sArgs.m_fSibs );
				if( iTotal <= c_iSizeCutoff )
					for( k = 0; k < pOnto->GetGenes( vecsTerms[ j ].m_iID, sArgs.m_fSibs ); ++k )
						vecpGenes.push_back( &pOnto->GetGene( vecsTerms[ j ].m_iID, k ) );
				else
					for( k = 0; k < Genes.GetGenes( ); ++k )
						if( pOnto->IsAnnotated( vecsTerms[ j ].m_iID, Genes.GetGene( k ) ) )
							vecpGenes.push_back( &Genes.GetGene( k ) ); }
			vecstrGenes.clear( );
			iWidth = FormatGenes( vecpGenes, vecstrGenes, &Genes ); }

		for( j = l; j < veciOnto[ i ]; ++j ) {
			PrintAnnotation( pOnto, vecsTerms[ j ].m_iID, sArgs, &vecsTerms[ j ] );
			if( !sArgs.m_fGenes ) {
				vecpGenes.clear( );
				iTotal = pOnto->GetGenes( vecsTerms[ j ].m_iID, sArgs.m_fSibs );
				if( iTotal <= c_iSizeCutoff )
					for( k = 0; k < pOnto->GetGenes( vecsTerms[ j ].m_iID, sArgs.m_fSibs ); ++k )
						vecpGenes.push_back( &pOnto->GetGene( vecsTerms[ j ].m_iID, k ) );
				else
					for( k = 0; k < Genes.GetGenes( ); ++k )
						if( pOnto->IsAnnotated( vecsTerms[ j ].m_iID, Genes.GetGene( k ),
							sArgs.m_fSibs ) )
							vecpGenes.push_back( &Genes.GetGene( k ) );
				PrintGenes( vecpGenes, iWidth, &Genes ); } } }

	return true; }

bool CParserConsole::ParseHelp( const vector<string>& vecstrLine ) {
	size_t	i;

	if( vecstrLine.size( ) > 1 )
		for( i = 0; c_aszParsers[ i ]; ++i )
			if( vecstrLine[ 1 ] == c_aszParsers[ i ] ) {
				cout << c_aszHelps[ i ] << endl;
				return true; }

	cout << c_szHelpHelp << endl;

	return true; }

bool CParserConsole::ParseLs( const vector<string>& vecstrLine ) {
	SLocation			sLoc;
	size_t				i;
	string				strLoc;
	SArgs				sArgs;
	vector<SLocation>	vecVisited;

	for( i = 1; i < vecstrLine.size( ); ++i )
		if( !( sArgs.Parse( vecstrLine[ i ] ) || strLoc.size( ) ) )
			strLoc = vecstrLine[ i ];
	sLoc = strLoc.size( ) ? GetLocation( strLoc ) : m_sLocation;
	if( !Recurse( sLoc, sArgs.m_fRecursive, sArgs.m_fZeroes, vecVisited ) ) {
		cout << "ls, illegal location: " << strLoc << endl;
		return false; }

	PrintLocations( vecVisited, sArgs );
	return true; }

void CParserConsole::PrintLocations( const vector<SLocation>& vecVisited,
	const SArgs& sArgs ) const {
	const IOntology*		pOnto;
	size_t					i, j;
	vector<const CGene*>	vecpGenes;
	string					strLoc;

	for( i = 0; i < vecVisited.size( ); ++i ) {
		const SLocation&	sLoc	= vecVisited[ i ];

		if( pOnto = sLoc.m_pOnto ) {
			if( sLoc.m_iNode == -1 ) {
				PrintOntology( pOnto, '-' );
				if( sArgs.m_fSibs ) {
					for( j = 0; j < pOnto->GetNodes( ); ++j )
						if( !pOnto->GetParents( j ) )
							PrintLink( pOnto, j, 'C', sArgs ); } }
			else {
				PrintLink( pOnto, sLoc.m_iNode, '-', sArgs );
				if( sArgs.m_fSibs ) {
					for( j = 0; j < pOnto->GetParents( sLoc.m_iNode ); ++j )
						PrintLink( pOnto, pOnto->GetParent( sLoc.m_iNode, j ), 'P', sArgs );
					for( j = 0; j < pOnto->GetChildren( sLoc.m_iNode ); ++j )
						PrintLink( pOnto, pOnto->GetChild( sLoc.m_iNode, j ), 'C', sArgs ); }
				if( sArgs.m_fGenes ) {
					if( sArgs.m_fSibs )
						vecpGenes.clear( );
					for( j = 0; j < pOnto->GetGenes( sLoc.m_iNode ); ++j )
						vecpGenes.push_back( &pOnto->GetGene( sLoc.m_iNode, j ) );
					if( sArgs.m_fSibs )
						PrintGenes( vecpGenes ); } } }
		else {
			PrintOntology( NULL, '-' );
			if( sArgs.m_fSibs )
				for( j = 0; j < m_vecpOntologies.size( ); ++j )
					PrintOntology( m_vecpOntologies[ j ], 'O' ); } }

	if( sArgs.m_fGenes && !sArgs.m_fSibs ) {
		set<const CGene*>	setpGenes;

		for( i = 0; i < vecpGenes.size( ); ++i )
			setpGenes.insert( vecpGenes[ i ] );
		vecpGenes.resize( setpGenes.size( ) );
		copy( setpGenes.begin( ), setpGenes.end( ), vecpGenes.begin( ) );
		PrintGenes( vecpGenes ); } }

void CParserConsole::PrintOntology( const IOntology* pOnto, char cType ) const {
	string	strLoc;

	strLoc = pOnto ? pOnto->GetID( ) : "ROOT";
	cout << cType << ' ' << strLoc;
	PrintSpaces( c_iWidthOnto - strLoc.length( ) );
	if( pOnto )
		cout << (unsigned int)m_Genome.CountGenes( pOnto );
	cout << endl; }

bool CParserConsole::ParseShell( const string& strCmd ) const {
	size_t	i;

	i = strCmd.find( c_cShell );
	system( strCmd.substr( i + 1 ).c_str( ) );
	return true; }

bool CParserConsole::ParseParentage( const vector<string>& vecstrLine ) {
	string				strOnto, strFile;
	CSlim				Slim;
	const IOntology*	pOnto;
	size_t				i, j;
	ifstream			ifsm;
	SArgs				sArgs;
	vector<bool>		vecfTerms;

	if( vecstrLine.size( ) < 2 )
		return false;
	for( i = 1; i < vecstrLine.size( ); ++i ) {
		if( sArgs.Parse( vecstrLine[ i ] ) )
			continue;
		if( strOnto.empty( ) )
			strOnto = vecstrLine[ i ];
		else if( strFile.empty( ) )
			strFile = vecstrLine[ i ]; }

	pOnto = NULL;
	for( i = 0; i < m_vecpOntologies.size( ); ++i )
		if( strOnto == m_vecpOntologies[ i ]->GetID( ) ) {
			pOnto = m_vecpOntologies[ i ];
			break; }
	if( !pOnto ) {
		cout << "parentage, can't find ontology: " << strOnto << endl;
		return false; }

	ifsm.open( strFile.c_str( ) );
	if( !( ifsm.is_open( ) && Slim.Open( ifsm, pOnto ) ) ) {
		cout << "parentage, can't open file: " << strFile << endl;
		return false; }
	ifsm.close( );

	vecfTerms.resize( pOnto->GetNodes( ) );
	for( i = 0; i < Slim.GetSlims( ); ++i )
		for( j = 0; j < Slim.GetNodes( i ); ++j )
			vecfTerms[ Slim.GetNode( i, j ) ] = true;
	for( i = 0; i < pOnto->GetNodes( ); ++i ) {
		set<size_t>					setiParents;
		set<size_t>::const_iterator	iterParent;
		vector<size_t>				veciIntersection;

		pOnto->GetParents( i, setiParents );
		for( iterParent = setiParents.begin( ); iterParent != setiParents.end( ); ++iterParent )
			if( vecfTerms[ *iterParent ] )
				veciIntersection.push_back( *iterParent );
		if( veciIntersection.empty( ) && !sArgs.m_fZeroes )
			continue;
		cout << pOnto->GetID( i );
		for( j = 0; j < veciIntersection.size( ); ++j )
			cout << '\t' << pOnto->GetID( veciIntersection[ j ] );
		cout << endl; }

	return true; }

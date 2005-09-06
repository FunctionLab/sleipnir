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
	"cat <gene>+          Displays information on individual genes.\n"
	"cd [path]            Display or change current term.\n"
	"find <filename> [p]  Runs term finder on the given gene list.\n"
	"help [command]       Provides help on command syntax.\n"
	"ls [path]            List parents, children, and annotations.";
const CParserConsole::TPFnParser	CParserConsole::c_apfnParsers[]	=
	{ &CParserConsole::ParseCat, &CParserConsole::ParseCd, &CParserConsole::ParseFind,
	&CParserConsole::ParseHelp, &CParserConsole::ParseLs, NULL };
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
	"find <filename> [p=0.05] [-l] [-b]\n\n"
	"Performs a hypergeometric test over each ontology using the given gene list.\n"
	"Only terms with probability less than p are displayed, with a default p-value\n"
	"of 0.05.  The total number of possible genes is assumed to be the entire\n"
	"genome.  The optional flags are:\n"
	"-l  Long listings; deactivates term gloss abbreviation.\n"
	"-b  Bonferroni correction; deactivates Bonferroni correction.\n"
	"-g  Genes; display genes associated with each ontology term.\n"
	"-a  All listings; display additional information for ontology terms.\n"
	"-s  Siblings; deactivates child annotations during analysis.\n"
	"-k  Background; uses whole genome background in place of ontology background.",
	CParserConsole::c_szHelpHelp,
	"ls [-l] [-a] [-g] [-s] [path]\n\n"
	"With no arguments, the ls command displays the parents, children, and gene\n"
	"annotations of the current term.  Given a path, it displays the same\n"
	"information for that target instead.  The four optional flags are:\n"
	"-l  Long listings; deactivates term gloss abbreviation.\n"
	"-a  All listings; includes terms with zero gene annotations.\n"
	"-g  Genes; deactives gene listings.\n"
	"-s  Siblings; deactivates parent and child listings.\n"
	"-r  Recursive; descend into child nodes.",
	NULL };

CParserConsole::SArgs::SArgs( ) : m_fGenes(m_afFlags[ 0 ]), m_fLong(m_afFlags[ 1 ]),
	m_fSibs(m_afFlags[ 2 ]), m_fZeroes(m_afFlags[ 3 ]), m_fBonferroni(m_afFlags[ 4 ]),
	m_fRecursive(m_afFlags[ 5 ]), m_fBackground(m_afFlags[ 6 ]) {

	m_fGenes = true;
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
	const SArgs& sArgs, const CGenes* pGenes, double dP ) {
	char	szBuf[ 128 ];
	string	strID;
	size_t	iWidth;

	iWidth = c_iWidthGloss + c_iWidthGenes;
	cout << ( strID = pOnto->GetID( iNode ) );
	PrintSpaces( c_iWidthID - strID.size( ) );
	if( dP >= 0 ) {
		iWidth -= c_iWidthGenes;
		sprintf( szBuf, "%g", dP );
		cout << szBuf;
		PrintSpaces( c_iWidthP - strlen( szBuf ) ); }
	if( pGenes ) {
		sprintf( szBuf, "%-4d %-4d %-4d %-4d ", pGenes->CountAnnotations( pOnto, iNode,
			sArgs.m_fSibs ), pGenes->GetGenes( ), pOnto->GetGenes( iNode, sArgs.m_fSibs ),
			sArgs.m_fBackground ? pGenes->GetGenome( ).GetGenes( ) :
			pGenes->GetGenome( ).CountGenes( pOnto ) );
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
		cout << pOnto->GetID( ) << ':';
		PrintSpaces( c_iWidthOnto - pOnto->GetID( ).size( ) - 1 );
		PrintAnnotation( pOnto, Gene.GetAnnotation( i, 0 ), sArgs );
		for( j = 1; j < Gene.GetAnnotations( i ); ++j ) {
			PrintSpaces( c_iWidthOnto );
			PrintAnnotation( pOnto, Gene.GetAnnotation( i, j ), sArgs ); } } }

void CParserConsole::PrintGenes( const vector<const CGene*>& vecpGenes, size_t iWidth ) {
	size_t			i, iCol, iCols, iSpaces;
	vector<string>	vecstrGenes;

	iSpaces = 1;
	i = FormatGenes( vecpGenes, vecstrGenes );
	if( !iWidth )
		iWidth = i;
	iCols = c_iWidthScreen / iWidth;
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
	vector<string>& vecstrGenes ) {
	size_t	i, j, iRet;

	vecstrGenes.resize( vecpGenes.size( ) );
	for( iRet = i = 0; i < vecpGenes.size( ); ++i ) {
		vecstrGenes[ i ] = vecpGenes[ i ]->GetName( );
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
	CGenes					Genes( (CGenome&)m_Genome );
	ifstream				ifsm;
	size_t					i, j, k, l, iWidth;
	vector<TPrID>			vecprTerms;
	vector<size_t>			veciOnto;
	string					strFile, strP;
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
			strP = vecstrLine[ i ]; }
	ifsm.open( strFile.c_str( ) );
	if( !( ifsm.is_open( ) && Genes.Open( ifsm ) ) ) {
		cout << "find, can't open file: " << strFile << endl;
		return false; }
	ifsm.close( );
	if( strP.length( ) )
		dP = (float)atof( strP.c_str( ) );

	CParser::TermFinder( Genes, dP, sArgs.m_fBonferroni, sArgs.m_fSibs, sArgs.m_fBackground,
		veciOnto, vecprTerms );

	for( i = j = 0; i < m_vecpOntologies.size( ); ++i ) {
		pOnto = m_vecpOntologies[ i ];
		if( j >= veciOnto[ i ] )
			continue;
		cout << pOnto->GetID( ) << ':' << endl;

		l = j;
		if( !sArgs.m_fGenes ) {
			vecpGenes.clear( );
			for( ; j < veciOnto[ i ]; ++j )
				for( k = 0; k < Genes.GetGenes( ); ++k )
					if( pOnto->IsAnnotated( vecprTerms[ j ].first, Genes.GetGene( k ) ) )
						vecpGenes.push_back( &Genes.GetGene( k ) );
			vecstrGenes.clear( );
			iWidth = FormatGenes( vecpGenes, vecstrGenes ); }

		for( j = l; j < veciOnto[ i ]; ++j ) {
			PrintAnnotation( pOnto, vecprTerms[ j ].first, sArgs, sArgs.m_fZeroes ?
				&Genes : NULL, vecprTerms[ j ].second );
			if( !sArgs.m_fGenes ) {
				vecpGenes.clear( );
				for( k = 0; k < Genes.GetGenes( ); ++k )
					if( pOnto->IsAnnotated( vecprTerms[ j ].first, Genes.GetGene( k ),
						sArgs.m_fSibs ) )
						vecpGenes.push_back( &Genes.GetGene( k ) );
				PrintGenes( vecpGenes, iWidth ); } } }

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

	for( i = 1; i < vecstrLine.size( ); ++i ) {
		if( !( sArgs.Parse( vecstrLine[ i ] ) || strLoc.size( ) ) )
			strLoc = vecstrLine[ i ]; }
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
				if( sArgs.m_fSibs ) {
					PrintOntology( pOnto, '-' );
					for( j = 0; j < pOnto->GetNodes( ); ++j )
						if( !pOnto->GetParents( j ) )
							PrintLink( pOnto, j, 'C', sArgs ); } }
			else {
				if( sArgs.m_fSibs ) {
					PrintLink( pOnto, sLoc.m_iNode, '-', sArgs );
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
		else if( sArgs.m_fSibs ) {
			PrintOntology( NULL, '-' );
			for( j = 0; j < m_vecpOntologies.size( ); ++j )
				PrintOntology( m_vecpOntologies[ j ], 'O' ); } }

	if( sArgs.m_fGenes && !sArgs.m_fSibs )
		PrintGenes( vecpGenes ); }

void CParserConsole::PrintOntology( const IOntology* pOnto, char cType ) const {
	string	strLoc;

	strLoc = pOnto ? pOnto->GetID( ) : "ROOT";
	cout << cType << ' ' << strLoc;
	PrintSpaces( c_iWidthOnto - strLoc.length( ) );
	if( pOnto )
		cout << (unsigned int)m_Genome.CountGenes( pOnto );
	cout << endl; }
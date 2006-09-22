#ifndef PARSERCONSOLE_H
#define PARSERCONSOLE_H

#include "parser.h"

class CParserConsole : public CParser {
public:
	CParserConsole( const IOntology**, const CGenome& );

	bool ProcessLine( const char* );
	SLocation GetLocation( const string& = c_szDot, bool = true ) const;

protected:
	typedef bool (CParserConsole::*TPFnParser)( const vector<string>& );

	struct SArgs {
		static const char*	c_aszFlags[];

		bool	m_afFlags[ 7 ];
		bool&	m_fGenes;
		bool&	m_fLong;
		bool&	m_fSibs;
		bool&	m_fZeroes;
		bool&	m_fBonferroni;
		bool&	m_fRecursive;
		bool&	m_fBackground;

		SArgs( );

		bool Parse( const string& );
	};

	static const size_t		c_iWidthGenes		= 5;
	static const size_t		c_iWidthID			= 19;
	static const size_t		c_iWidthGloss		= 43;
	static const size_t		c_iWidthOnto		= 6;
	static const size_t		c_iWidthP			= 14;
	static const size_t		c_iWidthScreen		= 80;
	static const size_t		c_iSizeCutoff		= 40;
	static const char		c_cSemicolon		= ';';
	static const char		c_cShell			= '!';
	static const char		c_szDotDotDot[];
	static const char		c_szBackground[];
	static const char		c_szBonferroni[];
	static const char		c_szGenes[];
	static const char		c_szLong[];
	static const char		c_szRecursive[];
	static const char		c_szSibs[];
	static const char		c_szStar[];
	static const char		c_szZeroes[];
	static const char		c_szHelpHelp[];
	static const TPFnParser	c_apfnParsers[];
	static const char*		c_aszHelps[];

	static void PrintLink( const IOntology*, size_t, char, const SArgs& );
	static void PrintNumber( size_t, size_t );
	static void PrintSpaces( size_t );
	static void PrintAnnotation( const IOntology*, size_t, const SArgs&,
		const STermFound* = NULL );
	static void PrintGloss( string, size_t, bool );
	static void PrintGene( const CGene&, const SArgs& );
	static void PrintGenes( const vector<const CGene*>&, size_t = 0, const CGenes* = NULL );
	static size_t FormatGenes( const vector<const CGene*>&, vector<string>&, const CGenes* = NULL );

	bool ParseCat( const vector<string>& );
	bool ParseCd( const vector<string>& );
	bool ParseHelp( const vector<string>& );
	bool ParseLs( const vector<string>& );
	bool ParseFind( const vector<string>& );
	bool ParseShell( const string& ) const;
	void PrintOntology( const IOntology*, char ) const;
	void PrintLocations( const vector<SLocation>&, const SArgs& ) const;
	void PrintGenes( const vector<SLocation>&, const SArgs& ) const;

	SLocation	m_sLocation;
};

#endif // PARSERCONSOLE_H

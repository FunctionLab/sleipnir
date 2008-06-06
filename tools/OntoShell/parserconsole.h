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
#ifndef PARSERCONSOLE_H
#define PARSERCONSOLE_H

#include "parser.h"

class CParserConsole : public CParser {
public:
	CParserConsole( const Sleipnir::IOntology**, const Sleipnir::CGenome& );

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
	static const size_t		c_iWidthOnto		= 32;
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

	static void PrintLink( const Sleipnir::IOntology*, size_t, char, const SArgs& );
	static void PrintNumber( size_t, size_t );
	static void PrintSpaces( size_t );
	static void PrintAnnotation( const Sleipnir::IOntology*, size_t, const SArgs&,
		const Sleipnir::STermFound* = NULL );
	static void PrintGloss( string, size_t, bool );
	static void PrintGene( const Sleipnir::CGene&, const SArgs& );
	static void PrintGenes( const std::vector<const Sleipnir::CGene*>&, size_t = 0,
		const Sleipnir::CGenes* = NULL );
	static size_t FormatGenes( const std::vector<const Sleipnir::CGene*>&, std::vector<string>&,
		const Sleipnir::CGenes* = NULL );

	bool ParseCat( const std::vector<string>& );
	bool ParseCd( const std::vector<string>& );
	bool ParseHelp( const std::vector<string>& );
	bool ParseLs( const std::vector<string>& );
	bool ParseFind( const std::vector<string>& );
	bool ParseShell( const string& ) const;
	void PrintOntology( const Sleipnir::IOntology*, char ) const;
	void PrintLocations( const std::vector<SLocation>&, const SArgs& ) const;
	void PrintGenes( const std::vector<SLocation>&, const SArgs& ) const;

	SLocation	m_sLocation;
};

#endif // PARSERCONSOLE_H

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
#ifndef PARSER_H
#define PARSER_H

class CParser {
public:
	static const char	c_cSep	= '/';
	static const char	c_szDot[];
	static const char	c_szDotDot[];
	static const char*	c_aszParsers[];

	struct SLocation {
		static const char	c_szRoot[];

		const IOntology*	m_pOnto;
		size_t				m_iNode;

		SLocation( );

		bool operator==( const SLocation& ) const;

		string ToString( bool ) const;
		bool IsValid( ) const;
		void Invalidate( );
	};

	static const char* GetCommand( size_t );

	CParser( const Sleipnir::IOntology**, const Sleipnir::CGenome& );

	size_t GetOntologies( ) const;
	const IOntology* GetOntology( size_t ) const;
	const CGenome& GetGenome( ) const;

protected:
	typedef set<const CGene*>	TSetPGenes;

	static bool SplitLocation( const string&, std::vector<string>& );
	static bool IsRooted( const string& );
	static SLocation GetLocation( const std::vector<const Sleipnir::IOntology*>&, const string& = c_szDot,
		bool = true, const SLocation* = NULL );
	static bool MoveLocation( SLocation&, const string&, const std::vector<const Sleipnir::IOntology*>& );
	static void CollectGenes( const std::vector<SLocation>&, TSetPGenes& );

	bool Recurse( SLocation, bool, bool, std::vector<SLocation>& ) const;
	void TermFinder( const Sleipnir::CGenes&, float, const Sleipnir::CGenes&, bool, bool, bool,
		std::vector<size_t>&, std::vector<Sleipnir::STermFound>& ) const;

	const CGenome&				m_Genome;
	vector<const IOntology*>	m_vecpOntologies;
};

#endif // PARSER_H

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

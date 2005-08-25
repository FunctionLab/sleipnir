#ifndef GENOME_H
#define GENOME_H

#include <string>

#include "genomei.h"

namespace libBioUtils {

class CGene : CGeneImpl {
public:
	CGene( const std::string& );

	bool AddAnnotation( const IOntology*, size_t );
	bool AddSynonym( const std::string& );
	size_t GetOntologies( ) const;
	size_t GetAnnotations( size_t ) const;
	size_t GetAnnotation( size_t, size_t ) const;
	const std::string& GetName( ) const;
	const IOntology* GetOntology( size_t ) const;
	size_t GetSynonyms( ) const;
	const std::string& GetSynonym( size_t ) const;
	bool IsAnnotated( const IOntology* ) const;
	bool IsAnnotated( const IOntology*, size_t ) const;
	void SetRNA( bool );
	bool GetRNA( ) const;
	void SetDubious( bool );
	bool GetDubious( ) const;
	void SetGloss( const std::string& );
	const std::string& GetGloss( ) const;
};

class CGenome : CGenomeImpl {
public:
	bool Open( std::istream& );
	CGene& AddGene( const std::string& );
	CGene& GetGene( size_t ) const;
	size_t GetGene( const std::string& ) const;
	size_t GetGenes( ) const;
	size_t CountGenes( const IOntology* ) const;
	bool AddSynonym( CGene&, const std::string& );
};

class CGenes : CGenesImpl {
public:
	CGenes( CGenome& );

	bool Open( istream& );
	bool Open( const std::vector<std::string>& );
	void Filter( const CGenes& );
	size_t GetGenes( ) const;
	const CGene& GetGene( size_t ) const;
	const CGenome& GetGenome( ) const;
	size_t CountAnnotations( const IOntology*, size_t, bool = true ) const;
	bool IsGene( const std::string& ) const;
};

}

#endif // GENOME_H

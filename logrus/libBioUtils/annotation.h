#ifndef ANNOTATION_H
#define ANNOTATION_H

#include <string>

namespace libBioUtils {

typedef std::pair<size_t,double>	TPrID;

}

#include "annotationi.h"

namespace libBioUtils {

class IOntology {
public:
	virtual const std::string& GetID( ) const = 0;
	virtual size_t GetNodes( ) const = 0;
	virtual const std::string& GetID( size_t ) const = 0;
	virtual const std::string& GetGloss( size_t ) const = 0;
	virtual size_t GetParents( size_t ) const = 0;
	virtual size_t GetParent( size_t, size_t ) const = 0;
	virtual size_t GetChildren( size_t ) const = 0;
	virtual size_t GetChild( size_t, size_t ) const = 0;
	virtual size_t GetGenes( size_t, bool = false ) const = 0;
	virtual const CGene& GetGene( size_t, size_t ) const = 0;
	virtual bool IsAnnotated( size_t, const CGene&, bool = true ) const = 0;
	virtual size_t GetNode( const std::string& ) const = 0;
	virtual void GetGeneNames( std::vector<std::string>& ) const = 0;
	virtual void TermFinder( const CGenes&, std::vector<TPrID>&, bool = true,
		bool = true, bool = false ) const = 0;
};

class COntologyKEGG : COntologyKEGGImpl, public IOntology {
public:
	COntologyKEGG( );
	bool Open( std::istream&, CGenome& );

	const std::string& GetID( ) const;
	size_t GetNodes( ) const;
	const std::string& GetID( size_t ) const;
	const std::string& GetGloss( size_t ) const;
	size_t GetParents( size_t ) const;
	size_t GetParent( size_t, size_t ) const;
	size_t GetChildren( size_t ) const;
	size_t GetChild( size_t, size_t ) const;
	size_t GetGenes( size_t, bool ) const;
	const CGene& GetGene( size_t, size_t ) const;
	bool IsAnnotated( size_t, const CGene&, bool ) const;
	size_t GetNode( const std::string& ) const;
	void GetGeneNames( std::vector<std::string>& ) const;
	void TermFinder( const CGenes&, std::vector<TPrID>&, bool, bool, bool ) const;
};

class COntologyGO : COntologyGOImpl, public IOntology {
public:
	enum ENamespace {
		ENamespaceCC,
		ENamespaceBP,
		ENamespaceMF
	};

	static bool Open( std::istream&, std::istream&, CGenome&, COntologyGO&, COntologyGO&,
		COntologyGO& );

	COntologyGO( );
	bool Open( std::istream&, std::istream&, CGenome&, ENamespace );

	const std::string& GetID( ) const;
	size_t GetNodes( ) const;
	const std::string& GetID( size_t ) const;
	const std::string& GetGloss( size_t ) const;
	size_t GetParents( size_t ) const;
	size_t GetParent( size_t, size_t ) const;
	size_t GetChildren( size_t ) const;
	size_t GetChild( size_t, size_t ) const;
	size_t GetGenes( size_t, bool ) const;
	const CGene& GetGene( size_t, size_t ) const;
	bool IsAnnotated( size_t, const CGene&, bool ) const;
	size_t GetNode( const std::string& ) const;
	void GetGeneNames( std::vector<std::string>& ) const;
	void TermFinder( const CGenes&, std::vector<TPrID>&, bool, bool, bool ) const;
};

class COntologyMIPS : protected COntologyMIPSImpl, public IOntology {
public:
	COntologyMIPS( );
	bool Open( std::istream&, std::istream&, CGenome& );

	const std::string& GetID( ) const;
	size_t GetNodes( ) const;
	const std::string& GetID( size_t ) const;
	const std::string& GetGloss( size_t ) const;
	size_t GetParents( size_t ) const;
	size_t GetParent( size_t, size_t ) const;
	size_t GetChildren( size_t ) const;
	size_t GetChild( size_t, size_t ) const;
	size_t GetGenes( size_t, bool ) const;
	const CGene& GetGene( size_t, size_t ) const;
	bool IsAnnotated( size_t, const CGene&, bool ) const;
	size_t GetNode( const std::string& ) const;
	void GetGeneNames( std::vector<std::string>& ) const;
	void TermFinder( const CGenes&, std::vector<TPrID>&, bool, bool, bool ) const;
};

class COntologyMIPSPhenotypes : public COntologyMIPS {
public:
	COntologyMIPSPhenotypes( );

protected:
	static const char	c_szMIPSPhen[];
};

class CSlim : CSlimImpl {
public:
	bool Open( std::istream&, const IOntology* );
	size_t GetGenes( size_t ) const;
	size_t GetSlims( ) const;
	const std::string& GetSlim( size_t ) const;
	const CGene& GetGene( size_t, size_t ) const;
	void GetGeneNames( std::vector<std::string>& ) const;
};

}

#endif // ANNOTATION_H

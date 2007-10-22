#ifndef DATABASEI_H
#define DATABASEI_H

#include <vector>

namespace libBioUtils {

class CDatabaseImpl {
protected:
	static const char	c_acDAB[];
	static const char	c_acExtension[];

	struct SHeader {
		static const size_t	c_iName	= 128;

		char	m_acName[ c_iName ];
		size_t	m_iGenes;
		size_t	m_iDatasets;
	};

	size_t GetOffset( size_t iGene ) const {

		return ( sizeof(SHeader) + ( m_iDatasets * iGene ) ); }

	size_t GetOffset( size_t iGene, size_t iDataset ) const {

		return ( GetOffset( iGene ) + iDataset ); }

	std::vector<std::string>	m_vecstrGenes;
	size_t						m_iDatasets;
	std::vector<FILE*>			m_vecfile;
};

}

#endif // DATABASEI_H

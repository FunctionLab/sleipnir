#ifndef MATHBI_H
#define MATHBI_H

namespace libBioUtils {

class CMathImpl {
protected:
	static const size_t	c_iFactCutoff	= 10000;

	static std::vector<double>	s_vecdLogFact;

	static double LogFactStirling( size_t );
	static double LogFactRec( size_t );
};

}

#endif // MATHBI_H

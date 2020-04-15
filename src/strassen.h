#ifndef STRASSEN_H
#define STRASSEN_H

//#include "seekbasic.h"
#include <sstream>
#include <fstream>

namespace Sleipnir {

class CStrassen { 
public:
// Set LEAF_SIZE to 1 if you want to the pure strassen algorithm
// otherwise, the ikj-algorithm will be applied when the split
// matrices are as small as LEAF_SIZE x LEAF_SIZE
/*
 * Implementation of the strassen algorithm, similar to
 * http://en.wikipedia.org/w/index.php?title=Strassen_algorithm&oldid=498910018#Source_code_of_the_Strassen_algorithm_in_C_language
 */
static void strassen(vector<vector<float> > &A,
vector<vector<float> > &B,
vector<vector<float> > &C, unsigned int tam);

private:

//int leafsize = 64;
static void ikjalgorithm(vector<vector<float> > A,
vector<vector<float> > B, vector<vector<float> > &C, int n);

static unsigned int nextPowerOfTwo(int n);
static void strassenR(vector<vector<float> > &A,
vector<vector<float> > &B,
vector<vector<float> > &C, int tam);
static void sum(vector<vector<float> > &A,
vector<vector<float> > &B,
vector<vector<float> > &C, int tam);
static void subtract(vector<vector<float> > &A,
vector<vector<float> > &B,
vector<vector<float> > &C, int tam);
//void printMatrix(vector<vector<float> > matrix, int n);
//void read(string filename, vector<vector<float> > &A, vector<vector<float> > &B);
};

}
#endif

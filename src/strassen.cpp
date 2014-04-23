#include "stdafx.h"
#include "strassen.h"


namespace Sleipnir {

void CStrassen::ikjalgorithm(vector<vector<float> > A,
vector<vector<float> > B, vector<vector<float> > &C, int n) {
	for (int i = 0; i < n; i++)
		for (int k = 0; k < n; k++)
			for (int j = 0; j < n; j++)
				C[i][j] += A[i][k] * B[k][j];
}
 
void CStrassen::strassenR(vector< vector<float> > &A,
vector< vector<float> > &B,
vector< vector<float> > &C, int tam) {
	int leafsize = 64;

	if (tam <= leafsize){
		CStrassen::ikjalgorithm(A, B, C, tam);
		return;
	}
	// other cases are treated here:
	else {
		//fprintf(stderr, "TAM: %d\n", tam);
		int newTam = tam/2;
		vector<float> inner (newTam);
		vector<vector<float> >
			a11(newTam,inner), a12(newTam,inner), a21(newTam,inner), a22(newTam,inner),
			b11(newTam,inner), b12(newTam,inner), b21(newTam,inner), b22(newTam,inner),
			c11(newTam,inner), c12(newTam,inner), c21(newTam,inner), c22(newTam,inner),
			p1(newTam,inner), p2(newTam,inner), p3(newTam,inner), p4(newTam,inner),
			p5(newTam,inner), p6(newTam,inner), p7(newTam,inner),
			aResult(newTam,inner), bResult(newTam,inner);
  
		int i, j;
		//dividing the matrices in 4 sub-matrices:
		for (i = 0; i < newTam; i++) {
			for (j = 0; j < newTam; j++) {
				a11[i][j] = A[i][j];
				a12[i][j] = A[i][j + newTam];
				a21[i][j] = A[i + newTam][j];
				a22[i][j] = A[i + newTam][j + newTam];

				b11[i][j] = B[i][j];
				b12[i][j] = B[i][j + newTam];
				b21[i][j] = B[i + newTam][j];
				b22[i][j] = B[i + newTam][j + newTam];
			}
		}
  
        // Calculating p1 to p7:
		CStrassen::sum(a11, a22, aResult, newTam); // a11 + a22
		CStrassen::sum(b11, b22, bResult, newTam); // b11 + b22
		CStrassen::strassenR(aResult, bResult, p1, newTam); // p1 = (a11+a22) * (b11+b22)

		CStrassen::sum(a21, a22, aResult, newTam); // a21 + a22
		CStrassen::strassenR(aResult, b11, p2, newTam); // p2 = (a21+a22) * (b11)
  
		CStrassen::subtract(b12, b22, bResult, newTam); // b12 - b22
		CStrassen::strassenR(a11, bResult, p3, newTam); // p3 = (a11) * (b12 - b22)
  
		CStrassen::subtract(b21, b11, bResult, newTam); // b21 - b11
		CStrassen::strassenR(a22, bResult, p4, newTam); // p4 = (a22) * (b21 - b11)
  
		CStrassen::sum(a11, a12, aResult, newTam); // a11 + a12
		CStrassen::strassenR(aResult, b22, p5, newTam); // p5 = (a11+a12) * (b22)  
  
		CStrassen::subtract(a21, a11, aResult, newTam); // a21 - a11
		CStrassen::sum(b11, b12, bResult, newTam); // b11 + b12
		CStrassen::strassenR(aResult, bResult, p6, newTam); // p6 = (a21-a11) * (b11+b12)
  
		CStrassen::subtract(a12, a22, aResult, newTam); // a12 - a22
		CStrassen::sum(b21, b22, bResult, newTam); // b21 + b22
		CStrassen::strassenR(aResult, bResult, p7, newTam); // p7 = (a12-a22) * (b21+b22)

        // calculating c21, c21, c11 e c22:
		CStrassen::sum(p3, p5, c12, newTam); // c12 = p3 + p5
		CStrassen::sum(p2, p4, c21, newTam); // c21 = p2 + p4
  
		CStrassen::sum(p1, p4, aResult, newTam); // p1 + p4
		CStrassen::sum(aResult, p7, bResult, newTam); // p1 + p4 + p7
		CStrassen::subtract(bResult, p5, c11, newTam); // c11 = p1 + p4 - p5 + p7

		CStrassen::sum(p1, p3, aResult, newTam); // p1 + p3
		CStrassen::sum(aResult, p6, bResult, newTam); // p1 + p3 + p6
		CStrassen::subtract(bResult, p2, c22, newTam); // c22 = p1 + p3 - p2 + p6

		// Grouping the results obtained in a single matrix:
		for (i = 0; i < newTam ; i++) {
			for (j = 0 ; j < newTam ; j++) {
				C[i][j] = c11[i][j];
				C[i][j + newTam] = c12[i][j];
				C[i + newTam][j] = c21[i][j];
				C[i + newTam][j + newTam] = c22[i][j];
			}
		}
	}
}
 
unsigned int CStrassen::nextPowerOfTwo(int n) {
	return (unsigned int)pow(2, int(ceil(log2(n))));
}
 
void CStrassen::strassen(vector< vector<float> > &A,
vector< vector<float> > &B,
vector< vector<float> > &C, unsigned int n) {

	//unsigned int n = tam;
	unsigned int m = CStrassen::nextPowerOfTwo(n);
	vector<float> inner(m);
	vector< vector<float> > APrep(m, inner), BPrep(m, inner), CPrep(m, inner);
 
	for(unsigned int i=0; i<n; i++) {
		for (unsigned int j=0; j<n; j++) {
			APrep[i][j] = A[i][j];
			BPrep[i][j] = B[i][j];
		}
	}
 
	CStrassen::strassenR(APrep, BPrep, CPrep, m);
	for(unsigned int i=0; i<n; i++)
		for (unsigned int j=0; j<n; j++)
			C[i][j] = CPrep[i][j];
}
 
void CStrassen::sum(vector< vector<float> > &A,
vector< vector<float> > &B,
vector< vector<float> > &C, int tam) {
	int i, j;
	for (i = 0; i < tam; i++)
		for (j = 0; j < tam; j++)
			C[i][j] = A[i][j] + B[i][j];
}
 
void CStrassen::subtract(vector< vector<float> > &A,
vector< vector<float> > &B,
vector< vector<float> > &C, int tam) {
	int i, j;
	for (i = 0; i < tam; i++)
		for (j = 0; j < tam; j++) 
			C[i][j] = A[i][j] - B[i][j]; 
}

/* 
int CStrassen::getMatrixSize(string filename) {
	string line;
	ifstream infile;
	infile.open(filename.c_str());
	getline(infile, line);
	return count(line.begin(), line.end(), '\t') + 1;
}*/
 
/*
void read(string filename, vector< vector<float> > &A, vector< vector<float> > &B) {
	string line;
	FILE* matrixfile = freopen(filename.c_str(), "r", stdin);
     
	if (matrixfile == 0) {
        cerr << "Could not read file " << filename << endl;
        return;
    }
 
    int i = 0, j, a;
    while (getline(cin, line) && !line.empty()) {
        istringstream iss(line);
        j = 0;
        while (iss >> a) {
            A[i][j] = a;
            j++;
        }
        i++;
    }
 
    i = 0;
    while (getline(cin, line)) {
        istringstream iss(line);
        j = 0;
        while (iss >> a) {
            B[i][j] = a;
            j++;
        }
        i++;
    }
 
    fclose (matrixfile);
}
 
void printMatrix(vector< vector<float> > matrix, int n) {
    for (int i=0; i < n; i++) {
        for (int j=0; j < n; j++) {
            if (j != 0) {
                cout << "\t";
            }
            cout << matrix[i][j];
        }
        cout << endl;
    }
}
 
int main (int argc, char* argv[]) {
    string filename;
    if (argc < 3) {
        filename = "2000.in";
    } else {
        filename = argv[2];
    }
 
    if (argc < 5) {
        leafsize = 16;
    } else {
        leafsize = atoi(argv[4]);
    }
 
    int n = getMatrixSize(filename);
    vector<float> inner (n);
    vector< vector<float> > A(n, inner), B(n, inner), C(n, inner);
    read (filename, A, B);
    strassen(A, B, C, n);
    printMatrix(C, n);
    return 0;
}*/

}

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <iomanip>
#include <utility>
#include <math.h>
#include <cassert>

using namespace std;

class Matrix {
    public:
        int size;
        Matrix() {}
        virtual ~Matrix() {}
        virtual void multWithVec(double* vec, double* into) = 0;
};

class FullMatrix : public Matrix {
    private:
        double* data;
    public:
        FullMatrix(const string& filename);
        ~FullMatrix();
        void multWithVec(double* vec, double* into);
};

class CrMatrix : public Matrix {
    private:
        double* data;
        int* ci;
        int* addr;
    public:
        CrMatrix(const string& filename);
        ~CrMatrix();
        void multWithVec(double* vec, double* into);
};

/**
 * Method for reading matrix from input.
 * Input is expected in format:
 *     size numer_of_nonzero_values
 *     x_index  y_index  value
 *     ...      ...      ...
 * Optional parameter oneIndexed says if indexes starts at 1 or not
 * @param input input stream of matrix
 * @param oneIndexed whether indexes starts at 1
 * @return size of the matrix and the matrix itself in form of 1D array
 */
pair <int, double*> readMatrix(istream& input, bool oneIndexed = false) {
    int size = 0;
    int nonZeros = 0;
    
    input >> size >> nonZeros;

    double* matrix = new double[size * size]();

    int i, j = 0;
    double num = 0.0;
    for (int x = 0; x < nonZeros; x++) {
        input >> i >> j >> num;
        if (oneIndexed) { i--; j--;}
        matrix[i * size + j] = num;
    }

    cout << "Loaded matrix of size [" << size << ", " << size << "] with " << nonZeros << " non zero numbers." << endl;

    return make_pair(size, matrix);
}

/**
 * Method for reading matrix from input.
 * Input is expected in format:
 *     size
 *     value  value  ...
 *     ...    ...    ...
 * @param input input stream of matrix
 * @return size of the matrix and the matrix itself in form of 1D array
 */
pair <int, double*> readMatrixFull(istream& input) {
    int size = 0;
    
    input >> size;

    double* matrix = new double[size * size]();

    double num = 0.0;
    for (int x = 0; x < size * size; x++) {
        input >> num;
        matrix[x] = num;
    }

    cout << "Loaded matrix of size [" << size << ", " << size << "] with " << size * size << " non zero numbers." << endl;

    return make_pair(size, matrix);
}

/**
 * Method for reading matrix from input.
 * Based on the filename, it calls one of the methods above
 * @param input input stream of matrix
 * @return size of the matrix and the matrix itself in form of 1D array
 */
pair <int, double*> readMatrix(const string& filename) {
    cout << "Loading matrix from " << filename << endl;

    ifstream ifs(filename);

    if (filename.find("mat_cct") != string::npos || filename.find("mat_ps") != string::npos) {
        return readMatrix(ifs);
    } else if (filename.find("mat_cr") != string::npos) {
        return readMatrix(ifs, true);
    } else {
        return readMatrixFull(ifs);
    }
}

FullMatrix::FullMatrix(const string& filename) {
    pair<int, double*> matrix = readMatrix(filename);
    size = matrix.first;
    data = matrix.second;
}

FullMatrix::~FullMatrix() {
    delete[] data;
}

void FullMatrix::multWithVec(double* vec, double* into) {
    double acc = 0.0;
    for (int i = 0; i < size; i++) {
        acc = 0.0;
        for (int j = 0; j < size; j++) {
            acc += data[i * size + j] * vec[j];
        }
        into[i] = acc;
    }
}

CrMatrix::CrMatrix(const string& filename) {    
    ifstream input(filename);

    bool oneIndexed;
    if (filename.find("mat_cct") != string::npos || filename.find("mat_ps") != string::npos) {
        oneIndexed = false;
    } else if (filename.find("mat_cr") != string::npos) {
        oneIndexed = true;
    }

    int nonZeros;
    input >> size >> nonZeros;

    data = new double[nonZeros];
    ci = new int[nonZeros];
    addr = new int[size + 1];
    addr[size] = nonZeros;

    int i, j = 0;
    int previ = -1;
    double num = 0.0;
    for (int x = 0; x < nonZeros; x++) {
        input >> i >> j >> num;
        if (oneIndexed) { i--; j--;}
        
        data[x] = num;
        ci[x] = j;
        if (previ < i) {
            previ = i;
            addr[i] = x;
        }
    }
}

CrMatrix::~CrMatrix() {
    delete[] data;
    delete[] ci;
    delete[] addr;
}

void CrMatrix::multWithVec(double* vec, double* into) {
    double acc = 0.0;
    for (int i = 0; i < size; i++) {
        acc = 0.0;
        int start = addr[i];
        int end = addr[i + 1];
        for (int jX = start; jX < end; jX++) {
            int j = ci[jX];
            acc += data[jX] * vec[j];
        }
        into[i] = acc;
    }
}

void printVec(int size, double* vec) {
    cout << "[";
    for (int i = 0; i < size; i++) {
        cout << vec[i] << ", ";
    }
    cout << "]" << endl;
}

pair <int, double*> readVector(int size, istream& input) {
    double* vector = new double[size]();

    double num = 0.0;
    for (int i = 0; i < size; i++) {
        input >> num;
        vector[i] = num;
    }

    cout << "Loaded vector of size [" << size << "]" << endl;

    return make_pair(size, vector);
}

pair <int, double*> readVector(const string& filename) {
    cout << "Loading vector from " << filename << endl;

    ifstream ifs(filename);

    int size;

    size_t found = filename.find("vektor");
    if (found != string::npos) {
        string num = filename.substr(found + 6, filename.length() - found - 6 - 4);
        size = stoi(num);
    } else {
        ifs >> size;
    }

    return readVector(size, ifs);
}

double multVecVec(int size, double* a, double* b) {
    double acc = 0.0;
    for (int i = 0; i < size; i++) {
        acc += a[i] * b[i];
    }
    return acc;
}

void multNumVec(int size, double num, double* vec, double* into) {
    for (int i = 0; i < size; i++) {
        into[i] = vec[i] * num;
    }
}

double norm(int size, double* vec) {
    return sqrt(multVecVec(size, vec, vec));
}

void vecPlusNumVec(int size, double* a, double num, double* b, double* into) {
    for (int i = 0; i < size; i++) {
        into[i] = a[i] + num * b[i];
    }
}

void vecMinusNumVec(int size, double* a, double num, double* b, double* into) {
    for (int i = 0; i < size; i++) {
        into[i] = a[i] - num * b[i];
    }
}

void gradientDescent(int size, Matrix* A, double* b, double* x) {
    double* r = new double[size];
    double* h = new double[size];
    double alpha;

    // INITIALIZE
    for (int i = 0; i < size; i++) {
        r[i] = b[i];
        x[i] = 0.0;
    }

    int iter = 0;
    while (true) {
        A->multWithVec(r, h);
        alpha = multVecVec(size, r, r) / multVecVec(size, r, h);
        vecPlusNumVec(size, x, alpha, r, x);
        vecMinusNumVec(size, r, alpha, h, r);
        
        alpha = norm(size, r);
        iter++;
        if (iter % 1 == 0) {
            cout << "||R" << iter << "|| = " << alpha << endl;
        }
        if (alpha < 0.17) {
            break;
        }        
    }

    delete[] r;
    delete[] h;
}

void sdruGrad(int size, Matrix* A, double *b, double* x) {
    double* rO = new double[size];
    double* rN = new double[size];
    double* s = new double[size];
    double* h = new double[size];
    double* tmp;
    double alpha, beta;

    // INITIALIZE
    for (int i = 0; i < size; i++) {
        rO[i] = b[i];
        s[i] = b[i];
        x[i] = 0.0;
    }

    int iter = 0;
    while (true) {
        A->multWithVec(s, h);
        alpha = multVecVec(size, rO, rO) / multVecVec(size, s, h);
        vecPlusNumVec(size, x, alpha, s, x);

        vecMinusNumVec(size, rO, alpha, h, rN);
        
        alpha = norm(size, rN);
        iter++;
        if (iter % 1 == 0) {
            cout << "||R" << iter << "|| = " << alpha << endl;
        }
        if (alpha < 0.01) {
            break;
        }

        beta = multVecVec(size, rN, rN) / multVecVec(size, rO, rO);
        vecPlusNumVec(size, rN, beta, s, s);
        tmp = rO;
        rO = rN;
        rN = tmp;
    }

    delete[] rO;
    delete[] rN;
    delete[] h;
    delete[] s;

}

int main(int argc, char** argv) {

    if (argc < 3) {
        cout << "Provide filename of matrix and vector, eventually matrix storing method, like this:" << endl;
        cout << "\t" << argv[0] << " matrix.txt vector.txt [full/cr]" << endl;
        return 1;
    }

    string matrixFN = argv[1];
    string vectorFN = argv[2];

    bool CR = false;
    if (argc >= 4 && string(argv[3]).compare("cr") == 0) {
        CR = true;
    }
    

    Matrix* matrix;
    if (CR) {
        matrix = new CrMatrix(matrixFN);
    } else {
        matrix = new FullMatrix(matrixFN);
    }
    
    pair<int, double*> vectorR = readVector(vectorFN);

    int size = matrix->size;
    double* vec = vectorR.second;

    double* result = new double[size];
    double* result2 = new double[size];

    gradientDescent(size, matrix, vec, result);

    // cout << endl << endl;
    // printVec(size, result);
    matrix->multWithVec(result, result2);

    int x = 0;
    cout << "-> " << vec[x] << " - " << result2[x] << endl;
    x = 1;
    cout << "-> " << vec[x] << " - " << result2[x] << endl;
    x = 3;
    cout << "-> " << vec[x] << " - " << result2[x] << endl;
    x = 8;
    cout << "-> " << vec[x] << " - " << result2[x] << endl;
    x = 9;
    cout << "-> " << vec[x] << " - " << result2[x] << endl;
    // x = 99;
    // cout << "-> " << vec[x] << " - " << result[x] << endl;
    // x = 333;
    // cout << "-> " << vec[x] << " - " << result[x] << endl;


    delete matrix;
    delete[] vec;
    delete[] result;
    delete[] result2;
    
    return 0;
}
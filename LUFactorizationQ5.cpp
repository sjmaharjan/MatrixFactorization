/*
 * LUFactorization.cpp
 *
 *  Created on: Nov 30, 2013
 *      Author: sjmaharjan
 */

#include <iostream>
#include<string.h>
#include <stdlib.h>
#include<string>
#include <sys/time.h>
#include <time.h>
#include<fstream>
#include<math.h>
#include<cmath>
using namespace std;

/*
 *
 * Helper functions
 */
int min(int a, int b) {
  if (a < b)
    return a;
  else
    return b;
}

//function that returns the seconds
double get_seconds() {
  struct timeval tval;
  gettimeofday(&tval, NULL);
  return ((double) tval.tv_sec + (double) tval.tv_usec / 1000000.0);
}

/*
 *
 * Vector class
 */
class Vector {
private:
  int N;
  double *array1D;
  void copy(const Vector& other);
  void free();
public:
  Vector() {
    N = 0;
    array1D = NULL;
//    N = 0;
//    array1D = new double[N];
//    for (int i = 0; i < N; i++)
//      array1D[i] = 0;
  }
  Vector(int n) {
    N = n;
    array1D = new double[N];
    //initialization to zero
    for (int i = 0; i < N; i++) {
      array1D[i] = 0;
    }
  }

  Vector(const Vector& other);
  Vector& operator=(const Vector& other);

  void display() const;
  double &get(int n);
  double get(int n) const;
  void set(int n, double value);
  int getElementCount() const;
  int maxIndex(int from);
  void swap(int i, int j);

};
Vector::Vector(const Vector& other) {
  //cout << "called copy constructor" << endl;
  copy(other);
}

Vector& Vector::operator=(const Vector& other) {
  // cout << "Called equal operator" << endl;
  if (this != &other) {
    free();
    copy(other);
  }
  return *this;
}

void Vector::copy(const Vector& other) {
  N = other.getElementCount();
  array1D = new double[N];
  for (int i = 0; i < N; i++) {
    array1D[i] = other.get(i);
  }
}

void Vector::free() {
  delete[] array1D;
}

int Vector::getElementCount() const {
  return N;
}

double &Vector::get(int n) {
  return array1D[n];

}
double Vector::get(int n) const {
  return array1D[n];

}

void Vector::set(int n, double value) {
  array1D[n] = value;
}

int Vector::maxIndex(int from) {
  int maxIndex = from;
  double max = abs(get(from));
  for (int i = from; i < N; i++) {
    if (max < abs(get(i))) {
      max = abs(get(i));
      maxIndex = i;
    }
  }
  return maxIndex;
}

void Vector::swap(int i, int j) {
  double temp = get(i);
  set(i, get(j));
  set(j, temp);
}

//display vector

void Vector::display() const {

  cout << "\n";
  cout << "Vector:\n----------------\n";

  for (int i = 0; i < N; i++) {
    cout << get(i) << "  ";
  }
  cout << "\n---------------------\n";
}

/************************************************************************************************/
/*
 *
 * Matrix class
 */
class Matrix {

private:
  int rows;
  int columns;
  double **array2D;
  void copy(const Matrix& other);
  void free();

public:

  Matrix() {
//    rows = 1;
//    columns = 1;
//    array2D = new double*[rows];
//    for (int i = 0; i < rows; i++)
//      array2D[i] = new double[columns];
    rows = 0;
    columns = 0;
    array2D = NULL;

  }
  Matrix(int m, int n) {
    rows = m;
    columns = n;
    array2D = new double*[rows];
    for (int i = 0; i < rows; i++)
      array2D[i] = new double[columns];
    //initialization to zero
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < columns; j++) {
        array2D[i][j] = 0;
      }
    }

  }

  Matrix(char* filename) {
    fstream dataFile;
    //open file
    dataFile.open(filename, ios::in);
    dataFile >> rows;
    dataFile >> columns;
    //dynamically allocate memory
    array2D = new double*[rows];
    for (int i = 0; i < rows; i++)
      array2D[i] = new double[columns];

    //populate matrix from file;
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < columns; j++) {
        dataFile >> array2D[i][j];
      }
    }
    //finally close the file
    dataFile.close();
  }

  Matrix(const Matrix& other);
  Matrix& operator=(const Matrix& other);
  void populateMatrix();
  void display() const;
  double &get(int r, int c);
  double get(int r, int c) const;
  void set(int r, int c, double value);
  int getRows() const;
  int getColumns() const;
  void getColumnsByIndex(int columnIndex, double *result);
  Vector getColumnVectorByIndex(int colIndex);
  Vector getColumnVectorArrangedBy(int colIndex, Vector p);
  Matrix operator +(const Matrix &);
  Matrix operator -(const Matrix &);
  Matrix multiplyIJK(const Matrix& other);
  Matrix subMatrixFrom(int rowIndex, int columnIndex, Vector p);
  Matrix getU12(int r, Vector p);
  Matrix getL21(int r, Vector p);

  ~Matrix() {
    free();
  }
};

//copy constructor
Matrix::Matrix(const Matrix& other) {
  //cout << "called copy constructor" << endl;
  copy(other);
}

//assignment operator overloading

Matrix& Matrix::operator=(const Matrix& other) {
  // cout << "Called equal operator" << endl;
  if (this != &other) {
    free();
    copy(other);
  }
  return *this;
}

void Matrix::copy(const Matrix& other) {

  rows = other.getRows();
  columns = other.getColumns();

  array2D = new double*[rows];
  for (int k = 0; k < rows; k++)
    array2D[k] = new double[columns];

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      array2D[i][j] = other.get(i, j);
    }
  }

}

void Matrix::free() {
  for (int i = 0; i < rows; ++i) {
    delete[] array2D[i];
  }
  delete[] array2D;

}

// get the row count
int Matrix::getRows() const {
  return rows;
}

//get the column count
int Matrix::getColumns() const {
  return columns;
}

//get element stored at given row and column index

double &Matrix::get(int r, int c) {
  return array2D[r][c];

}

double Matrix::get(int r, int c) const {
  return array2D[r][c];

}

//set the value at given row and column index

void Matrix::set(int r, int c, double value) {
  array2D[r][c] = value;
}

//populate the content of matrix

void Matrix::populateMatrix() {
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++)
      set(i, j, (double) (i + j * 2 +3*(10*i+5*j+11)));
  }
}

//display matrix

void Matrix::display() const {

  cout << "\n";
  cout << "Matrix:\n----------------\n";
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++)
      cout << get(i, j) << "  ";
    cout << "\n";
  }
  cout << "---------------------\n";
}

Matrix Matrix::operator +(const Matrix &m) {
  Matrix final(m.columns, m.columns);
  for (int i = 0; i < m.columns; i++) {
    for (int j = 0; j < m.columns; j++) {
      final.array2D[i][j] = array2D[i][j] + m.array2D[i][j];
    }
  }
  return final;
}

Matrix Matrix::operator -(const Matrix &m) {
  Matrix final(m.columns, m.columns);
  for (int i = 0; i < m.columns; i++) {
    for (int j = 0; j < m.columns; j++) {
      final.array2D[i][j] = array2D[i][j] - m.array2D[i][j];
    }
  }
  return final;
}

//multiply ijk -method 1
Matrix Matrix::multiplyIJK(const Matrix& other) {

  Matrix result(rows, other.getColumns());
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < other.getColumns(); j++) {

      for (int k = 0; k < columns; k++) {
        result.set(i, j, result.get(i, j) + get(i, k) * other.get(k, j));
      }
    }
  }

  return result;
}

void Matrix::getColumnsByIndex(int columnIndex, double *result) {
  for (int i = 0; i < rows; i++)
    result[i] = array2D[i][columnIndex];
}

Vector Matrix::getColumnVectorByIndex(int colIndex) {
  Vector v(rows);
  for (int i = 0; i < rows; i++)
    v.set(i, array2D[i][colIndex]);
  return v;

}
Vector Matrix::getColumnVectorArrangedBy(int colIndex, Vector p) {
  Vector v(rows);
  for (int i = 0; i < rows; i++)
    v.set(i, array2D[(int) p.get(i)][colIndex]);
  return v;
}

Matrix Matrix::subMatrixFrom(int rowIndex, int columnIndex, Vector p) {
  Matrix sub(rows - columnIndex, columns - columnIndex);
  for (int i = 0; i < sub.getRows(); i++) {
    for (int j = 0; j < sub.getColumns(); j++) {
      sub.set(i, j, get(p.get(i + rowIndex), j + columnIndex));
    }
  }
  return sub;
}

Matrix Matrix::getU12(int r, Vector p) {
  Matrix U12(r, columns - r);
  for (int i = 0; i < U12.getRows(); i++) {
    for (int j = 0; j < U12.getColumns(); j++) {
      U12.set(i, j, get(p.get(i), j + r));
    }
  }
  return U12;
}
Matrix Matrix::getL21(int r, Vector p) {
  Matrix L21(rows - r, r);
  for (int i = 0; i < L21.getRows(); i++) {
    for (int j = 0; j < L21.getColumns(); j++) {
      L21.set(i, j, get(p.get(i + r), j));
    }
  }
  return L21;
}

/*******************************************************************************************************/
/*
 *
 *
 * LU Factorization class
 */

class LUFactorization {
private:
  Matrix A;
  Matrix LU;
  Vector P;

public:
  LUFactorization(Matrix& m) {
    A = m;
    LU = Matrix(m.getRows(), m.getColumns());
    LU = A;
    P = Vector(m.getRows());
    //initialize p with normal order
    for (int i = 0; i < P.getElementCount(); i++) {
      P.set(i, i);
    }
  }
  Vector getP();
  Matrix getLU();
  Matrix getFinalLU();

  void factorize(int rows, int columns);
  void factorize();
  void setLU(Matrix lu);
  void setP(Vector p);
  LUFactorization rightLookingBlockLU(Matrix& A, int n, int r);
  LUFactorization rightLookingBlockLU(int r);
  void display();
  void displayP();
  void displayLU();

};

void LUFactorization::setLU(Matrix lu) {
  LU = lu;
}

void LUFactorization::setP(Vector p) {
  P = p;
}
Vector LUFactorization::getP() {
  return P;
}
Matrix LUFactorization::getLU() {
  return LU;
}
Matrix LUFactorization::getFinalLU() {
  Matrix result(LU.getRows(), LU.getColumns());
  for (int i = 0; i < LU.getRows(); i++) {
    for (int j = 0; j < LU.getColumns(); j++)
      result.set(i, j, LU.get(P.get(i), j));
  }
  return result;
}

void LUFactorization::factorize(int rows, int columns) {
  for (int i = 0; i < columns; i++) {
    //pivoting starts
    Vector colT = LU.getColumnVectorArrangedBy(i, P);
    int maxIndex = colT.maxIndex(i);
    if (maxIndex > i) {
      P.swap(i, maxIndex);
    }
    //pivot ends
    //P.display();
    for (int j = i + 1; j < rows; j++) {
      LU.set(P.get(j), i, LU.get(P.get(j), i) / LU.get(P.get(i), i));
      for (int k = i + 1; k < columns; k++) {
        double partial = LU.get(P.get(j), k)
            - LU.get(P.get(j), i) * LU.get(P.get(i), k);
        LU.set(P.get(j), k, partial);
      }
    }
   // LU.display();
  }
}
void LUFactorization::factorize() {
  factorize(A.getRows(), A.getColumns());
}

LUFactorization LUFactorization::rightLookingBlockLU(Matrix& A, int n, int r) {
  if (n <= r) {
    LUFactorization lu(A);
    lu.factorize();
    return lu;
  } else {
    LUFactorization lu(A);
    lu.factorize(A.getRows(), r);
    Vector p = lu.getP();
    Matrix LU = lu.getLU();
    // p.display();
    // LU.display();
    //solve for u12 using L11*U12=A12
    for (int i = 0; i < r; i++) {
      for (int j = r; j < n; j++) {
        double sum = 0, result = 0;
        for (int k = 0; k < i; k++) {
          sum = sum + LU.get(p.get(i), k) * LU.get(p.get(k), j);
        }
        result = LU.get(p.get(i), j) - sum;
        LU.set(p.get(i), j, result);
      }
    }
    // LU.display();
    //get the submatrix and comput schur's complement
    Matrix L21 = LU.getL21(r, p);
    Matrix U12 = LU.getU12(r, p);
   // L21.display();
    // U12.display();
    Matrix L21U12 = L21.multiplyIJK(U12);
     //L21U12.display();
    Matrix ASub = LU.subMatrixFrom(r, r, p);
   // ASub.display();
    Matrix Abar = ASub - L21U12;
   // Abar.display();
    LUFactorization result = rightLookingBlockLU(Abar, n - r, r);
    //update P
    Vector PDash = result.getP();
    //PDash.display();
    Vector temp(PDash.getElementCount());
    for (int i = 0; i < PDash.getElementCount(); i++) {
      temp.set(i, p.get(r  + PDash.get(i)));
    }
   // cout << "Temp" << endl;
   // temp.display();
    for (int i = 0; i < PDash.getElementCount(); i++) {
      p.set(r + i , temp.get(i));
    }
   // cout << "Vector P" << endl;
   // p.display();
    //fill up the result
    Matrix LUDash = result.getLU();
    for (int i = 0; i < LUDash.getRows(); i++) {
      for (int j = 0; j < LUDash.getColumns(); j++) {
        LU.set(p.get(i + r), j + r, LUDash.get(PDash.get(i), j));
      }
    }
    lu.setP(p);
    lu.setLU(LU);
    return lu;
  }
}

LUFactorization LUFactorization::rightLookingBlockLU(int r) {
  return rightLookingBlockLU(A, A.getColumns(), r);
}
void LUFactorization::displayP() {
  for (int i = 0; i < P.getElementCount(); i++) {
    for (int j = 0; j < P.getElementCount(); j++) {
      if (j == P.get(i))
        cout << "1" << "\t";
      else
        cout << "0" << "\t";
    }
    cout << endl;
  }
  cout << endl;

}
void LUFactorization::displayLU() {
  cout << "L Matrix:" << endl;
  for (int i = 0; i < LU.getRows(); i++) {
    for (int j = 0; j < min(LU.getRows(), LU.getColumns()); j++) {
      if (j == i)
        cout << "1" << "\t";
      else if (j > i)
        cout << "0" << "\t";
      else
        cout << LU.get(P.get(i), j) << "\t";
    }
    cout << endl;
  }

  cout << "U Matrix:" << endl;
  for (int i = 0; i < min(LU.getRows(), LU.getColumns()); i++) {
    for (int j = 0; j < LU.getColumns(); j++) {
      if (j < i)
        cout << "0" << "\t";
      else
        cout << LU.get(P.get(i), j) << "\t";
    }
    cout << endl;
  }

}

void LUFactorization::display() {
  cout << "Permutation Matrix P:" << endl;
  displayP();
  cout << "LU Matrix:" << endl;
  displayLU();

}


int main(int argc, char *argv[]) {


  //Question no 5 N=128,256,512,1024
 fstream dataFile;
  double sec0, sec;
  int blockSizes[] = { 2, 4, 8, 16,20, 30,32,40,50,60, 64,70,80,100, 128,200,300,400,500,600,700,800,900,1000 };
 dataFile.open("output_time.csv", ios::out);
  dataFile << "Block Size,128,256,512,1024\n";
  //create matrix
  Matrix N128(128, 128), N256(256, 256), N512(512, 512), N1024(1024,1024);
  N128.populateMatrix();
  N256.populateMatrix();
  N512.populateMatrix();
  N1024.populateMatrix();
  for (int count = 0; count < 24; count++) {
    dataFile << blockSizes[count];
    if(blockSizes[count] <= 128){
    LUFactorization luN128(N128);
    sec0 = get_seconds();
    LUFactorization resultN128 = luN128.rightLookingBlockLU(blockSizes[count]);
    sec = get_seconds();
    dataFile << "," << (sec - sec0);
    }

    if(blockSizes[count] <= 256){
    LUFactorization luN256(N256);
    sec0 = get_seconds();
    LUFactorization resultN256 = luN256.rightLookingBlockLU(blockSizes[count]);
    sec = get_seconds();
    dataFile << "," << (sec - sec0);
    }

    if(blockSizes[count] <= 512){
    LUFactorization luN512(N512);
    sec0 = get_seconds();
    LUFactorization resultN512 = luN512.rightLookingBlockLU(blockSizes[count]);
    sec = get_seconds();
    dataFile << "," << (sec - sec0);
    }

    if(blockSizes[count] <= 1024){
    LUFactorization luN1024(N1024);
    sec0 = get_seconds();
    LUFactorization resultN1024 = luN1024.rightLookingBlockLU(blockSizes[count]);
    sec = get_seconds();
    dataFile << "," << (sec - sec0);
    }
    dataFile << "\n";
    cout<<"Loop "<< count <<" done"<<endl;
  }
  dataFile.close();
  return 0;
}

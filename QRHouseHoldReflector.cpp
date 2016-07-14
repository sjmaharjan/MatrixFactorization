//============================================================================
// Name        : MatrixFactorization.cpp
// Author      : suraj
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

//============================================================================

#include <iostream>
#include<string.h>
#include <stdlib.h>
#include<string>
#include <time.h>
#include<fstream>
#include<math.h>
using namespace std;

class Matrix {

private:
  int rows;
  int columns;

  double **array2D;
  void copy(const Matrix& other);
  void free();

public:

  Matrix() {
    rows = 5;
    columns = 5;
    array2D = new double*[rows];
    for (int i = 0; i < rows; i++)
      array2D[i] = new double[columns];
  }
  Matrix(int m, int n) {
    rows = m;
    columns = n;
    array2D = new double*[rows];
    for (int i = 0; i < rows; i++)
      array2D[i] = new double[columns];
    //initialization
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
  Matrix subMatrix(int index);
  Matrix operator +(const Matrix &);
  Matrix operator -(const Matrix &);

  Matrix multiplyIJK(const Matrix& other);

  ~Matrix() {

    free();
  }
};

//copy constructor
Matrix::Matrix(const Matrix& other) {
  cout << "called copy constructor" << endl;
  copy(other);
}

//assignment operator overloading

Matrix& Matrix::operator=(const Matrix& other) {
  cout << "Called equal operator" << endl;
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
      set(i, j, (double) (i + j * 2));
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

Matrix Matrix::subMatrix(int index) {
  Matrix sub = Matrix(rows, columns);
  for (int i = 0; i < index; i++)
    sub.get(i, i) = 1;
  for (int i = index; i < sub.getRows(); i++)
    for (int j = index; j < sub.getColumns(); j++)
      sub.get(i, j) = array2D[i][j];
  return sub;

}

class QRHouseHolderReflector {

private:
  Matrix A;
  Matrix R;

public:
  QRHouseHolderReflector(Matrix& m) {
    A = m;
  }

  void factorize();
  double vectorNorm(double v[], int n);
  Matrix getA() {
    return A;
  }
  Matrix getR() {
    return R;
  }

  Matrix houseHoldMatrix(double v[], int n);

};

/* m = I -2 v v^T/vT.v */
Matrix QRHouseHolderReflector::houseHoldMatrix(double v[], int n) {
  Matrix H = Matrix(n, n);
  //V^T.V
  double vTv = 0;
  for (int k = 0; k < n; k++) {
    vTv += v[k] * v[k];
  }
  if (vTv != 0) {
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        H.get(i, j) = -2 * v[i] * v[j] / vTv;
  }
  for (int i = 0; i < n; i++)
    H.get(i, i) += 1;
  return H;
}
;

void QRHouseHolderReflector::factorize() {
  Matrix Z = A, temp;
  for (int i = 0; i < A.getColumns() && i < (A.getRows() - 1); i++) {
    //first find vector
    temp = Z.subMatrix(i);
    double vector[temp.getRows()], a[temp.getRows()];
    temp.getColumnsByIndex(i, a);
    double normOfa = vectorNorm(a, temp.getRows());
    int signOfa = 1;
    if (a[i] < 0)
      signOfa = -1;

    for (int j = 0; j < temp.getRows(); j++) {
      if (j < i)
        vector[j] = 0;
      else if (j == i)
        vector[j] = temp.get(j, i) - signOfa * normOfa * 1;
      else
        vector[j] = temp.get(j, i);

    }
    //household matrix
    Matrix household = houseHoldMatrix(vector, temp.getRows());
    Z = household.multiplyIJK(Z);

  }
  R = Z;
  Z.display();

}
;

/* ||x|| */
double QRHouseHolderReflector::vectorNorm(double v[], int n) {
  double sum = 0;
  for (int i = 0; i < n; i++)
    sum += v[i] * v[i];
  return sqrt(sum);
}

int main(int argc, char *argv[]) {
  //Usage
  if (argc == 1) {
    cout << "Usage: QR foo.txt" << endl;
    exit(1);
  }

  char* matrixFile = argv[argc - 1];
  Matrix A(matrixFile);
  A.display();
  QRHouseHolderReflector qr(A);

  qr.factorize();
  qr.getR().display();

  return 0;
}

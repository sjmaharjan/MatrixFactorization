/*
 * QRFactorization.cpp
 *
 *  Created on: Nov 11, 2013
 *      Author: sjmaharjan
 */

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
#include<cmath>
using namespace std;

//vector functions

double vectorDotProduct(double v[], int n) {
  double vect = 0;
  for (int k = 0; k < n; k++) {
    vect += v[k] * v[k];
  }
  return vect;
}

void swapVector(double* v, int i, int j) {
  double temp = v[i];
  v[i] = v[j];
  v[j] = temp;
}

class Matrix {

private:
  int rows;
  int columns;

  double **array2D;
  void copy(const Matrix& other);
  void free();

public:

  Matrix() {
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
  void swapColumns(int i, int j);

  void getColumnsByIndex(int columnIndex, double *result);
  Matrix subMatrix(int index);
  Matrix operator +(const Matrix &);
  Matrix operator -(const Matrix &);
  Matrix transpose();

  Matrix multiplyIJK(const Matrix& other);

  ~Matrix() {

    free();
  }
};

//copy constructor
Matrix::Matrix(const Matrix& other) {
// cout << "called copy constructor" << endl;
  copy(other);
}

//assignment operator overloading

Matrix& Matrix::operator=(const Matrix& other) {
//  cout << "Called equal operator" << endl;
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
//  cout << "====================================" << endl;
//  cout << "H: " << endl;
//  display();
//  cout << "A: " << endl;
//  other.display();

  Matrix result(rows, other.getColumns());
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < other.getColumns(); j++) {

      for (int k = 0; k < columns; k++) {
        result.set(i, j, result.get(i, j) + get(i, k) * other.get(k, j));
      }
    }
  }

//  cout << "updated A: " << endl;
//  result.display();
//  cout << "====================================" << endl;
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

Matrix Matrix::transpose() {
  Matrix transpose(getRows(), getColumns());
  for (int i = 0; i < getRows(); i++) {
    for (int j = 0; j < getColumns(); j++) {
      transpose.set(i, j, get(j, i));
    }
  }
  return transpose;
}

void Matrix::swapColumns(int i, int j) {

  for (int k = 0; k < rows; k++) {
    double tmp = get(k, i);
    set(k, i, get(k, j));
    set(k, j, tmp);

  }

}

class QRHouseHolderReflector {

private:
  Matrix A;
  Matrix R;
  Matrix Q;

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
  Matrix getQ() {
    return Q;
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
        H.get(i, j) = -(2 * v[i] * v[j] / vTv);
  }
  for (int k = 0; k < n; k++)
    H.get(k, k) += 1;
  return H;
}

void QRHouseHolderReflector::factorize() {
//permutation vector
  double perm[A.getColumns()];
  for (int p = 0; p < A.getColumns(); p++) {
    perm[p] = p;
  }

  double columnNorms[A.getColumns()];
  for (int p = 0; p < A.getColumns(); p++) {
    double col[A.getRows()];
    A.getColumnsByIndex(p, col);
    columnNorms[p] = vectorDotProduct(col, A.getRows());
    // cout << columnNorms[p] << endl;
  }

  Matrix Z = A, temp;
  Matrix Identity(A.getRows(), A.getRows());
  for (int m = 0; m < Identity.getRows(); m++) {
    Identity.set(m, m, 1);
  }
  Q = Identity;

  for (int i = 0; i < A.getColumns() && i < (A.getRows() - 1); i++) {

//pivoting choosing p such that colnorm(p)= max(colmorms(j:n))
    double max = columnNorms[i];
    int colp = i;
    for (int p = i; p < A.getColumns(); p++) {

      if (max < columnNorms[p]) {
        max = columnNorms[p];
        colp = p;
      }
    }

    if (max == 0)
      return;

    if (i != colp) {
//swap columns
//swap p
      swapVector(perm, i, colp);
      swapVector(columnNorms, i, colp);

//swap Z
      Z.swapColumns(i, colp);
      // Z.display();

    }

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
    // Z.display();

    Q = Q.multiplyIJK(household);

//update the column norms
    for (int q = (i + 1); q < A.getColumns(); q++) {

      columnNorms[q] = columnNorms[q] - (Z.get(i, q) * Z.get(i, q));

    }

  }
  R = Z;
  //cout << "Q matrix" << endl;
  //Q.display();
//Z.display();

}
;

/* ||x|| */
double QRHouseHolderReflector::vectorNorm(double v[], int n) {
  double sum = 0;
  for (int i = 0; i < n; i++)
    sum += v[i] * v[i];
  return sqrt(sum);
}

/**
 * QR - Givens
 */

class QRGivens {

private:
  Matrix A;
  Matrix R;
  Matrix Q;

public:
  QRGivens(Matrix& m) {
    A = m;
  }

  void factorize();

  Matrix getA() {
    return A;
  }
  Matrix getR() {
    return R;
  }
  Matrix getQ() {
    return Q;
  }

  Matrix givens(double a, double b);
  Matrix G(int i, int j, Matrix& rotation);

};

/* m = I -2 v v^T/vT.v */
Matrix QRGivens::givens(double a, double b) {
  Matrix rotation(2, 2);
  double cos = 0, sine = 0, r = 0;
  if (b == 0) {
    cos = 1;
    sine = 0;
  } else {
    if (std::abs(b) > std::abs(a)) {
      r = -a / b;
      sine = 1.0 / sqrt(1 + r * r);
      cos = sine * r;
    } else {
      r = -b / a;
      cos = 1.0 / sqrt(1 + r * r);
      sine = cos * r;
    }
  }

  rotation.set(0, 0, cos);
  rotation.set(1, 1, cos);
  rotation.set(0, 1, -sine);
  rotation.set(1, 0, sine);

  return rotation;

}
;

Matrix QRGivens::G(int i, int j, Matrix& rotation) {
  Matrix G(A.getRows(), A.getRows());
  for (int k = 0; k < G.getRows(); k++) {
    G.set(k, k, 1);
  }

  G.set(i, i, rotation.get(1, 1));
  G.set(i - 1, i - 1, rotation.get(0, 0));
  G.set(i - 1, i, rotation.get(0, 1));
  G.set(i, i - 1, rotation.get(1, 0));

  return G;
}

void QRGivens::factorize() {
  R = A;
  Matrix Identity(A.getRows(), A.getRows());
  for (int m = 0; m < Identity.getRows(); m++) {
    Identity.set(m, m, 1);
  }
  Q = Identity;

//permutation vector
  double perm[A.getColumns()];
  for (int p = 0; p < A.getColumns(); p++) {
    perm[p] = p;
  }

  double columnNorms[A.getColumns()];
  for (int p = 0; p < A.getColumns(); p++) {
    double col[A.getRows()];
    A.getColumnsByIndex(p, col);
    columnNorms[p] = vectorDotProduct(col, A.getRows());
  }
  for (int j = 0; j < A.getColumns(); j++) {

//pivoting choosing p such that colnorm(p)= max(colmorms(j:n))
    double max = columnNorms[j];
    int colp = j;
    for (int p = j; p < A.getColumns(); p++) {

      if (max < columnNorms[p]) {
        max = columnNorms[p];
        colp = p;
      }
    }

    if (max == 0)
      return;

    if (j != colp) {
//swap columns
//swap p
      swapVector(perm, j, colp);
      swapVector(columnNorms, j, colp);

//swap Z
      R.swapColumns(j, colp);

    }

    for (int i = A.getRows() - 1; i >= (j + 1); i--) {
      Matrix rotation = givens(R.get(i - 1, j), R.get(i, j));
      Matrix g = G(i, j, rotation);
      // cout << "G matrix";
      //g.display();
      R = g.multiplyIJK(R);
      //cout<<"R matrix";
      //R.display();
      Q = Q.multiplyIJK(g.transpose());
    }

//update the column norms
    for (int q = j + 1; q < A.getColumns(); q++) {
      columnNorms[q] = columnNorms[q] - R.get(j, q) * R.get(j, q);

    }

  }
//R.display();
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

  cout << "QR- Household Reflector" << endl;
  QRHouseHolderReflector qrH(A);
  qrH.factorize();
  cout << "Q matrix" << endl;
  qrH.getQ().display();
  cout << "R matrix" << endl;
  qrH.getR().display();
  cout << "QR- Givens Rotation" << endl;
  QRGivens qrG(A);
  qrG.factorize();
  cout << "Q matrix" << endl;
  qrG.getQ().display();
  cout << "R matrix" << endl;
  qrG.getR().display();

  return 0;
}

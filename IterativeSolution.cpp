/*
 * GaussSeidel.cpp
 *
 *  Created on: Dec 7, 2013
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

int min(int a, int b) {
  if (a < b)
    return a;
  else
    return b;
}

int max(int a, int b) {
  if (a > b)
    return a;
  else
    return b;
}

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
  }
  Vector(int n) {
    N = n;
    array1D = new double[N];
    //initialization to zero
    for (int i = 0; i < N; i++) {
      array1D[i] = 0;
    }
  }

  Vector(char* filename) {
    int temp;
    fstream dataFile;
    //open file
    dataFile.open(filename, ios::in);
    dataFile >> N;
    dataFile >> temp;

    array1D = new double[N];
    for (int i = 0; i < N; i++) {
      dataFile >> array1D[i];
    }

    //finally close the file
    dataFile.close();
  }

  Vector(const Vector& other);
  Vector& operator=(const Vector& other);

  void display() const;
  double &get(int n);
  double get(int n) const;
  void set(int n, double value);
  int getElementCount() const;
  Vector subtract(Vector V);
  double absMax();
  void setInitialGuess();

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

Vector Vector::subtract(Vector V) {
  Vector result(V.getElementCount());
  for (int i = 0; i < V.getElementCount(); i++) {
    result.set(i, get(i) - V.get(i));
  }
  return result;
}

double Vector::absMax() {
  double max = abs(get(0));
  for (int i = 0; i < getElementCount(); i++) {
    if (max < abs(get(i))) {
      max = abs(get(i));
    }
  }
  return max;
}

void Vector::setInitialGuess() {
  for (int i = 0; i < getElementCount(); i++) {
    set(i, 1);
  }
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
  Matrix operator +(const Matrix &);
  Matrix operator -(const Matrix &);
  Matrix multiplyIJK(const Matrix& other);
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
      set(i, j, (double) (i + j * 2 + 3 * (10 * i + 5 * j + 11)));
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

/**********************************************************************************************/

class Jacobi {
private:
  Matrix A;
  Vector B, X;
  double e;
  int N;

public:
  Jacobi(Matrix& a, Vector& b, Vector& x) {
    A = a;
    B = b;
    X = x;
    e = 1E-5;
    N = 100;
  }

  void solve();
  void display();

};

void Jacobi::solve() {
  int n = A.getRows();
  bool solFound = false;
  for (int m = 0; m < N; m++) {
    Vector O = X;
    for (int j = 0; j < n; j++) {
      double l = 0;
      for (int k = 0; k < n; k++) {
        if (k != j)
          l += A.get(j, k) * O.get(k);
      }
      X.set(j, (B.get(j) - l) / A.get(j, j));
    }
    if (X.subtract(O).absMax() < e) {
     // X.display();
      solFound=true;
      break;
    }
  }
  if (!solFound)
  cout << "No solution satisfying the tolerance condition obtained after " << N
      << "iterations" << endl;
}

void Jacobi::display() {
  cout << "Input Matrix A:" << endl;
  A.display();
  cout << "Input Vector B:" << endl;
  B.display();

  cout << "Solution vector X:" << endl;
  X.display();

}

/**********************************************************************************************/

class GaussSeidel {
private:
  Matrix A;
  Vector B, X;
  double e;
  int N;

public:
  GaussSeidel(Matrix& a, Vector& b, Vector& x) {
    A = a;
    B = b;
    X = x;
    e = 1E-5;
    N = 100;
  }

  void solve();
  void display();

};

void GaussSeidel::solve() {
  int n = A.getRows();
  bool solFound = false;
  for (int m = 0; m < N; m++) {
    Vector O = X;
    for (int j = 0; j < n; j++) {
      double l = 0, u = 0;
      for (int k = 0; k < j; k++) {
        l += A.get(j, k) * X.get(k);
      }
      for (int k = j + 1; k < n; k++) {
        u += A.get(j, k) * O.get(k);
      }
      X.set(j, (B.get(j) - l - u) / A.get(j, j));
    }
    if (X.subtract(O).absMax() < e) {
      //X.display();
      solFound = true;
      break;
    }
  }
  if (!solFound)
    cout << "No solution satisfying the tolerance condition obtained after "
        << N << "iterations" << endl;
}

void GaussSeidel::display() {
  cout << "Input Matrix A:" << endl;
  A.display();
  cout << "Input Vector B:" << endl;
  B.display();

  cout << "Solution vector X:" << endl;
  X.display();

}

int main(int argc, char *argv[]) {
  //Usage
  if (argc < 3) {
    cout << "Usage: GS foo.txt vector.txt" << endl;
    exit(1);
  }

  char* matrixFile = argv[1];

  char* vectorFile = argv[2];
  Matrix A(matrixFile);

  Vector B(vectorFile);
  Vector I(B.getElementCount());
  I.setInitialGuess();
  cout << "Jacobi solution" << endl;
  Jacobi JSol(A, B, I);
  JSol.solve();
  JSol.display();
  cout << "Gauss Seidel solution" << endl;
  GaussSeidel sol(A, B, I);
  sol.solve();
  sol.display();

  return 0;
}

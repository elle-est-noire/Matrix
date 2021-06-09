#include <iostream>
#include <cmath>
#include <complex>
#include <random>
#include <set>
#include "Matrix.h"
#include "prob3-3.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#define THRESHOLD 1e-12
using namespace std;
using float100 = boost::multiprecision::cpp_dec_float_100;
using cp = boost::multiprecision::cpp_complex_100;


template <typename T>
bool compare(const vector<T> &v) {
  int n = v.size();
  // Vandermonde 行列
  Matrix<T> A(n, n), p, l, u;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      A.at(i, j) = pow(v[j], i);
    }
  }

  // LU分解
  A.lu(&p, &l, &u);
  T lu_val = T(1);
  for (int i = 0; i < n; i++)
    lu_val *= u.at(i, i);

  // 定義通りの行列式
  T det_val = A.det();

  // Vandermonde 行列式
  T poly_val = T(1);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i >= j) continue;
      poly_val *= v[j] - v[i];
    }
  }

  if (abs(lu_val - poly_val) >= THRESHOLD || abs(det_val - poly_val) >= THRESHOLD) {
    cout << "not equal!" << endl;
    for (const T &c : v)cout << c << ", ";
    cout << endl;
    cout << setprecision(30) << "lu=" << lu_val << ", det=" << det_val << ", poly=" << poly_val << endl;
    cout << lu_val - poly_val << " " << det_val - poly_val << endl;
    return false;
  }
  return true;
}

void prob33() {
  // PLU 分解の確認
  Matrix<double> A = { {0,2,3},{3,4,5},{5,6,9} }, p, l, u;
  A.lu(&p, &l, &u);
  cout << "P = ";
  p.print();
  cout << "\nL = ";
  l.print();
  cout << "\nU = ";
  u.print();
  cout << "\nP * L * U = ";
  (p * l * u).print();
  cout << endl;

  // Vandermonde 行列式
  int n = 5;
  random_device seed_gen;
  mt19937 mt(seed_gen());
  uniform_int_distribution<> rand(-1e5, 1e5);
  uniform_int_distribution<> crand(-1e5, 1e5);
  vector<cp> cv(n);
  vector<float100> v(n);

  // 実 Vandermonde 行列
  bool flg = true;
  for (int rep = 0; rep < 100; rep++) {
    for (int i = 0; i < n; i++)
      v[i] = rand(mt);
    flg = flg && compare(v);
  }
  cout << flg << endl;

  // 複素 Vandermonde 行列
  flg = true;
  for (int rep = 0; rep < 100; rep++) {
    for (int i = 0; i < n; i++)
      cv[i] = cp(crand(mt), crand(mt));
    flg = flg && compare(cv);
  }
  cout << flg << endl;
}
#include <iostream>
#include <cmath>
#include <complex>
#include <random>
#include <set>
#include "Matrix.h"
#include "prob3-3_noboost.h"
#include "prob3-3.h"
#define THRESHOLD 1e-12
using namespace std;
using cp = complex<double>;

template <typename T>
bool compare(const vector<T> &v, int &big, bool disp) {
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
    if (disp) {
      cout << "not equal!" << endl;
      for (const T &c : v)cout << c << ", ";
      cout << endl;
      cout << setprecision(30) << "lu=" << lu_val << ", det=" << det_val << ", poly=" << poly_val << endl;
      cout << lu_val - poly_val << " " << det_val - poly_val << endl;
    }
    big = abs(lu_val - poly_val) > abs(det_val - poly_val);
    return false;
  }
  return true;
}

void prob33_noboost() {
  // Vandermonde 行列式
  int n = 4;
  random_device seed_gen;
  mt19937 mt(seed_gen());
  uniform_int_distribution<> rand(-1e2, 1e2);
  uniform_int_distribution<> crand(-1e2, 1e2);
  vector<cp> cv(n);
  vector<long> v(n);

  // 実 Vandermonde 行列
  bool flg = true; int dummy;
  for (int rep = 0; rep < 100; rep++) {
    for (int i = 0; i < n; i++)
      v[i] = rand(mt);
    flg = flg && compare(v, dummy, true);
  }
  cout << flg << endl;

  // 複素 Vandermonde 行列
  flg = true;
  int sum = 0;
  for (int rep = 0; rep < 100; rep++) {
    for (int i = 0; i < n; i++)
      cv[i] = cp(crand(mt), crand(mt));
    int big = 0;
    flg = compare(cv, big, false) && flg;
    sum += big;
  }
  cout << flg << endl;
  cout << "sum = " << sum;
}
#pragma once
#include <vector>
#include <initializer_list>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <complex>

template <class T>
class Matrix {
private:
  std::vector<std::vector<T>> table;
  //up以上down未満の行，left以上right未満の列を切り出す
  int top;
  int bottom;
  int left;
  int right;

  Matrix<T> &resize(int n, int m);
  //行を入れ替える
  Matrix<T> &swapRows(int i, int j);
  //i 行に j 行の c 倍を加える
  Matrix<T> &addMultpliedRow(int i, int j, T c);
  //i 行を c 倍する
  Matrix<T> &multplyRow(int i, T c);
  //i 列に j 列の c 倍を加える
  Matrix<T> &addMultpliedColumn(int i, int j, T c);
  //i 列を c 倍する
  Matrix<T> &multplyColumn(int i, T c);
public:
  Matrix() {}
  Matrix(int n, int m) :top(0), bottom(n), left(0), right(m) { resize(n, m); }
  Matrix(std::initializer_list<std::initializer_list<T>> init) :table(init.begin(), init.end()), top(0), bottom(init.size()), left(0), right((*init.begin()).size()) {}
  Matrix(const Matrix<T> &mat) :table(mat.table), top(mat.top), bottom(mat.bottom), left(mat.left), right(mat.right) {}
  inline T at(int i, int j, bool absolutePos = false) const { return absolutePos ? table[i][j] : table[i + top][j + left]; }
  inline T &at(int i, int j, bool absolutePos = false) { return absolutePos ? table[i][j] : table[i + top][j + left]; }
  inline Matrix<T> & restrict(int up, int down, int left, int right);
  inline Matrix<T> &unrestrict();
  inline int height() const { return bottom - top; }
  inline int width() const { return right - left; }
  void print() const;
  T innerProduct(bool isVertical, int i, Matrix<T> &mat, bool isVertical2, int j) const;
  //1-ノルム Σ|a_ij|
  T norm1() const;
  //2-ノルムの二乗 Σ|a_ij|^2
  T sqNorm2() const;
  //max{|x_ij|}
  T absmax() const;
  //行列式
  T det() const;
  Matrix<T> operator+(const Matrix<T> &mat) const;
  Matrix<T> operator-(const Matrix<T> &mat) const;
  Matrix<T> operator*(const Matrix<T> &mat) const;
  Matrix<T> operator*(const T &c) const;
  Matrix<T> &operator+=(const Matrix<T> &mat);
  Matrix<T> &operator-=(const Matrix<T> &mat);
  Matrix<T> &operator*=(const T &c);
  //逆行列
  Matrix<T> inverse() const;
  //mat に逆行列を代入
  void inverse(Matrix<T> &mat) const;
  //転置行列
  Matrix<T> trans() const;
  //mat に転置行列を代入
  void trans(Matrix<T> &mat) const;
  //エルミート行列
  Matrix<T> hermite() const;
  //mat にエルミート行列を代入
  void hermite(Matrix<T> &mat) const;
  //Householder 変換によるQR分解
  //縦長行列かつ実行列を想定
  Matrix<T> qr(Matrix<T> *q) const;
  //GramSchmidt の直交化法による QR 分解
  //エルミート行列まで可
  Matrix<T> qrGramSchmidt(Matrix<T> *q) const;
  //A を QR 分解 → A = RQ を反復し，A (上三角行列)を返す．
  //Q は単位行列に近づき，R(A) には固有値が絶対値の大きい順に並ぶ
  //gs が true なら GramSchmidt の直交化法で QR 分解
  Matrix<T> qrAsymp(int rep, bool gs) const;
  //Wilkinson シフト
  //gs : trueならGramSchmidt
  Matrix<T> qrAsympShift(int rep, bool gs) const;
  Matrix<T> qrAsympShiftCp(int rep, bool gs) const;
  //LU 分解
  void lu(Matrix<T> *p, Matrix<T> *l, Matrix<T> *u) const;
};


template<class T>
inline Matrix<T> &Matrix<T>::restrict(int up, int down, int left, int right) {
  up = up, down = down, left = left, right = right;
  return *this;
}

template<class T>
inline Matrix<T> &Matrix<T>::unrestrict() {
  top = 0, bottom = table.size(), left = 0, right = table[0].size();
  return *this;
}

template<class T>
inline void Matrix<T>::print() const {
  std::cout << std::fixed << std::setprecision(10) << std::endl;
  for (int i = 0; i < height(); i++) {
    for (int j = 0; j < width(); j++)
      std::cout << at(i, j) << " ";
    std::cout << std::endl;
  }
}

template<class T>
inline T Matrix<T>::innerProduct(bool isVertical, int l1, Matrix<T> &mat, bool isVertical2, int l2) const {
  T ret(0);
  if (isVertical && isVertical2)
    for (int i = 0; i < height(); i++)
      ret += std::conj(at(i, l1)) * mat.at(i, l2);
  else if (isVertical && !isVertical2)
    for (int i = 0; i < height(); i++)
      ret += std::conj(at(i, l1)) * mat.at(l2, i);
  else if (!isVertical && isVertical2)
    for (int i = 0; i < width(); i++)
      ret += std::conj(at(l1, i)) * mat.at(i, l2);
  else if (!isVertical && !isVertical2)
    for (int i = 0; i < width(); i++)
      ret += std::conj(at(l1, i)) * mat.at(l2, i);
  return ret;
}

template<class T>
inline T Matrix<T>::norm1() const {
  T sum = T(0);
  for (int i = 0; i < height(); i++)
    for (int j = 0; j < width(); j++)
      sum += std::abs(at(i, j));
  return sum;
}

template<class T>
inline T Matrix<T>::sqNorm2() const {
  T ret = 0;
  for (int i = 0; i < height(); i++)
    for (int j = 0; j < width(); j++)
      ret += std::norm(at(i, j));
  return ret;
}

template<class T>
inline T Matrix<T>::absmax() const {
  T ret = abs(at(0, 0));
  for (int i = 0; i < height(); i++)
    for (int j = 0; j < width(); j++)
      if (abs(at(i, j)) > ret) ret = abs(at(i, j));
  return ret;
}

inline bool nextPermutation(std::vector<int> &vs, int &sgn) {
  int a = vs.size() - 2;
  while (0 <= a && vs[a + 1] <= vs[a])--a;
  if (a < 0)return false;
  int b = vs.size() - 1;
  while (vs[b] <= vs[a])--b;
  int t = vs[b];
  vs[b] = vs[a], vs[a] = t;
  sgn *= -1;
  for (int i = a + 1, j = vs.size() - 1; i < j; ++i, --j) {
    t = vs[i], vs[i] = vs[j], vs[j] = t;
    sgn *= -1;
  }
  return true;
}

template<class T>
T Matrix<T>::det() const {
  T sum = T(0);
  const int n = height();
  std::vector<int> perm(n);
  for (int i = 0; i < n; i++) perm[i] = i;
  int sgn = 1;
  do {
    T prod = T(sgn);
    for (int i = 0; i < n; i++) prod *= at(i, perm[i]);
    sum += prod;
  } while (nextPermutation(perm, sgn));
  return sum;
}

template<class T>
inline Matrix<T> Matrix<T>::operator+(const Matrix<T> &mat) const {
  Matrix<T> ret(*this);
  ret += mat;
  return ret;
}

template<class T>
inline Matrix<T> Matrix<T>::operator-(const Matrix<T> &mat) const {
  Matrix<T> ret(*this);
  ret -= mat;
  return ret;
}

template<class T>
inline Matrix<T> Matrix<T>::operator*(const Matrix<T> &mat) const {
  Matrix<T> ret(height(), mat.width());
  for (int i = 0; i < ret.height(); i++)
    for (int j = 0; j < ret.width(); j++) {
      ret.at(i, j) = T(0);
      for (int k = 0; k < table[i].size(); k++)
        ret.at(i, j) += at(i, k) * mat.at(k, j);
    }
  return ret;
}

template<class T>
inline Matrix<T> Matrix<T>::operator*(const T &c) const {
  Matrix<T> ret(*this);
  ret *= c;
  return ret;
}

template<class T>
inline Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &mat) {
  for (int i = 0; i < height(); i++)
    for (int j = 0; j < width(); j++)
      at(i, j) += mat.at(i, j);
  return *this;
}

template<class T>
inline Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &mat) {
  for (int i = 0; i < height(); i++)
    for (int j = 0; j < width(); j++)
      at(i, j) -= mat.at(i, j);
  return *this;
}

template<class T>
inline Matrix<T> &Matrix<T>::operator*=(const T &c) {
  for (int i = 0; i < height(); i++)
    for (int j = 0; j < width(); j++)
      at(i, j) *= c;
  return *this;
}

template<class T>
inline Matrix<T> Matrix<T>::inverse() const {
  Matrix<T> ret;
  inverse(ret);
  return ret;
}

template<class T>
inline void Matrix<T>::inverse(Matrix<T> &mat) const {
  const int n = height();
  mat.resize(n, n);
  for (int i = 0; i < n; i++)
    mat.at(i, i) = T(1);
  Matrix<T> org(*this);
  int rowIdx = 0;
  for (int colIdx = 0; colIdx < n; ++colIdx) {
    int pivotRowIdx = -1;
    for (int i = rowIdx; i < n; ++i) {
      if (at(i, colIdx) != T(0)) {
        pivotRowIdx = i;
        break;
      }
    }
    if (pivotRowIdx < 0) continue;
    if (pivotRowIdx != rowIdx) {
      org.swapRows(pivotRowIdx, rowIdx);
      mat.swapRows(pivotRowIdx, rowIdx);
    }
    for (int i = 0; i < n; i++) {
      if (i == rowIdx) continue;
      if (org.at(i, colIdx) == T(0)) continue;
      T c1 = T(1) / org.at(rowIdx, colIdx);
      org.multplyRow(rowIdx, c1);
      mat.multplyRow(rowIdx, c1);
      T c2 = -org.at(i, colIdx);
      org.addMultpliedRow(i, rowIdx, c2);
      mat.addMultpliedRow(i, rowIdx, c2);
    }
    ++rowIdx;
  }
}

template<class T>
inline Matrix<T> Matrix<T>::trans() const {
  Matrix<T> ret;
  trans(ret);
  return ret;
}

template<class T>
inline void Matrix<T>::trans(Matrix<T> &mat) const {
  mat.resize(width(), height());
  for (int i = 0; i < height(); i++)
    for (int j = 0; j < width(); j++)
      mat.at(j, i) = at(i, j);
}

template<class T>
inline Matrix<T> Matrix<T>::hermite() const {
  Matrix<T> ret;
  hermite(ret);
  return ret;
}

template<class T>
inline void Matrix<T>::hermite(Matrix<T> &mat) const {
  mat.resize(width(), height());
  for (int i = 0; i < height(); i++)
    for (int j = 0; j < width(); j++)
      mat.at(j, i) = conj(at(i, j));
}

template<class T>
inline Matrix<T> Matrix<T>::qr(Matrix<T> *q) const {
  Matrix<T> R(*this);
  Matrix<T> Q;
  int n = R.width();
  if (q != nullptr) {
    Q.resize(n, n);
    for (int i = 0; i < n; i++)
      Q.at(i, i) = T(1);
  }
  for (int diagIdx = 0; diagIdx < n - 1; diagIdx++) {
    Matrix<T> x(n - diagIdx, 1);
    for (int i = diagIdx; i < n; i++)
      x.at(i - diagIdx, 0) = R.at(i, diagIdx);
    //xの第１成分と同符号を選んだ方が，桁落ちが少ない
    x.at(0, 0) += std::real(x.at(0, 0)) > 0 ? sqrt(x.sqNorm2()) : -sqrt(x.sqNorm2());
    T nrm = -T(2) / x.sqNorm2();
    x = x * x.trans();
    x *= nrm;
    Matrix<T> H(n, n);
    for (int i = 0; i < n; i++)
      H.at(i, i) = T(1);
    for (int i = diagIdx; i < n; i++)
      for (int j = diagIdx; j < n; j++)
        H.at(i, j) += x.at(i - diagIdx, j - diagIdx);
    R = H * R;
    if (q != nullptr) Q = H * Q;
  }
  if (q != nullptr) Q.trans(*q);
  return R;
}

template<class T>
inline Matrix<T> Matrix<T>::qrGramSchmidt(Matrix<T> *Q) const {
  *Q = *this;
  //横の長さ
  int n = Q->width();
  Matrix<T> R(n, n);
  for (int i = 0; i < n; i++)
    R.at(i, i) = T(1);
  std::vector<T> sqNorms(n);
  sqNorms[0] = Q->innerProduct(true, 0, *Q, true, 0);
  for (int i = 1; i < n; i++) {
    for (int j = 0; j < i; j++) {
      T x = this->innerProduct(true, i, *Q, true, j) / sqNorms[j];
      for (int k = 0; k < n; k++)
        Q->at(k, i) -= x * Q->at(k, j);
      R.at(j, i) = x;
    }
    sqNorms[i] = Q->innerProduct(true, i, *Q, true, i);
  }
  for (int i = 0; i < n; i++) {
    sqNorms[i] = std::sqrt(sqNorms[i]);
    Q->multplyColumn(i, T(1) / sqNorms[i]);
    R.multplyRow(i, sqNorms[i]);
  }
  return R;
}

template<class T>
inline Matrix<T> Matrix<T>::qrAsymp(int rep, bool gs) const {
  Matrix<T> A(*this);
  Matrix<T> Q, R;
  for (int i = 0; i < rep; i++) {
    if (gs) R = A.qrGramSchmidt(&Q);
    else R = A.qr(&Q);
    A = R * Q;
  }
  return A;
}

template<class T>
inline Matrix<T> Matrix<T>::qrAsympShift(int rep, bool gs) const {
  Matrix<T> A(*this);
  const int n = A.height();
  Matrix<T> Q(n, n), R(n, n);
  for (int i = 0; i < rep; i++) {
    T a = A.at(n - 2, n - 2), b = A.at(n - 2, n - 1), c = A.at(n - 1, n - 2), d = A.at(n - 1, n - 1);
    T delta = (a - d) / T(2);
    T mu = delta * delta + b * c;
    if (mu < 0) mu = d;
    else {
      mu = b * c / (abs(delta) + sqrt(mu));
      if (delta < 0) mu = -mu;
      mu = d - mu;
    }
    for (int j = 0; j < n; j++)
      A.at(j, j) -= mu;
    if (gs) R = A.qrGramSchmidt(&Q);
    else R = A.qr(&Q);
    A = R * Q;
    for (int j = 0; j < n; j++)
      A.at(j, j) += mu;
  }
  return A;
}

template<class T>
inline Matrix<T> Matrix<T>::qrAsympShiftCp(int rep, bool gs) const {
  Matrix<T> A(*this);
  const int n = A.height();
  Matrix<T> Q(n, n), R(n, n);
  for (int i = 0; i < rep; i++) {
    T a = A.at(n - 2, n - 2), b = A.at(n - 2, n - 1), c = A.at(n - 1, n - 2), d = A.at(n - 1, n - 1);
    T delta = (a - d) / T(2);
    T mu = sqrt(delta * delta + b * c);
    if (std::norm(delta + mu) < std::norm(delta - mu)) mu = d + delta + mu;
    else mu = d + delta - mu;
    for (int j = 0; j < n; j++)
      A.at(j, j) -= mu;
    if (gs) R = A.qrGramSchmidt(&Q);
    else R = A.qr(&Q);
    bool nan = false;
    for (int j = R.height() - 1; !nan && j >= 0; j--)
      for (int k = R.width() - 1; !nan && k >= 0; k--)
        if (std::isnan(R.at(j, k).real()) || std::isnan(R.at(j, k).imag())
          || std::isnan(Q.at(j, k).real()) || std::isnan(Q.at(j, k).imag()))
          nan = true;
    if (nan) {
      for (int j = 0; j < n; j++)
        A.at(j, j) += mu;
      break;
    }
    A = R * Q;
    for (int j = 0; j < n; j++)
      A.at(j, j) += mu;
    //A.print();
  }
  return A;
}

template<class T>
inline Matrix<T> &Matrix<T>::resize(int n, int m) {
  if (table.size() != n) table.resize(n);
  if (table[0].size() != m)
    for (int i = 0; i < n; i++)
      table[i].resize(m);
  top = 0, bottom = n, left = 0, right = m;
  return *this;
}

template<class T>
inline Matrix<T> &Matrix<T>::swapRows(int i, int j) {
  std::vector<T> tmp = table[i + top];
  table[i + top] = table[j + top], table[j + top] = tmp;
  return *this;
}

template<class T>
inline Matrix<T> &Matrix<T>::addMultpliedRow(int i, int j, T c) {
  for (int k = 0; k < width(); k++)
    at(i, k) += at(j, k) * c;
  return *this;
}

template<class T>
inline Matrix<T> &Matrix<T>::multplyRow(int i, T c) {
  for (int k = 0; k < width(); k++)
    at(i, k) *= c;
  return *this;
}

template<class T>
inline Matrix<T> &Matrix<T>::addMultpliedColumn(int i, int j, T c) {
  for (int k = 0; k < width(); k++)
    at(k, i) += at(k, j) * c;
  return *this;
}

template<class T>
inline Matrix<T> &Matrix<T>::multplyColumn(int i, T c) {
  for (int k = 0; k < height(); k++)
    at(k, i) *= c;
  return *this;
}

template<class T>
inline void Matrix<T>::lu(Matrix<T> *p, Matrix<T> *l, Matrix<T> *u) const {
  int n = width();
  p->resize(n, n);
  l->resize(n, n);
  u->resize(n, n);
  *l = *this;
  for (int i = 0; i < n; i++)
    p->at(i, i) = 1;
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      if (l->at(i, i) == T(0)) {
        for (int k = i + 1; k < n; k++)
          if (l->at(k, i) != T(0)) {
            l->swapRows(i, k);
            p->at(i, i) = p->at(k, k) = 0;
            p->at(i, k) = p->at(k, i) = 1;
            break;
          }
      } else
        l->at(j, i) /= l->at(i, i);
      for (int k = i + 1; k < n; k++)
        l->at(j, k) -= l->at(j, i) * l->at(i, k);
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i > j) { u->at(i, j) = 0; continue; }
      u->at(i, j) = l->at(i, j);
      if (i == j)l->at(i, j) = 1;
      else l->at(i, j) = 0;
    }
  }
}

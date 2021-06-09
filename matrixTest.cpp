#include <iostream>
#include <complex>
#include "Matrix.h"
using namespace std;
using cp = complex<double>;

void test() {

  Matrix<double> M0 = { {-1,2,3},{1,3,5},{2,3,4} };
  Matrix<double> M = { {1,1,2},{1,2,1},{2,1,1} };
  Matrix<double> M2 = { {1,-1,1},{1,1,3},{1,0,0} };
  Matrix<double> M3 = { {5,4,3,3,2},{4,3,2,1,1},{3,2,1,1,1},{2,1,1,1,1},{1,1,1,1,1} };
  Matrix<double> M4 = { {4,3,2,1},{3,2,1,1},{2,1,1,1},{1,1,1,1} };
  Matrix<double> E = { {1,0,0},{0,1,0},{0,0,1} };
  Matrix<cp> M5 = { {cp(1,0),cp(-1,0)},{cp(1,0),cp(0,0)} };
  Matrix<cp> M_ = {
    {cp(0,0),cp(-1,0),cp(1,0)},
    {cp(0,0),cp(0,0),cp(-1,0)},
    {cp(-1,0),cp(0,0),cp(0,0)}
  };
  // Matrix<double> N = { {3,1,1},{16,-3,2},{-20,10,5} };
  Matrix<double> ode2 = { {-4,-12,8},{2,6,-4},{1,3,-2} };
  Matrix<double> A = { {1,1,1},{6,0,2},{-10,5,3} };
  // Matrix<double> S = A - N;
  Matrix<double> P1 = { {0,0,0},{-2,1,0},{2,-1,0} };
  Matrix<double> P2 = E - P1;
  Matrix<double> P1_ = (A - E * 3) * (A - E * 3) * (1.0 / 25.0);
  Matrix<double> P2_ = (E * 8 - A) * (A + E * 2) * (1.0 / 25.0);
  Matrix<double> S_ = { {0,1,0},{-1,2,0},{1,0,1} };
  Matrix<double> J = { {-2,0,0},{0,3,0},{0,0,3} };
  Matrix<double> S_1 = { {2,-1,0},{1,0,0},{-2,1,1} };
  auto S = S_ * J * S_1;
  auto N = A - S;
  //S.print();
  (P1_ * N).print();
  (P2_ * N).print();
  //P1_.print();
  //P2_.print();

  //for (int i = 0; i < 5; i++) {
  //  N = (N * N);
  //  N.print();
 // }
  /*
  cout << "M : " << endl;
  M_.print();
  cout << "\nR : " << endl;
  auto R = M_.qrAsympShiftCp(40, true);
  R.print();
  /*R.at(0, 0) += cp(2);
  R.at(1, 1) += cp(2);
  R.at(2, 2) += cp(2);
  Matrix<cp> R2,Q;
  R2=R.qrGramSchmidt(&Q);
  R2.print();
  /*cout << "M : " << endl;
  M.print();
  cout << "\nR : " << endl;
  M.qrAsympShift(36, true).print();
  //*/
  /*Matrix<cp> Q;
  Matrix<cp> R = M_.qrGramSchmidt(&Q);
  cout << "\nQ : " << endl;
  Q.print();
  cout << Q.det() << endl;
  cout << "\nR : " << endl;
  R.print();
  (Q.trans() * Q).print();
  (Q * R).print();
  //*/
}
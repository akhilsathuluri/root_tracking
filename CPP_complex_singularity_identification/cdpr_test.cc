#include<iostream>
#include<eigen3/Eigen/Dense>
#include<math.h>
#include"cdpr.hh"
#include"cdpr_inputs.hh"

using namespace Eigen;

int main(int argc, char const *argv[]) {
  std::cout << initsols << '\n';
  std::cout << " " << '\n';

  VectorXd l(3), q(9);
  l << 7.5, 10.0, 9.5;
  q << (initsols.row(4)).transpose(), l;
  std::cout << q << '\n';
  std::cout << Jetaphi(q) << '\n';
  std::cout << eta(q) << '\n';
  return 0;
}

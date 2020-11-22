#include<iostream>
#include<eigen3/Eigen/Dense>
#include "root_tracker.h"

using namespace Eigen;

/*!
@todo See how to pass a function into another function and verify NRTracker
*/
int main(int argc, char const *argv[]) {
  RootTracker rt;
  rt.methods();
  // VectorXd a(2);
  // a = rt.NRTracker();
  // std::cout << a << std::endl;
  return 0;
}

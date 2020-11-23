#ifndef RootTracker_H
#define RootTracker_H

#include<eigen3/Eigen/Dense>
#include<functional>

using namespace Eigen;

/*!
The RootTracker class consists of implementations of the following:
- Newton-Raphson method based tracking or NRTracker
- Davidenkos method based tracking or DMTracker
- Nearest neighbour based tracking or NNTracker
*/
class RootTracker
{
public:
  VectorXd NRTracker(VectorXd x, VectorXd y, std::function<VectorXd (VectorXd)> f,\
   std::function<MatrixXd (VectorXd)> Jfy, double eps = pow(10, -10));
  // void DMTracker();
  // void NNTracker();
  int methods();
};

#endif

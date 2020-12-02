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
@todo Add chrono to profile function evaluation timings
@todo See how things change when openMP is enabled with Eigen
@todo See how things change when BLAS and LAPACK are used with Eigen
@todo Implement an event identification method that can handle or atleast
alert when the system moves close to a singularity
*/
class RootTracker
{
public:
  VectorXd NRTracker(VectorXd x, VectorXd y, std::function<VectorXd (VectorXd)> f,\
   std::function<MatrixXd (VectorXd)> Jfy, double eps = pow(10, -10));
  VectorXd DMTracker(VectorXd xprev, VectorXd x, VectorXd y, std::function<MatrixXd (VectorXd)> Jfx, \
     std::function<MatrixXd (VectorXd)> Jfy, double eps = 0, std::function<VectorXd (VectorXd)> f = NULL);
  VectorXd NNTracker(VectorXd ys, MatrixXd ysols, int index);
  int Methods();
};

#endif

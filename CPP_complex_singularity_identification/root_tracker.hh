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

@todo See how things change when openMP is enabled with Eigen
@todo See how things change when BLAS and LAPACK are used with Eigen
@todo Implement an event identification method that can handle or atleast
alert when the system moves close to a singularity

@todo modify the paper with a definition of distance as all the individual
elements of the vector being within the eps radius ball. So once we define the
eps, we take 2 steps, which gives us 4 solution values, i.e. two branches two sols.
Then use a multivariable interpolation scheme to interpolate from both sides. Then
find their intersection  to find the location of singularity.
*/
class RootTracker
{
public:
  int Methods();
  VectorXd NRTracker(VectorXd x, VectorXd y, std::function<VectorXd (VectorXd)> f,\
   std::function<MatrixXd (VectorXd)> Jfy, double eps = pow(10, -10));
  VectorXd DMTracker(VectorXd xprev, VectorXd x, VectorXd y, std::function<MatrixXd (VectorXd)> Jfx, \
     std::function<MatrixXd (VectorXd)> Jfy, double eps = 0, std::function<VectorXd (VectorXd)> f = NULL);
  VectorXd NNTracker(VectorXd ys, MatrixXd ysols, int index);

  MatrixXd trackAllBranches(VectorXd x, MatrixXd y, std::function<VectorXd (VectorXd)> f, \
    std::function<MatrixXd (VectorXd)> Jfy);

  VectorXcd NRCTracker(VectorXcd x, VectorXcd y, std::function<VectorXcd (VectorXcd)> f,\
   std::function<MatrixXcd (VectorXcd)> Jfy, double eps = pow(10, -10));
  VectorXcd DMCTracker(VectorXcd xprev, VectorXcd x, VectorXcd y, std::function<MatrixXcd (VectorXcd)> Jfx, \
   std::function<MatrixXcd (VectorXcd)> Jfy, double eps = 0, std::function<VectorXcd (VectorXcd)> f = NULL);
   VectorXcd NNCTracker(VectorXcd ys, MatrixXcd ysols, int index);

  VectorXd SingularityEventIdentifier(VectorXd ys, MatrixXd ysols, int index, double eps = pow(10, -2));
};

#endif

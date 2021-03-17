#include<iostream>
#include<eigen3/Eigen/Dense>
#include<math.h>
#include"manipulator.hh"
#include"utils.hh"
#include<cassert>

using namespace Eigen;

int SEI(MatrixXd allroots, double alpha, int selectedroot, function computefromalpha, Ref<VectorXd> alphahist, Ref<MatrixXd> disthist, \
  std::function<VectorXd (VectorXd)> f, std::function<MatrixXd (VectorXd)> Jfy){
  VectorXd x(6), dist(allroots.rows()), graddist(allroots.rows()), preddist(allroots.rows());
  double stepsize;
  assert(pushHist(alphahist, alpha) && "Pushing alpha to history unsuccessful");
  // Edit the required path parametrisation in the computeXfromParam function
  x = computeXfromParam(alpha);
  allroots = RootTracker::trackAllBranches(x, allroots, f, Jfy);
  // Compute the distance between the selected root and the rest
  // Returns the distance with itself too, size = nuber of branches
  dist = computeDist(allroots, selectedroot);
  assert(pushHist(disthist, dist) && "Pushing alpha to history unsuccessful");
  stepsize = (alphahist(0)-alphahist(1));
  graddist = (disthist.row(0)-disthist.row(1))/stepsize;
  preddist = ((disthist.row(0)).transpose())+graddist*stepsize;
  for (size_t i = 0; i < preddist.size(); i++) {
    if (preddist(i) < 0 &&  i!= selectedroot){
      std::cout << "Singularity approaching. Branches " << selectedroot << " and " << i << " are going to merge!" << '\n';
      // quad = findExtrapCoeffs(alpha, dist.col(i));
      // x = computealphasandxfromcoeff(quad);
      // std::cout << "The estimated singular configuration is "<< x << '\n';
      return 1;
    }
  }
  return 0;
}

int main(int argc, char const *argv[]) {
  MatrixXd distvec(3, 2);
  VectorXd alpha(3), grad(2);
  distvec << 1, 2, 3, 4, 5, 6;
  alpha << 1, 2, 3;
  double stepsize = alpha(0) - alpha(1);
  grad = (distvec.row(0) - distvec.row(1))/stepsize;
  std::cout << ((distvec.row(0)).transpose())+grad*stepsize << '\n';
  // std::cout << (distvec.row(0) - distvec.row(1))/alpha(1) << '\n';

  return 0;
}

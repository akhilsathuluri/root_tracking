#include<iostream>
#include<eigen3/Eigen/Dense>
#include<math.h>
#include"inputs.hh"

using namespace Eigen;

VectorXd computeDist(MatrixXd allroots, int selectedroot, int qx=1){
  VectorXd dist((allroots.rows())-1);
  int j = 0;
  for (size_t i = 0; i < allroots.rows(); i++) {
    if (i!=selectedroot) {
      if (qx == 1) {
        dist(j) = (((allroots.row(selectedroot)).head(6)) - ((allroots.row(i)).head(6))).norm();
        j++;
      }
    }
  }
  return dist;
}

VectorXd findExtrapCoeffs(VectorXd xlist, VectorXd ylist, int degree = 2){
  MatrixXd Amat(xlist.size(), degree+1);
  VectorXd coeff(degree+1);
  for (size_t i = 0; i < xlist.size(); i++) {
    for (size_t j = 0; j <= degree; j++) {
      Amat(i, j) = pow(xlist(i), degree-j);
    }
  }
  if (xlist.size()==degree+1) {
    coeff = Amat.inverse()*ylist;
  }
  else {
    coeff = (Amat.transpose()*Amat).inverse()*Amat.transpose()*ylist;
  }
  return coeff;
}

int main(int argc, char const *argv[]) {
  VectorXd xlist(4), ylist(4), coeffs(3);
  xlist << -1, 0, 1, 1;
  ylist << 0, -1, 0, 1;

  coeffs = findExtrapCoeffs(xlist, ylist);
  // VectorXd dist(1);
  // dist(0) = (((initsols.row(1)).head(6)) - ((initsols.row(2)).head(6))).norm();
  std::cout << coeffs << '\n';
  return 0;
}

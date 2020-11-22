#include<iostream>
#include<math.h>
#include<eigen3/Eigen/Dense>
#include "root_tracker.h"
#include "utils.h"

using namespace Eigen;

/*!
The methods function can be called to view names of all the
implemented root tracking methods.
*/
int RootTracker::methods(){
  std::cout << "NRTracker, DMTracker, NNTracker" << std::endl;
  return 0;
}

/*!
The NRTracker uses the Newton-Raphson method iteratively to find
the solutions satisfied by the constrain equations. The output is the
values of the unknown variables at each tracking step.
@param x The set of input/known variables at each tracking step
@param y The set of output/unknown variables at each step
@param f The set of expressions of the non-linear functions relating x, y
@param Jfy The expression of the Jacobian matrix of f with respect to y, the unknown variables
@param eps The required tolerance to which the computed solutions are to satisfy the non-linear equations.
Default value is set to 10^-10
@todo Decide to return q or y in the NRTracker method
*/
VectorXd RootTracker::NRTracker(VectorXd x, VectorXd y, VectorXd f, MatrixXd Jfy, float eps = pow(10, -10)){
  VectorXd q, tempy, dy, fval;
  MatrixXd Jval;
  int loopcounter = 0;
  q << y, x;
  fval = f(q);
  tempy = y;
  while ((fval.cwiseAbs()).maxCoeff()>=eps) {
		if (loopcounter>=100){
      std::cout << "Warning::Max iterations reached. Convergence \
      not achived for the input tolerance, eps" << std::endl;
      return tempy;
    }
		loopcounter++;
    Jval = Jfy(q);
		dy = LinearSolve(Jval, fval);
		tempy -= dy;
		q << tempy, x;
		fval = f(q);
	}
  // Decide to return q or only y
	return tempy;
  }

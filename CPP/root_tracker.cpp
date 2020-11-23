#include<iostream>
#include<math.h>
#include<eigen3/Eigen/Dense>
#include "root_tracker.h"
#include "utils.h"
#include<functional>

using namespace Eigen;

/*!
The methods function can be called to view names of all the
implemented root tracking methods.
@todo Give provision for an FK solver along with root trackers
to use when trackers fail.
*/
int RootTracker::methods(){
  std::cout << "NRTracker, DMTracker, NNTracker" << std::endl;
  return 0;
}

/*!
The NRTracker uses the Newton-Raphson method iteratively to find
the solutions satisfied by the constrain equations. The output is the
values of the unknown variables at each tracking step.
The function accepts functions as arguments using the C++11 style using
the functional library.
@param x The set of input/known variables at each tracking step
@param y The set of output/unknown variables at each step
@param f The set of expressions of the non-linear functions relating x, y.
Takes in a single argument of type `VectorXd` consisting y, x and outputs the evaluation of f.
@param Jfy The expression of the Jacobian matrix of f with respect to y, the unknown variables.
Takes in a single argument of type `VectorXd` consisting y, x and outputs the evaluation of Jfy.
@param eps The required tolerance to which the computed solutions are to satisfy the non-linear equations.
Default value is set to 10^-10
@todo Decide to return q or y in the NRTracker method
@todo fval returning nan is not being handled in the while loop
*/
VectorXd RootTracker::NRTracker(VectorXd x, VectorXd y, std::function<VectorXd (VectorXd)> f, \
  std::function<MatrixXd (VectorXd)> Jfy, double eps /*= pow(10, -10)*/){
  VectorXd q(x.size()+y.size()), tempy(y.size()), dy(y.size()), fval;
  MatrixXd Jval;
  int loopcounter = 0;
  q << y, x;
  fval = f(q);
  tempy = y;
  while ((fval.cwiseAbs()).maxCoeff()>=eps) {
		if (loopcounter>=100){
      std::cout << "Warning::Max iterations reached. Convergence not achived for the input tolerance, eps" << std::endl;
      return tempy;
    }
		loopcounter++;
    Jval = Jfy(q);
		dy = LinearSolve(Jval, fval);
    // std::cout << dy << std::endl;
		tempy = tempy - dy;
		q << tempy, x;
		fval = f(q);
	}
	return tempy;
  }

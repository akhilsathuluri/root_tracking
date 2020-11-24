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
This function accepts functions as arguments using the `C++11` style functional library.
@param x The set of input/known variables at current tracking step
@param y The set of output/unknown variables at current step
@param f The set of expressions of the non-linear functions relating x, y.
Takes in a single argument of type `VectorXd` consisting y, x and outputs the evaluation of f.
@param Jfy The expression of the Jacobian matrix of f with respect to y, the unknown variables.
Takes in a single argument of type `VectorXd` consisting y, x and outputs the evaluation of Jfy.
@param eps The required tolerance to which the computed solutions are to satisfy the non-linear equations.
Default value is set to 10^-10
@todo Decide to return q or y in the NRTracker method
@todo fval returning nan is not being handled in the while loop
@todo Input argument format and other relevant checks
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
		dy = linearSolve(Jval, fval);
		tempy = tempy - dy;
		q << tempy, x;
		fval = f(q);
	}
	return tempy;
  }

  /*!
  The DMTracker uses the Davindenkos' integration method to find
  the solutions satisfied by the constrain equations. The output is the
  values of the unknown variables at each tracking step. The problem of tracking
  is solved as an initial value problem using the first order derivative form of
  the constraint equations. This function accepts functions as arguments using the `C++11` style
  functional library. This method by default uses the Explicit Euler integration scheme and is the only
  supported scheme currently.
  @param x The set of input/known variables at the current tracking step
  @param xnext The set of input/known variables at the next tracking step
  @param y The set of output/unknown variables at current step
  @param Jfx The expression of the Jacobian matrix of f with respect to x, the known variables.
  Takes in a single argument of type `VectorXd` consisting y, x and outputs the evaluation of Jfx
  @param Jfy The expression of the Jacobian matrix of f with respect to y, the unknown variables.
  Takes in a single argument of type `VectorXd` consisting y, x and outputs the evaluation of Jfy.
  @param eps The tolerance of `drift` to which the computed solutions are to satisfy the non-linear equations.
  If the drift exceeds the given tolerance, a Newton-Raphson (NR) step is used to bring the
  variables back to the constraint manifold. The default value of `eps` is set to 0, meaning the NR step correction
  is `off` by default.
  @param f The set of expressions of the non-linear functions relating x, y.
  This parameter is only required if `eps` is defined.
  @todo Provide support for different integration methods
  @todo This isnt working. Implement it using prevtheta like in Mathematica rather than using xnext. That should fix it.
  @todo The same scheme estimating dx as xnext-x is unstable
  */
  VectorXd RootTracker::DMTracker(VectorXd xprev, VectorXd x, VectorXd y, std::function<MatrixXd (VectorXd)> Jfx, \
     std::function<MatrixXd (VectorXd)> Jfy, double eps /*= 0*/, std::function<VectorXd (VectorXd)> f /*= NULL*/){
    VectorXd q(x.size()+y.size()), tempy(y.size()), fval;
    tempy = y;
    q << y, x;
    tempy += -(Jfy(q).inverse()*Jfx(q))*(x-xprev);
    // Check to apply an NR step
    if (eps > 0){
      q << tempy, x;
      VectorXd fval;
      fval = f(q);
      if ((fval.cwiseAbs()).maxCoeff()>=eps){
        tempy = RootTracker::NRTracker(x, tempy, f, Jfy);
      }
    }
    return tempy;
  }

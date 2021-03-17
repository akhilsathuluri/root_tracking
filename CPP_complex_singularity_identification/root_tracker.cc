/*
This file `root_tracker.h` implements three different root-tracking methods.
Copyright (C) 2020  Akhil Sathuluri

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include<iostream>
#include<math.h>
#include<eigen3/Eigen/Dense>
#include "root_tracker.hh"
#include "utils.hh"
#include<functional>
#include<cassert>

using namespace Eigen;

/*!
The methods function can be called to view names of all the
implemented root tracking methods.
@todo Give provision for an FK solver along with root trackers
to use when trackers fail.
*/
int RootTracker::Methods(){
  std::cout << "Listing available trackers:" << std::endl;
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
  // Jval = Jfy(q);
  // dy = linearSolve(Jval, fval);
  // std::cout << tempy - dy << '\n';
  while (fval.norm()>=eps) {
		if (loopcounter>=100){
      std::cout << "NRTracker::Warning::Max iterations reached. Convergence not achived for the input tolerance, eps" << std::endl;
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
  The NRCTracker uses the Newton-Raphson method iteratively to find
  the solutions satisfied by the constrain equations.
  @todo Write overloaded function for saveDate
  @todo Write overloaded function for eta and J
  @todo Make the tracker itself overloaded
  */
  VectorXcd RootTracker::NRCTracker(VectorXcd x, VectorXcd y, std::function<VectorXcd (VectorXcd)> f, \
    std::function<MatrixXcd (VectorXcd)> Jfy, double eps /*= pow(10, -10)*/){
    VectorXcd q(x.size()+y.size()), tempy(y.size()), dy(y.size()), fval;
    MatrixXcd Jval;
    int loopcounter = 0;
    q << y, x;
    fval = f(q);
    tempy = y;
    while (fval.norm()>=eps) {
  		if (loopcounter>=100){
        std::cout << "NRCTracker::Warning::Max iterations reached. Convergence not achived for the input tolerance, eps" << std::endl;
        return tempy;
      }
  		loopcounter++;
      Jval = Jfy(q);
  		// dy = linearSolve(Jval, fval);
      // Using Eigens linear solver so that it handles
      // the complex elements
      dy = Jval.householderQr().solve(fval);
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
  @param xprev The set of input/known variables at the previous tracking step
  @param x The set of input/known variables at the current tracking step
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
      if (fval.norm()>=eps){
        tempy = RootTracker::NRTracker(x, tempy, f, Jfy);
      }
    }
    return tempy;
  }


  /*!
  The DMCTracker uses the Davindenkos' integration method to find
  the solutions satisfied by the constrain equations.
  */
  VectorXcd RootTracker::DMCTracker(VectorXcd xprev, VectorXcd x, VectorXcd y, std::function<MatrixXcd (VectorXcd)> Jfx, \
     std::function<MatrixXcd (VectorXcd)> Jfy, double eps /*= 0*/, std::function<VectorXcd (VectorXcd)> f /*= NULL*/){
    VectorXcd q(x.size()+y.size()), tempy(y.size()), fval;
    tempy = y;
    q << y, x;
    tempy += -(Jfy(q).inverse()*Jfx(q))*(x-xprev);
    // Check to apply an NR step
    if (eps > 0){
      q << tempy, x;
      VectorXcd fval;
      fval = f(q);
      if (fval.norm()>=eps){
        tempy = RootTracker::NRCTracker(x, tempy, f, Jfy);
      }
    }
    return tempy;
  }



  /*!
  The NNTracker uses the nearest neighbour method to identify the roots
  belonging to a required branch. The output is the
  selected root at each tracking step. This method assumes the existance of a
  solver method, `Solve`, which computes all the roots for given input
  variables, `x`. The function expects the solutions to be ordered with all
  the reals together and the variables belonging to S1 together.

  @param ys The initiation of the known root of the required branch
  @param ysols All the solutions obtained by the `Solve` method used
  @param index The index upto which the reals are present.

  @todo The method currently deals only with variables belonging to R and S1.
  @todo Treating the Rodriques parameters as locally belonging to R.
  @todo Add assertions for cols of ys and ysols to be same.
  */
VectorXd RootTracker::NNTracker(VectorXd ys, MatrixXd ysols, int index){
  VectorXd rys, sys, rysols, sysols, dist, dmaxlist, rdist, s1dist;
  double dmax;
  rys = ys.head(index);
  sys = ys.tail(ys.size()-index);
  dmaxlist = VectorXd::Constant(ysols.rows(), std::numeric_limits<double>::infinity());

  for (int i = 0; i < ysols.rows(); i++) {
    rysols = ysols.row(i).head(index);
    sysols = ysols.row(i).tail(ys.size()-index);

    dist = (rys-rysols).cwiseAbs(), s1Dist(sys, sysols);

    dmax = dist.maxCoeff();
    dmaxlist(i) = dmax;
  }

  MatrixXf::Index minIndex;
  dmaxlist.minCoeff(&minIndex);
  return ysols.row(minIndex);
}



  /*!
  The NNCTracker uses the nearest neighbour method to identify the roots
  belonging to a required branch.
  */
VectorXcd RootTracker::NNCTracker(VectorXcd ys, MatrixXcd ysols, int index){
  VectorXcd rys, sys, rysols, sysols, rdist, s1dist;
  VectorXd dist, dmaxlist;
  double dmax;
  rys = ys.head(index);
  sys = ys.tail(ys.size()-index);
  dmaxlist = VectorXd::Constant(ysols.rows(), std::numeric_limits<double>::infinity());

  for (int i = 0; i < ysols.rows(); i++) {
    rysols = ysols.row(i).head(index);
    sysols = ysols.row(i).tail(ys.size()-index);

    dist = (rys-rysols).cwiseAbs(), s1CDist(sys, sysols);

    dmax = dist.maxCoeff();
    dmaxlist(i) = dmax;
  }

  MatrixXf::Index minIndex;
  dmaxlist.minCoeff(&minIndex);
  return ysols.row(minIndex);
}


/*!
The trackAllBranches routine tracks one step in all the branches given the next steps input variable
and an initial guess of the unknown variables using the relating function f and its Jacobian Jfy.
*/

MatrixXd RootTracker::trackAllBranches(VectorXd x, MatrixXd y, std::function<VectorXd (VectorXd)> f, \
  std::function<MatrixXd (VectorXd)> Jfy){
    MatrixXd roots(y.rows(), y.cols());
    // roots = MatrixXd::Zero(y.rows(), y.cols());
    roots = y;
    // std::cout << "Input y to trackAllBranches" << y << '\n';
    for (size_t i = 0; i < y.rows(); i++) {
      roots.row(i) = RootTracker::NRTracker(x, y.row(i), f, Jfy);
    }
    // std::cout << roots << '\n';
    // std::cout << " " << '\n';
    return roots;
  }

/*!
The SEI uses the distance between different branches of the solutions to identify
and estimate the singular configuration of the given set of non-linear equations.
Further, it uses a quadratic extrapolation scheme to estimate the singular configuration.
This function needs the computation of all the roots (currently handles only real roots).

@param allroots Solutions of all the branches of the non-linear equations
@param alpha Parameter value
@param selectedroot The index of the selected branch (Use index starting with 1)
@param computeXfromParam Compute the values of the input variables based on
the path parametrisation
@param alphahist Input vector (global to main) to store history of alphas for
estimating the singular configuration
@param disthist Input vector (global to main) to store history of distances for
estimating the singular configuration
@param f Set of non-linear equations
@param Jft Jacobian matrix of the set of equations
@param computeqExtfromParam Function which produces the extended configuration values
given the path parameter

@todo Optionally Bertini can be used to compute all the roots. Note that the method is
limited to single parameter paths.
@todo Check the alpha selection with Aditya
@todo Assumes a quadratic interpolation, but can be made to accept an nth degree polynomial
*/

int RootTracker::SEI(MatrixXd allroots, double alpha, int selectedroot, \
  std::function<VectorXd (double)> computeXfromParam, Ref<VectorXd> alphahist, \
  Ref<MatrixXd> disthist, std::function<VectorXd (VectorXd)> f, \
  std::function<MatrixXd (VectorXd)> Jfy, std::function<VectorXd (double)> computeqExtfromParam){
  // Convert branch number to solution index
  selectedroot = selectedroot-1;
  VectorXd alphaest(2), x(6), dist(allroots.rows()), graddist(allroots.rows()), preddist(allroots.rows()), coeff(3);
  MatrixXd currentroots(allroots.rows(), allroots.cols());
  currentroots = allroots;
  double stepsize, aa, bb, cc;
  assert(pushHist(alphahist, alpha) && "Pushing alpha to history unsuccessful");
  dist = computeDist(currentroots, selectedroot);
  assert(pushHist(disthist, dist) && "Pushing alpha to history unsuccessful");
  stepsize = (alphahist(0)-alphahist(1));
  graddist = (disthist.row(0)-disthist.row(1))/stepsize;
  preddist = ((disthist.row(0)).transpose())+graddist*stepsize;
  for (size_t i = 0; i < preddist.size(); i++) {
    if (preddist(i) < 0 &&  i!= selectedroot){
      std::cout << "Approaching singularity. Branches " << selectedroot+1 << " and " << i+1 << " are going to merge!" << '\n';
      std::cout << "The corresponding alpha is " << alphahist(0) << '\n';
      // Estimating the singular configuration
      // This assumes a quadratic extrapolation
      coeff = findExtrapCoeffs(alphahist, disthist.col(i));
      aa = coeff(0);
      bb = coeff(1);
      cc = coeff(2);
      alphaest(0) = (-bb+pow(pow(bb, 2)-4*aa*cc, 0.5))/(2*aa);
      alphaest(1) = (-bb-pow(pow(bb, 2)-4*aa*cc, 0.5))/(2*aa);
      // Returning the configuration for both the possibilities
      if(abs(alphaest(0)-alphahist(0))<abs(alphaest(1)-alphahist(0))){
        std::cout << "The estimated singular configuration is " << '\n' << \
        computeqExtfromParam(alphaest(0)) << '\n';
      }
      else{
        std::cout << "The estimated singular configuration is " << '\n' << \
        computeqExtfromParam(alphaest(1)) << '\n';
      }
      return 1;
    }
  }
  return 0;
}

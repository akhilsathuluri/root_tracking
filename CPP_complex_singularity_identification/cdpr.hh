#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include<iostream>
#include<eigen3/Eigen/Dense>

using namespace Eigen;

/////////////////////////////////////Variables as Real///////////////////////////////////////

/*!
The funcion etaext is the set of constraint equations of an SRSPM formulated in
the extended configuration space of the manipulator. This is for the example
problem demosntrating the use of the root-trackers.
The formulation of the symbolic equations can be found in the accompanying
Mathematica notebook named, `SRSPM_root_tracking.nb`.

@param q Takes in the 24-dimensional extended-configuration-space variables as input
*/
VectorXd eta(VectorXd q){
  VectorXd eta(6);
    double x = q(0), y = q(1), z = q(2), c1 = q(3), c2 = q(4), c3 = q(5), l1 = q(6), l2 = q(7), l3 = q(8);

  eta << 1 + 2*x + 4*c3*y - 4*c2*z + 4*c1*(c2*y + c3*z) + pow(c3,2) - 2*x*pow(c3,2) - pow(l1,2) - pow(c3,2)*pow(l1,2) + pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) +
    pow(c3,2)*pow(y,2) + pow(z,2) + pow(c3,2)*pow(z,2) + pow(c2,2)*(1 - 2*x - pow(l1,2) + pow(x,2) + pow(y,2) + pow(z,2)) +
    pow(c1,2)*(1 + 2*x - pow(l1,2) + pow(x,2) + pow(y,2) + pow(z,2)),
   100 + 40*c3 - 22*x - 4*c3*x + 2*y - 4*c3*y + 4*c2*(1 + c3)*z + 4*c1*(c2*(-10 + x - y) + z - c3*z) + 100*pow(c3,2) - 18*x*pow(c3,2) - 2*y*pow(c3,2) + pow(l1,2) +
    pow(c3,2)*pow(l1,2) + pow(c2,2)*(2*(50 - 9*x + y) + pow(l1,2) - pow(l2,2)) + pow(c1,2)*(-2*(-50 + 11*x + y) + pow(l1,2) - pow(l2,2)) - pow(l2,2) -
    pow(c3,2)*pow(l2,2),144 - 2*x - 24*y - 4*c3*y + 2*z + 4*c2*(x + c3*(-12 + y) + z) - 4*c1*(-12 + y + c2*y + c3*(-x + z)) + 144*pow(c3,2) + 2*x*pow(c3,2) -
    24*y*pow(c3,2) + 2*z*pow(c3,2) + pow(l1,2) + pow(c3,2)*pow(l1,2) + pow(c2,2)*(2*(72 + x - 12*y - z) + pow(l1,2) - pow(l3,2)) +
    pow(c1,2)*(-2*(-72 + x + 12*y + z) + pow(l1,2) - pow(l3,2)) - pow(l3,2) - pow(c3,2)*pow(l3,2),
   2*(60 + 54*x - 110*c3*x + 55*y - 6*x*y - 2*c3*x*y - 120*pow(c3,2) - 12*x*pow(c3,2) + 10*y*pow(c3,2) - 130*x*pow(c3,3) - 2*x*y*pow(c3,3) + 60*pow(c3,4) -
      66*x*pow(c3,4) - 65*y*pow(c3,4) + 6*x*y*pow(c3,4) - 6*pow(x,2) + 12*c3*pow(x,2) + 12*pow(c3,3)*pow(x,2) + 6*pow(c3,4)*pow(x,2) +
      2*pow(c2,2)*(5*(-12 + y)*y*pow(c3,2) + x*(60 - c3*(55 + y) + 6*y*pow(c3,2)) + 6*(-1 + c3)*pow(x,2)) -
      2*pow(c2,3)*(-6 + 6*x + 49*y + 5*c3*(-1 + x)*(1 + y) - 5*pow(y,2)) -
      2*c2*(6 + 6*x + 61*y + c3*(115 + 7*y + x*(-7 + 5*y)) + 5*(-1 + x)*(-1 + y)*pow(c3,3) + pow(c3,2)*(6 - 6*x + 59*y - 5*pow(y,2)) - 5*pow(y,2)) - 5*pow(y,2) +
      5*pow(c3,4)*pow(y,2) + pow(c1,4)*(-6*x*(9 + y) + 6*pow(x,2) - 5*(12 - 13*y + pow(y,2))) + pow(c2,4)*(6*x*(11 + y) - 6*pow(x,2) + 5*(-12 - 11*y + pow(y,2))) +
      2*pow(c1,3)*(5*(-1 + y) + x*(-5 - 65*c2 + 5*y + 11*c2*y) + 6*c2*pow(x,2) + c3*(6 + 6*x - 71*y + 5*pow(y,2))) +
      2*pow(c1,2)*(c3*x*(-65 + 6*x - y) + (60 - 6*x - 5*y)*y + (60 + 6*x - 5*y)*pow(c2,2) + 6*(-10 + x)*x*pow(c3,2) +
         c2*(6 + 6*x - 71*y + c3*(-125 + 17*x + 17*y - 5*x*y) + 5*pow(y,2))) +
      2*c1*(5*(1 + x)*(1 + y) + x*(6*x + 11*(-5 + y))*pow(c2,3) + (-125 + 7*x - 7*y + 5*x*y)*pow(c3,2) + c2*x*(-55 + 6*x + 11*y + (-65 + 6*x + 11*y)*pow(c3,2)) -
         c3*(6 + 6*x + 61*y - 5*pow(y,2)) + pow(c3,3)*(-6 + 6*x - 59*y + 5*pow(y,2)) + pow(c2,2)*(-115 + 17*x + 17*y + 5*x*y + c3*(6 - 6*x - 49*y + 5*pow(y,2))))),
   pow(c1,4)*(10*(-12 + y)*z + x*(1 - 11*y + 11*z) + pow(x,2)) - pow(c2,4)*(10*(-12 + y)*z + x*(-1 + 9*y + 13*z) + pow(x,2)) +
    2*pow(c2,3)*(x*(-120 + 10*y - z) - 10*y*z + c3*(-10 + x)*(-12 + 12*x + y + z) + 13*pow(x,2)) +
    2*pow(c1,3)*(-10*(12 + z + y*(-1 + c3*z)) + x*(-108 - y + z + c3*(120 - 10*y + z) + c2*(y - 2*(5 + 6*z))) + (12 + c2 - 11*c3)*pow(x,2)) -
    (1 + pow(c3,2))*(10*(-12 + y)*z*(-1 + pow(c3,2)) - x*(1 + 11*y + 13*z + 2*c3*(10 + y + 12*z) + (1 + 9*y - 11*z)*pow(c3,2)) + (-1 + 2*c3 + pow(c3,2))*pow(x,2)) +
    2*pow(c2,2)*(x*(1 + y) + c3*(-20*y + x*(10 + y + 12*z) - pow(x,2)) - pow(c3,2)*(20 + 10*(-12 + y)*z + 3*x*(-7 + 4*z) + pow(x,2))) +
    2*c2*(x*(-120 + 10*y - z) - 10*y*z + (-10 + x)*(-12 + 12*x + y + z)*pow(c3,3) + 13*pow(x,2) + c3*(10*(-12 + y - z) + x*(-132 + y + z) + 12*pow(x,2)) +
       pow(c3,2)*(x*(-140 + 10*y - z) - 10*y*z + 13*pow(x,2))) + 2*c1*
     ((-10 + x)*(12 + 12*x - y + z) + x*(x + y - 2*(5 + 6*z))*pow(c2,3) + c3*(-10*y*z + x*(100 - 10*y + z) - 11*pow(x,2)) +
       pow(c3,3)*(-10*y*z + x*(120 - 10*y + z) - 11*pow(x,2)) + pow(c3,2)*(x*(-108 - y + z) - 10*(-12 + y + z) + 12*pow(x,2)) +
       pow(c2,2)*(-10*(-12 + y + z + c3*y*z) + x*(-108 - y + z + c3*(100 - 10*y + z)) + (12 - 11*c3)*pow(x,2)) +
       c2*(-40*c3 - 20*y + x*(-10 + y - 12*z) + pow(x,2) + pow(c3,2)*(-20*y + x*(-10 + y - 12*z) + pow(x,2)))) -
    2*pow(c1,2)*(10*(2 + 2*c3*y - (-12 + y)*z) + x*(-1 + 10*y + z)*pow(c2,2) + x*(19 - 12*z - c3*(10 + y + 12*z) + (-1 + y)*pow(c3,2)) + (-1 + c3)*pow(x,2) +
       c2*(10*y*z + x*(140 - 10*y + z) - 13*pow(x,2) - c3*(10*(-12 + y - z) + x*(-132 + y + z) + 12*pow(x,2)))),
   -2*(y*(12 + 12*x - y + 11*z + c3*(120 - 11*x - 10*y + z)) + c2*((-12 + y)*(-10 + y + 10*z) + x*(-12 + y + 12*z)))*pow(c1,3) +
    (-(x*(y + 12*(-1 + z))) + (-12 + y)*(-1 + 11*y + z))*pow(c1,4) - 2*y*(-120 + 13*x + 10*y + c3*(-12 + 12*x + y - 9*z) - z)*pow(c2,3) +
    (x*(y + 12*(-1 + z)) + (-12 + y)*(-1 + 9*y + z))*pow(c2,4) - 2*c2*y*(-120 + 13*x + 10*y + c3*(-12 + 12*x + y - 9*z) - z)*(1 + pow(c3,2)) +
    (1 + pow(c3,2))*(-((-12 + y)*(1 + 11*y + z)) + x*(12 - y + 12*z) + 2*c3*(x*(y - 12*(1 + z)) - (-12 + y)*(y + 10*(1 + z))) +
       (-((-12 + y)*(1 + 9*y + z)) + x*(y - 12*(1 + z)))*pow(c3,2)) +
    2*pow(c2,2)*(-12 - 13*y + 12*x*z + (-1 + x)*(12 + y)*pow(c3,2) + c3*(x*(12 + y - 12*z) + 120*(-1 + z) - 2*y*(11 + 5*z) - pow(y,2)) - pow(y,2)) +
    2*pow(c1,2)*(-((1 + x)*(12 + y)) + c2*y*(120 - 13*x - 10*y - c3*(-12 + 12*x + y - 9*z) + z) + (-12 + y)*(-1 + 10*y + z)*pow(c2,2) +
       c3*(x*(12 + y - 12*z) + 120*(-1 + z) - 2*y*(11 + 5*z) - pow(y,2)) + pow(c3,2)*(11*y - 12*(1 + x*z) + pow(y,2))) -
    2*c1*(y*(12 + 12*x - y + 11*z + c3*(120 - 11*x - 10*y + z))*pow(c2,2) + ((-12 + y)*(y + 10*(-1 + z)) + x*(y + 12*(-1 + z)))*pow(c2,3) -
       y*(-12 - 12*x + y + c3*(-120 + 11*x + 10*y - z) - 11*z)*(1 + pow(c3,2)) + c2*(1 + pow(c3,2))*(-120*(1 + z) + 2*y*(1 + 5*z) + x*(y + 12*(1 + z)) + pow(y,2)));
  return eta;
}


/*!
The funcion Jetaextphi is the constraint Jacobian matrix of an SRSPM formulated in
the extended configuration space of the manipulator. This is for the example
problem demonstrating the use of the root-trackers.
The formulation of the symbolic equations can be found in the accompanying
Mathematica notebook named, `SRSPM_root_tracking.nb`.

@param q Takes in the 24-dimensional extended-configuration-space variables as input
*/
MatrixXd Jetaphi(VectorXd q){
  MatrixXd Jetaphi(6, 6);
    double x = q(0), y = q(1), z = q(2), c1 = q(3), c2 = q(4), c3 = q(5), l1 = q(6), l2 = q(7), l3 = q(8);

    Jetaphi << 4*(c2*y + c3*z) + 2*c1*(1 + 2*x - pow(l1,2) + pow(x,2) + pow(y,2) + pow(z,2)),4*c1*y - 4*z + 2*c2*(1 - 2*x - pow(l1,2) + pow(x,2) + pow(y,2) + pow(z,2)),
   4*(y + c1*z) + 2*c3*(1 - 2*x - pow(l1,2) + pow(x,2) + pow(y,2) + pow(z,2)),2*(1 + x + (1 + x)*pow(c1,2) + (-1 + x)*pow(c2,2) - pow(c3,2) + x*pow(c3,2)),
   2*(2*c1*c2 + 2*c3 + y + y*pow(c1,2) + y*pow(c2,2) + y*pow(c3,2)),2*(-2*c2 + 2*c1*c3 + z + z*pow(c1,2) + z*pow(c2,2) + z*pow(c3,2)),
   4*(c2*(-10 + x - y) + z - c3*z) + 2*c1*(100 - 22*x - 2*y + pow(l1,2) - pow(l2,2)),
   4*c1*(-10 + x - y) + 4*(1 + c3)*z + 2*c2*(2*(50 - 9*x + y) + pow(l1,2) - pow(l2,2)),
   -2*(2*(-10 + x + y + c1*z - c2*z) + c3*(2*(-50 + 9*x + y) - pow(l1,2) + pow(l2,2))),-2*(11 - 2*c1*c2 + 2*c3 + 11*pow(c1,2) + 9*pow(c2,2) + 9*pow(c3,2)),
   -2*(-1 + 2*c1*c2 + 2*c3 + pow(c1,2) - pow(c2,2) + pow(c3,2)),4*(c1 + c2 - c1*c3 + c2*c3),
   -4*(-12 + y + c2*y + c3*(-x + z)) + 2*c1*(-2*(-72 + x + 12*y + z) + pow(l1,2) - pow(l3,2)),
   4*(x + c3*(-12 + y) - c1*y + z) + 2*c2*(144 + 2*x - 24*y - 2*z + pow(l1,2) - pow(l3,2)),
   4*c1*x + 4*c2*(-12 + y) - 4*y - 4*c1*z + 2*c3*(2*(72 + x - 12*y + z) + pow(l1,2) - pow(l3,2)),2*(-1 + 2*c2 + 2*c1*c3 - pow(c1,2) + pow(c2,2) + pow(c3,2)),
   -4*(6 + c1 + c1*c2 + c3 - c2*c3 + 6*pow(c1,2) + 6*pow(c2,2) + 6*pow(c3,2)),-2*(-1 - 2*c2 + 2*c1*c3 + pow(c1,2) + pow(c2,2) - pow(c3,2)),
   4*(5 - 6*c3 + 5*x - 6*c3*x + 5*y - 61*c3*y + 5*x*y + x*(6*x + 11*(-5 + y))*pow(c2,3) - 125*pow(c3,2) + 7*x*pow(c3,2) - 7*y*pow(c3,2) + 5*x*y*pow(c3,2) +
      c2*x*(6*x + 11*(-5 + y) + (-65 + 6*x + 11*y)*pow(c3,2)) - 6*pow(c3,3) + 6*x*pow(c3,3) - 59*y*pow(c3,3) + 5*c3*pow(y,2) + 5*pow(c3,3)*pow(y,2) +
      2*pow(c1,3)*(-6*x*(9 + y) + 6*pow(x,2) - 5*(12 - 13*y + pow(y,2))) +
      3*pow(c1,2)*(5*(-1 + y) + x*(-5 - 65*c2 + 5*y + 11*c2*y) + 6*c2*pow(x,2) + c3*(6 + 6*x - 71*y + 5*pow(y,2))) +
      pow(c2,2)*(-115 + 17*y + x*(17 + 5*y) + c3*(6 - 6*x - 49*y + 5*pow(y,2))) +
      2*c1*(c3*x*(-65 + 6*x - y) + (60 - 6*x - 5*y)*y + (60 + 6*x - 5*y)*pow(c2,2) + 6*(-10 + x)*x*pow(c3,2) +
         c2*(6 + 6*x - 71*y + c3*(-125 + 17*x + 17*y - 5*x*y) + 5*pow(y,2)))),
   4*(-6 - 115*c3 - 6*x + 7*c3*x - 61*y - 7*c3*y - 5*c3*x*y + x*(-65 + 6*x + 11*y)*pow(c1,3) - 6*pow(c3,2) + 6*x*pow(c3,2) - 59*y*pow(c3,2) - 5*pow(c3,3) +
      5*x*pow(c3,3) + 5*y*pow(c3,3) - 5*x*y*pow(c3,3) + 2*c2*(5*(-12 + y)*y*pow(c3,2) + x*(60 - c3*(55 + y) + 6*y*pow(c3,2)) + 6*(-1 + c3)*pow(x,2)) -
      3*pow(c2,2)*(-6 + 6*x + 49*y + 5*c3*(-1 + x)*(1 + y) - 5*pow(y,2)) - 2*pow(c2,3)*(60 + 55*y - 6*x*(11 + y) + 6*pow(x,2) - 5*pow(y,2)) + 5*pow(y,2) +
      5*pow(c3,2)*pow(y,2) + pow(c1,2)*(6 + 6*x + 2*c2*(60 + 6*x - 5*y) - 71*y + c3*(-125 + 17*x + 17*y - 5*x*y) + 5*pow(y,2)) +
      c1*(3*x*(6*x + 11*(-5 + y))*pow(c2,2) + x*(-55 + 6*x + 11*y + (-65 + 6*x + 11*y)*pow(c3,2)) +
         2*c2*(-115 + 17*x + 17*y + 5*x*y + c3*(6 - 6*x - 49*y + 5*pow(y,2))))),
   4*(-120*c3 - 55*x - 12*c3*x + 10*c3*y - x*y - (x*(65 - 12*c3*(-10 + x) - 6*x + y) + c2*(125 - 17*y + x*(-17 + 5*y)))*pow(c1,2) - 5*(-1 + x)*(1 + y)*pow(c2,3) -
      195*x*pow(c3,2) - 3*x*y*pow(c3,2) + 120*pow(c3,3) - 132*x*pow(c3,3) - 130*y*pow(c3,3) + 12*x*y*pow(c3,3) + 6*pow(x,2) + 18*pow(c3,2)*pow(x,2) +
      12*pow(c3,3)*pow(x,2) + pow(c2,2)*(10*c3*(-12 + y)*y + x*(-55 + (-1 + 12*c3)*y) + 6*pow(x,2)) + 10*pow(c3,3)*pow(y,2) +
      pow(c1,3)*(6 + 6*x - 71*y + 5*pow(y,2)) - c2*(115 + 7*y + x*(-7 + 5*y) + 15*(-1 + x)*(-1 + y)*pow(c3,2) - 2*c3*(-6 + 6*x - 59*y + 5*pow(y,2))) +
      c1*(-6 - 6*x - 61*y + 2*c2*c3*x*(-65 + 6*x + 11*y) + 2*c3*(-125 + 7*x - 7*y + 5*x*y) + 5*pow(y,2) + 3*pow(c3,2)*(-6 + 6*x - 59*y + 5*pow(y,2)) +
         pow(c2,2)*(6 - 6*x - 49*y + 5*pow(y,2)))),4*(27 - 55*c3 - 6*x + 12*c3*x - 3*y - c3*y + (-5 + 6*c3 + 5*y + c2*(-65 + 12*x + 11*y))*pow(c1,3) +
      3*(-9 + 2*x - y)*pow(c1,4) - (6 + 5*c3*(1 + y))*pow(c2,3) + (33 - 6*x + 3*y)*pow(c2,4) - 6*pow(c3,2) +
      pow(c1,2)*(c2*(6 + c3*(17 - 5*y)) + c3*(-65 + 12*x - y) - 6*y + 6*pow(c2,2) + 12*(-5 + x)*pow(c3,2)) +
      pow(c2,2)*(-12*(-5 + x) + c3*(-55 + 12*x - y) + 6*y*pow(c3,2)) - 65*pow(c3,3) + 12*x*pow(c3,3) - y*pow(c3,3) +
      c1*(-6*c3 + 5*(1 + y) + (17 - 6*c3 + 5*y)*pow(c2,2) + (12*x + 11*(-5 + y))*pow(c2,3) + (7 + 5*y)*pow(c3,2) +
         c2*(-55 + 12*x + 11*y + (-65 + 12*x + 11*y)*pow(c3,2)) + 6*pow(c3,3)) + c2*(-6 + c3*(7 - 5*y) + 6*pow(c3,2) - 5*(-1 + y)*pow(c3,3)) - 33*pow(c3,4) +
      6*x*pow(c3,4) + 3*y*pow(c3,4)),-2*(-55 + 6*x + 2*c3*x + 10*y - 2*(5 + (5 + 11*c2)*x + c3*(-71 + 10*y))*pow(c1,3) + (-65 + 6*x + 10*y)*pow(c1,4) -
      2*c3*(-x + 6*c3*x + 10*c3*(-6 + y))*pow(c2,2) + 2*pow(c1,2)*((6 + c3)*x + c2*(71 + c3*(-17 + 5*x) - 10*y) + 10*(-6 + y) + 5*pow(c2,2)) +
      2*(49 + 5*c3*(-1 + x) - 10*y)*pow(c2,3) + (55 - 6*x - 10*y)*pow(c2,4) - 10*pow(c3,2) + 2*x*pow(c3,3) +
      2*c2*(61 + c3*(7 + 5*x) - 10*y + (59 - 10*y)*pow(c3,2) + 5*(-1 + x)*pow(c3,3)) -
      2*c1*(5*(1 + x) + c3*(-61 + 10*y) + (17 + 5*x + c3*(-49 + 10*y))*pow(c2,2) + 11*x*pow(c2,3) + (-7 + 5*x)*pow(c3,2) + 11*c2*x*(1 + pow(c3,2)) +
         (-59 + 10*y)*pow(c3,3)) + 65*pow(c3,4) - 6*x*pow(c3,4) - 10*y*pow(c3,4)),0,
   2*(-120 - 108*x + 100*c3*x + 10*y - x*y - 10*c3*x*y - 10*z + x*z + c3*x*z - 10*c3*y*z + x*(x + y - 2*(5 + 6*z))*pow(c2,3) + 120*pow(c3,2) - 108*x*pow(c3,2) -
      10*y*pow(c3,2) - x*y*pow(c3,2) - 10*z*pow(c3,2) + x*z*pow(c3,2) + 120*x*pow(c3,3) - 10*x*y*pow(c3,3) + x*z*pow(c3,3) - 10*y*z*pow(c3,3) + 12*pow(x,2) -
      11*c3*pow(x,2) + 12*pow(c3,2)*pow(x,2) - 11*pow(c3,3)*pow(x,2) + 2*pow(c1,3)*(10*(-12 + y)*z + x*(1 - 11*y + 11*z) + pow(x,2)) +
      pow(c2,2)*(-10*(-12 + y + z + c3*y*z) + x*(-108 - y + z + c3*(100 - 10*y + z)) + (12 - 11*c3)*pow(x,2)) +
      3*pow(c1,2)*(-10*(12 + z + y*(-1 + c3*z)) + x*(-108 - y + z + c3*(120 - 10*y + z) + c2*(y - 2*(5 + 6*z))) + (12 + c2 - 11*c3)*pow(x,2)) +
      c2*(-40*c3 - 20*y + x*(-10 + y - 12*z) + pow(x,2) + pow(c3,2)*(-20*y + x*(-10 + y - 12*z) + pow(x,2))) -
      2*c1*(10*(2 + 2*c3*y - (-12 + y)*z) + x*(-1 + 10*y + z)*pow(c2,2) + x*(19 - 12*z - c3*(10 + y + 12*z) + (-1 + y)*pow(c3,2)) + (-1 + c3)*pow(x,2) +
         c2*(10*y*z + x*(140 - 10*y + z) - 13*pow(x,2) - c3*(10*(-12 + y - z) + x*(-132 + y + z) + 12*pow(x,2))))),
   2*(-120*x + 2*c2*x + 10*x*y + 2*c2*x*y - x*z - 10*y*z + x*(x + y - 2*(5 + 6*z))*pow(c1,3) - 360*x*pow(c2,2) + 30*x*y*pow(c2,2) - 3*x*z*pow(c2,2) -
      30*y*z*pow(c2,2) + 2*x*pow(c2,3) - 18*x*y*pow(c2,3) + 240*z*pow(c2,3) - 26*x*z*pow(c2,3) - 20*y*z*pow(c2,3) + (-10 + x)*(-12 + 12*x + y + z)*pow(c3,3) +
      13*pow(x,2) + 39*pow(c2,2)*pow(x,2) - 2*pow(c2,3)*pow(x,2) +
      c3*(10*(-12 + y - z) + x*(-132 + y + z) + 3*(-10 + x)*(-12 + 12*x + y + z)*pow(c2,2) + c2*(-40*y + 2*x*(10 + y + 12*z) - 2*pow(x,2)) + 12*pow(x,2)) -
      pow(c3,2)*(10*y*z + x*(140 - 10*y + z) - 13*pow(x,2) + 2*c2*(20 + 10*(-12 + y)*z + 3*x*(-7 + 4*z) + pow(x,2))) +
      pow(c1,2)*(-10*y*z - x*(140 - 10*y + z + 2*c2*(-1 + 10*y + z)) + 13*pow(x,2) + c3*(10*(-12 + y - z) + x*(-132 + y + z) + 12*pow(x,2))) +
      c1*(-40*c3 - 10*x - 20*y + x*y - 12*x*z + 3*x*(-10 + x + y - 12*z)*pow(c2,2) + pow(x,2) + pow(c3,2)*(-20*y + x*(-10 + y - 12*z) + pow(x,2)) -
         2*c2*(x*(108 + y + c3*(-100 + 10*y - z) - z) + 10*(-12 + y + z + c3*y*z) + (-12 + 11*c3)*pow(x,2)))),
   2*(10*x + 2*c3*x + x*y + 20*c3*x*y + 12*x*z + 2*c3*x*z + (-10 + x)*(-12 + 12*x + y + z)*pow(c2,3) + 30*x*pow(c3,2) + 3*x*y*pow(c3,2) + 36*x*z*pow(c3,2) +
      2*x*pow(c3,3) + 18*x*y*pow(c3,3) + 240*z*pow(c3,3) - 22*x*z*pow(c3,3) - 20*y*z*pow(c3,3) + pow(c1,3)*(-10*y*z + x*(120 - 10*y + z) - 11*pow(x,2)) +
      pow(c1,2)*(-20*y + x*(10 - 2*c3*(-1 + y) + y + 12*z) - pow(x,2)) - pow(x,2) - 3*pow(c3,2)*pow(x,2) - 2*pow(c3,3)*pow(x,2) -
      pow(c2,2)*(x*(-10 - y + c1*(-100 + 10*y - z) - 12*z) + 10*y*(2 + c1*z) + (1 + 11*c1)*pow(x,2) + 2*c3*(20 + 10*(-12 + y)*z + 3*x*(-7 + 4*z) + pow(x,2))) -
      c1*(x*(-100 + 10*y - z) + 10*y*z + c3*(2*x*(108 + y - z) + 20*(-12 + y + z) - 24*pow(x,2)) + 11*pow(x,2) +
         pow(c3,2)*(30*x*y + 30*y*z - 3*x*(120 + z) + 33*pow(x,2))) +
      c2*(-120 - 132*x + 10*y + x*y - 10*z + x*z + 3*(-10 + x)*(-12 + 12*x + y + z)*pow(c3,2) + 12*pow(x,2) +
         pow(c1,2)*(10*(-12 + y - z) + x*(-132 + y + z) + 12*pow(x,2)) + c3*(20*x*y - 20*y*z - 2*x*(140 + z) + 26*pow(x,2)) +
         2*c1*(-20 + c3*(-20*y + x*(-10 + y - 12*z) + pow(x,2))))),
   2*(-108 + 24*x - y + c2*(-10 + 2*x + y - 12*z) + z + c3*(120 - 22*x - 10*y + z))*pow(c1,3) + (1 + 2*x - 11*y + 11*z)*pow(c1,4) +
    2*(-120 + 26*x + 10*y - z + c3*(-132 + 24*x + y + z))*pow(c2,3) - (-1 + 2*x + 9*y + 13*z)*pow(c2,4) +
    2*pow(c1,2)*(-19 + 2*x + 12*z + c3*(10 - 2*x + y + 12*z) + c2*(-140 + 26*x + 10*y - z + c3*(-132 + 24*x + y + z)) - (-1 + 10*y + z)*pow(c2,2) -
       (-1 + y)*pow(c3,2)) + 2*pow(c2,2)*(1 + y + c3*(10 - 2*x + y + 12*z) + (21 - 2*x - 12*z)*pow(c3,2)) -
    (1 + pow(c3,2))*(-1 - 2*x + 4*c3*x - 11*y - 13*z - 2*c3*(10 + y + 12*z) + (-1 + 2*x - 9*y + 11*z)*pow(c3,2)) +
    2*c1*(-108 + 24*x - y + z + c3*(100 - 22*x - 10*y + z) + (-108 + 24*x - y + z + c3*(100 - 22*x - 10*y + z))*pow(c2,2) + (-10 + 2*x + y - 12*z)*pow(c2,3) +
       (-108 + 24*x - y + z)*pow(c3,2) + c2*(-10 + 2*x + y - 12*z)*(1 + pow(c3,2)) + (120 - 22*x - 10*y + z)*pow(c3,3)) +
    2*c2*(-120 + 26*x + 10*y - z + c3*(-132 + 24*x + y + z) + (-140 + 26*x + 10*y - z)*pow(c3,2) + (-132 + 24*x + y + z)*pow(c3,3)),
   2*(10 + (-1 + c2 - 10*c3)*x - 10*c3*z)*pow(c1,3) + (-11*x + 10*z)*pow(c1,4) + 2*(c3*(-10 + x) + 10*(x - z))*pow(c2,3) - (9*x + 10*z)*pow(c2,4) +
    pow(c1,2)*(2*c3*(-20 + x + c2*(10 + x)) - 20*(-1 + c2)*(c2*x + z) - 2*x*pow(c3,2)) + 2*pow(c2,2)*(c3*(-20 + x) + x - 10*z*pow(c3,2)) +
    (1 + pow(c3,2))*(-10*z*(-1 + pow(c3,2)) + x*(11 + 2*c3 + 9*pow(c3,2))) + 2*c2*(c3*(10 + x) + 10*(x - z) + 10*(x - z)*pow(c3,2) + (-10 + x)*pow(c3,3)) +
    2*c1*(10 - x - 10*c3*(x + z) - (10 + x + 10*c3*x + 10*c3*z)*pow(c2,2) + x*pow(c2,3) - (10 + x)*pow(c3,2) + c2*(-20 + x)*(1 + pow(c3,2)) - 10*(x + z)*pow(c3,3)),
   (1 + pow(c1,2) + pow(c2,2) + pow(c3,2))*(-120 + 13*x + 24*c3*x + 2*c2*(c3*(-10 + x) - x - 10*y) + 10*y + c1*(2*(1 - 12*c2 + c3)*x - 20*(1 + c3*y)) +
      (11*x + 10*(-12 + y))*pow(c1,2) + (120 - 13*x - 10*y)*pow(c2,2) + 120*pow(c3,2) - 11*x*pow(c3,2) - 10*y*pow(c3,2)),
   -2*(3*(y*(12 + 12*x - y + 11*z + c3*(120 - 11*x - 10*y + z)) + c2*((-12 + y)*(-10 + y + 10*z) + x*(-12 + y + 12*z)))*pow(c1,2) +
      2*(x*(y + 12*(-1 + z)) - (-12 + y)*(-1 + 11*y + z))*pow(c1,3) + y*(12 + 12*x - y + 11*z + c3*(120 - 11*x - 10*y + z))*pow(c2,2) +
      ((-12 + y)*(y + 10*(-1 + z)) + x*(y + 12*(-1 + z)))*pow(c2,3) - y*(-12 - 12*x + y + c3*(-120 + 11*x + 10*y - z) - 11*z)*(1 + pow(c3,2)) +
      c2*(1 + pow(c3,2))*(-120*(1 + z) + 2*y*(1 + 5*z) + x*(y + 12*(1 + z)) + pow(y,2)) +
      2*c1*((1 + x)*(12 + y) + c2*y*(-120 + 13*x + 10*y + c3*(-12 + 12*x + y - 9*z) - z) - (-12 + y)*(-1 + 10*y + z)*pow(c2,2) +
         c3*(120 + 22*y - x*(12 + y - 12*z) - 120*z + 10*y*z + pow(y,2)) - pow(c3,2)*(11*y - 12*(1 + x*z) + pow(y,2)))),
   -2*((y*(-120 + 13*x + 10*y + c3*(-12 + 12*x + y - 9*z) - z) - 2*c2*(-12 + y)*(-1 + 10*y + z))*pow(c1,2) +
      ((-12 + y)*(y + 10*(-1 + z)) + x*(y + 12*(-1 + z)))*pow(c1,3) + 3*y*(-120 + 13*x + 10*y + c3*(-12 + 12*x + y - 9*z) - z)*pow(c2,2) -
      2*(x*(y + 12*(-1 + z)) + (-12 + y)*(-1 + 9*y + z))*pow(c2,3) + y*(-120 + 13*x + 10*y + c3*(-12 + 12*x + y - 9*z) - z)*(1 + pow(c3,2)) +
      2*c2*(12 + 13*y - 12*x*z - (-1 + x)*(12 + y)*pow(c3,2) + pow(y,2) + c3*(120 + 22*y - x*(12 + y - 12*z) - 120*z + 10*y*z + pow(y,2))) +
      c1*(2*c2*y*(12 + 12*x - y + 11*z + c3*(120 - 11*x - 10*y + z)) + 3*((-12 + y)*(y + 10*(-1 + z)) + x*(y + 12*(-1 + z)))*pow(c2,2) +
         (1 + pow(c3,2))*(-120*(1 + z) + 2*y*(1 + 5*z) + x*(y + 12*(1 + z)) + pow(y,2)))),
   2*(120 + 24*c3 - 12*x + 2*y + 238*c3*y + x*y + 120*z + 24*c3*z - 12*x*z - 10*y*z - 2*c3*y*z + y*(-120 + 11*x + 10*y - z)*pow(c1,3) -
      y*(-12 + 12*x + y - 9*z)*pow(c2,3) + 360*pow(c3,2) - 36*x*pow(c3,2) + 6*y*pow(c3,2) + 3*x*y*pow(c3,2) + 360*z*pow(c3,2) - 36*x*z*pow(c3,2) -
      30*y*z*pow(c3,2) - c2*y*(-12 + 12*x + y + c3*(-240 + 26*x + 20*y - 2*z) - 9*z + 3*(-12 + 12*x + y - 9*z)*pow(c3,2)) + 24*pow(c3,3) - 24*x*pow(c3,3) +
      214*y*pow(c3,3) + 2*x*y*pow(c3,3) + 24*z*pow(c3,3) - 24*x*z*pow(c3,3) - 2*y*z*pow(c3,3) +
      pow(c2,2)*(-120 - 22*y + 2*c3*(-1 + x)*(12 + y) + x*(12 + y - 12*z) + 120*z - 10*y*z - pow(y,2)) - pow(y,2) - 20*c3*pow(y,2) - 3*pow(c3,2)*pow(y,2) -
      18*pow(c3,3)*pow(y,2) - pow(c1,2)*(120 + 22*y - 12*c2*y - 120*z + 10*y*z - 9*c2*y*z + x*(-12 - y + 12*c2*y + 12*z) + pow(y,2) + c2*pow(y,2) -
         2*c3*(11*y - 12*(1 + x*z) + pow(y,2))) + c1*(y*(-120 + 11*x + 10*y - z)*pow(c2,2) +
         y*(-120 + 11*x + 10*y - z - 2*c3*(12 + 12*x - y + 11*z) + (-360 + 33*x + 30*y - 3*z)*pow(c3,2)) -
         2*c2*c3*(-120*(1 + z) + 2*y*(1 + 5*z) + x*(y + 12*(1 + z)) + pow(y,2)))),
   -2*((12 - 11*c3)*y + c2*(y + 12*(-1 + z)))*pow(c1,3) - (y + 12*(-1 + z))*pow(c1,4) - 2*(13 + 12*c3)*y*pow(c2,3) + (y + 12*(-1 + z))*pow(c2,4) +
    2*pow(c2,2)*(c3*(12 + y - 12*z) + 12*z + (12 + y)*pow(c3,2)) - 2*pow(c1,2)*(12 + y + 13*c2*y + (-1 + 12*c2)*c3*y + 12*c3*(-1 + z) + 12*z*pow(c3,2)) -
    2*c1*((12 - 11*c3)*y*pow(c2,2) + (y + 12*(-1 + z))*pow(c2,3) + c2*(y + 12*(1 + z))*(1 + pow(c3,2)) + y*(12 - 11*c3 + 12*pow(c3,2) - 11*pow(c3,3))) -
    2*c2*y*(13 + 12*c3 + 13*pow(c3,2) + 12*pow(c3,3)) + (y - 12*(1 + z))*(-1 + 2*c3 + 2*pow(c3,3) + pow(c3,4)),
   -2*(12 + 12*x - 2*y + 11*z + c3*(120 - 11*x - 20*y + z) + c2*(x + 2*(-11 + y + 5*z)))*pow(c1,3) + (-133 - x + 22*y + z)*pow(c1,4) +
    2*(120 - 13*x - 20*y + z + c3*(12 - 12*x - 2*y + 9*z))*pow(c2,3) + (-109 + x + 18*y + z)*pow(c2,4) -
    2*c2*(-120 + 13*x + 20*y + c3*(-12 + 12*x + 2*y - 9*z) - z)*(1 + pow(c3,2)) + 2*pow(c2,2)*(-13 - 2*y + c3*(x - 2*(11 + y + 5*z)) + (-1 + x)*pow(c3,2)) +
    2*pow(c1,2)*(-1 - x + c3*(x - 2*(11 + y + 5*z)) + c2*(120 - 13*x - 20*y + z + c3*(12 - 12*x - 2*y + 9*z)) + (-121 + 20*y + z)*pow(c2,2) +
       (11 + 2*y)*pow(c3,2)) + (1 + pow(c3,2))*(131 - x - 22*y - z + 2*c3*(x - 2*(-1 + y + 5*z)) + (107 + x - 18*y - z)*pow(c3,2)) -
    2*c1*((12 + 12*x - 2*y + 11*z + c3*(120 - 11*x - 20*y + z))*pow(c2,2) + (x + 2*(-11 + y + 5*z))*pow(c2,3) -
       (-12 - 12*x + 2*y + c3*(-120 + 11*x + 20*y - z) - 11*z)*(1 + pow(c3,2)) + c2*(x + 2*(1 + y + 5*z))*(1 + pow(c3,2))),
   -((1 + pow(c1,2) + pow(c2,2) + pow(c3,2))*(-12 - 240*c3 - 12*x + 24*c3*x + y + 20*c3*y + 2*c1*(11 + c3)*y + 4*c1*c2*(-60 + 6*x + 5*y) - 2*c2*(y + 9*c3*y) +
        (12 + 12*x - y)*pow(c1,2) - (-12 + 12*x + y)*pow(c2,2) - 12*pow(c3,2) + 12*x*pow(c3,2) + y*pow(c3,2)));

    return Jetaphi;
}

/*!
The computeXfromParam returns input values x from the path parametrisation.
Edit this function to suite the required path paramterisation and the mapping
between the parametrisation and the actual input values. In the case of an SRSPM,
we use the IK mapping to find the leg lengths from the task space parametrisation
*/

VectorXd computeXfromParam(double alpha){
  // The first example is for the SRSPM so using the parametrisation of the SRSPM path
  VectorXd l(3);
  l << 7.5+2.5*alpha, 10.0+alpha, 9.5-1.5*alpha;
  return l;
}


/*!
The computeqExtfromParam returns the extended space coordinates for the given
configuration
*/

// VectorXd computeqExtfromParam(double alpha){
//   // The first example is for the SRSPM so using the parametrisation of the SRSPM path
//   VectorXd qx(6), qext(24);
//   double scale = 0.2;
//   qx << scale-alpha, alpha, 1.28-alpha, scale, 0, 0;
//   qext << qx, computeIK(qx);
//   return qext;
// }



// #endif CONSTRAINTS_H
#endif

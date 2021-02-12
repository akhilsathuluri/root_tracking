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
VectorXd etaext(VectorXd q){
  VectorXd eta(18);
    double x = q(0), y = q(1), z = q(2), c1 = q(3), c2 = q(4), c3 = q(5),\
    p1 = q(6), p2 = q(7), p3 = q(8), p4 = q(9), p5=q(10), p6=q(11),\
    s1=q(12), s2=q(13), s3=q(14), s4=q(15), s5=q(16), s6=q(17),\
    l1=q(18), l2=q(19), l3=q(20), l4=q(21), l5=q(22), l6=q(23);

  eta << 0.7325760476016683 - x - 0.22056834933218009*(2*c1*c2 - 2*c3)*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) -
       0.5367473272154016*(1 + pow (c1, 2) - pow (c2, 2) -
      pow (c3, 2))*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) +
   l1*cos (s1)*sin (p1),
     0.6806851948443702 - y - 1.0734946544308033*(c1*c2 + c3)*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) -
       0.22056834933218009*(1 - pow (c1, 2) + pow (c2, 2) -
      pow (c3, 2))*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) -
   l1*sin (s1),
     -z + l1*cos (p1)*cos (s1) - 1.0734946544308033*(-c2 + c1*c3)*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) -
       0.44113669866436017*(c1 + c2*c3)*
    pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1),
     0.22320264691435115 - x - 0.5751209954480265*(2*c1*c2 - 2*c3)*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) +
       0.07735586981523229*(1 + pow (c1, 2) - pow (c2, 2) -
      pow (c3, 2))*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) +
   l2*cos (s2)*sin (p2),
     0.974772064849228 - y + 0.15471173963046458*(c1*c2 + c3)*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) -
       0.5751209954480265*(1 - pow (c1, 2) + pow (c2, 2) -
      pow (c3, 2))*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) -
   l2*sin (s2),
     -z + l2*cos (p2)*cos (s2) + 0.15471173963046458*(-c2 + c1*c3)*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) -
       1.150241990896053*(c1 + c2*c3)*
    pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1),
     -0.955778694516019 - x - 0.35455264611584636*(2*c1*c2 - 2*c3)*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) +
       0.4593914574001692*(1 + pow (c1, 2) - pow (c2, 2) -
      pow (c3, 2))*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) +
   l3*cos (s3)*sin (p3),
     0.2940868700048578 - y + 0.9187829148003384*(c1*c2 + c3)*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) -
       0.35455264611584636*(1 - pow (c1, 2) + pow (c2, 2) -
      pow (c3, 2))*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) -
   l3*sin (s3),
     -z + l3*cos (p3)*cos (s3) + 0.9187829148003384*(-c2 + c1*c3)*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) -
       0.7091052922316927*(c1 + c2*c3)*
    pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1),
     -0.9557786945160192 - x + 0.35455264611584636*(2*c1*c2 - 2*c3)*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) +
       0.4593914574001694*(1 + pow (c1, 2) - pow (c2, 2) -
      pow (c3, 2))*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) +
   l4*cos (s4)*sin (p4),
     -0.2940868700048576 - y + 0.9187829148003388*(c1*c2 + c3)*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) +
       0.35455264611584636*(1 - pow (c1, 2) + pow (c2, 2) -
      pow (c3, 2))*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) -
   l4*sin (s4),
     -z + l4*cos (p4)*cos (s4) + 0.9187829148003388*(-c2 + c1*c3)*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) +
       0.7091052922316927*(c1 + c2*c3)*
    pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1),
     0.22320264691435077 - x + 0.5751209954480264*(2*c1*c2 - 2*c3)*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) +
       0.07735586981523246*(1 + pow (c1, 2) - pow (c2, 2) -
      pow (c3, 2))*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) +
   l5*cos (s5)*sin (p5),
     -0.974772064849228 - y + 0.15471173963046492*(c1*c2 + c3)*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) +
       0.5751209954480264*(1 - pow (c1, 2) + pow (c2, 2) -
      pow (c3, 2))*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) -
   l5*sin (s5),
     -z + l5*cos (p5)*cos (s5) + 0.15471173963046492*(-c2 + c1*c3)*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) +
       1.1502419908960528*(c1 + c2*c3)*
    pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1),
     0.7325760476016681 - x + 0.22056834933218022*(2*c1*c2 - 2*c3)*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) -
       0.5367473272154015*(1 + pow (c1, 2) - pow (c2, 2) -
      pow (c3, 2))*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) +
   l6*cos (s6)*sin (p6),
     -0.6806851948443704 - y - 1.073494654430803*(c1*c2 + c3)*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) +
       0.22056834933218022*(1 - pow (c1, 2) + pow (c2, 2) -
      pow (c3, 2))*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) -
   l6*sin (s6),
     -z + l6*cos (p6)*cos (s6) - 1.073494654430803*(-c2 + c1*c3)*
         pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1) +
       0.44113669866436045*(c1 + c2*c3)*
    pow (1 + pow (c1, 2) + pow (c2, 2) + pow (c3, 2), -1);
  return eta;
}


/*!
The funcion Jetaexttheta is the Jacobian matrix of an SRSPM formulated in
the extended configuration space of the manipulator. This is for the example
problem demonstrating the use of the root-trackers.
The formulation of the symbolic equations can be found in the accompanying
Mathematica notebook named, `SRSPM_root_tracking.nb`.

@param q Takes in the 24-dimensional extended-configuration-space variables as input
*/
MatrixXd Jetaexttheta(VectorXd q){
  MatrixXd Jetatheta(18, 6);
  double x = q(0), y = q(1), z = q(2), c1 = q(3), c2 = q(4), c3 = q(5),\
    p1 = q(6), p2 = q(7), p3 = q(8), p4 = q(9), p5=q(10), p6=q(11),\
    s1=q(12), s2=q(13), s3=q(14), s4=q(15), s5=q(16), s6=q(17),\
    l1=q(18), l2=q(19), l3=q(20), l4=q(21), l5=q(22), l6=q(23);

  Jetatheta << cos(s1)*sin(p1),0,0,0,0,0,-sin(s1),0,0,0,0,0,cos(p1)*cos(s1),0,0,0,0,0,0,
   cos(s2)*sin(p2),0,0,0,0,0,-sin(s2),0,0,0,0,0,cos(p2)*cos(s2),0,0,0,0,0,0,
   cos(s3)*sin(p3),0,0,0,0,0,-sin(s3),0,0,0,0,0,cos(p3)*cos(s3),0,0,0,0,0,0,
   cos(s4)*sin(p4),0,0,0,0,0,-sin(s4),0,0,0,0,0,cos(p4)*cos(s4),0,0,0,0,0,0,
   cos(s5)*sin(p5),0,0,0,0,0,-sin(s5),0,0,0,0,0,cos(p5)*cos(s5),0,0,0,0,0,0,
   cos(s6)*sin(p6),0,0,0,0,0,-sin(s6),0,0,0,0,0,cos(p6)*cos(s6);

   return Jetatheta;
 }


/*!
The funcion Jetaextphi is the constraint Jacobian matrix of an SRSPM formulated in
the extended configuration space of the manipulator. This is for the example
problem demonstrating the use of the root-trackers.
The formulation of the symbolic equations can be found in the accompanying
Mathematica notebook named, `SRSPM_root_tracking.nb`.

@param q Takes in the 24-dimensional extended-configuration-space variables as input
*/
MatrixXd Jetaextphi(VectorXd q){
  MatrixXd Jetaphi(18, 18);
    double x = q(0), y = q(1), z = q(2), c1 = q(3), c2 = q(4), c3 = q(5),\
    p1 = q(6), p2 = q(7), p3 = q(8), p4 = q(9), p5=q(10), p6=q(11),\
    s1=q(12), s2=q(13), s3=q(14), s4=q(15), s5=q(16), s6=q(17),\
    l1=q(18), l2=q(19), l3=q(20), l4=q(21), l5=q(22), l6=q(23);

    Jetaphi << -1,0,0,0.44113669866436017*c1*(2*c1*c2 - 2*c3)*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.0734946544308033*c1*(1 + pow(c1,2) - pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.0734946544308033*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.44113669866436017*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   0.44113669866436017*c2*(2*c1*c2 - 2*c3)*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.0734946544308033*c2*(1 + pow(c1,2) - pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.44113669866436017*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    1.0734946544308033*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   0.44113669866436017*(2*c1*c2 - 2*c3)*c3*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.0734946544308033*c3*(1 + pow(c1,2) - pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.44113669866436017*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    1.0734946544308033*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   l1*cos(p1)*cos(s1),0,0,0,0,0,-(l1*sin(p1)*sin(s1)),0,0,0,0,0,0,-1,0,
   2.1469893088616065*c1*(c1*c2 + c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.44113669866436017*c1*(1 - pow(c1,2) + pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.44113669866436017*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    1.0734946544308033*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   2.1469893088616065*c2*(c1*c2 + c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.44113669866436017*c2*(1 - pow(c1,2) + pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.0734946544308033*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.44113669866436017*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   2.1469893088616065*c3*(c1*c2 + c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.44113669866436017*c3*(1 - pow(c1,2) + pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.0734946544308033*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.44113669866436017*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0,0,0,0,0,0,
   -(l1*cos(s1)),0,0,0,0,0,0,0,-1,2.1469893088616065*c1*(-c2 + c1*c3)*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.8822733973287203*c1*(c1 + c2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.44113669866436017*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    1.0734946544308033*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   2.1469893088616065*c2*(-c2 + c1*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.8822733973287203*c2*(c1 + c2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.0734946544308033*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.44113669866436017*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   2.1469893088616065*c3*(-c2 + c1*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.8822733973287203*c3*(c1 + c2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.0734946544308033*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.44113669866436017*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -(l1*cos(s1)*sin(p1)),0,0,0,0,0,-(l1*cos(p1)*sin(s1)),0,0,0,0,0,-1,0,0,
   1.150241990896053*c1*(2*c1*c2 - 2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.15471173963046458*c1*(1 + pow(c1,2) - pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.15471173963046458*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    1.150241990896053*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   1.150241990896053*c2*(2*c1*c2 - 2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.15471173963046458*c2*(1 + pow(c1,2) - pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.150241990896053*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.15471173963046458*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   1.150241990896053*(2*c1*c2 - 2*c3)*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.15471173963046458*c3*(1 + pow(c1,2) - pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.150241990896053*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.15471173963046458*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0,
   l2*cos(p2)*cos(s2),0,0,0,0,0,-(l2*sin(p2)*sin(s2)),0,0,0,0,0,-1,0,
   -0.30942347926092917*c1*(c1*c2 + c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.150241990896053*c1*(1 - pow(c1,2) + pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.150241990896053*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.15471173963046458*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.30942347926092917*c2*(c1*c2 + c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.150241990896053*c2*(1 - pow(c1,2) + pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.15471173963046458*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    1.150241990896053*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.30942347926092917*c3*(c1*c2 + c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.150241990896053*c3*(1 - pow(c1,2) + pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.15471173963046458*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    1.150241990896053*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0,0,0,0,0,0,0,
   -(l2*cos(s2)),0,0,0,0,0,0,-1,-0.30942347926092917*c1*(-c2 + c1*c3)*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    2.300483981792106*c1*(c1 + c2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.150241990896053*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.15471173963046458*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.30942347926092917*c2*(-c2 + c1*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    2.300483981792106*c2*(c1 + c2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.15471173963046458*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    1.150241990896053*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.30942347926092917*c3*(-c2 + c1*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    2.300483981792106*c3*(c1 + c2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.15471173963046458*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    1.150241990896053*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0,
   -(l2*cos(s2)*sin(p2)),0,0,0,0,0,-(l2*cos(p2)*sin(s2)),0,0,0,0,-1,0,0,
   0.7091052922316927*c1*(2*c1*c2 - 2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.9187829148003384*c1*(1 + pow(c1,2) - pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.9187829148003384*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.7091052922316927*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   0.7091052922316927*c2*(2*c1*c2 - 2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.9187829148003384*c2*(1 + pow(c1,2) - pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.7091052922316927*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.9187829148003384*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   0.7091052922316927*(2*c1*c2 - 2*c3)*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.9187829148003384*c3*(1 + pow(c1,2) - pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.7091052922316927*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.9187829148003384*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0,0,
   l3*cos(p3)*cos(s3),0,0,0,0,0,-(l3*sin(p3)*sin(s3)),0,0,0,0,-1,0,
   -1.8375658296006767*c1*(c1*c2 + c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.7091052922316927*c1*(1 - pow(c1,2) + pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.7091052922316927*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.9187829148003384*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.8375658296006767*c2*(c1*c2 + c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.7091052922316927*c2*(1 - pow(c1,2) + pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.9187829148003384*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.7091052922316927*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.8375658296006767*c3*(c1*c2 + c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.7091052922316927*c3*(1 - pow(c1,2) + pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.9187829148003384*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.7091052922316927*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0,0,0,0,0,0,0,0,
   -(l3*cos(s3)),0,0,0,0,0,-1,-1.8375658296006767*c1*(-c2 + c1*c3)*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.4182105844633854*c1*(c1 + c2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.7091052922316927*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.9187829148003384*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.8375658296006767*c2*(-c2 + c1*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.4182105844633854*c2*(c1 + c2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.9187829148003384*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.7091052922316927*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.8375658296006767*c3*(-c2 + c1*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.4182105844633854*c3*(c1 + c2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.9187829148003384*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.7091052922316927*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0,0,
   -(l3*cos(s3)*sin(p3)),0,0,0,0,0,-(l3*cos(p3)*sin(s3)),0,0,0,-1,0,0,
   -0.7091052922316927*c1*(2*c1*c2 - 2*c3)*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.9187829148003388*c1*(1 + pow(c1,2) - pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.9187829148003388*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.7091052922316927*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.7091052922316927*c2*(2*c1*c2 - 2*c3)*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.9187829148003388*c2*(1 + pow(c1,2) - pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.7091052922316927*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.9187829148003388*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.7091052922316927*(2*c1*c2 - 2*c3)*c3*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.9187829148003388*c3*(1 + pow(c1,2) - pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.7091052922316927*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.9187829148003388*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0,0,0,
   l4*cos(p4)*cos(s4),0,0,0,0,0,-(l4*sin(p4)*sin(s4)),0,0,0,-1,0,
   -1.8375658296006776*c1*(c1*c2 + c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.7091052922316927*c1*(1 - pow(c1,2) + pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.7091052922316927*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.9187829148003388*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.8375658296006776*c2*(c1*c2 + c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.7091052922316927*c2*(1 - pow(c1,2) + pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.9187829148003388*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.7091052922316927*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.8375658296006776*c3*(c1*c2 + c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.7091052922316927*c3*(1 - pow(c1,2) + pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.9187829148003388*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.7091052922316927*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0,0,0,0,0,0,0,0,0,
   -(l4*cos(s4)),0,0,0,0,-1,-1.8375658296006776*c1*(-c2 + c1*c3)*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.4182105844633854*c1*(c1 + c2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.7091052922316927*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.9187829148003388*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.8375658296006776*c2*(-c2 + c1*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.4182105844633854*c2*(c1 + c2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.9187829148003388*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.7091052922316927*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.8375658296006776*c3*(-c2 + c1*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.4182105844633854*c3*(c1 + c2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.9187829148003388*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.7091052922316927*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0,0,0,
   -(l4*cos(s4)*sin(p4)),0,0,0,0,0,-(l4*cos(p4)*sin(s4)),0,0,-1,0,0,
   -1.1502419908960528*c1*(2*c1*c2 - 2*c3)*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.15471173963046492*c1*(1 + pow(c1,2) - pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.15471173963046492*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    1.1502419908960528*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.1502419908960528*c2*(2*c1*c2 - 2*c3)*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.15471173963046492*c2*(1 + pow(c1,2) - pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.1502419908960528*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.15471173963046492*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.1502419908960528*(2*c1*c2 - 2*c3)*c3*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.15471173963046492*c3*(1 + pow(c1,2) - pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.1502419908960528*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.15471173963046492*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0,0,0,0,
   l5*cos(p5)*cos(s5),0,0,0,0,0,-(l5*sin(p5)*sin(s5)),0,0,-1,0,
   -0.30942347926092983*c1*(c1*c2 + c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.1502419908960528*c1*(1 - pow(c1,2) + pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.1502419908960528*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.15471173963046492*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.30942347926092983*c2*(c1*c2 + c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.1502419908960528*c2*(1 - pow(c1,2) + pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.15471173963046492*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    1.1502419908960528*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.30942347926092983*c3*(c1*c2 + c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.1502419908960528*c3*(1 - pow(c1,2) + pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.15471173963046492*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    1.1502419908960528*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0,0,0,0,0,0,0,0,0,
   0,-(l5*cos(s5)),0,0,0,-1,-0.30942347926092983*c1*(-c2 + c1*c3)*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    2.3004839817921057*c1*(c1 + c2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.1502419908960528*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.15471173963046492*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.30942347926092983*c2*(-c2 + c1*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    2.3004839817921057*c2*(c1 + c2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.15471173963046492*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    1.1502419908960528*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.30942347926092983*c3*(-c2 + c1*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    2.3004839817921057*c3*(c1 + c2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.15471173963046492*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    1.1502419908960528*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0,0,0,0,
   -(l5*cos(s5)*sin(p5)),0,0,0,0,0,-(l5*cos(p5)*sin(s5)),0,-1,0,0,
   -0.44113669866436045*c1*(2*c1*c2 - 2*c3)*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.073494654430803*c1*(1 + pow(c1,2) - pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.073494654430803*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.44113669866436045*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.44113669866436045*c2*(2*c1*c2 - 2*c3)*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.073494654430803*c2*(1 + pow(c1,2) - pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.44113669866436045*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    1.073494654430803*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.44113669866436045*(2*c1*c2 - 2*c3)*c3*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.073494654430803*c3*(1 + pow(c1,2) - pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.44113669866436045*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    1.073494654430803*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0,0,0,0,0,
   l6*cos(p6)*cos(s6),0,0,0,0,0,-(l6*sin(p6)*sin(s6)),0,-1,0,
   2.146989308861606*c1*(c1*c2 + c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.44113669866436045*c1*(1 - pow(c1,2) + pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.44113669866436045*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    1.073494654430803*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   2.146989308861606*c2*(c1*c2 + c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.44113669866436045*c2*(1 - pow(c1,2) + pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.073494654430803*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.44113669866436045*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   2.146989308861606*c3*(c1*c2 + c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.44113669866436045*c3*(1 - pow(c1,2) + pow(c2,2) - pow(c3,2))*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.073494654430803*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.44113669866436045*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0,0,0,0,0,0,0,0,
   0,0,0,-(l6*cos(s6)),0,0,-1,2.146989308861606*c1*(-c2 + c1*c3)*
     pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.8822733973287209*c1*(c1 + c2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.44113669866436045*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    1.073494654430803*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   2.146989308861606*c2*(-c2 + c1*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.8822733973287209*c2*(c1 + c2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.073494654430803*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.44113669866436045*c3*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   2.146989308861606*c3*(-c2 + c1*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.8822733973287209*c3*(c1 + c2*c3)*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.073494654430803*c1*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.44113669866436045*c2*pow(1 + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0,0,0,0,0,
   -(l6*cos(s6)*sin(p6)),0,0,0,0,0,-(l6*cos(p6)*sin(s6));

    return Jetaphi;
}

//////////////////////////////////////Variables as Complex/////////////////////////////////

/*!
The funcion etaext is the set of constraint equations of an SRSPM formulated in
the extended configuration space of the manipulator. This is for the example
problem demosntrating the use of the root-trackers.
The formulation of the symbolic equations can be found in the accompanying
Mathematica notebook named, `SRSPM_root_tracking.nb`.

@param q Takes in the 24-dimensional extended-configuration-space variables as input
*/
VectorXcd etaextC(VectorXcd q){
  VectorXcd eta(18);
    std::complex<double> x = q(0), y = q(1), z = q(2), c1 = q(3), c2 = q(4), c3 = q(5),\
    p1 = q(6), p2 = q(7), p3 = q(8), p4 = q(9), p5=q(10), p6=q(11),\
    s1=q(12), s2=q(13), s3=q(14), s4=q(15), s5=q(16), s6=q(17),\
    l1=q(18), l2=q(19), l3=q(20), l4=q(21), l5=q(22), l6=q(23);

  eta << 0.7325760476016683 - 1.*x - 0.22056834933218009*(2.*c1*c2 - 2.*c3)*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.5367473272154016*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) + l1*cos(s1)*sin(p1),
   0.6806851948443702 - 1.*y - 1.0734946544308033*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.22056834933218009*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) - 1.*l1*sin(s1),
   -1.*z + l1*cos(p1)*cos(s1) - 1.0734946544308033*(-1.*c2 + c1*c3)*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.44113669866436017*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   0.22320264691435115 - 1.*x - 0.5751209954480265*(2.*c1*c2 - 2.*c3)*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.07735586981523229*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) + l2*cos(s2)*sin(p2),
   0.974772064849228 - 1.*y + 0.15471173963046458*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.5751209954480265*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) - 1.*l2*sin(s2),
   -1.*z + l2*cos(p2)*cos(s2) + 0.15471173963046458*(-1.*c2 + c1*c3)*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    1.150241990896053*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.955778694516019 - 1.*x - 0.35455264611584636*(2.*c1*c2 - 2.*c3)*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.4593914574001692*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) + l3*cos(s3)*sin(p3),
   0.2940868700048578 - 1.*y + 0.9187829148003384*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.35455264611584636*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) - 1.*l3*sin(s3),
   -1.*z + l3*cos(p3)*cos(s3) + 0.9187829148003384*(-1.*c2 + c1*c3)*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.7091052922316927*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.9557786945160192 - 1.*x + 0.35455264611584636*(2.*c1*c2 - 2.*c3)*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.4593914574001694*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) + l4*cos(s4)*sin(p4),
   -0.2940868700048576 - 1.*y + 0.9187829148003388*(c1*c2 + c3)*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.35455264611584636*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) - 1.*l4*sin(s4),
   -1.*z + l4*cos(p4)*cos(s4) + 0.9187829148003388*(-1.*c2 + c1*c3)*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.7091052922316927*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   0.22320264691435077 - 1.*x + 0.5751209954480264*(2.*c1*c2 - 2.*c3)*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.07735586981523246*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) + l5*cos(s5)*sin(p5),
   -0.974772064849228 - 1.*y + 0.15471173963046492*(c1*c2 + c3)*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.5751209954480264*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) - 1.*l5*sin(s5),
   -1.*z + l5*cos(p5)*cos(s5) + 0.15471173963046492*(-1.*c2 + c1*c3)*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    1.1502419908960528*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   0.7325760476016681 - 1.*x + 0.22056834933218022*(2.*c1*c2 - 2.*c3)*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.5367473272154015*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) + l6*cos(s6)*sin(p6),
   -0.6806851948443704 - 1.*y - 1.073494654430803*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.22056834933218022*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) - 1.*l6*sin(s6),
   -1.*z + l6*cos(p6)*cos(s6) - 1.073494654430803*(-1.*c2 + c1*c3)*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.44113669866436045*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1);
  return eta;
}


/*!
The funcion Jetaexttheta is the Jacobian matrix of an SRSPM formulated in
the extended configuration space of the manipulator. This is for the example
problem demonstrating the use of the root-trackers.
The formulation of the symbolic equations can be found in the accompanying
Mathematica notebook named, `SRSPM_root_tracking.nb`.

@param q Takes in the 24-dimensional extended-configuration-space variables as input
*/
MatrixXcd JetaextthetaC(VectorXcd q){
  MatrixXcd Jetatheta(18, 6);
  std::complex<double> x = q(0), y = q(1), z = q(2), c1 = q(3), c2 = q(4), c3 = q(5),\
    p1 = q(6), p2 = q(7), p3 = q(8), p4 = q(9), p5=q(10), p6=q(11),\
    s1=q(12), s2=q(13), s3=q(14), s4=q(15), s5=q(16), s6=q(17),\
    l1=q(18), l2=q(19), l3=q(20), l4=q(21), l5=q(22), l6=q(23);

    Jetatheta << cos(s1)*sin(p1),0.0,0.0,0.0,0.0,0.0,-sin(s1),0.0,0.0,0.0,0.0,0.0,cos(p1)*cos(s1),0.0,0.0,0.0,0.0,0.0,0.0,
     cos(s2)*sin(p2),0.0,0.0,0.0,0.0,0.0,-sin(s2),0.0,0.0,0.0,0.0,0.0,cos(p2)*cos(s2),0.0,0.0,0.0,0.0,0.0,0.0,
     cos(s3)*sin(p3),0.0,0.0,0.0,0.0,0.0,-sin(s3),0.0,0.0,0.0,0.0,0.0,cos(p3)*cos(s3),0.0,0.0,0.0,0.0,0.0,0.0,
     cos(s4)*sin(p4),0.0,0.0,0.0,0.0,0.0,-sin(s4),0.0,0.0,0.0,0.0,0.0,cos(p4)*cos(s4),0.0,0.0,0.0,0.0,0.0,0.0,
     cos(s5)*sin(p5),0.0,0.0,0.0,0.0,0.0,-sin(s5),0.0,0.0,0.0,0.0,0.0,cos(p5)*cos(s5),0.0,0.0,0.0,0.0,0.0,0.0,
     cos(s6)*sin(p6),0.0,0.0,0.0,0.0,0.0,-sin(s6),0.0,0.0,0.0,0.0,0.0,cos(p6)*cos(s6);

   return Jetatheta;
 }


/*!
The funcion Jetaextphi is the constraint Jacobian matrix of an SRSPM formulated in
the extended configuration space of the manipulator. This is for the example
problem demonstrating the use of the root-trackers.
The formulation of the symbolic equations can be found in the accompanying
Mathematica notebook named, `SRSPM_root_tracking.nb`.

@param q Takes in the 24-dimensional extended-configuration-space variables as input
*/
MatrixXcd JetaextphiC(VectorXcd q){
  MatrixXcd Jetaphi(18, 18);
    std::complex<double> x = q(0), y = q(1), z = q(2), c1 = q(3), c2 = q(4), c3 = q(5),\
    p1 = q(6), p2 = q(7), p3 = q(8), p4 = q(9), p5=q(10), p6=q(11),\
    s1=q(12), s2=q(13), s3=q(14), s4=q(15), s5=q(16), s6=q(17),\
    l1=q(18), l2=q(19), l3=q(20), l4=q(21), l5=q(22), l6=q(23);

    Jetaphi << -1.,0.,0.,0.44113669866436017*c1*(2.*c1*c2 - 2.*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.0734946544308033*c1*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.0734946544308033*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.44113669866436017*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   0.44113669866436017*c2*(2.*c1*c2 - 2.*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.0734946544308033*c2*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.44113669866436017*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    1.0734946544308033*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   0.44113669866436017*(2.*c1*c2 - 2.*c3)*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.0734946544308033*c3*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.44113669866436017*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    1.0734946544308033*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),l1*cos(p1)*cos(s1),0.,0.,0.,0.,0.,
   -1.*l1*sin(p1)*sin(s1),0.,0.,0.,0.,0.,0.,-1.,0.,
   2.1469893088616065*c1*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.44113669866436017*c1*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.44113669866436017*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    1.0734946544308033*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   2.1469893088616065*c2*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.44113669866436017*c2*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.0734946544308033*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.44113669866436017*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   2.1469893088616065*c3*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.44113669866436017*c3*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.0734946544308033*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.44113669866436017*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0.,0.,0.,0.,0.,0.,-1.*l1*cos(s1),0.,
   0.,0.,0.,0.,0.,0.,-1.,2.1469893088616065*c1*(-1.*c2 + c1*c3)*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.8822733973287203*c1*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.44113669866436017*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    1.0734946544308033*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   2.1469893088616065*c2*(-1.*c2 + c1*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.8822733973287203*c2*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.0734946544308033*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.44113669866436017*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   2.1469893088616065*c3*(-1.*c2 + c1*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.8822733973287203*c3*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.0734946544308033*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.44113669866436017*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),-1.*l1*cos(s1)*sin(p1),0.,0.,0.,0.,0.,
   -1.*l1*cos(p1)*sin(s1),0.,0.,0.,0.,0.,-1.,0.,0.,
   1.150241990896053*c1*(2.*c1*c2 - 2.*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.15471173963046458*c1*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.15471173963046458*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    1.150241990896053*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   1.150241990896053*c2*(2.*c1*c2 - 2.*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.15471173963046458*c2*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.150241990896053*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.15471173963046458*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   1.150241990896053*(2.*c1*c2 - 2.*c3)*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.15471173963046458*c3*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.150241990896053*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.15471173963046458*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0.,l2*cos(p2)*cos(s2),0.,0.,0.,0.,0.,
   -1.*l2*sin(p2)*sin(s2),0.,0.,0.,0.,0.,-1.,0.,
   -0.30942347926092917*c1*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.150241990896053*c1*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.150241990896053*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.15471173963046458*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.30942347926092917*c2*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.150241990896053*c2*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.15471173963046458*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    1.150241990896053*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.30942347926092917*c3*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.150241990896053*c3*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.15471173963046458*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    1.150241990896053*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0.,0.,0.,0.,0.,0.,0.,-1.*l2*cos(s2),0.,
   0.,0.,0.,0.,0.,-1.,-0.30942347926092917*c1*(-1.*c2 + c1*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    2.300483981792106*c1*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.150241990896053*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.15471173963046458*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.30942347926092917*c2*(-1.*c2 + c1*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    2.300483981792106*c2*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.15471173963046458*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    1.150241990896053*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.30942347926092917*c3*(-1.*c2 + c1*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    2.300483981792106*c3*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.15471173963046458*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    1.150241990896053*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0.,-1.*l2*cos(s2)*sin(p2),0.,0.,0.,0.,
   0.,-1.*l2*cos(p2)*sin(s2),0.,0.,0.,0.,-1.,0.,0.,
   0.7091052922316927*c1*(2.*c1*c2 - 2.*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.9187829148003384*c1*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.9187829148003384*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.7091052922316927*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   0.7091052922316927*c2*(2.*c1*c2 - 2.*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.9187829148003384*c2*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.7091052922316927*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.9187829148003384*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   0.7091052922316927*(2.*c1*c2 - 2.*c3)*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.9187829148003384*c3*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.7091052922316927*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.9187829148003384*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0.,0.,l3*cos(p3)*cos(s3),0.,0.,0.,0.,
   0.,-1.*l3*sin(p3)*sin(s3),0.,0.,0.,0.,-1.,0.,
   -1.8375658296006767*c1*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.7091052922316927*c1*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.7091052922316927*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.9187829148003384*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.8375658296006767*c2*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.7091052922316927*c2*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.9187829148003384*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.7091052922316927*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.8375658296006767*c3*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.7091052922316927*c3*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.9187829148003384*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.7091052922316927*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0.,0.,0.,0.,0.,0.,0.,0.,-1.*l3*cos(s3),
   0.,0.,0.,0.,0.,-1.,-1.8375658296006767*c1*(-1.*c2 + c1*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.4182105844633854*c1*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.7091052922316927*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.9187829148003384*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.8375658296006767*c2*(-1.*c2 + c1*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.4182105844633854*c2*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.9187829148003384*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.7091052922316927*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.8375658296006767*c3*(-1.*c2 + c1*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.4182105844633854*c3*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.9187829148003384*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.7091052922316927*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0.,0.,-1.*l3*cos(s3)*sin(p3),0.,0.,0.,
   0.,0.,-1.*l3*cos(p3)*sin(s3),0.,0.,0.,-1.,0.,0.,
   -0.7091052922316927*c1*(2.*c1*c2 - 2.*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.9187829148003388*c1*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.9187829148003388*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.7091052922316927*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.7091052922316927*c2*(2.*c1*c2 - 2.*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.9187829148003388*c2*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.7091052922316927*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.9187829148003388*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.7091052922316927*(2.*c1*c2 - 2.*c3)*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.9187829148003388*c3*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.7091052922316927*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.9187829148003388*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0.,0.,0.,l4*cos(p4)*cos(s4),0.,0.,0.,
   0.,0.,-1.*l4*sin(p4)*sin(s4),0.,0.,0.,-1.,0.,
   -1.8375658296006776*c1*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.7091052922316927*c1*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.7091052922316927*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.9187829148003388*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.8375658296006776*c2*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.7091052922316927*c2*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.9187829148003388*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.7091052922316927*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.8375658296006776*c3*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.7091052922316927*c3*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.9187829148003388*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.7091052922316927*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0.,0.,0.,0.,0.,0.,0.,0.,0.,
   -1.*l4*cos(s4),0.,0.,0.,0.,-1.,-1.8375658296006776*c1*(-1.*c2 + c1*c3)*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.4182105844633854*c1*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.7091052922316927*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.9187829148003388*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.8375658296006776*c2*(-1.*c2 + c1*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.4182105844633854*c2*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.9187829148003388*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.7091052922316927*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.8375658296006776*c3*(-1.*c2 + c1*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.4182105844633854*c3*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.9187829148003388*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.7091052922316927*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0.,0.,0.,-1.*l4*cos(s4)*sin(p4),0.,0.,
   0.,0.,0.,-1.*l4*cos(p4)*sin(s4),0.,0.,-1.,0.,0.,
   -1.1502419908960528*c1*(2.*c1*c2 - 2.*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.15471173963046492*c1*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.15471173963046492*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    1.1502419908960528*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.1502419908960528*c2*(2.*c1*c2 - 2.*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.15471173963046492*c2*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.1502419908960528*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.15471173963046492*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -1.1502419908960528*(2.*c1*c2 - 2.*c3)*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.15471173963046492*c3*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.1502419908960528*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.15471173963046492*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0.,0.,0.,0.,l5*cos(p5)*cos(s5),0.,0.,
   0.,0.,0.,-1.*l5*sin(p5)*sin(s5),0.,0.,-1.,0.,
   -0.30942347926092983*c1*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.1502419908960528*c1*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.1502419908960528*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.15471173963046492*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.30942347926092983*c2*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.1502419908960528*c2*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.15471173963046492*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    1.1502419908960528*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.30942347926092983*c3*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.1502419908960528*c3*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.15471173963046492*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    1.1502419908960528*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
   -1.*l5*cos(s5),0.,0.,0.,-1.,-0.30942347926092983*c1*(-1.*c2 + c1*c3)*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    2.3004839817921057*c1*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.1502419908960528*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.15471173963046492*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.30942347926092983*c2*(-1.*c2 + c1*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    2.3004839817921057*c2*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.15471173963046492*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    1.1502419908960528*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.30942347926092983*c3*(-1.*c2 + c1*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    2.3004839817921057*c3*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.15471173963046492*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    1.1502419908960528*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0.,0.,0.,0.,-1.*l5*cos(s5)*sin(p5),0.,
   0.,0.,0.,0.,-1.*l5*cos(p5)*sin(s5),0.,-1.,0.,0.,
   -0.44113669866436045*c1*(2.*c1*c2 - 2.*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.073494654430803*c1*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.073494654430803*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.44113669866436045*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.44113669866436045*c2*(2.*c1*c2 - 2.*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.073494654430803*c2*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.44113669866436045*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    1.073494654430803*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   -0.44113669866436045*(2.*c1*c2 - 2.*c3)*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.073494654430803*c3*(1. + pow(c1,2) - 1.*pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.44113669866436045*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    1.073494654430803*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0.,0.,0.,0.,0.,l6*cos(p6)*cos(s6),0.,0.,
   0.,0.,0.,-1.*l6*sin(p6)*sin(s6),0.,-1.,0.,
   2.146989308861606*c1*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.44113669866436045*c1*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.44113669866436045*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    1.073494654430803*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   2.146989308861606*c2*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.44113669866436045*c2*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.073494654430803*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.44113669866436045*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   2.146989308861606*c3*(c1*c2 + c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.44113669866436045*c3*(1. - 1.*pow(c1,2) + pow(c2,2) - 1.*pow(c3,2))*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.073494654430803*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    0.44113669866436045*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
   -1.*l6*cos(s6),0.,0.,-1.,2.146989308861606*c1*(-1.*c2 + c1*c3)*
     pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.8822733973287209*c1*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    0.44113669866436045*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) -
    1.073494654430803*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   2.146989308861606*c2*(-1.*c2 + c1*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.8822733973287209*c2*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) +
    1.073494654430803*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.44113669866436045*c3*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
   2.146989308861606*c3*(-1.*c2 + c1*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    0.8822733973287209*c3*(c1 + c2*c3)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2) -
    1.073494654430803*c1*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1) +
    0.44113669866436045*c2*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),0.,0.,0.,0.,0.,-1.*l6*cos(s6)*sin(p6),
   0.,0.,0.,0.,0.,-1.*l6*cos(p6)*sin(s6);

    return Jetaphi;
}

#endif

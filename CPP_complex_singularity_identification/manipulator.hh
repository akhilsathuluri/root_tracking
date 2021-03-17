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

// IK handles only reals for now
VectorXd computeLvals(VectorXd qx){
  double x = qx(0), y = qx(1), z = qx(2), c1 = qx(3), c2 = qx(4), c3 = qx(5);
  VectorXd lvals(6);
  lvals << pow((0.2500563992522112 - 0.8150914777120721*c3 - 0.3916574407725333*x - 0.8822733973287203*c3*x - 0.9202336910243802*y + 2.1469893088616065*c3*y +
       c1*c2*(-2.1077561943504124 + 0.8822733973287203*x + 2.1469893088616065*y) + c2*(-2.1469893088616065 + 0.8822733973287203*c3)*z +
       c1*(0.8822733973287203 + 2.1469893088616065*c3)*z + 2.4234397807477883*pow(c3,2) - 2.5386467496341396*x*pow(c3,2) - 1.8025070883531005*y*pow(c3,2) + pow(x,2) +
       pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) + pow(c3,2)*pow(z,2) +
       pow(c1,2)*(0.8506068386189155 - 0.3916574407725333*x - 1.8025070883531005*y + pow(x,2) + pow(y,2) + pow(z,2)) +
       pow(c2,2)*(1.8228893413810843 - 2.5386467496341396*x - 0.9202336910243802*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
    0.5),pow((0.25005639925221135 + 0.8150914777120721*c3 - 0.6011170334591669*x - 2.3004839817921057*c3*x +
       c1*c2*(-0.21185675012805605 + 2.3004839817921057*x - 0.30942347926092917*y) - 0.7993021388024031*y - 0.30942347926092917*c3*y +
       c1*(2.3004839817921057 - 0.30942347926092917*c3)*z + c2*(0.30942347926092917 + 2.3004839817921057*c3)*z + 2.423439780747789*pow(c3,2) -
       0.2916935541982377*x*pow(c3,2) - 3.0997861205945085*y*pow(c3,2) + pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) +
       pow(c3,2)*pow(z,2) + pow(c1,2)*(2.4925039203362758 - 0.6011170334591669*x - 3.0997861205945085*y + pow(x,2) + pow(y,2) + pow(z,2)) +
       pow(c2,2)*(0.1809922596637239 - 0.2916935541982377*x - 0.7993021388024031*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
    0.5),pow((0.2500563992522111 - 0.8150914777120724*c3 + 0.9927744742316997*x - 1.4182105844633854*c3*x +
       c1*c2*(1.895899444222358 + 1.4182105844633854*x - 1.8375658296006767*y) + 0.12093155222197716*y - 1.8375658296006767*c3*y +
       c1*(1.4182105844633854 - 1.8375658296006767*c3)*z + c2*(1.8375658296006767 + 1.4182105844633854*c3)*z + 2.4234397807477888*pow(c3,2) +
       2.8303403038323767*x*pow(c3,2) - 1.2972790322414083*y*pow(c3,2) + pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) +
       pow(c3,2)*pow(z,2) + pow(c1,2)*(0.6671335110448082 + 0.9927744742316997*x - 1.2972790322414083*y + pow(x,2) + pow(y,2) + pow(z,2)) +
       pow(c2,2)*(2.0063626689551914 + 2.8303403038323767*x + 0.12093155222197716*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
    0.5),pow((0.25005639925221135 + 0.815091477712072*c3 + 0.9927744742317001*x + 1.418210584463385*c3*x +
       c1*c2*(-1.8958994442223576 - 1.4182105844633848*x - 1.8375658296006767*y) - 0.12093155222197716*y - 1.8375658296006767*c3*y +
       c1*(-1.4182105844633848 - 1.8375658296006767*c3)*z + c2*(1.837565829600677 - 1.4182105844633848*c3)*z + 2.4234397807477888*pow(c3,2) +
       2.8303403038323767*x*pow(c3,2) + 1.2972790322414078*y*pow(c3,2) + pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) +
       pow(c3,2)*pow(z,2) + pow(c2,2)*(2.0063626689551923 + 2.8303403038323767*x - 0.12093155222197716*y + pow(x,2) + pow(y,2) + pow(z,2)) +
       pow(c1,2)*(0.6671335110448082 + 0.9927744742317001*x + 1.2972790322414078*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
    0.5),pow((0.2500563992522112 - 0.815091477712072*c3 - 0.6011170334591664*x + 2.3004839817921057*c3*x +
       c1*c2*(0.21185675012805427 - 2.3004839817921057*x - 0.30942347926092983*y) + 0.7993021388024031*y - 0.30942347926092983*c3*y +
       c2*(0.30942347926092983 - 2.3004839817921057*c3)*z + c1*(-2.3004839817921057 - 0.30942347926092983*c3)*z + 2.4234397807477888*pow(c3,2) -
       0.2916935541982366*x*pow(c3,2) + 3.0997861205945085*y*pow(c3,2) + pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) +
       pow(c3,2)*pow(z,2) + pow(c2,2)*(0.180992259663724 - 0.2916935541982366*x + 0.7993021388024031*y + pow(x,2) + pow(y,2) + pow(z,2)) +
       pow(c1,2)*(2.4925039203362758 - 0.6011170334591664*x + 3.0997861205945085*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
    0.5),pow((0.25005639925221124 + 0.8150914777120724*c3 - 0.39165744077253295*x + 0.8822733973287207*c3*x + 0.9202336910243805*y + 2.146989308861606*c3*y +
       c1*c2*(2.1077561943504133 - 0.8822733973287208*x + 2.146989308861606*y) + c2*(-2.146989308861606 - 0.8822733973287208*c3)*z +
       c1*(-0.8822733973287208 + 2.146989308861606*c3)*z + 2.4234397807477888*pow(c3,2) - 2.5386467496341396*x*pow(c3,2) + 1.8025070883531011*y*pow(c3,2) + pow(x,2) +
       pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) + pow(c3,2)*pow(z,2) +
       pow(c2,2)*(1.8228893413810836 - 2.5386467496341396*x + 0.9202336910243805*y + pow(x,2) + pow(y,2) + pow(z,2)) +
       pow(c1,2)*(0.8506068386189165 - 0.39165744077253295*x + 1.8025070883531011*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
    0.5);

  return lvals;
}


/*!
The computeXfromParam returns input values x from the path parametrisation.
Edit this function to suite the required path paramterisation and the mapping
between the parametrisation and the actual input values. In the case of an SRSPM,
we use the IK mapping to find the leg lengths from the task space parametrisation
*/

VectorXd computeXfromParam(double alpha){
  // The first example is for the SRSPM so using the parametrisation of the SRSPM path
  VectorXd qx(6), l(6);
  double scale = 0.2;
  qx << scale-alpha, alpha, 1.28-alpha, scale, 0, 0;
  // std::cout << qx << '\n';
  // Now since the NRTracker needs the legvals for the simulation we use the
  // IKsolver to compute legvals from qx
  l = computeLvals(qx);
  return l;
}


VectorXd computeIK(VectorXd qx){
  VectorXd iksol(18), lvals(6), phivals(6), psivals(6);
  double x = qx(0), y = qx(1), z = qx(2), c1 = qx(3), c2 = qx(4), c3 = qx(5);
  lvals << pow((0.2500563992522112 - 0.8150914777120721*c3 - 0.3916574407725333*x - 0.8822733973287203*c3*x - 0.9202336910243802*y + 2.1469893088616065*c3*y +
       c1*c2*(-2.1077561943504124 + 0.8822733973287203*x + 2.1469893088616065*y) + c2*(-2.1469893088616065 + 0.8822733973287203*c3)*z +
       c1*(0.8822733973287203 + 2.1469893088616065*c3)*z + 2.4234397807477883*pow(c3,2) - 2.5386467496341396*x*pow(c3,2) - 1.8025070883531005*y*pow(c3,2) + pow(x,2) +
       pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) + pow(c3,2)*pow(z,2) +
       pow(c1,2)*(0.8506068386189155 - 0.3916574407725333*x - 1.8025070883531005*y + pow(x,2) + pow(y,2) + pow(z,2)) +
       pow(c2,2)*(1.8228893413810843 - 2.5386467496341396*x - 0.9202336910243802*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
    0.5),pow((0.25005639925221135 + 0.8150914777120721*c3 - 0.6011170334591669*x - 2.3004839817921057*c3*x +
       c1*c2*(-0.21185675012805605 + 2.3004839817921057*x - 0.30942347926092917*y) - 0.7993021388024031*y - 0.30942347926092917*c3*y +
       c1*(2.3004839817921057 - 0.30942347926092917*c3)*z + c2*(0.30942347926092917 + 2.3004839817921057*c3)*z + 2.423439780747789*pow(c3,2) -
       0.2916935541982377*x*pow(c3,2) - 3.0997861205945085*y*pow(c3,2) + pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) +
       pow(c3,2)*pow(z,2) + pow(c1,2)*(2.4925039203362758 - 0.6011170334591669*x - 3.0997861205945085*y + pow(x,2) + pow(y,2) + pow(z,2)) +
       pow(c2,2)*(0.1809922596637239 - 0.2916935541982377*x - 0.7993021388024031*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
    0.5),pow((0.2500563992522111 - 0.8150914777120724*c3 + 0.9927744742316997*x - 1.4182105844633854*c3*x +
       c1*c2*(1.895899444222358 + 1.4182105844633854*x - 1.8375658296006767*y) + 0.12093155222197716*y - 1.8375658296006767*c3*y +
       c1*(1.4182105844633854 - 1.8375658296006767*c3)*z + c2*(1.8375658296006767 + 1.4182105844633854*c3)*z + 2.4234397807477888*pow(c3,2) +
       2.8303403038323767*x*pow(c3,2) - 1.2972790322414083*y*pow(c3,2) + pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) +
       pow(c3,2)*pow(z,2) + pow(c1,2)*(0.6671335110448082 + 0.9927744742316997*x - 1.2972790322414083*y + pow(x,2) + pow(y,2) + pow(z,2)) +
       pow(c2,2)*(2.0063626689551914 + 2.8303403038323767*x + 0.12093155222197716*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
    0.5),pow((0.25005639925221135 + 0.815091477712072*c3 + 0.9927744742317001*x + 1.418210584463385*c3*x +
       c1*c2*(-1.8958994442223576 - 1.4182105844633848*x - 1.8375658296006767*y) - 0.12093155222197716*y - 1.8375658296006767*c3*y +
       c1*(-1.4182105844633848 - 1.8375658296006767*c3)*z + c2*(1.837565829600677 - 1.4182105844633848*c3)*z + 2.4234397807477888*pow(c3,2) +
       2.8303403038323767*x*pow(c3,2) + 1.2972790322414078*y*pow(c3,2) + pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) +
       pow(c3,2)*pow(z,2) + pow(c2,2)*(2.0063626689551923 + 2.8303403038323767*x - 0.12093155222197716*y + pow(x,2) + pow(y,2) + pow(z,2)) +
       pow(c1,2)*(0.6671335110448082 + 0.9927744742317001*x + 1.2972790322414078*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
    0.5),pow((0.2500563992522112 - 0.815091477712072*c3 - 0.6011170334591664*x + 2.3004839817921057*c3*x +
       c1*c2*(0.21185675012805427 - 2.3004839817921057*x - 0.30942347926092983*y) + 0.7993021388024031*y - 0.30942347926092983*c3*y +
       c2*(0.30942347926092983 - 2.3004839817921057*c3)*z + c1*(-2.3004839817921057 - 0.30942347926092983*c3)*z + 2.4234397807477888*pow(c3,2) -
       0.2916935541982366*x*pow(c3,2) + 3.0997861205945085*y*pow(c3,2) + pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) +
       pow(c3,2)*pow(z,2) + pow(c2,2)*(0.180992259663724 - 0.2916935541982366*x + 0.7993021388024031*y + pow(x,2) + pow(y,2) + pow(z,2)) +
       pow(c1,2)*(2.4925039203362758 - 0.6011170334591664*x + 3.0997861205945085*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
    0.5),pow((0.25005639925221124 + 0.8150914777120724*c3 - 0.39165744077253295*x + 0.8822733973287207*c3*x + 0.9202336910243805*y + 2.146989308861606*c3*y +
       c1*c2*(2.1077561943504133 - 0.8822733973287208*x + 2.146989308861606*y) + c2*(-2.146989308861606 - 0.8822733973287208*c3)*z +
       c1*(-0.8822733973287208 + 2.146989308861606*c3)*z + 2.4234397807477888*pow(c3,2) - 2.5386467496341396*x*pow(c3,2) + 1.8025070883531011*y*pow(c3,2) + pow(x,2) +
       pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) + pow(c3,2)*pow(z,2) +
       pow(c2,2)*(1.8228893413810836 - 2.5386467496341396*x + 0.9202336910243805*y + pow(x,2) + pow(y,2) + pow(z,2)) +
       pow(c1,2)*(0.8506068386189165 - 0.39165744077253295*x + 1.8025070883531011*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
    0.5);

    double l1 = lvals(0), l2 = lvals(1), l3 = lvals(2), l4 = lvals(3), l5 = lvals(4), l6 = lvals(5);

    psivals << atan2((0.4601168455121902 - 1.0734946544308033*c1*c2 - 1.0734946544308033*c3 - 1.*y + (0.9012535441765503 - 1.*y)*pow(c1,2) +
       (0.4601168455121902 - 1.*y)*pow(c2,2) + 0.9012535441765503*pow(c3,2) - 1.*y*pow(c3,2))*pow(l1,-1)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
    pow(pow(l1,-2)*(0.03834888772812256 + 0.17277447042972754*c3 - 0.3916574407725333*x - 0.8822733973287203*c3*x +
        (c2*(-0.17277447042972754 + 0.8822733973287203*x) + (0.8822733973287203 + 2.1469893088616065*c3)*z)*pow(c1,3) +
        (-2.1469893088616065 + 0.8822733973287203*c3)*z*pow(c2,3) + 0.691741531402099*pow(c3,2) - 2.930304190406673*x*pow(c3,2) +
        0.8822733973287202*c2*(-2.433473926973312 + 1.*c3)*z*(1. + 1.*pow(c3,2)) +
        c1*((0.8822733973287203 + 2.1469893088616065*c3)*z*pow(c2,2) + (-1.1198902462086129 + 0.8822733973287203*x)*pow(c2,3) +
           2.1469893088616065*(0.4109351610122951 + 1.*c3)*z*(1. + 1.*pow(c3,2)) +
           c2*(-1.1198902462086129 - 2.3047815461830194*c3 + 0.8822733973287203*x + (-0.17277447042972754 + 0.8822733973287203*x)*pow(c3,2))) +
        1.1198902462086129*pow(c3,3) - 0.8822733973287203*x*pow(c3,3) + 1.6111818298569958*pow(c3,4) - 2.5386467496341396*x*pow(c3,4) + pow(x,2) +
        2.*pow(c3,2)*pow(x,2) + pow(c3,4)*pow(x,2) + pow(z,2) + 2.*pow(c3,2)*pow(z,2) + pow(c3,4)*pow(z,2) +
        pow(c2,4)*(1.6111818298569958 - 2.5386467496341396*x + pow(x,2) + pow(z,2)) + pow(c1,4)*(0.03834888772812256 - 0.3916574407725333*x + pow(x,2) + pow(z,2)) +
        pow(c2,2)*(1.6495307175851184 + c3*(0.17277447042972754 - 0.8822733973287203*x) - 2.930304190406673*x + 2.*pow(x,2) + 2.*pow(z,2) +
           pow(c3,2)*(3.416965246622482 - 5.077293499268279*x + 2.*pow(x,2) + 2.*pow(z,2))) +
        pow(c1,2)*(0.2712993623647356 + c3*(1.1198902462086129 - 0.8822733973287203*x) - 0.7833148815450666*x + c2*(-2.1469893088616065 + 0.8822733973287203*c3)*z +
           2.*pow(x,2) + 2.*pow(z,2) + pow(c2,2)*(0.691741531402099 - 2.930304190406673*x + 2.*pow(x,2) + 2.*pow(z,2)) +
           pow(c3,2)*(1.6495307175851184 - 2.930304190406673*x + 2.*pow(x,2) + 2.*pow(z,2))))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2),0.5)),
   atan2((0.39965106940120154 + 0.15471173963046464*c1*c2 + 0.15471173963046464*c3 - 1.*y + (1.5498930602972543 - 1.*y)*pow(c1,2) +
       (0.39965106940120154 - 1.*y)*pow(c2,2) + 1.5498930602972543*pow(c3,2) - 1.*y*pow(c3,2))*pow(l2,-1)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
    pow(pow(l2,-2)*(0.0903354219786873 + 0.6914300533276011*c3 - 0.6011170334591668*x - 2.3004839817921057*c3*x +
        (c2*(-0.6914300533276011 + 2.3004839817921057*x) + (2.3004839817921057 - 0.30942347926092917*c3)*z)*pow(c1,3) +
        (0.30942347926092917 + 2.3004839817921057*c3)*z*pow(c2,3) + 1.4107276196099183*pow(c3,2) - 0.8928105876574046*x*pow(c3,2) +
        2.3004839817921052*c2*(0.13450364432439318 + 1.*c3)*z*(1. + 1.*pow(c3,2)) +
        c1*((2.3004839817921057 - 0.30942347926092917*c3)*z*pow(c2,2) + (-0.3355181745125263 + 2.3004839817921057*x)*pow(c2,3) -
           0.3094234792609292*(-7.434742790970185 + 1.*c3)*z*(1. + 1.*pow(c3,2)) +
           c2*(-0.3355181745125263 - 0.04787144475896937*c3 + 2.3004839817921057*x + (-0.6914300533276011 + 2.3004839817921057*x)*pow(c3,2))) +
        0.3355181745125263*pow(c3,3) - 2.3004839817921057*x*pow(c3,3) + 0.02127128239019993*pow(c3,4) - 0.2916935541982376*x*pow(c3,4) + pow(x,2) +
        2.*pow(c3,2)*pow(x,2) + pow(c3,4)*pow(x,2) + pow(z,2) + 2.*pow(c3,2)*pow(z,2) + pow(c3,4)*pow(z,2) +
        pow(c1,4)*(0.09033542197868727 - 0.6011170334591668*x + pow(x,2) + pow(z,2)) + pow(c2,4)*(0.02127128239019993 - 0.2916935541982376*x + pow(x,2) + pow(z,2)) +
        pow(c1,2)*(1.5037274815778905 + c3*(0.3355181745125263 - 2.3004839817921057*x) - 1.2022340669183336*x + c2*(0.30942347926092917 + 2.3004839817921057*c3)*z +
           2.*pow(x,2) + 2.*pow(z,2) + pow(c3,2)*(0.11160670436888731 - 0.8928105876574046*x + 2.*pow(x,2) + 2.*pow(z,2)) +
           pow(c2,2)*(1.4107276196099183 - 0.8928105876574046*x + 2.*pow(x,2) + 2.*pow(z,2))) +
        pow(c2,2)*(0.11160670436888731 + c3*(0.6914300533276011 - 2.3004839817921057*x) - 0.8928105876574046*x + 2.*pow(x,2) + 2.*pow(z,2) +
           pow(c3,2)*(1.3655992024009158 - 0.5833871083964752*x + 2.*pow(x,2) + 2.*pow(z,2))))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2),0.5)),
   atan2((-0.06046577611098858 + 0.9187829148003385*c1*c2 + 0.9187829148003385*c3 - 1.*y + (0.648639516120704 - 1.*y)*pow(c1,2) +
       (-0.06046577611098858 - 1.*y)*pow(c2,2) + 0.648639516120704*pow(c3,2) - 1.*y*pow(c3,2))*pow(l3,-1)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
    1.*pow(pow(l3,-2)*(0.24640028917150705 - 0.7039816336702347*c3 + 0.9927744742316997*x - 1.4182105844633854*c3*x +
         (c2*(0.7039816336702347 + 1.4182105844633854*x) + (1.4182105844633854 - 1.8375658296006767*c3)*z)*pow(c1,3) +
         (1.8375658296006767 + 1.4182105844633854*c3)*z*pow(c2,3) + 1.907775118987983*pow(c3,2) + 3.823114778064076*x*pow(c3,2) +
         1.418210584463385*c2*(1.2956932134983075 + 1.*c3)*z*(1. + 1.*pow(c3,2)) +
         c1*((1.4182105844633854 - 1.8375658296006767*c3)*z*pow(c2,2) + (2.0070092882641957 + 1.4182105844633854*x)*pow(c2,3) -
            1.8375658296006767*(-0.7717876342811502 + 1.*c3)*z*(1. + 1.*pow(c3,2)) +
            c2*(2.0070092882641957 - 1.6883240890580118*c3 + 1.4182105844633854*x + (0.7039816336702347 + 1.4182105844633854*x)*pow(c3,2))) -
         2.0070092882641952*pow(c3,3) - 1.4182105844633854*x*pow(c3,3) + 2.0027065588744875*pow(c3,4) + 2.8303403038323767*x*pow(c3,4) + pow(x,2) +
         2.*pow(c3,2)*pow(x,2) + pow(c3,4)*pow(x,2) + pow(z,2) + 2.*pow(c3,2)*pow(z,2) + pow(c3,4)*pow(z,2) +
         pow(c1,4)*(0.24640028917150703 + 0.9927744742316997*x + pow(x,2) + pow(z,2)) + pow(c2,4)*(2.0027065588744875 + 2.8303403038323767*x + pow(x,2) + pow(z,2)) +
         pow(c1,2)*(0.9956308938140086 + c3*(-2.0070092882641952 - 1.4182105844633854*x) + 1.9855489484633995*x + c2*(1.8375658296006767 + 1.4182105844633854*c3)*z +
            2.*pow(x,2) + 2.*pow(z,2) + pow(c2,2)*(1.907775118987983 + 3.823114778064076*x + 2.*pow(x,2) + 2.*pow(z,2)) +
            pow(c3,2)*(2.2491068480459946 + 3.823114778064076*x + 2.*pow(x,2) + 2.*pow(z,2))) +
         pow(c2,2)*(2.2491068480459946 + c3*(-0.7039816336702347 - 1.4182105844633854*x) + 3.823114778064076*x + 2.*pow(x,2) + 2.*pow(z,2) +
            pow(c3,2)*(4.5082434332199695 + 5.660680607664753*x + 2.*pow(x,2) + 2.*pow(z,2))))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2),0.5)),
   atan2((0.06046577611098858 + 0.9187829148003385*c1*c2 + 0.9187829148003385*c3 - 1.*y + (-0.648639516120704 - 1.*y)*pow(c1,2) +
       (0.06046577611098858 - 1.*y)*pow(c2,2) - 0.648639516120704*pow(c3,2) - 1.*y*pow(c3,2))*pow(l4,-1)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
    1.*pow(pow(l4,-2)*(0.24640028917150714 + 0.7039816336702344*c3 + 0.9927744742316998*x + 1.418210584463385*c3*x +
         (c2*(-0.7039816336702341 - 1.4182105844633848*x) + (-1.4182105844633848 - 1.8375658296006767*c3)*z)*pow(c1,3) +
         (1.837565829600677 - 1.4182105844633848*c3)*z*pow(c2,3) + 1.907775118987983*pow(c3,2) + 3.8231147780640766*x*pow(c3,2) -
         1.418210584463385*c2*(-1.2956932134983083 + 1.*c3)*z*(1. + 1.*pow(c3,2)) +
         c1*((-1.4182105844633848 - 1.8375658296006767*c3)*z*pow(c2,2) + (-2.007009288264195 - 1.4182105844633848*x)*pow(c2,3) -
            1.8375658296006767*(0.7717876342811498 + 1.*c3)*z*(1. + 1.*pow(c3,2)) +
            c2*(-2.007009288264195 - 1.6883240890580127*c3 - 1.4182105844633848*x + (-0.7039816336702341 - 1.4182105844633848*x)*pow(c3,2))) +
         2.007009288264195*pow(c3,3) + 1.418210584463385*x*pow(c3,3) + 2.0027065588744883*pow(c3,4) + 2.8303403038323767*x*pow(c3,4) + pow(x,2) +
         2.*pow(c3,2)*pow(x,2) + pow(c3,4)*pow(x,2) + pow(z,2) + 2.*pow(c3,2)*pow(z,2) + pow(c3,4)*pow(z,2) +
         pow(c1,4)*(0.24640028917150708 + 0.9927744742316998*x + pow(x,2) + pow(z,2)) + pow(c2,4)*(2.0027065588744883 + 2.8303403038323767*x + pow(x,2) + pow(z,2)) +
         pow(c1,2)*(0.9956308938140075 + 1.9855489484633997*x + c3*(2.007009288264195 + 1.418210584463385*x) + c2*(1.837565829600677 - 1.4182105844633848*c3)*z +
            2.*pow(x,2) + 2.*pow(z,2) + pow(c2,2)*(1.907775118987983 + 3.8231147780640766*x + 2.*pow(x,2) + 2.*pow(z,2)) +
            pow(c3,2)*(2.249106848045995 + 3.8231147780640766*x + 2.*pow(x,2) + 2.*pow(z,2))) +
         pow(c2,2)*(2.249106848045995 + 3.8231147780640766*x + c3*(0.7039816336702344 + 1.418210584463385*x) + 2.*pow(x,2) + 2.*pow(z,2) +
            pow(c3,2)*(4.5082434332199695 + 5.660680607664753*x + 2.*pow(x,2) + 2.*pow(z,2))))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2),0.5)),
   atan2((-0.39965106940120154 + 0.15471173963046464*c1*c2 + 0.15471173963046464*c3 - 1.*y + (-1.5498930602972543 - 1.*y)*pow(c1,2) +
       (-0.3996510694012015 - 1.*y)*pow(c2,2) - 1.5498930602972543*pow(c3,2) - 1.*y*pow(c3,2))*pow(l5,-1)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
    1.*pow(pow(l5,-2)*(0.09033542197868717 - 0.691430053327601*c3 - 0.6011170334591664*x + 2.3004839817921057*c3*x +
         (c2*(0.691430053327601 - 2.3004839817921057*x) + (-2.3004839817921057 - 0.30942347926092983*c3)*z)*pow(c1,3) +
         (0.30942347926092983 - 2.3004839817921057*c3)*z*pow(c2,3) + 1.4107276196099177*pow(c3,2) - 0.8928105876574031*x*pow(c3,2) -
         2.3004839817921052*c2*(-0.13450364432439346 + 1.*c3)*z*(1. + 1.*pow(c3,2)) +
         c1*((-2.3004839817921057 - 0.30942347926092983*c3)*z*pow(c2,2) + (0.33551817451252575 - 2.3004839817921057*x)*pow(c2,3) -
            0.30942347926092983*(7.434742790970169 + 1.*c3)*z*(1. + 1.*pow(c3,2)) +
            c2*(0.33551817451252575 - 0.04787144475896954*c3 - 2.3004839817921057*x + (0.691430053327601 - 2.3004839817921057*x)*pow(c3,2))) -
         0.33551817451252575*pow(c3,3) + 2.3004839817921057*x*pow(c3,3) + 0.021271282390199764*pow(c3,4) - 0.2916935541982366*x*pow(c3,4) + pow(x,2) +
         2.*pow(c3,2)*pow(x,2) + pow(c3,4)*pow(x,2) + pow(z,2) + 2.*pow(c3,2)*pow(z,2) + pow(c3,4)*pow(z,2) +
         pow(c1,4)*(0.09033542197868716 - 0.6011170334591664*x + pow(x,2) + pow(z,2)) +
         pow(c2,4)*(0.021271282390199764 - 0.2916935541982366*x + pow(x,2) + pow(z,2)) +
         pow(c1,2)*(1.50372748157789 - 1.202234066918333*x + c3*(-0.33551817451252575 + 2.3004839817921057*x) + c2*(0.30942347926092983 - 2.3004839817921057*c3)*z +
            2.*pow(x,2) + 2.*pow(z,2) + pow(c3,2)*(0.1116067043688872 - 0.8928105876574031*x + 2.*pow(x,2) + 2.*pow(z,2)) +
            pow(c2,2)*(1.4107276196099177 - 0.8928105876574031*x + 2.*pow(x,2) + 2.*pow(z,2))) +
         pow(c2,2)*(0.1116067043688872 - 0.8928105876574031*x + c3*(-0.691430053327601 + 2.3004839817921057*x) + 2.*pow(x,2) + 2.*pow(z,2) +
            pow(c3,2)*(1.3655992024009151 - 0.5833871083964732*x + 2.*pow(x,2) + 2.*pow(z,2))))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2),0.5)),
   atan2((-0.46011684551219006 - 1.073494654430803*c1*c2 - 1.073494654430803*c3 - 1.*y + (-0.9012535441765505 - 1.*y)*pow(c1,2) +
       (-0.46011684551219006 - 1.*y)*pow(c2,2) - 0.9012535441765505*pow(c3,2) - 1.*y*pow(c3,2))*pow(l6,-1)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
    1.*pow(pow(l6,-2)*(0.038348887728122616 - 0.1727744704297275*c3 - 0.39165744077253284*x + 0.8822733973287207*c3*x +
         (c2*(0.1727744704297277 - 0.8822733973287208*x) + (-0.8822733973287208 + 2.146989308861606*c3)*z)*pow(c1,3) +
         (-2.146989308861606 - 0.8822733973287208*c3)*z*pow(c2,3) + 0.6917415314020992*pow(c3,2) - 2.9303041904066722*x*pow(c3,2) -
         0.8822733973287209*c2*(2.4334739269733108 + 1.*c3)*z*(1. + 1.*pow(c3,2)) +
         c1*((-0.8822733973287208 + 2.146989308861606*c3)*z*pow(c2,2) + (1.119890246208613 - 0.8822733973287208*x)*pow(c2,3) +
            2.1469893088616057*(-0.4109351610122953 + 1.*c3)*z*(1. + 1.*pow(c3,2)) +
            c2*(1.119890246208613 - 2.30478154618302*c3 - 0.8822733973287208*x + (0.1727744704297277 - 0.8822733973287208*x)*pow(c3,2))) -
         1.1198902462086138*pow(c3,3) + 0.8822733973287207*x*pow(c3,3) + 1.6111818298569953*pow(c3,4) - 2.5386467496341396*x*pow(c3,4) + pow(x,2) +
         2.*pow(c3,2)*pow(x,2) + pow(c3,4)*pow(x,2) + pow(z,2) + 2.*pow(c3,2)*pow(z,2) + pow(c3,4)*pow(z,2) +
         pow(c2,4)*(1.6111818298569953 - 2.5386467496341396*x + pow(x,2) + pow(z,2)) +
         pow(c1,4)*(0.038348887728122505 - 0.39165744077253284*x + pow(x,2) + pow(z,2)) +
         pow(c2,2)*(1.6495307175851175 + c3*(-0.1727744704297275 + 0.8822733973287207*x) - 2.9303041904066722*x + 2.*pow(x,2) + 2.*pow(z,2) +
            pow(c3,2)*(3.416965246622481 - 5.077293499268279*x + 2.*pow(x,2) + 2.*pow(z,2))) +
         pow(c1,2)*(0.2712993623647352 + c3*(-1.1198902462086138 + 0.8822733973287207*x) - 0.7833148815450657*x + c2*(-2.146989308861606 - 0.8822733973287208*c3)*z +
            2.*pow(x,2) + 2.*pow(z,2) + pow(c2,2)*(0.6917415314020992 - 2.9303041904066722*x + 2.*pow(x,2) + 2.*pow(z,2)) +
            pow(c3,2)*(1.6495307175851175 - 2.9303041904066722*x + 2.*pow(x,2) + 2.*pow(z,2))))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-2),0.5));

    double p1 = psivals(0), p2 = psivals(1), p3 = psivals(2), p4 = psivals(3), p5 = psivals(4), p6 = psivals(5);

    phivals << atan2((-0.19582872038626664 + 0.44113669866436017*c1*c2 - 0.44113669866436017*c3 + x + (-0.19582872038626664 + x)*pow(c1,2) +
       (-1.2693233748170698 + x)*pow(c2,2) - 1.2693233748170698*pow(c3,2) + x*pow(c3,2))*pow(cos(p1),-1)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1)*
     pow((0.2500563992522112 - 0.8150914777120721*c3 - 0.3916574407725333*x - 0.8822733973287203*c3*x - 0.9202336910243802*y + 2.1469893088616065*c3*y +
         c1*c2*(-2.1077561943504124 + 0.8822733973287203*x + 2.1469893088616065*y) + c2*(-2.1469893088616065 + 0.8822733973287203*c3)*z +
         c1*(0.8822733973287203 + 2.1469893088616065*c3)*z + 2.4234397807477883*pow(c3,2) - 2.5386467496341396*x*pow(c3,2) - 1.8025070883531005*y*pow(c3,2) +
         pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) + pow(c3,2)*pow(z,2) +
         pow(c1,2)*(0.8506068386189155 - 0.3916574407725333*x - 1.8025070883531005*y + pow(x,2) + pow(y,2) + pow(z,2)) +
         pow(c2,2)*(1.8228893413810843 - 2.5386467496341396*x - 0.9202336910243802*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
      -0.5),(-1.0734946544308033*c2 + 0.44113669866436017*c2*c3 + c1*(0.44113669866436017 + 1.0734946544308033*c3) + z + z*pow(c1,2) + z*pow(c2,2) + z*pow(c3,2))*
     pow(cos(p1),-1)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1)*
     pow((0.2500563992522112 - 0.8150914777120721*c3 - 0.3916574407725333*x - 0.8822733973287203*c3*x - 0.9202336910243802*y + 2.1469893088616065*c3*y +
         c1*c2*(-2.1077561943504124 + 0.8822733973287203*x + 2.1469893088616065*y) + c2*(-2.1469893088616065 + 0.8822733973287203*c3)*z +
         c1*(0.8822733973287203 + 2.1469893088616065*c3)*z + 2.4234397807477883*pow(c3,2) - 2.5386467496341396*x*pow(c3,2) - 1.8025070883531005*y*pow(c3,2) +
         pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) + pow(c3,2)*pow(z,2) +
         pow(c1,2)*(0.8506068386189155 - 0.3916574407725333*x - 1.8025070883531005*y + pow(x,2) + pow(y,2) + pow(z,2)) +
         pow(c2,2)*(1.8228893413810843 - 2.5386467496341396*x - 0.9202336910243802*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
      -0.5)),atan2((-0.30055851672958345 + 1.1502419908960528*c1*c2 - 1.1502419908960528*c3 + x + (-0.3005585167295834 + x)*pow(c1,2) +
       (-0.1458467770991188 + x)*pow(c2,2) - 0.1458467770991188*pow(c3,2) + x*pow(c3,2))*pow(cos(p2),-1)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1)*
     pow((0.25005639925221135 + 0.8150914777120721*c3 - 0.6011170334591669*x - 2.3004839817921057*c3*x +
         c1*c2*(-0.21185675012805605 + 2.3004839817921057*x - 0.30942347926092917*y) - 0.7993021388024031*y - 0.30942347926092917*c3*y +
         c1*(2.3004839817921057 - 0.30942347926092917*c3)*z + c2*(0.30942347926092917 + 2.3004839817921057*c3)*z + 2.423439780747789*pow(c3,2) -
         0.2916935541982377*x*pow(c3,2) - 3.0997861205945085*y*pow(c3,2) + pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) +
         pow(c3,2)*pow(z,2) + pow(c1,2)*(2.4925039203362758 - 0.6011170334591669*x - 3.0997861205945085*y + pow(x,2) + pow(y,2) + pow(z,2)) +
         pow(c2,2)*(0.1809922596637239 - 0.2916935541982377*x - 0.7993021388024031*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
      -0.5),(0.15471173963046458*c2 + c1*(1.150241990896053 - 0.15471173963046458*c3) + 1.150241990896053*c2*c3 + z + z*pow(c1,2) + z*pow(c2,2) + z*pow(c3,2))*
     pow(cos(p2),-1)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1)*
     pow((0.25005639925221135 + 0.8150914777120721*c3 - 0.6011170334591669*x - 2.3004839817921057*c3*x +
         c1*c2*(-0.21185675012805605 + 2.3004839817921057*x - 0.30942347926092917*y) - 0.7993021388024031*y - 0.30942347926092917*c3*y +
         c1*(2.3004839817921057 - 0.30942347926092917*c3)*z + c2*(0.30942347926092917 + 2.3004839817921057*c3)*z + 2.423439780747789*pow(c3,2) -
         0.2916935541982377*x*pow(c3,2) - 3.0997861205945085*y*pow(c3,2) + pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) +
         pow(c3,2)*pow(z,2) + pow(c1,2)*(2.4925039203362758 - 0.6011170334591669*x - 3.0997861205945085*y + pow(x,2) + pow(y,2) + pow(z,2)) +
         pow(c2,2)*(0.1809922596637239 - 0.2916935541982377*x - 0.7993021388024031*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
      -0.5)),atan2((0.49638723711584987 + 0.7091052922316927*c1*c2 - 0.7091052922316927*c3 + x + (0.49638723711584987 + x)*pow(c1,2) +
       (1.4151701519161883 + x)*pow(c2,2) + 1.4151701519161883*pow(c3,2) + x*pow(c3,2))*pow(cos(p3),-1)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1)*
     pow((0.2500563992522111 - 0.8150914777120724*c3 + 0.9927744742316997*x - 1.4182105844633854*c3*x +
         c1*c2*(1.895899444222358 + 1.4182105844633854*x - 1.8375658296006767*y) + 0.12093155222197716*y - 1.8375658296006767*c3*y +
         c1*(1.4182105844633854 - 1.8375658296006767*c3)*z + c2*(1.8375658296006767 + 1.4182105844633854*c3)*z + 2.4234397807477888*pow(c3,2) +
         2.8303403038323767*x*pow(c3,2) - 1.2972790322414083*y*pow(c3,2) + pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) +
         pow(c3,2)*pow(z,2) + pow(c1,2)*(0.6671335110448082 + 0.9927744742316997*x - 1.2972790322414083*y + pow(x,2) + pow(y,2) + pow(z,2)) +
         pow(c2,2)*(2.0063626689551914 + 2.8303403038323767*x + 0.12093155222197716*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1)
       ,-0.5),(0.9187829148003384*c2 + c1*(0.7091052922316927 - 0.9187829148003384*c3) + 0.7091052922316927*c2*c3 + z + z*pow(c1,2) + z*pow(c2,2) + z*pow(c3,2))*
     pow(cos(p3),-1)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1)*
     pow((0.2500563992522111 - 0.8150914777120724*c3 + 0.9927744742316997*x - 1.4182105844633854*c3*x +
         c1*c2*(1.895899444222358 + 1.4182105844633854*x - 1.8375658296006767*y) + 0.12093155222197716*y - 1.8375658296006767*c3*y +
         c1*(1.4182105844633854 - 1.8375658296006767*c3)*z + c2*(1.8375658296006767 + 1.4182105844633854*c3)*z + 2.4234397807477888*pow(c3,2) +
         2.8303403038323767*x*pow(c3,2) - 1.2972790322414083*y*pow(c3,2) + pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) +
         pow(c3,2)*pow(z,2) + pow(c1,2)*(0.6671335110448082 + 0.9927744742316997*x - 1.2972790322414083*y + pow(x,2) + pow(y,2) + pow(z,2)) +
         pow(c2,2)*(2.0063626689551914 + 2.8303403038323767*x + 0.12093155222197716*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1)
       ,-0.5)),atan2((0.49638723711585 - 0.7091052922316924*c1*c2 + 0.7091052922316925*c3 + x + (0.4963872371158499 + x)*pow(c1,2) +
       (1.4151701519161883 + x)*pow(c2,2) + 1.4151701519161883*pow(c3,2) + x*pow(c3,2))*pow(cos(p4),-1)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1)*
     pow((0.25005639925221135 + 0.815091477712072*c3 + 0.9927744742317001*x + 1.418210584463385*c3*x +
         c1*c2*(-1.8958994442223576 - 1.4182105844633848*x - 1.8375658296006767*y) - 0.12093155222197716*y - 1.8375658296006767*c3*y +
         c1*(-1.4182105844633848 - 1.8375658296006767*c3)*z + c2*(1.837565829600677 - 1.4182105844633848*c3)*z + 2.4234397807477888*pow(c3,2) +
         2.8303403038323767*x*pow(c3,2) + 1.2972790322414078*y*pow(c3,2) + pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) +
         pow(c3,2)*pow(z,2) + pow(c2,2)*(2.0063626689551923 + 2.8303403038323767*x - 0.12093155222197716*y + pow(x,2) + pow(y,2) + pow(z,2)) +
         pow(c1,2)*(0.6671335110448082 + 0.9927744742317001*x + 1.2972790322414078*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
      -0.5),(0.9187829148003388*c2 + c1*(-0.7091052922316927 - 0.9187829148003388*c3) - 0.7091052922316927*c2*c3 + z + z*pow(c1,2) + z*pow(c2,2) + z*pow(c3,2))*
     pow(cos(p4),-1)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1)*
     pow((0.25005639925221135 + 0.815091477712072*c3 + 0.9927744742317001*x + 1.418210584463385*c3*x +
         c1*c2*(-1.8958994442223576 - 1.4182105844633848*x - 1.8375658296006767*y) - 0.12093155222197716*y - 1.8375658296006767*c3*y +
         c1*(-1.4182105844633848 - 1.8375658296006767*c3)*z + c2*(1.837565829600677 - 1.4182105844633848*c3)*z + 2.4234397807477888*pow(c3,2) +
         2.8303403038323767*x*pow(c3,2) + 1.2972790322414078*y*pow(c3,2) + pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) +
         pow(c3,2)*pow(z,2) + pow(c2,2)*(2.0063626689551923 + 2.8303403038323767*x - 0.12093155222197716*y + pow(x,2) + pow(y,2) + pow(z,2)) +
         pow(c1,2)*(0.6671335110448082 + 0.9927744742317001*x + 1.2972790322414078*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
      -0.5)),atan2((-0.3005585167295832 - 1.1502419908960528*c1*c2 + 1.1502419908960528*c3 + x + (-0.3005585167295832 + x)*pow(c1,2) +
       (-0.1458467770991183 + x)*pow(c2,2) - 0.1458467770991183*pow(c3,2) + x*pow(c3,2))*pow(cos(p5),-1)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1)*
     pow((0.2500563992522112 - 0.815091477712072*c3 - 0.6011170334591664*x + 2.3004839817921057*c3*x +
         c1*c2*(0.21185675012805427 - 2.3004839817921057*x - 0.30942347926092983*y) + 0.7993021388024031*y - 0.30942347926092983*c3*y +
         c2*(0.30942347926092983 - 2.3004839817921057*c3)*z + c1*(-2.3004839817921057 - 0.30942347926092983*c3)*z + 2.4234397807477888*pow(c3,2) -
         0.2916935541982366*x*pow(c3,2) + 3.0997861205945085*y*pow(c3,2) + pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) +
         pow(c3,2)*pow(z,2) + pow(c2,2)*(0.180992259663724 - 0.2916935541982366*x + 0.7993021388024031*y + pow(x,2) + pow(y,2) + pow(z,2)) +
         pow(c1,2)*(2.4925039203362758 - 0.6011170334591664*x + 3.0997861205945085*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
      -0.5),(0.15471173963046492*c2 + c1*(-1.1502419908960528 - 0.15471173963046492*c3) - 1.1502419908960528*c2*c3 + z + z*pow(c1,2) + z*pow(c2,2) + z*pow(c3,2))*
     pow(cos(p5),-1)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1)*
     pow((0.2500563992522112 - 0.815091477712072*c3 - 0.6011170334591664*x + 2.3004839817921057*c3*x +
         c1*c2*(0.21185675012805427 - 2.3004839817921057*x - 0.30942347926092983*y) + 0.7993021388024031*y - 0.30942347926092983*c3*y +
         c2*(0.30942347926092983 - 2.3004839817921057*c3)*z + c1*(-2.3004839817921057 - 0.30942347926092983*c3)*z + 2.4234397807477888*pow(c3,2) -
         0.2916935541982366*x*pow(c3,2) + 3.0997861205945085*y*pow(c3,2) + pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) +
         pow(c3,2)*pow(z,2) + pow(c2,2)*(0.180992259663724 - 0.2916935541982366*x + 0.7993021388024031*y + pow(x,2) + pow(y,2) + pow(z,2)) +
         pow(c1,2)*(2.4925039203362758 - 0.6011170334591664*x + 3.0997861205945085*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1),
      -0.5)),atan2((-0.19582872038626636 - 0.4411366986643604*c1*c2 + 0.44113669866436034*c3 + x + (-0.19582872038626636 + x)*pow(c1,2) +
       (-1.2693233748170698 + x)*pow(c2,2) - 1.2693233748170698*pow(c3,2) + x*pow(c3,2))*pow(cos(p6),-1)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1)*
     pow((0.25005639925221124 + 0.8150914777120724*c3 - 0.39165744077253295*x + 0.8822733973287207*c3*x + 0.9202336910243805*y + 2.146989308861606*c3*y +
         c1*c2*(2.1077561943504133 - 0.8822733973287208*x + 2.146989308861606*y) + c2*(-2.146989308861606 - 0.8822733973287208*c3)*z +
         c1*(-0.8822733973287208 + 2.146989308861606*c3)*z + 2.4234397807477888*pow(c3,2) - 2.5386467496341396*x*pow(c3,2) + 1.8025070883531011*y*pow(c3,2) +
         pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) + pow(c3,2)*pow(z,2) +
         pow(c2,2)*(1.8228893413810836 - 2.5386467496341396*x + 0.9202336910243805*y + pow(x,2) + pow(y,2) + pow(z,2)) +
         pow(c1,2)*(0.8506068386189165 - 0.39165744077253295*x + 1.8025070883531011*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1)
       ,-0.5),(-1.073494654430803*c2 - 0.44113669866436045*c2*c3 + c1*(-0.44113669866436045 + 1.073494654430803*c3) + z + z*pow(c1,2) + z*pow(c2,2) + z*pow(c3,2))*
     pow(cos(p6),-1)*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1)*
     pow((0.25005639925221124 + 0.8150914777120724*c3 - 0.39165744077253295*x + 0.8822733973287207*c3*x + 0.9202336910243805*y + 2.146989308861606*c3*y +
         c1*c2*(2.1077561943504133 - 0.8822733973287208*x + 2.146989308861606*y) + c2*(-2.146989308861606 - 0.8822733973287208*c3)*z +
         c1*(-0.8822733973287208 + 2.146989308861606*c3)*z + 2.4234397807477888*pow(c3,2) - 2.5386467496341396*x*pow(c3,2) + 1.8025070883531011*y*pow(c3,2) +
         pow(x,2) + pow(c3,2)*pow(x,2) + pow(y,2) + pow(c3,2)*pow(y,2) + pow(z,2) + pow(c3,2)*pow(z,2) +
         pow(c2,2)*(1.8228893413810836 - 2.5386467496341396*x + 0.9202336910243805*y + pow(x,2) + pow(y,2) + pow(z,2)) +
         pow(c1,2)*(0.8506068386189165 - 0.39165744077253295*x + 1.8025070883531011*y + pow(x,2) + pow(y,2) + pow(z,2)))*pow(1. + pow(c1,2) + pow(c2,2) + pow(c3,2),-1)
       ,-0.5));


    iksol << phivals, psivals, lvals;

    return iksol;
}


// #endif CONSTRAINTS_H
#endif

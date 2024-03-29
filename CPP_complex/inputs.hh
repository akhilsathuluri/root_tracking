#ifndef INPUTS_H
#define INPUTS_H

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

/*!
An array of leg values corresponding to the path to be tracked in the task-space
*/
const double legvals_data[] = {1.4458208723874426,
   1.4060252316390276,
   1.3708899598114825,
   1.3403156514115564,
   1.3132546895462032,
   1.2881484000129075,
   1.263443138859881,
   1.2380553048623144,
   1.211745268649149,
   1.185388299791834,
   1.1610732176018594,
   1.1418722234161967,
   1.1311320840917536,
   1.1313661511986288,
   1.1432107771774895,
   1.1650585693451805,
   1.1935872752440149,
   1.2248021573719514,
   1.2550126767829055,
   1.2814337246761802,
   1.3024334715438375,
   1.3175791912107189,
   1.3275829446473018,
   1.3341408835158517,
   1.3395927442126638,
   1.3463695596368885,
   1.356351124710815,
   1.3704022458423257,
   1.38832306204365,
   1.4092141706965846,
   1.4320259156883506,
   1.4560198698916078,
   1.4809871670462658,
   1.5071977387443252,
   1.5351244027325475,
   1.5650237495550017,
   1.5965058230797666,
   1.6282735799507622,
   1.658188646424637,
   1.6836755810855188,
   1.7022879555614236,
   1.7121803954479813,
   1.7123255090072766,
   1.7024866163420196,
   1.683069475196287,
   1.6549821733884242,
   1.6195708145003058,
   1.5786204210694783,
   1.5343440888187254,
   1.489251214647897,
   1.4458208723874426,
   1.568679047566447,
   1.5698603026795732,
   1.5717537721590564,
   1.5731650379454716,
   1.5726066008225281,
   1.568597529116758,
   1.5598431815903437,
   1.5453140163790533,
   1.5242991496985907,
   1.4965057537674302,
   1.462227134595034,
   1.4225297323549373,
   1.3793398277837687,
   1.335295291530446,
   1.2933155479545482,
   1.2560119082445755,
   1.2251895676141775,
   1.2016615000138375,
   1.1854197232441477,
   1.176031384240781,
   1.1730577072601807,
   1.176326536784823,
   1.185958884090406,
   1.2021273711896059,
   1.2246271468186232,
   1.2524595032937231,
   1.2836758266949604,
   1.3156109115465828,
   1.3454062626863312,
   1.3705781982568874,
   1.3894317600928279,
   1.4012736122201563,
   1.4064860381329416,
   1.406526321840845,
   1.4038398024840395,
   1.4015872298272118,
   1.4030714768923358,
   1.4108853342135725,
   1.4260672795356202,
   1.4477183790961154,
   1.473342398908298,
   1.499720315554311,
   1.5238463443849297,
   1.5435557607830577,
   1.5577628253534497,
   1.5664206258378706,
   1.5703389080206165,
   1.5709291224985105,
   1.5698819456773074,
   1.5687764807078677,
   1.568679047566447,1.57865471573012,
   1.5772578276277807,
   1.5778550434070688,
   1.5804763044388688,
   1.5844641972877689,
   1.5888630480863817,
   1.5927561560814267,
   1.595446977387423,
   1.5964815651629782,
   1.5955644938018765,
   1.5924330814940793,
   1.5867498657306927,
   1.5780605576281015,
   1.5658375014045909,
   1.5495866105943508,
   1.528961666381093,1.50383320479944,
   1.4743000080842352,
   1.4406756438352373,
   1.40349862085779,
   1.3635976630869306,
   1.3222045854202875,
   1.2810559418775362,
   1.24237401197219,
   1.2086054101042907,
   1.1818819491644472,
   1.1633656544579651,
   1.1528157801790948,
   1.1486603317592439,
   1.148555947915099,
   1.1501524549925461,
   1.1517679103259142,
   1.1528440961679616,
   1.1541758422111095,
   1.1578953753598673,
   1.1671073725633816,
   1.185063064303903,
   1.2139773945843437,
   1.2539686839009128,
   1.3027423014606274,
   1.3562322726387697,
   1.4097848293457194,
   1.4592416713122303,
   1.5015573360610353,
   1.5349660686605489,
   1.5588967728694139,
   1.5738149545365978,
   1.5810640655554495,
   1.5826765171647856,
   1.5810785225308954,
   1.57865471573012,
   1.3393859531342835,
   1.368575630574579,
   1.3990402579634633,
   1.4288222696814576,
   1.4560204177960676,
   1.4792536849334956,
   1.4979004662831201,
   1.5121226533330083,
   1.5227379402845689,
   1.5309852696873718,
   1.538200251250867,
   1.5454304682575488,
   1.5530814556663117,
   1.5607387986983303,
   1.567277222847329,
   1.5712281052866268,
   1.5712373719657666,
   1.5664278762145638,
   1.5565814524800494,
   1.5421687593469104,
   1.5242967920395236,
   1.5046133362009289,
   1.4851492680829543,
   1.468043545906806,
   1.4551244591464778,
   1.4474256187466994,
   1.4448311580040185,
   1.4460444323070634,
   1.448908150047061,
   1.4509042076548107,
   1.449609006197701,
   1.4429911308797956,
   1.4295823954903513,
   1.4086219983903405,
   1.3802505848254434,
   1.3457450507661326,
   1.3076700025348424,
   1.269733252724591,
   1.2361677233858674,
   1.2106988517041266,
   1.1955019779514906,
   1.1906976421825644,
   1.1946459606372457,
   1.2048048595183738,
   1.2186701808472291,
   1.2344455213081487,
   1.2513377907280816,
   1.2695107554340237,
   1.2897514909301053,
   1.3129260105501703,
   1.3393859531342835,
   1.152479064455009,
   1.1595753765506795,
   1.1721049992320172,
   1.18887065543626,
   1.2082237126188393,
   1.2287137112670097,
   1.2496153967682926,
   1.271205099575883,
   1.294755760043515,1.32223440782205,
   1.3556932835577098,
   1.39644398567447,
   1.4442926386071655,
   1.4972251643060046,
   1.551765418849317,
   1.6038422568459851,
   1.6497250895566526,
   1.6866458219427636,
   1.7129988860704202,
   1.728242084930226,
   1.7326875702832472,
   1.7273157444390375,
   1.7136498271758644,
   1.6936489673279618,
   1.6695400137505436,
   1.643533013863014,
   1.6174536287772314,
   1.5924224975496022,
   1.568726664307539,
   1.5459313915737352,
   1.5231548106995787,
   1.4993786035995285,
   1.4737102502578532,
   1.445581392985175,
   1.4149016475530765,
   1.3821720077434474,
   1.348516891487555,
   1.3155607836244794,
   1.2851041180302696,
   1.2586598005680696,
   1.237030707426415,
   1.2201240883552567,
   1.2070750031070991,
   1.1965923286574918,
   1.187374323174641,1.17847893833186,
   1.1695966986532988,
   1.1611935747797215,
   1.154467481718555,
   1.1510521884154628,
   1.152479064455009,
   1.286875291086054,
   1.3058097061451694,
   1.3282124617412332,
   1.3504796245654223,
   1.3692233542469587,
   1.3817554335147053,
   1.386320639522505,
   1.3821931794665863,
   1.3697650289970196,
   1.3506710908047719,
   1.3278655903434193,
   1.3054421811294783,
   1.2879748728406808,
   1.27937368575113,
   1.2816708299306843,
   1.294416757562172,
   1.3150661385423303,
   1.3400693303583646,
   1.3660265964756986,
   1.3904549898083338,
   1.4120921979464427,
   1.430861208929243,1.44761804576645,
   1.463739689938029,
   1.4805963056788591,
   1.4990200036491172,
   1.5189598125080592,
   1.5394789171783798,
   1.559077833967764,1.57615099355277,
   1.5893554869769675,
   1.5977829431774786,
   1.6009558121152316,
   1.5987316031269319,
   1.5911944867441836,
   1.5785812388458356,
   1.561254321091324,
   1.5397076612593057,
   1.5145750608784607,
   1.4866145773003316,
   1.456665922936953,
   1.4256063286923566,
   1.3943406361440123,
   1.3638457627512215,
   1.335259281176652,1.30996952184303,
   1.289635087646386,
   1.2760442639628007,
   1.2707508478933789,
   1.2745312491735887,
   1.286875291086054};

   /*!
   Leg values corresponding to the path to be tracked in the task-space typecasted into a matrix
   */
Matrix<double, 51, 6> legvals(legvals_data);
// MatrixXcd Clegvals(51,6);
// Clegvals = legvals.cast<std::complex<double>>();

#endif

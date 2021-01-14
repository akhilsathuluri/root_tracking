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
  1.437740549618771,
  1.4301750144554164,
  1.4232526924372804,
  1.4170929249620312,
  1.4118034164194326,
  1.4074778134381123,
  1.4041935080226569,
   1.4020097561628253,1.4009661964384608,1.401081839037638,1.4023545751649704,1.4047612316468712,1.4082581679548272,1.4127823855547674,1.4182530950961945,
   1.4245736676741485,1.431633883666936,1.4393123870021307,1.4474792537816439,1.455998590946067,1.4647310915770808,1.4735364868612368,1.4822758490907673,
   1.4908137140387256,1.4990200036491172,1.5067717406273915,1.5139545549318294,1.5204639883408673,1.5262066073746174,1.5311009371552886,1.5350782296192125,
   1.5380830791600333,1.5400738975762298,1.5410232583605405,1.5409181181004703,1.5397599202100503,1.537564583499942,1.5343623753077038,1.5301976661325702,
   1.5251285600307438,1.5192263925237772,1.5125750855769478,1.5052703474770286,1.4974187043853568,1.4891363502159405,1.4805478025848398,1.471784355216493,
   1.4629823216886901,1.4542810720096184,1.4458208723874426,1.568679047566447,1.5618733068576758,1.5557552821919862,1.5504301250135062,1.5459903978985194,
   1.5425140239685817,1.5400625018083371,1.5386794581333458,1.5383896002329995,1.5391981149754592,1.5410905419262735,1.544033126463039,1.5479736365720547,
   1.552842606272007,1.5585549510989358,1.5650118881316748,1.5721030853622748,1.5797089629052978,1.5877030711033164,1.595954477106505,1.6043301008241522,
   1.6126969520519938,1.620924231955427,1.6288852730191665,1.6364593013933115,1.643533013863014,1.6500019682585494,1.655771790997585,1.6607592087302185,1.6648929129338004,
   1.6681142670068714,1.6703778651716694,1.6716519515288135,1.671918706105823,1.671174402863929,1.6694294425093272,1.6667082607092314,1.6630491100414042,1.658503711807115,
   1.6531367718179315,1.6470253525486935,1.6402580927787045,1.6329342651951075,1.6251626626097202,1.6170603046687282,1.6087509594468574,1.600363478330416,
   1.5920299482708633,1.583883672898011,1.5760570030347674,1.568679047566447,1.57865471573012,1.5786858762809939,1.5777245005249374,1.5757839566643683,1.5728912995025544,
   1.5690870354022246,1.564424762187769,1.5589706773664527,1.5528029462807174,1.546010920502448,1.5386941961230856,1.5309615017953166,1.5229294076831317,
   1.5147208491273614,1.5064634630656981,1.4982877412350304,1.4903250119970213,1.4827052721671505,1.4755549011725202,1.468994301628968,1.4631355221486841,
   1.4580799287529902,1.4539159993990467,1.4507173205490451,1.4485408643110766,1.4474256187466994,1.447391632327658,1.4484395167712556,1.450550431808698,1.453686552613863,
   1.4577919977259683,1.4627941744329427,1.468605481541178,1.4751252974995497,1.4822421755530997,1.4898361668714928,1.4977811967446515,1.505947426876223,
   1.5142035472551558,1.5224189527686902,1.5304657715107652,1.5382207227564222,1.5455667922121474,1.5523947200782944,1.5586043035622712,1.5641055198188865,
   1.5688194780371918,1.5726792107714411,1.575630314877259,1.577631451804296,1.57865471573012,1.3393859531342835,1.338179740423019,1.335819938882744,1.3323377649639636,
   1.3277795655248554,1.3222065258817914,1.3156942635123556,1.3083322922268081,1.3002233381790358,1.2914824860886722,1.2822361318113487,1.2726207163588856,
   1.2627812171855364,1.252869375661966,1.2430416458342528,1.2334568594299642,1.2242736160453973,1.2156474255636882,1.2077276515710331,1.2006543285689855,
   1.1945549500087291,1.1895413457324064,1.1857067829682377,1.1831234313042052,1.1818403264775703,1.1818819491644474,1.1832475039211232,1.1859109427345245,
   1.1898217317165696,1.1949063137287428,1.2010701795189946,1.208200429636932,1.2161686915599879,1.224834251671936,1.2340472686815562,1.2436519511067623,1.253489603270958,
   1.2634014685513082,1.2732313225579168,1.282827790444165,1.2920463804133953,1.3007512391957632,1.3088166448874967,1.316128258508609,1.3225841585709188,
   1.3280956835348232,1.332588105907944,1.336001159444594,1.338289437869533,1.3394226800848106,1.3393859531342835,1.152479064455009,1.1625014909874525,1.1730905885761809,
   1.1840664660079032,1.1952477996243895,1.2064550092126523,1.2175130671448406,1.228253931267404,1.2385186127714847,1.2481589044104724,1.257038802770728,
   1.2650356619754752,1.2720411163995822,1.2779618077793342,1.2827199484102165,1.2862537476136349,1.2885177237899088,1.2894829194524782,1.2891370318097328,
   1.2874844667756504,1.2845463197273148,1.280360281823172,1.2749804661701052,1.2684771435057347,1.260936372300444,1.252459503293723,1.2431625335609993,1.2331752804819927,
   1.22264034183132,1.2117118052123232,1.2005536689975802,1.1893379388253849,1.1782423696863136,1.1674478349016066,1.1571353208413877,1.1474825705545975,
   1.1386604302109224,1.1308289877725795,1.1241336305151752,1.1187011823308062,1.1146363075675099,1.1120183797388983,1.1108990061378963,1.1113003710134555,
   1.1132145117672834,1.1166035795540257,1.1214010656886053,1.1275139078000385,1.1348253336092669,1.1431982620734027,1.152479064455009,1.286875291086054,
   1.2963731932518796,1.3061268289680983,1.315978422329431,1.3257721582215505,1.3353565401461447,1.3445864777022527,1.353325104610952,1.361445338338907,1.3688311990590512,
   1.3753789093911828,1.380997797686222,1.3856110271416466,1.3891561712816816,1.3915856537281097,1.3928670670480097,1.3929833820203195,1.391933055061295,1.389730037876076,
   1.386403689704945,1.3819985888322872,1.3765742363456357,1.3702046415109121,1.3629777746447167,1.354994870148935,1.3463695596368885,1.337226813125795,1.3277016655010963,
   1.3179377063761346,1.3080853146523803,1.298299625130333,1.2887382239636156,1.2795585829317864,1.2709152594252315,1.2629569091648396,1.2558231808366036,
   1.2496415840993966,1.2445244422622304,1.2405660553612874,1.237840205463837,1.236398131451888,1.2362670841447967,1.2374495448915785,1.239923153969473,1.2436413530826158,
   1.248534703718808,1.2545128049382657,1.261466704398732,1.2692716776755009,1.2777902440520963,1.286875291086054};

   /*!
   Leg values corresponding to the path to be tracked in the task-space typecasted into a matrix
   */
Matrix<double, 51, 6> legvals(legvals_data);
// MatrixXcd Clegvals(51,6);
// Clegvals = legvals.cast<std::complex<double>>();

#endif

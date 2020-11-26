#ifndef UTILS_H
#define UTILS_H

#include<iostream>
#include<eigen3/Eigen/Dense>
#include<gsl/gsl_linalg.h>
#include<math.h>

using namespace Eigen;

/*!
The linearSolve function is a simple, usable wraper around the GSLs
linear solver using the Eigen library.
@param Amat Input matrix
@param bvec Input vector
@todo Add assertions for square matrix verification and verification of
sizes of A and b
*/
VectorXd linearSolve(MatrixXd Amat, VectorXd bvec){
  // Transposing so that a_dat gets the right sequence
  int bsize = bvec.size();
  Amat.transposeInPlace();
  double *a_data = Amat.data(), *b_data=bvec.data();
  gsl_matrix_view m
    = gsl_matrix_view_array (a_data, Amat.rows(), Amat.cols());

  gsl_vector_view b
    = gsl_vector_view_array (b_data, bsize);

  gsl_vector *x = gsl_vector_alloc (bsize);

  int s;

  gsl_permutation * p = gsl_permutation_alloc (bsize);

  gsl_linalg_LU_decomp (&m.matrix, p, &s);

  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

  VectorXd xval(bsize);
  for (size_t i = 0; i < bsize; i++) {
    xval(i) = gsl_vector_get(x, i);
  }

  gsl_permutation_free (p);
  gsl_vector_free (x);
  return xval;
}

/*!
The s1Dist returns the distance between two elements in the S1 space.
@param theta1 First element
@param theta2 Second element

@todo Can do more effecient computation if the angles mapping
is saved in the memory
@todo Assert sizes of theta1 and theta2 to be same
*/
VectorXd s1Dist(VectorXd theta1, VectorXd theta2){
  for (int i = 0; i < theta1.size(); i++) {
    theta1(i) = fmod(theta1(i), 2*M_PI);
    theta2(i) = fmod(theta2(i), 2*M_PI);

    if (theta1(i) >= M_PI){
      theta1(i) = 2*M_PI - theta1(i);
    }
    if (theta2(i) >= M_PI){
      theta2(i) = 2*M_PI - theta2(i);
    }
  }
  return (theta1 - theta2).cwiseAbs();
}

#endif

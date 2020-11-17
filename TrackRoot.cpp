# include <iostream>
# include <eigen3/Eigen/Dense>
# include "test.h"
#include <stdio.h>
#include <gsl/gsl_linalg.h>

using namespace Eigen;

# define etalength 12
# define philength 12

VectorXd TrackRoot(VectorXd thetai, VectorXd phii);
VectorXd LinearSolve(MatrixXd Amat, VectorXd bvec);

int main(int argc, char const *argv[]) {
	VectorXd etavalnew(12), etaval(12), q(18);
	MatrixXd lvals(51,6);
	lvals << 1.00575, 1.20209, 1.47276, 1.47276, 1.20209, 1.00575, \
	1.00509, 1.2049, 1.47547, 1.46993, 1.1993, 1.00652, 1.00454, 1.20774, \
	1.47808, 1.467, 1.19655, 1.00739, 1.0041, 1.21061, 1.48057, 1.46396, \
	1.19384, 1.00836, 1.00377, 1.2135, 1.48296, 1.46082, 1.19116, \
	1.00944, 1.00354, 1.21641, 1.48522, 1.45758, 1.18851, 1.01062, \
	1.00343, 1.21935, 1.48738, 1.45424, 1.18591, 1.0119, 1.00344, 1.2223, \
	1.48942, 1.45079, 1.18335, 1.01328, 1.00355, 1.22527, 1.49134, \
	1.44725, 1.18083, 1.01476, 1.00378, 1.22825, 1.49314, 1.44362, \
	1.17836, 1.01633, 1.00412, 1.23124, 1.49482, 1.43989, 1.17593, \
	1.01799, 1.00457, 1.23424, 1.49639, 1.43608, 1.17355, 1.01974, \
	1.00514, 1.23726, 1.49784, 1.43217, 1.17122, 1.02159, 1.00581, \
	1.24027, 1.49916, 1.42818, 1.16894, 1.02352, 1.00661, 1.24329, \
	1.50037, 1.4241, 1.16671, 1.02553, 1.00751, 1.24632, 1.50145, \
	1.41994, 1.16453, 1.02762, 1.00853, 1.24934, 1.50241, 1.41571, \
	1.16241, 1.0298, 1.00965, 1.25236, 1.50325, 1.41139, 1.16034, \
	1.03205, 1.01089, 1.25538, 1.50397, 1.407, 1.15833, 1.03438, 1.01224, \
	1.25839, 1.50456, 1.40254, 1.15638, 1.03678, 1.0137, 1.26139, \
	1.50504, 1.398, 1.15448, 1.03925, 1.01527, 1.26438, 1.50539, 1.3934, \
	1.15265, 1.04179, 1.01694, 1.26736, 1.50562, 1.38874, 1.15087, \
	1.04439, 1.01872, 1.27033, 1.50573, 1.38401, 1.14915, 1.04706, \
	1.02061, 1.27328, 1.50572, 1.37922, 1.14749, 1.04979, 1.02259, \
	1.27621, 1.50558, 1.37438, 1.1459, 1.05258, 1.02469, 1.27913, \
	1.50533, 1.36948, 1.14436, 1.05542, 1.02688, 1.28203, 1.50496, \
	1.36453, 1.14288, 1.05831, 1.02917, 1.2849, 1.50446, 1.35953, \
	1.14147, 1.06126, 1.03156, 1.28775, 1.50385, 1.35448, 1.14012, \
	1.06426, 1.03405, 1.29058, 1.50312, 1.34939, 1.13882, 1.0673, \
	1.03664, 1.29338, 1.50227, 1.34425, 1.13759, 1.07038, 1.03931, \
	1.29615, 1.50131, 1.33908, 1.13642, 1.07351, 1.04208, 1.29889, \
	1.50023, 1.33387, 1.13531, 1.07667, 1.04494, 1.3016, 1.49904, \
	1.32863, 1.13425, 1.07988, 1.04789, 1.30428, 1.49774, 1.32335, \
	1.13326, 1.08311, 1.05092, 1.30693, 1.49632, 1.31805, 1.13232, \
	1.08638, 1.05404, 1.30954, 1.49479, 1.31272, 1.13144, 1.08968, \
	1.05724, 1.31212, 1.49315, 1.30737, 1.13062, 1.09301, 1.06052, \
	1.31466, 1.4914, 1.302, 1.12986, 1.09637, 1.06388, 1.31716, 1.48955, \
	1.29661, 1.12915, 1.09975, 1.06732, 1.31963, 1.48758, 1.2912, \
	1.12849, 1.10315, 1.07083, 1.32205, 1.48552, 1.28578, 1.12788, \
	1.10657, 1.07441, 1.32443, 1.48335, 1.28035, 1.12733, 1.11002, \
	1.07807, 1.32677, 1.48108, 1.27491, 1.12683, 1.11348, 1.08179, \
	1.32907, 1.47871, 1.26947, 1.12637, 1.11695, 1.08558, 1.33133, \
	1.47624, 1.26402, 1.12596, 1.12044, 1.08944, 1.33353, 1.47367, \
	1.25857, 1.1256, 1.12395, 1.09336, 1.3357, 1.471, 1.25312, 1.12529, \
	1.12746, 1.09733, 1.33782, 1.46825, 1.24767, 1.12502, 1.13098, \
	1.10137, 1.33989, 1.4654, 1.24223, 1.12479, 1.13452;

	q << -0.0415145, -0.0835474, 0.520437, 0.520437, -0.0835474, -0.0415145, 0.475166, 0.338916, -0.0410677, 0.0410677, -0.338916, -0.475166, 1.00575, 1.20209, 1.47276, 1.47276, 1.20209, 1.00575;

	VectorXd temp(12);
	MatrixXd sols(12,1);
	sols << q.head(12);
	temp = q.head(12);
	for (int i = 0; i < lvals.rows(); i++) {
	temp = TrackRoot(lvals.row(i), temp);
	sols.conservativeResize(sols.rows(), sols.cols()+1);
	sols.col(sols.cols()-1) = temp;
}
	std::cout<<sols<<std::endl;
	return 0;
}

//Root Tracking via Newton-Raphson
VectorXd TrackRoot(VectorXd thetai, VectorXd phii){
	VectorXd qval(18), tempphi(12), dphi(12), nvec(12);
	MatrixXd Jmat(12,12);
	int loopcounter=0;
	qval << phii, thetai;
	nvec = eta(qval);
	loopcounter = 0;
	tempphi = phii;
	while ((nvec.cwiseAbs()).maxCoeff()>=pow(10, -10)) {
		if (loopcounter>=100)
			break;
		loopcounter++;
		// qval << tempphi, thetai;
		Jmat = Jetaphi(qval);
		dphi = LinearSolve(Jmat, nvec);
		tempphi = tempphi - dphi;
		qval << tempphi, thetai;
		nvec = eta(qval);
	}
	return tempphi;
}

//Fixed Dimension Linear Solve
VectorXd LinearSolve(MatrixXd Amat, VectorXd bvec){
  // Transposing so that a_dat gets the right sequence
  Amat.transposeInPlace();
  double *a_data = Amat.data(), *b_data=bvec.data();
  gsl_matrix_view m
    = gsl_matrix_view_array (a_data, etalength, philength);

  gsl_vector_view b
    = gsl_vector_view_array (b_data, philength);

  gsl_vector *x = gsl_vector_alloc (philength);

  int s;

  gsl_permutation * p = gsl_permutation_alloc (philength);

  gsl_linalg_LU_decomp (&m.matrix, p, &s);

  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

  VectorXd xval(philength);
  for (size_t i = 0; i < philength; i++) {
    xval(i) = gsl_vector_get(x, i);
  }

  gsl_permutation_free (p);
  gsl_vector_free (x);
  return xval;
}

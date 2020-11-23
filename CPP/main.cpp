#include<iostream>
#include<eigen3/Eigen/Dense>
#include "root_tracker.h"
#include "inputs.h"
#include<functional>

using namespace Eigen;

/*!
@todo Write an example root-tracker for SRSPM
*/
int main(int argc, char const *argv[]) {
  RootTracker rt;
  rt.methods();

  // VectorXd q(24);
//   q << 0.20000000000000112, 0, 1.2799999999999998, 0.19999999999999996, 0, 0,
//   0.003056244705087313, -0.06688548576654774, 0.45696277640159066, \
// 0.5469562635040984,
//   -0.09469005512919006, 0.0034901115067618653, 0.3362762851850709, \
// 0.2868906143657115,
//   -0.021027360477035077, 0.024784412575612605, -0.3953833637520812, \
// -0.37979521343323663,
//   1.4458208723874424, 1.568679047566447, 1.57865471573012, \
// 1.339385953134284,
//   1.1524790644550087, 1.2868752910860537;

  // Initial solution
  VectorXd phii(18), temp(18);
  phii << 0.2,0,1.28,0.2,0,0,
   0.003056244705087313,-0.06688548576654774,0.45696277640159066,0.5469562635040984,
   -0.09469005512919006,0.0034901115067618653,0.3362762851850709,0.2868906143657115,
   -0.021027360477035077,0.024784412575612605,-0.3953833637520812,-0.37979521343323663;

  MatrixXd sols(18,1);
  sols << phii;
  temp = phii;

  for (int i = 1; i < legvals.rows(); i++){
    temp = rt.NRTracker(legvals.row(i), temp, etaext, Jetaextphi);
    sols.conservativeResize(sols.rows(), sols.cols()+1);
  	sols.col(sols.cols()-1) = temp;
  }

  std::cout << sols.rows() << " " << sols.cols() << std::endl;
  std::cout << sols.col(50) << std::endl;

  return 0;
}

#include<iostream>
#include<eigen3/Eigen/Dense>
#include "root_tracker.hh"
#include "inputs.hh"
#include <fstream>
#include<functional>

using namespace Eigen;

/*!
The saveData function takes in an input `MatrixXd` and saves it in a .txt file

@param Amat Input matrix
@param file_name Name of the file to save the input matrix in
*/
void saveData(MatrixXd A, const char file_name[]){
  std::ofstream file(file_name);
  if(file.is_open()){
    file << A;
  }
  std::cout << "Saved data to file " << file_name << std::endl;
}

/*!
Example problem demonstrating the usage of the root trackers
for solving a task space path following problem of an SRSPM.
*/
int main(int argc, char const *argv[]) {
  RootTracker rt;
  rt.Methods();

  // Initial solution
  VectorXd phii(18), tempNR(18), tempDM(18);
  phii << 0.2,0,1.28,0.2,0,0,
   0.003056244705087313,-0.06688548576654774,0.45696277640159066,0.5469562635040984,
   -0.09469005512919006,0.0034901115067618653,0.3362762851850709,0.2868906143657115,
   -0.021027360477035077,0.024784412575612605,-0.3953833637520812,-0.37979521343323663;

  MatrixXd solsNR(18,1);
  solsNR << phii;
  tempNR = phii;

  // Tracking using NRTracker
  for (int i = 1; i < legvals.rows(); i++){
    tempNR = rt.NRTracker(legvals.row(i), tempNR, etaext, Jetaextphi);
    solsNR.conservativeResize(solsNR.rows(), solsNR.cols()+1);
  	solsNR.col(solsNR.cols()-1) = tempNR;
  }

  // Display results
  // std::cout << solsNR.rows() << " " << solsNR.cols() << std::endl;
  // std::cout << solsNR.col(50) << std::endl;

  MatrixXd solsDM(18,1);
  solsDM << phii;
  tempDM = phii;

  // Tracking using DMTracker
  for (int i = 1; i < legvals.rows(); i++){
    // Without NR step correction
    // tempDM = rt.DMTracker(legvals.row(i-1), legvals.row(i), tempDM, Jetaexttheta, Jetaextphi);

    // With NR step correction
    /*!
    Choosing small enough `eps` makes DMTracker work like an NRTracker
    */
    tempDM = rt.DMTracker(legvals.row(i-1), legvals.row(i), tempDM, Jetaexttheta, Jetaextphi, 0.05, etaext);

    solsDM.conservativeResize(solsDM.rows(), solsDM.cols()+1);
  	solsDM.col(solsDM.cols()-1) = tempDM;
  }

  // Display results
  // std::cout << solsDM.rows() << " " << solsDM.cols() << std::endl;
  // std::cout << solsDM.col(50) << std::endl;

  saveData(solsDM, "DMTracker_NR.txt");

  return 0;
}

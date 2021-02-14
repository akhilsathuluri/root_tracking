#include<iostream>
#include<eigen3/Eigen/Dense>
#include "root_tracker.hh"
#include "inputs.hh"
#include "manipulator.hh"
#include<fstream>
#include<functional>
#include<chrono>

using namespace Eigen;
// #define constants like 18 and 6 here
/*!
The saveData function takes in an input `MatrixXd` and saves it in a .txt file

@param A Input matrix
@param file_name Name of the file to save the input matrix in
*/

// Define all teh constants so that for a new problem everything can be setup from here
#define ySize 18
#define iters 51
#define xSize 6
#define simIters 18
#define branchesSize 14
#define epsval 0.05

void saveData(MatrixXd A, const char file_name[]){
  std::ofstream file(file_name);
  if(file.is_open()){
    file << A;
  }
  std::cout << "Saved data to file " << file_name << std::endl;
}

void saveCData(MatrixXcd A, const char file_name[]){
  std::ofstream file(file_name);
  if(file.is_open()){
    file << A;
  }
  std::cout << "Saved data to file " << file_name << std::endl;
}

/*!
Example problem demonstrating the usage of the root trackers
for solving a task space path following problem of an SRSPM using
all the discussed root-trackers.
*/
int main(int argc, char const *argv[])
{
  RootTracker rt;
  rt.Methods();

  // Add the input dimensions of your input vector
  MatrixXcd Clegvals(iters,xSize);
  Clegvals = legvals.cast<std::complex<double>>();
  // Declare input Variables
  MatrixXcd inputVars;
  inputVars = Clegvals;
  // Initial solution
  VectorXd tempNR(ySize), tempDM(ySize), tempNN(ySize);

  int singular_iter = simIters;
  VectorXcd phiic(ySize), tempNRC(ySize), tempDMC(ySize), tempNNC(ySize);
  phiic = phii.cast<std::complex<double>>();
  MatrixXcd solsNNC(ySize,1), solsDMC(ySize,1), solsNRC(ySize,1), tempFKC(branchesSize, ySize);

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++NEWTON RAPHSON METHOD+++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  solsNRC << phii2;
  tempNRC = phii2;

  // Tracking using NRTracker
  for (int i = 1; i < singular_iter; i++){
    tempNRC = rt.NRCTracker(Clegvals.row(i), tempNRC, etaextC, JetaextphiC);
    std::cout << tempNRC << '\n';

    solsNRC.conservativeResize(solsNRC.rows(), solsNRC.cols()+1);
  	solsNRC.col(solsNRC.cols()-1) = tempNRC;
  }

  // ToDo: Make saveData an overloaded function
  saveCData(solsNRC.transpose(), "NRCTracker_2.txt");

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++DAVIDENKOS METHOD+++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  solsDMC << phii2;
  tempDMC = phii2;

  // Tracking using DMTracker
  for (int i = 1; i < singular_iter; i++){
    // Without NR step correction
    // tempDMC = rt.DMCTracker(Clegvals.row(i-1), Clegvals.row(i), tempDMC, Jetaexttheta, Jetaextphi);

    // With NR step correction
    /*!
    Choosing small enough `eps` makes DMTracker work like an NRTracker
    */
    tempDMC = rt.DMCTracker(Clegvals.row(i-1), Clegvals.row(i), tempDMC, JetaextthetaC, JetaextphiC, epsval, etaextC);

    solsDMC.conservativeResize(solsDMC.rows(), solsDMC.cols()+1);
  	solsDMC.col(solsDMC.cols()-1) = tempDMC;
  }


  saveCData(solsDMC.transpose(), "DMCTracker_NRC_2.txt");

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++NEAREST NEIGHBOUR METHOD++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  solsNNC << phii2;
  tempNNC = phii2;
  // Initialise all the known roots including the complex root
  tempFKC = initsols2;

  // Tracking using NNCTracker
  for (int i = 1; i < singular_iter; i++){
    for (int j = 0; j < tempFKC.rows(); j++){
      tempFKC.row(j) = rt.NRCTracker(Clegvals.row(i), tempFKC.row(j), etaextC, JetaextphiC);
    }

    tempNNC = rt.NNCTracker(tempNNC, tempFKC, 6);

    solsNNC.conservativeResize(solsNNC.rows(), solsNNC.cols()+1);
  	solsNNC.col(solsNNC.cols()-1) = tempNNC;
  }

  saveCData(solsNNC.transpose(), "NNCTracker_NRC.txt");

  // Verifying the tracked roots with the NR roots
  // std::cout << (solsNNC - solsNRC).maxCoeff() << std::endl;

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++SINGULARITY EVENT IDENTIFICATION++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // solsNRC << phii2;
  // tempNRC = phii2;
  //
  // // Tracking using NRTracker with SingularityEventIdentification
  // for (int i = 1; i < singular_iter; i++){
  //   tempNRC = rt.NRCTracker(Clegvals.row(i), tempNRC, etaextC, JetaextphiC);
  //   std::cout << tempNRC << '\n';
  //
  //   solsNRC.conservativeResize(solsNRC.rows(), solsNRC.cols()+1);
  // 	solsNRC.col(solsNRC.cols()-1) = tempNRC;
  // }
  //
  // // ToDo: Make saveData an overloaded function
  // saveCData(solsNRC.transpose(), "NRCTracker_2.txt");

  VectorXd qx(6);
  qx << 0.2, 0, 1.28, 0.2, 0, 0;
  std::cout << computeIK(qx) << std::endl;

  return 0;
}

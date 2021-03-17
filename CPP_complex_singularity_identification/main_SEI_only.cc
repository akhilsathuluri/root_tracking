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

  VectorXd tempNR(ySize);

  // int singular_iter = simIters;
  VectorXd tempNRC(ySize);
  MatrixXd solsNR(ySize,1);

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++ trackAllBranches check ++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   MatrixXd allRoots(initsols.rows(), initsols.cols());
//   allRoots = rt.trackAllBranches(legvals.row(2), initsols, etaext, Jetaextphi);
//
// std::cout << allRoots<< '\n';

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++ NEWTON RAPHSON METHOD +++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++ WITH SEI +++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  VectorXd ys(ySize);
  ys << 0.2, 0., 1.28, 0.2, 0., 0., 0.00305624, -0.0668855, 0.456963, \
0.546956, -0.0946901, 0.00349011, 0.336276, 0.286891, -0.0210274, \
0.0247844, -0.395383, -0.379795;
  // solsNR << ys;
  // tempNR = ys;
  tempNR = (initsols.row(3)).transpose();

  double stepsize = M_PI/150, alpha = 0;
  // Declare variables to save the distance and alpha histories
  MatrixXd disthist(3, initsols.rows()), currentroots(initsols.rows(), initsols.cols());
  VectorXd alphahist(3);
  // Initialise history Variables
  alphahist << 0,0,0;
  disthist = MatrixXd::Zero(3, initsols.rows());
  currentroots = initsols;
  for (alpha = stepsize; alpha <= M_PI/3; alpha = alpha+stepsize){
    currentroots = rt.trackAllBranches(computeXfromParam(alpha), currentroots, etaext, Jetaextphi);
    if(rt.SEI(currentroots, alpha, 2, computeXfromParam, alphahist, disthist, etaext, Jetaextphi)==1){
      std::cout << "Singularity predicted. Exiting simulation." << '\n';
      break;
    }
    else{
      tempNR = rt.NRTracker(computeXfromParam(alpha), tempNR, etaext, Jetaextphi);

      solsNR.conservativeResize(solsNR.rows(), solsNR.cols()+1);
      solsNR.col(solsNR.cols()-1) = tempNR;
    }
  }

  // ToDo: Make saveData an overloaded function
  saveCData(solsNR.transpose(), "NRCTracker_SEI.txt");

  return 0;
}

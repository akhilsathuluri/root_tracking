#include<iostream>
#include<eigen3/Eigen/Dense>
#include "root_tracker.hh"
#include "cdpr_inputs.hh"
#include "cdpr.hh"
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

// Define all the constants so that for a new problem everything can be setup from here
#define ySize 6
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
for solving a task space path following problem of a CDPR using
all the discussed root-trackers.
*/
int main(int argc, char const *argv[])
{
  RootTracker rt;
  rt.Methods();

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++ SINGULARITY EVENT IDENTIFICATION ++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // Initial solution
  VectorXd tempNR(ySize);
  MatrixXd solsNR(ySize,1);

  int selectedroot = 1;
  tempNR = initsols.row(selectedroot-1);
  std::cout << tempNR << '\n';

  double stepsize = 0.001, alpha = 0;
  // Declare variables to save the distance and alpha histories
  MatrixXd disthist(3, initsols.rows()), currentroots(initsols.rows(), initsols.cols());
  VectorXd alphahist(3);
  // Initialise history Variables
  alphahist << 0,0,0;
  disthist = MatrixXd::Zero(3, initsols.rows());
  currentroots = initsols;
  // std::cout << computeXfromParam(0) << '\n';
  // VectorXd q(9);
  // std::cout << computeXfromParam(0) << '\n';
  // q << tempNR, computeXfromParam(0);
  // std::cout << eta(q) << '\n';
  // std::cout << Jetaphi(q) << '\n';
  std::cout << rt.NRTracker(computeXfromParam(0.01), tempNR, eta, Jetaphi) << '\n';
  for (alpha = stepsize; alpha <= 1; alpha = alpha+stepsize){
    // if(rt.SEI(currentroots, alpha, selectedroot, computeXfromParam, alphahist, disthist, eta, Jetaphi, computeXfromParam)==1){
    //   std::cout << "Singularity predicted. Exiting simulation." << '\n';
    //   break;
    // }
    // else{
      tempNR = rt.NRTracker(computeXfromParam(alpha), tempNR, eta, Jetaphi);

      // solsNR.conservativeResize(solsNR.rows(), solsNR.cols()+1);
      // solsNR.col(solsNR.cols()-1) = tempNR;
      //
      // std::cout << currentroots << '\n';
      // currentroots = rt.trackAllBranches(computeXfromParam(alpha), currentroots, eta, Jetaphi);
      // std::cout << tempNR << '\n';
    // }
  }

  // ToDo: Make saveData an overloaded function
  // saveData(solsNR.transpose(), "NRTracker_SEI_CDPR.txt");
  // std::cout << "Simulation done" << '\n';

  return 0;
}

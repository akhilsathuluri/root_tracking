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

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++ trackAllBranches check ++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  MatrixXd allRoots(initsols.rows(), initsols.cols());
  allRoots = rt.trackAllBranches(legvals.row(2), initsols, etaext, Jetaextphi);

std::cout << allRoots<< '\n';

  return 0;
}

#include<iostream>
#include<eigen3/Eigen/Dense>
#include "root_tracker.hh"
#include "inputs.hh"
#include <fstream>
#include<functional>
#include<chrono>

using namespace Eigen;

/*!
The saveData function takes in an input `MatrixXd` and saves it in a .txt file

@param A Input matrix
@param file_name Name of the file to save the input matrix in
*/
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
int main(int argc, char const *argv[]) {
  RootTracker rt;
  rt.Methods();
  MatrixXcd Clegvals(51,6);
  Clegvals = legvals.cast<std::complex<double>>();
  // Initial solution
  VectorXd phii(18), tempNR(18), tempDM(18), tempNN(18);
  VectorXcd phii2(18);
  phii << 0.2,0,1.28,0.2,0,0,
   0.003056244705087313,-0.06688548576654774,0.45696277640159066,0.5469562635040984,
   -0.09469005512919006,0.0034901115067618653,0.3362762851850709,0.2868906143657115,
   -0.021027360477035077,0.024784412575612605,-0.3953833637520812,-0.37979521343323663;

   phii2 << std::complex<double>(1.8700438348728161,0),std::complex<double>(-0.5235870043873335,0),std::complex<double>(0.,-0.926898973258719),std::complex<double>(0.,-0.10569096861540776),
      std::complex<double>(0.,0.7154784037290103),0,std::complex<double>(1.5707963267948966,1.5262959577414468),
      std::complex<double>(1.5707963267948966,0.6913993396133373),std::complex<double>(1.5707963267948966,-0.1997919002830573),
      std::complex<double>(1.5707963267948966,-0.5230490884925065),std::complex<double>(1.5707963267948966,0.37718441615032156),
      std::complex<double>(1.5707963267948966,1.3985097660833001),std::complex<double>(0.5889695511742916,-1.7974114210661284e-16),
      std::complex<double>(0.6273584382908439,-9.456869853007833e-17),std::complex<double>(0.3845316700324741,8.803289786100487e-17),
      std::complex<double>(0.5910427650086764,-9.488696919310518e-17),std::complex<double>(0.1527945222318842,-1.6367096645383268e-17),
      std::complex<double>(-0.07498723590723597,2.5098272394581245e-17);

   MatrixXd initsols(8, 18);
   initsols << 0.21833808165665786,-0.7906213868033382,0.7240859710802774,-0.7154994439464201,-0.6452575011458942,0,
   -0.11224708177891103,0.7453128304748228,1.4299633878665543,0.8300450230008743,-0.2868399561208049,
   -0.24735531615022688,0.8761914454555555,1.3561766375736213,0.8054121712576314,0.7196393797208634,
   0.10661311496251118,-0.03391249606302923,0.20000000000000112,0,1.2799999999999998,0.19999999999999996,0,0,
   0.003056244705087313,-0.06688548576654774,0.45696277640159066,0.5469562635040984,-0.09469005512919006,
   0.0034901115067618653,0.3362762851850709,0.2868906143657115,-0.021027360477035077,0.024784412575612605,
   -0.3953833637520812,-0.37979521343323663,-0.5136951579345258,-0.5199956193725788,-0.7987229831575837,
   -0.3337590445754842,-0.9516994614991616,0,-1.885474722853147,-2.65366482193082,2.782283527303565,
   2.89251489256385,-2.2052227902999824,-1.742941092611359,0.6160601733585822,0.6974690501752805,
   0.4197629659432657,0.5377001757469011,0.07050254829194261,-0.10394753137509125,0.5244336967602403,
   0.26850949397918333,1.0307988965821262,0.6037494337274671,-0.18856114122522924,0,0.19200440495122106,
   0.08892583707114052,0.6828422266917203,1.0712524507882133,0.5587642531431387,0.32979660408779493,
   0.2758192040473777,0.26974569330836556,-0.13921668983029753,-0.3564442780964989,-1.0169800563674098,
   -0.6281098269592473,0.20000000000000112,0,-1.2799999999999998,-0.19999999999999996,0,0,3.138536408884706,
   -3.0747071678232456,2.6846298771882027,2.594636390085695,-3.046902598460603,3.1381025420830313,
   0.3362762851850709,0.2868906143657115,-0.021027360477035077,0.024784412575612605,-0.3953833637520812,
   -0.37979521343323663,-0.5136951579345258,-0.5199956193725788,0.7987229831575837,0.3337590445754842,
   0.9516994614991616,0,-1.2561179307366463,-0.48792783165897324,0.3593091262862285,0.24907776102594312,
   -0.936369863289811,-1.398651560978434,0.6160601733585822,0.6974690501752805,0.4197629659432657,
   0.5377001757469011,0.07050254829194261,-0.10394753137509125,0.21833808165665786,-0.7906213868033382,
   -0.7240859710802774,0.7154994439464201,0.6452575011458942,0,-3.0293455718108824,2.3962798231149707,
   1.7116292657232388,2.311547630588919,-2.8547526974689883,-2.8942373374395665,0.8761914454555555,
   1.3561766375736213,0.8054121712576314,0.7196393797208634,0.10661311496251118,-0.03391249606302923,
   0.5244336967602403,0.26850949397918333,-1.0307988965821262,-0.6037494337274671,0.18856114122522924,0,
   2.9495882486385723,3.0526668165186526,2.4587504268980727,2.07034020280158,2.5828284004466546,
   2.8117960495019982,0.2758192040473777,0.26974569330836556,-0.13921668983029753,-0.3564442780964989,
   -1.0169800563674098,-0.6281098269592473;

  VectorXcd phiic(18), tempNRC(18), tempDMC(18), tempNNC(18);
  phiic = phii.cast<std::complex<double>>();
  MatrixXcd solsNNC(18,1), solsDMC(18,1), solsNRC(18,1), initsolsc(8, 18), tempFKC(8, 18);
  initsolsc = initsols.cast<std::complex<double>>();

  // solsNRC << phiic;
  // tempNRC = phiic;

  solsNRC << phii2;
  tempNRC = phii2;

  // Tracking using NRTracker
  for (int i = 1; i < legvals.rows(); i++){
    tempNRC = rt.NRCTracker(Clegvals.row(i), tempNRC, etaext, Jetaextphi);
    std::cout << tempNRC << '\n';

    solsNRC.conservativeResize(solsNRC.rows(), solsNRC.cols()+1);
  	solsNRC.col(solsNRC.cols()-1) = tempNRC;
  }
  // // Display results
  // std::cout << solsNRC.rows() << " " << solsNRC.cols() << std::endl;
  // std::cout << solsNRC.col(50) << std::endl;

  // ToDo: Make saveData an overloaded function
  saveCData(solsNRC.transpose(), "NRCTracker_2.txt");

  // MatrixXd solsNR(18,1);
  // solsNR << phii;
  // tempNR = phii;
  //
  // // Tracking using NRTracker
  // for (int i = 1; i < legvals.rows(); i++){
  //   tempNR = rt.NRTracker(legvals.row(i), tempNR, etaext, Jetaextphi);
  //
  //   solsNR.conservativeResize(solsNR.rows(), solsNR.cols()+1);
  // 	solsNR.col(solsNR.cols()-1) = tempNR;
  // }
  // // // Display results
  // // std::cout << solsNR.rows() << " " << solsNR.cols() << std::endl;
  // // std::cout << solsNR.col(50) << std::endl;
  //
  // saveData(solsNR.transpose(), "NRTracker_new.txt");
  //

  // solsDMC << phiic;
  // tempDMC = phiic;

  solsDMC << phii2;
  tempDMC = phii2;

  // Tracking using DMTracker
  for (int i = 1; i < Clegvals.rows(); i++){
    // Without NR step correction
    // tempDMC = rt.DMCTracker(Clegvals.row(i-1), Clegvals.row(i), tempDMC, Jetaexttheta, Jetaextphi);

    // With NR step correction
    /*!
    Choosing small enough `eps` makes DMTracker work like an NRTracker
    */
    tempDMC = rt.DMCTracker(Clegvals.row(i-1), Clegvals.row(i), tempDMC, Jetaexttheta, Jetaextphi, 0.05, etaext);

    solsDMC.conservativeResize(solsDMC.rows(), solsDMC.cols()+1);
  	solsDMC.col(solsDMC.cols()-1) = tempDMC;
  }
  // // Display results
  // std::cout << solsDM.rows() << " " << solsDM.cols() << std::endl;
  // std::cout << solsDM.col(50) << std::endl;

  saveCData(solsDMC.transpose(), "DMCTracker_NRC_2.txt");

  // MatrixXd solsDM(18,1);
  // solsDM << phii;
  // tempDM = phii;
  //
  // // Tracking using DMTracker
  // for (int i = 1; i < legvals.rows(); i++){
  //   // Without NR step correction
  //   // tempDM = rt.DMTracker(legvals.row(i-1), legvals.row(i), tempDM, Jetaexttheta, Jetaextphi);
  //
  //   // With NR step correction
  //   /*!
  //   Choosing small enough `eps` makes DMTracker work like an NRTracker
  //   */
  //   tempDM = rt.DMTracker(legvals.row(i-1), legvals.row(i), tempDM, Jetaexttheta, Jetaextphi, 0.05, etaext);
  //
  //   // solsDM.conservativeResize(solsDM.rows(), solsDM.cols()+1);
  // 	// solsDM.col(solsDM.cols()-1) = tempDM;
  // }
  // // // Display results
  // // std::cout << solsDM.rows() << " " << solsDM.cols() << std::endl;
  // // std::cout << solsDM.col(50) << std::endl;
  //
  // saveData(solsDM.transpose(), "DMTracker_NR.txt");
  //
  // MatrixXd solsNN(18,1), tempFK(8, 18);
  // solsNN << phii;
  // tempNN = phii;
  // // Initialise all the known roots
  // tempFK = initsols;
  //
  // // Tracking using NNTracker
  // for (int i = 1; i < legvals.rows(); i++){
  //   for (int j = 0; j < tempFK.rows(); j++) {
  //     tempFK.row(j) = rt.NRTracker(legvals.row(i), tempFK.row(j), etaext, Jetaextphi);
  //   }
  //
  //   tempNN = rt.NNTracker(tempNN, tempFK, 6);
  //
  //   solsNN.conservativeResize(solsNN.rows(), solsNN.cols()+1);
  // 	solsNN.col(solsNN.cols()-1) = tempNN;
  // }
  // // // Display results
  // // std::cout << solsNN.rows() << " " << solsNN.cols() << std::endl;
  // // std::cout << solsNN.col(50) << std::endl;
  //
  // saveData(solsNN.transpose(), "NNTracker_NR.txt");
  //
  // // // Verifying the tracked roots with the NR roots
  // // std::cout << (solsNN - solsNR).maxCoeff() << std::endl;

  solsNNC << phiic;
  tempNNC = phiic;
  // Initialise all the known roots including the complex root
  tempFKC = initsolsc;

  // Tracking using NNCTracker
  for (int i = 1; i < Clegvals.rows(); i++){
    for (int j = 0; j < tempFKC.rows(); j++) {
      tempFKC.row(j) = rt.NRCTracker(Clegvals.row(i), tempFKC.row(j), etaext, Jetaextphi);
    }

    tempNNC = rt.NNCTracker(tempNNC, tempFKC, 6);

    solsNNC.conservativeResize(solsNNC.rows(), solsNNC.cols()+1);
  	solsNNC.col(solsNNC.cols()-1) = tempNNC;
  }
  // // Display results
  // std::cout << solsNN.rows() << " " << solsNN.cols() << std::endl;
  // std::cout << solsNN.col(50) << std::endl;

  saveCData(solsNNC.transpose(), "NNCTracker_NRC.txt");

  // Verifying the tracked roots with the NR roots
  // std::cout << (solsNNC - solsNRC).maxCoeff() << std::endl;



  return 0;
}
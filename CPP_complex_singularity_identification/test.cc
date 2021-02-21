#include<iostream>
#include<eigen3/Eigen/Dense>

using namespace Eigen;

int main(int argc, char const *argv[]) {
  const double phii_data[] = {
    0.2,0,1.28,0.2,0,0,
     0.003056244705087313,-0.06688548576654774,0.45696277640159066,0.5469562635040984,
     -0.09469005512919006,0.0034901115067618653,0.3362762851850709,0.2868906143657115,
     -0.021027360477035077,0.024784412575612605,-0.3953833637520812,-0.37979521343323663
  };


const std::complex<double> phii_data2[] = {
  std::complex<double>(1.8700438348728161,0),std::complex<double>(-0.5235870043873335,0),std::complex<double>(0.,-0.926898973258719),std::complex<double>(0.,-0.10569096861540776),
     std::complex<double>(0.,0.7154784037290103),0,std::complex<double>(1.5707963267948966,1.5262959577414468),
     std::complex<double>(1.5707963267948966,0.6913993396133373),std::complex<double>(1.5707963267948966,-0.1997919002830573),
     std::complex<double>(1.5707963267948966,-0.5230490884925065),std::complex<double>(1.5707963267948966,0.37718441615032156),
     std::complex<double>(1.5707963267948966,1.3985097660833001),std::complex<double>(0.5889695511742916,-1.7974114210661284e-16),
     std::complex<double>(0.6273584382908439,-9.456869853007833e-17),std::complex<double>(0.3845316700324741,8.803289786100487e-17),
     std::complex<double>(0.5910427650086764,-9.488696919310518e-17),std::complex<double>(0.1527945222318842,-1.6367096645383268e-17),
     std::complex<double>(-0.07498723590723597,2.5098272394581245e-17)
};


  Matrix<std::complex<double>, 18, 1>phii(phii_data2);
  std::cout << phii << std::endl;

  return 0;
}
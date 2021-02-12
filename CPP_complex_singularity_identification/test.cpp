#include<iostream>
#include<complex>
#include<eigen3/Eigen/Dense>
#include<math.h>

using namespace Eigen;

// Can input reals and cast them to complex
// Can add a complex double to a matrix after we created it
// via append
// Data type of elements of complex type is std::complex<double>
// Complex matrices and vectors as MatrixXcd
// Trigonometric functions of complex numbers work
// Exponential functions also work with complex numbers
// Scalar multiplication also works

int main(int argc, char const *argv[]) {

  std::complex<double> x(1, 2), y(3, 4);
  VectorXcd test(2);
  test << x, y;

  std::cout << test.norm()  << '\n';

  return 0;
}

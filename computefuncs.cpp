//Dont know how to input the th and ph values

#include <iostream>
#include <math.h>

using namespace std;

double** computeJ(double* th, double* ph);
double* computef(double* th, double* ph);

int main(){
  double th[] = {0.1, 1.12, 2.31};
  double ph[] = {0.8898, 2.567, -0.12};

  double** Jmat;
  Jmat = computeJ(th, ph);

  cout<<Jmat[0][0]<<" "<<Jmat[0][1]<<" "<<Jmat[0][2]<<"\n"
  <<Jmat[1][0]<<" "<<Jmat[1][1]<<" "<<Jmat[1][2]<<"\n"
  <<Jmat[2][0]<<" "<<Jmat[2][1]<<" "<<Jmat[2][2]<<"\n"<<endl;

  double* fvec;
  fvec = computef(th, ph);
  cout<<fvec[0]<<" "<<fvec[1]<<" "<<fvec[2]<<endl;

return 0;
}

double** computeJ(double* th, double* ph){
double** tempJ = new double*[3];
for(int i = 0; i < 3; i++) {
  tempJ[i] = new double[3];
}
tempJ[0][0]= (2*sin(ph[0]))/3. + (2*cos(ph[1])*sin(ph[0]))/9. - (cos(th[0])*sin(ph[0]))/3. + (cos(th[1])*sin(ph[0]))/3. - (2*cos(ph[0])*sin(ph[1]))/9. + (cos(ph[0])*sin(th[0]))/3. - (cos(ph[0])*sin(th[1]))/3.;
tempJ[0][1] = (-2*cos(ph[1])*sin(ph[0]))/9. - (2*sin(ph[1]))/3. + (2*cos(ph[0])*sin(ph[1]))/9. + (cos(th[0])*sin(ph[1]))/3. - (cos(th[1])*sin(ph[1]))/3. - (cos(ph[1])*sin(th[0]))/3. + (cos(ph[1])*sin(th[1]))/3.;
tempJ[0][2] = 0;
tempJ[1][0] = -(cos(ph[0])*pow(3,-0.5)) + sin(ph[0])/3. + (2*cos(ph[2])*sin(ph[0]))/9. - (cos(th[0])*sin(ph[0]))/3. + (cos(th[2])*sin(ph[0]))/3. - (2*cos(ph[0])*sin(ph[2]))/9. + (cos(ph[0])*sin(th[0]))/3. - (cos(ph[0])*sin(th[2]))/3.;
tempJ[1][1] = 0;
tempJ[1][2] = cos(ph[2])*pow(3,-0.5) - (2*cos(ph[2])*sin(ph[0]))/9. - sin(ph[2])/3. + (2*cos(ph[0])*sin(ph[2]))/9. + (cos(th[0])*sin(ph[2]))/3. - (cos(th[2])*sin(ph[2]))/3. - (cos(ph[2])*sin(th[0]))/3. + (cos(ph[2])*sin(th[2]))/3.;
tempJ[2][0] = 0;
tempJ[2][1] = -(cos(ph[1])*pow(3,-0.5)) - sin(ph[1])/3. + (2*cos(ph[2])*sin(ph[1]))/9. - (cos(th[1])*sin(ph[1]))/3. + (cos(th[2])*sin(ph[1]))/3. - (2*cos(ph[1])*sin(ph[2]))/9. + (cos(ph[1])*sin(th[1]))/3. - (cos(ph[1])*sin(th[2]))/3.;
tempJ[2][2] = cos(ph[2])*pow(3,-0.5) - (2*cos(ph[2])*sin(ph[1]))/9. + sin(ph[2])/3. + (2*cos(ph[1])*sin(ph[2]))/9. + (cos(th[1])*sin(ph[2]))/3. - (cos(th[2])*sin(ph[2]))/3. - (cos(ph[2])*sin(th[1]))/3. + (cos(ph[2])*sin(th[2]))/3.;

return tempJ;
}

double* computef(double* th, double* ph){
float l=0.5, r= 0.33336, a= 0.25, b= 1.00;
double* tempf = new double[3];
tempf[0] = -2*b*r*cos(ph[0]) + 2*b*r*cos(ph[1]) - 2*b*l*cos(th[0]) + 2*l*r*cos(ph[0])*cos(th[0]) - 2*l*r*cos(ph[1])*cos(th[0]) + 2*b*l*cos(th[1]) - 2*l*r*cos(ph[0])*cos(th[1]) + 2*l*r*cos(ph[1])*cos(th[1]) - pow(a,2) + pow(b,2) + 2*pow(l,2) - 2*cos(th[0])*cos(th[1])*pow(l,2) + 2*pow(r,2) - 2*cos(ph[0])*cos(ph[1])*pow(r,2) - 2*pow(r,2)*sin(ph[0])*sin(ph[1]) + 2*l*r*sin(ph[0])*sin(th[0]) - 2*l*r*sin(ph[1])*sin(th[0]) - 2*l*r*sin(ph[0])*sin(th[1]) + 2*l*r*sin(ph[1])*sin(th[1]) - 2*pow(l,2)*sin(th[0])*sin(th[1]);
tempf[1] = -(b*r*cos(ph[0])) + b*r*cos(ph[2]) - b*l*cos(th[0]) + 2*l*r*cos(ph[0])*cos(th[0]) - 2*l*r*cos(ph[2])*cos(th[0]) + b*l*cos(th[2]) - 2*l*r*cos(ph[0])*cos(th[2]) + 2*l*r*cos(ph[2])*cos(th[2]) - pow(a,2) + pow(b,2) + 2*pow(l,2) - 2*cos(th[0])*cos(th[2])*pow(l,2) + 2*pow(r,2) - 2*cos(ph[0])*cos(ph[2])*pow(r,2) - b*r*pow(3,0.5)*sin(ph[0]) + b*r*pow(3,0.5)*sin(ph[2]) - 2*pow(r,2)*sin(ph[0])*sin(ph[2]) - b*l*pow(3,0.5)*sin(th[0]) + 2*l*r*sin(ph[0])*sin(th[0]) - 2*l*r*sin(ph[2])*sin(th[0]) + b*l*pow(3,0.5)*sin(th[2]) - 2*l*r*sin(ph[0])*sin(th[2]) + 2*l*r*sin(ph[2])*sin(th[2]) - 2*pow(l,2)*sin(th[0])*sin(th[2]);
tempf[2] = b*r*cos(ph[1]) - b*r*cos(ph[2]) + b*l*cos(th[1]) + 2*l*r*cos(ph[1])*cos(th[1]) - 2*l*r*cos(ph[2])*cos(th[1]) - b*l*cos(th[2]) - 2*l*r*cos(ph[1])*cos(th[2]) + 2*l*r*cos(ph[2])*cos(th[2]) - pow(a,2) + pow(b,2) + 2*pow(l,2) - 2*cos(th[1])*cos(th[2])*pow(l,2) + 2*pow(r,2) - 2*cos(ph[1])*cos(ph[2])*pow(r,2) - b*r*pow(3,0.5)*sin(ph[1]) + b*r*pow(3,0.5)*sin(ph[2]) - 2*pow(r,2)*sin(ph[1])*sin(ph[2]) - b*l*pow(3,0.5)*sin(th[1]) + 2*l*r*sin(ph[1])*sin(th[1]) - 2*l*r*sin(ph[2])*sin(th[1]) + b*l*pow(3,0.5)*sin(th[2]) - 2*l*r*sin(ph[1])*sin(th[2]) + 2*l*r*sin(ph[2])*sin(th[2]) - 2*pow(l,2)*sin(th[1])*sin(th[2]);

return tempf;
}

double** retthvec(){
  double** thvec;
  thvec = new double*[1001];
  for(int i = 0; i < 3; i++) {
    tempJ[i] = new double[3];
}

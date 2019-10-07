#include <iostream>
#include <math.h>
#include <gsl/gsl_statistics.h>
// #include <stdlib.h>

using namespace std;

double* in2pi(double* th, int size);
double* polynn(double* curr, int sizec, double* prev, int sizep);

int main(){
  double prev[] = {0.1353*M_PI, -2.14124*M_PI, 5.135*M_PI};
  double curr[] = {-2.97252*M_PI, -1.435*M_PI, -0.123*M_PI, 0.2536*M_PI, 5.924*M_PI};

  // cout<<"in2pis are: "<<in2pi(prev, 3)[0]<<endl;
  // cout<<"in2pis are: "<<in2pi(prev, 3)[1]<<endl;
  // cout<<"in2pis are: "<<in2pi(prev, 3)[2]<<endl;
  // cout<<"\n"<<endl;
  // cout<<"in2pis are: "<<in2pi(curr, 5)[0]<<endl;
  // cout<<"in2pis are: "<<in2pi(curr, 5)[1]<<endl;
  // cout<<"in2pis are: "<<in2pi(curr, 5)[2]<<endl;
  // cout<<"in2pis are: "<<in2pi(curr, 5)[3]<<endl;
  // cout<<"in2pis are: "<<in2pi(curr, 5)[4]<<endl;

  double* newsol;
  newsol = polynn(curr, 5, prev, 3);
  cout<<newsol[0]<<endl;
  cout<<newsol[1]<<endl;
  cout<<newsol[2]<<endl;
return 0;
}

double* in2pi(double* th, int size){
  double ref  = 2*M_PI;
  for(int i=0;i<size;i++)
    th[i] = ((th[i]/ref) - floor(th[i]/ref))*ref;
  return th;
}

double* polynn(double* curr, int sizec, double* prev, int sizep){
  curr = in2pi(curr, sizec);
  prev = in2pi(prev, sizep);

  double* temp = new double[sizep];
  int k = 0, loop = 0;
  double tempprev = 15.00;
  for(int i=0;i<sizep;i++){
    for(int j=0;j<sizec;j++){
      if(fabs(prev[i]-curr[j])<tempprev){
        k = j;
      };
        tempprev=fabs(prev[i]-curr[k]);
        // cout<<"tempprev is "<<i<<" "<<j<<" "<<tempprev<<endl;
    };
    // cout<<"chosen k is "<<k<<endl;
    temp[loop] = curr[k];
    loop++;
    tempprev = 15.00;
  };
  return temp;
}

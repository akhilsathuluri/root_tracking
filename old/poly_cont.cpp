// Take one known set of roots as input and check if the values are correct

#include <iostream>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_linalg.h>
#include <math.h>

using namespace std;


double* polycont(int sof_rows_J, int sof_cols_J, double* th, double* ph);
double* absarray(double* data, int size);
double* LinearSolve(double* A_data, int sof_A_rows, int sof_A_cols, double* b_data);
double* computeJ(double* th, double* ph);
double* computef(double* th, double* ph);


int main(){
double* phval;
double th[3] = {0.1353, -2.14124, 5.135};
double ph[3] = {-0.134, 0.2452, 1.442652};
 phval = ph;
 phval = polycont(3, 3, th, ph);
 cout<<phval[0]<<" "<<phval[1]<<" "<<phval[2]<<" "<<endl;
return 0;
}

double* polycont(int sof_rows_J, int sof_cols_J, double* th, double* ph){
  double* J = new double[sof_cols_J*sof_rows_J];
  double* f = new double[sof_cols_J];
  do{
    J = computeJ(th, ph);
    f = computef(th, ph);
    for(int a=0;a<sof_cols_J;a++)
      ph[a]=ph[a]-(LinearSolve(J, sof_rows_J, sof_cols_J, f))[a];
  }while(gsl_stats_max(absarray(f, sof_cols_J), 1, sof_cols_J)<gsl_sf_pow_int(10.0, -12));
  return ph;
}

double* absarray(double* data, int size){
for(int i=0;i<size;i = i+1){
if(data[i]<0){
data[i]=-data[i];
}
}
return data;
}

double* LinearSolve(double *A_data, int sof_A_rows, int sof_A_cols, double *b_data){

  gsl_matrix_view m
    = gsl_matrix_view_array (A_data, sof_A_rows, sof_A_cols);

  gsl_vector_view b
    = gsl_vector_view_array (b_data, sof_A_rows);

  gsl_vector *x = gsl_vector_alloc (sof_A_cols);

  int s;

  gsl_permutation * p = gsl_permutation_alloc (sof_A_cols);

  gsl_linalg_LU_decomp (&m.matrix, p, &s);

  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

  return x->data;
}

double* computeJ(double* th, double* ph){
double* tempJ = new double[9];
tempJ[0]= (2*sin(ph[0]))/3. + (2*cos(ph[1])*sin(ph[0]))/9. - (cos(th[0])*sin(ph[0]))/3. + (cos(th[1])*sin(ph[0]))/3. - (2*cos(ph[0])*sin(ph[1]))/9. + (cos(ph[0])*sin(th[0]))/3. - (cos(ph[0])*sin(th[1]))/3.;
tempJ[1] = (-2*cos(ph[1])*sin(ph[0]))/9. - (2*sin(ph[1]))/3. + (2*cos(ph[0])*sin(ph[1]))/9. + (cos(th[0])*sin(ph[1]))/3. - (cos(th[1])*sin(ph[1]))/3. - (cos(ph[1])*sin(th[0]))/3. + (cos(ph[1])*sin(th[1]))/3.;
tempJ[2] = 0;
tempJ[3] = -(cos(ph[0])*pow(3,-0.5)) + sin(ph[0])/3. + (2*cos(ph[2])*sin(ph[0]))/9. - (cos(th[0])*sin(ph[0]))/3. + (cos(th[2])*sin(ph[0]))/3. - (2*cos(ph[0])*sin(ph[2]))/9. + (cos(ph[0])*sin(th[0]))/3. - (cos(ph[0])*sin(th[2]))/3.;
tempJ[4] = 0;
tempJ[5] = cos(ph[2])*pow(3,-0.5) - (2*cos(ph[2])*sin(ph[0]))/9. - sin(ph[2])/3. + (2*cos(ph[0])*sin(ph[2]))/9. + (cos(th[0])*sin(ph[2]))/3. - (cos(th[2])*sin(ph[2]))/3. - (cos(ph[2])*sin(th[0]))/3. + (cos(ph[2])*sin(th[2]))/3.;
tempJ[6] = 0;
tempJ[7] = -(cos(ph[1])*pow(3,-0.5)) - sin(ph[1])/3. + (2*cos(ph[2])*sin(ph[1]))/9. - (cos(th[1])*sin(ph[1]))/3. + (cos(th[2])*sin(ph[1]))/3. - (2*cos(ph[1])*sin(ph[2]))/9. + (cos(ph[1])*sin(th[1]))/3. - (cos(ph[1])*sin(th[2]))/3.;
tempJ[8] = cos(ph[2])*pow(3,-0.5) - (2*cos(ph[2])*sin(ph[1]))/9. + sin(ph[2])/3. + (2*cos(ph[1])*sin(ph[2]))/9. + (cos(th[1])*sin(ph[2]))/3. - (cos(th[2])*sin(ph[2]))/3. - (cos(ph[2])*sin(th[1]))/3. + (cos(ph[2])*sin(th[2]))/3.;

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

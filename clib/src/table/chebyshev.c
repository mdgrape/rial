#include <stdlib.h>
#include <math.h>

//
// Fit function by Chebyshev
// size of c: norder+1 
//
void mdg_cheb_fit(double (*f)(double), double xmin, double xmax, int norder,
	      double *c) {

  for(int j=0;j<=norder;++j) {
    c[j] = 0.0;
    for(int k=1;k<=norder+1;++k) {
      double z=M_PI*(k-0.5)/(norder+1);
      double x=xmin+(cos(z)+1.0)*(xmax-xmin)*0.5;
      c[j] += (*f)(x) * cos(j*z);
    }
    c[j] *= 2.0/(norder+1);
  }
  c[0] *= 0.5;
}

//
// Fit function by Chebyshev
// Divide xmin-xmax to ndiv region
// size of c: ndiv * (norder+1)
//
double* mdg_cheb_fit_div(double (*f)(double), double xmin, double xmax, 
			 int ndiv, int norder) {
  
  double w=(xmax-xmin)/ndiv;
  double *c = malloc(sizeof(double)*ndiv*(norder+1));
  for(int j=0;j<ndiv;++j) {
    double min = xmin +w*j;
    double max = min + w;
    mdg_cheb_fit(f,min,max,norder,c+j*(norder+1));
  }
  return c;
}

//
// Convert Chebyshev coefficients to polynomial coefficients
// T0 = 1, T1=x, Tn+1+Tn-1 = 2xTn
// by Clenshaw's recurrence
//    http://en.wikipedia.org/wiki/Clenshaw_algorithm
//
void mdg_cheb_conv(int norder, double c[], double c_p[]) {

  double *b=malloc(sizeof(double)*(norder+1));

  for(int j=0;j<=norder;++j) c_p[j]=b[j]=0.0;

  // b_N = c_n
  c_p[0] = c[norder];
  // c_p : coefficients 
  // b   : coefficients to keep previous step

  // To evaluate b_j
  // b_j = 2 x b_j+1 - b_j+2 + c_j
  for(int j=norder-1;j>=1;--j) {
    for(int k=norder-j;k>=1;--k) {

      double tmp = c_p[k];
      c_p[k] = 2.0 * c_p[k-1] - b[k];
      b[k] = tmp;
    }
    double tmp = c_p[0];
    c_p[0] = -b[0] + c[j];
    b[0] = tmp;
  }
  // Final step differs
  // b_0 = x b_1 - b_2 + c_0
  for(int k=norder;k>=1;--k) {
    c_p[k] = c_p[k-1] - b[k];
  }
  c_p[0] = - b[0] + c[0];

  free(b);

}


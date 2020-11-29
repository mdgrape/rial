#ifndef __CHEBYSHEV_H
#define __CHEBYSHEV_H

#ifdef __cplusplus
extern "C" {
#endif

/* chebyshev.c */
void mdg_cheb_fit(double (*f)(double), double xmin, double xmax, int norder, double *c);
double *mdg_cheb_fit_div(double (*f)(double), double xmin, double xmax, int ndiv, int norder);
void mdg_cheb_conv(int norder, double c[], double c_p[]);

#ifdef __cplusplus
}
#endif

#endif

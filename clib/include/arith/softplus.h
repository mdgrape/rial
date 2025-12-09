#ifndef RIAL_CLIB_SOFTPLUS_H
#define RIAL_CLIB_SOFTPLUS_H

#ifdef __cplusplus
extern "C" {
#endif

  double rial_softplus(double x);
  double rial_softplus_deriv(double x);
  double rial_softplus_parabolic(double x, double a);
  double rial_softplus_parabolic_deriv(double x, double a);
  double rial_softplus_quadratic(double x, double a);
  double rial_softplus_quadratic_deriv(double x, double a);

#ifdef __cplusplus
}
#endif

#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#include "softplus.h"

/** @file softplus.c */

/**
   @brief Softplus in definition
   @param[in] x
   @retval log(1+exp(x))
*/
double
rial_softplus(double x)
{
  return log(1.0+exp(x));
}

/**
   @brief Derivative of softplus in definition
   @param[in] x
   @retval 1.0/(1+exp(-x))
*/
double
rial_softplus_deriv(double x)
{
  return 1.0/(1.0+exp(-x));
}

/**
   @brief Softplus approximated by 2nd order 
   @param[in] x : input
   @param[in] a : parabolic in [-a, a]
   @retval approximated softplus
*/
double
rial_softplus_parabolic(double x, double a)
{
  if (x<=-a) return 0.0;
  if (x>=a) return x;
  double z = x/a;
  return 0.25*(z+1.0)*(z+1.0)*a;
}

/**
   @brief Derivative of softplus approximated by 2nd order 
   @param[in] x : input
   @param[in] a : parabolic in [-a, a]
   @retval approximated derivative of softplus
*/
double
rial_softplus_parabolic_deriv(double x, double a)
{
  if (x<=-a) return 0.0;
  if (x>=a) return 1.0;
  double z = x/a;
  return 0.5*(z+1.0);
}

/**
   @brief Softplus approximated by 4th order 
   @param[in] x : input
   @param[in] a : quadratic in [-a, a]
   @retval approximated softplus
*/
double
rial_softplus_quadratic(double x, double a)
{
  if (x<=-a) return 0.0;
  if (x>=a) return x;
  double z = x/a;
  return (-0.0625*z*z*z*z + 0.375*z*z + 0.5*z + 0.1875)*a;
}

/**
   @brief Derivative of softplus approximated by 4th order 
   @param[in] x : input
   @param[in] a : quadratic in [-a, a]
   @retval approximated derivative of softplus
*/
// This takes min=0 at x=-a, max=1 at x=a
// integral of -b(x^2-a^2) = -b/3 x^3 + b a^2 x + c
// b a^3 /3 - b a^3 + c = 0, c = 2/3 b a^3
// 4/3 b a^3 = 1, b=3/(4 a^3), c=0.5
double
rial_softplus_quadratic_deriv(double x, double a)
{
  if (x<=-a) return 0.0;
  if (x>=a) return 1.0;
  double z = x/a;
  return -0.25*z*z*z + 0.75*z + 0.5;
}


//
//  nonbond_functions.c
//  
//  Library of nonbond functions 
//  used in MDG4A chip
//
//  Copyright (c) 2016 by Makoto Taiji/RIKEN QBiC
//  taiji@riken.jp
//
//  Force returns absolute magnitude of force.
//  Function evaluator should return force value/r.
//
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "nonbond_functions.h"

////////////////////////////////////////////////////////////////////////
// Ewald normal


double mdg_ewald_potential_core(double r, double alpha) {
  return erfc(alpha*r);
}

double mdg_ewald_potential(double r, double alpha) {
  return erfc(alpha*r)/r;
}

double mdg_ewald_force_core(double r, double alpha) {
  double x=alpha*r;
  return erfc(x) + 2.0*x*exp(-x*x)/sqrt(M_PI);
}

double mdg_ewald_force(double r, double alpha) {
  return mdg_ewald_force_core(r,alpha)/(r*r);
}

double mdg_ewald_potential_excl(double r, double alpha) {
  return -erf(alpha*r)/r;
}

double mdg_ewald_force_excl(double r, double alpha) {
  double x=alpha*r;
  return -erf(x)/(r*r) + 2.0*alpha*exp(-x*x)/sqrt(M_PI)/r;
}

double mdg_ewald_potential_excl14(double r, double alpha) {
  return mdg_ewald_potential(r, alpha) - 0.5/r;
}

double mdg_ewald_force_excl14(double r, double alpha) {
  return mdg_ewald_force(r, alpha) - 0.5/(r*r);
}

////////////////////////////////////////////////////////////////////////
// Zero dipole

static double zd_u(double r, double rc, double alpha) {
  return mdg_ewald_potential(r, alpha)+0.5*mdg_ewald_force(rc, alpha)*r*r;
}

double mdg_zd_potential(double r, double rc, double alpha) {
  if (r>=rc) return 0.0;
  return zd_u(r, rc, alpha) - zd_u(rc, rc, alpha);
}

double mdg_zd_force(double r, double rc, double alpha) {
  if (r>=rc) return 0.0;
  return mdg_ewald_force(r,alpha) - mdg_ewald_force(rc,alpha);
}

////////////////////////////////////////////////////////////////////////
// Soft core (linear exploration from 1/r)

double mdg_coulomb_softcore_potential(double r, double r0) {
  if (r>=r0) return 1.0/r;
  double grad = 1.0/(r0*r0);
  double dfdr = -2.0*pow(r0,-3.0);
  return 1.0/r0 - grad*(r-r0) - 0.5*dfdr*(r-r0)*(r-r0);
}

double mdg_coulomb_softcore_force(double r, double r0) {
  if (r>=r0) return pow(r, -2.0);
  double dfdr = -2.0*pow(r0,-3.0);
  return 1.0/(r0*r0) + dfdr * (r-r0);
}

////////////////////////////////////////////////////////////////////////
// Soft core (linear exploration from Ewald)

double mdg_ewald_softcore_potential(double r, double r0, double alpha) {
  double y = mdg_ewald_potential(r,alpha);
  if (r>=r0) return y;
  y -= 1/r;
  y += mdg_coulomb_softcore_potential(r, r0);
  return y;
}

double mdg_ewald_softcore_force(double r, double r0, double alpha) {
  double y = mdg_ewald_force(r, alpha);
  if (r>=r0) return y;
  y -= pow(r,-2.0);
  y += mdg_coulomb_softcore_force(r, r0);
  return y;
}

double mdg_ewald_softcore_potential_excl(double r, double r0, double alpha) {
  return mdg_ewald_potential_excl(r, alpha);
}

double mdg_ewald_softcore_force_excl(double r, double r0, double alpha) {
  return mdg_ewald_force_excl(r, alpha);
}

double mdg_ewald_softcore_potential_excl14(double r, double r0, double alpha) {
  if (r>=r0) return mdg_ewald_potential_excl14(r, alpha);
  double y = mdg_ewald_potential_excl(r, alpha);
  y += 0.5*mdg_coulomb_softcore_potential(r, r0);
  return y;
}

double mdg_ewald_softcore_force_excl14(double r, double r0, double alpha) {
  if (r>=r0) return mdg_ewald_force_excl14(r, alpha);
  double y = mdg_ewald_force_excl(r, alpha);
  y += 0.5*mdg_coulomb_softcore_force(r, r0);
  return y;
}

////////////////////////////////////////////////////////////////////////
// Soft core (potential linear exploration from 1/r)

double mdg_coulomb_cf_softcore_potential(double r, double r0) {
  if (r>=r0) return 1.0/r;
  double grad = 1.0/(r0*r0);
  return 1.0/r0 - grad*(r-r0);
}

double mdg_coulomb_cf_softcore_force(double r, double r0) {
  if (r>=r0) return 1.0/(r*r);
  return 1.0/(r0*r0);
}

////////////////////////////////////////////////////////////////////////
// Soft core (force propotional to r below r0)

double mdg_lf_softcore_potential(double r, double r0, double f0r0, double phi0) {
  return phi0 + 0.5*f0r0*(r0*r0 - r*r);
}

double mdg_lf_softcore_force(double r, double r0, double f0r0) {
  return f0r0*r;
}

double mdg_coulomb_lf_softcore_potential(double r, double r0) {
  if (r>=r0) return 1.0/r;
  double pot = 1.0/r0;
  double grad = pot*pot*pot;
  return mdg_lf_softcore_potential(r, r0, grad, pot);
}

double mdg_coulomb_lf_softcore_force(double r, double r0) {
  if (r>=r0) return 1.0/(r*r);
  double grad = 1.0/(r0*r0*r0);
  return mdg_lf_softcore_force(r, r0, grad);
}

////////////////////////////////////////////////////////////////////////
// Soft core (Gromacs)

// Actually, rA = (alpha sigma^6 lambda^p + r^6)^(1/6)
// Here, alpha=lambda=1, so give sigma to include alpha/lambda
// Gromacs suggests alpha=0.7 for p=1 or alpha=0.5 for p=2.
// In the latter case, alpha * (lambda)^p = 1/8 for lambda=0.5.
// Thus, set sigma <- sigma / sqrt(2) in the case.

static double gromacs_softcore_r(double r, double sigma) {
  return pow(pow(sigma,6.0)+pow(r,6.0), 1.0/6.0);
}

double mdg_coulomb_gromacs_softcore_potential(double r, double sigma) {
  double ra = gromacs_softcore_r(r, sigma);
  return 1.0/ra;
}

double mdg_coulomb_gromacs_softcore_force(double r, double sigma) {
  double ra = gromacs_softcore_r(r, sigma);
  return 1.0/(ra*ra) * pow(r/ra, 5.0);
}


////////////////////////////////////////////////////////////////////////
// Cutoff (Gromacs)

////////////////////////////////////////////////////////////////////////
// van der Waals

double mdg_lj_potential(double r) {
  return pow(r, -12.0) - pow(r, -6.0);
}

double mdg_lj_force(double r) {
  return 12.0*pow(r, -13.0) - 6.0*pow(r, -7.0);
}

#ifdef COMPILE_MAIN

int main(int argc, char *argv[]) {
  char sep=' ';
  
  printf("\"r\"%c",sep); // 1
  printf("\"x(=r^2)\"%c",sep); // 2
  printf("\"ewald potential\"%c",sep); // 3
  printf("\"ewald force\"%c",sep);     // 4
  printf("\"ewald potential excl\"%c",sep); // 5
  printf("\"ewald force excl\"%c",sep);     // 6
  printf("\"ewald potential excl14\"%c",sep); // 7
  printf("\"ewald force excl14\"%c",sep);     // 8
  printf("\"zd potential\"%c",sep); // 9
  printf("\"zd force\"%c",sep);     // 10
  printf("\"coulomb softcore potential\"%c",sep); // 11
  printf("\"coulomb softcore force\"%c",sep);     // 12
  printf("\"ewald softcore potential\"%c",sep); // 13
  printf("\"ewald softcore force\"%c",sep);     // 14
  printf("\"ewald softcore potential excl\"%c",sep); // 15
  printf("\"ewald softcore force excl\"%c",sep);     // 16
  printf("\"ewald softcore potential excl14\"%c",sep); // 17
  printf("\"ewald softcore force excl14\"%c",sep);     // 18
  printf("\"coulomb cf softcore potential\"%c",sep); // 19
  printf("\"coulomb cf softcore force\"%c",sep);     // 20
  printf("\"coulomb lf softcore potential\"%c",sep); // 21
  printf("\"coulomb lf softcore force\"%c",sep);     // 22
  printf("\"coulomb gromacs softcore potential\"%c",sep); // 23
  printf("\"coulomb gromacs softcore force\"%c",sep);     // 24
  printf("\"lj potential\"%c",sep); // 25
  printf("\"lj force\"%c",sep);     // 26
  printf("\n");
  double rmax  = 2.0; // nm
  double alpha = 3.5; // 1/nm
  double zd_rc = 1.2; // nm
  double softcore_rq = 0.3; // nm
  double vdw_sigma = 0.2; // nm

  int n=1000;
  for(int i=1;i<=n;++i) {
    double r = rmax*i/n;
    double x = alpha*alpha*r*r;
    printf("%le%c%le%c", r,sep, x,sep);
    printf("%le%c", mdg_ewald_potential(r, alpha),sep);
    printf("%le%c", mdg_ewald_force(r, alpha),sep);
    printf("%le%c", mdg_ewald_potential_excl(r, alpha),sep);
    printf("%le%c", mdg_ewald_force_excl(r, alpha),sep);
    printf("%le%c", mdg_ewald_potential_excl14(r, alpha),sep);
    printf("%le%c", mdg_ewald_force_excl14(r, alpha),sep);
    printf("%le%c", mdg_zd_potential(r, zd_rc, alpha),sep);
    printf("%le%c", mdg_zd_force(r, zd_rc, alpha),sep);
    printf("%le%c", mdg_coulomb_softcore_potential(r, softcore_rq),sep);
    printf("%le%c", mdg_coulomb_softcore_force(r, softcore_rq),sep);
    printf("%le%c", mdg_ewald_softcore_potential(r, softcore_rq, alpha),sep);
    printf("%le%c", mdg_ewald_softcore_force(r, softcore_rq, alpha),sep);
    printf("%le%c", mdg_ewald_softcore_potential_excl(r, softcore_rq, alpha),sep);
    printf("%le%c", mdg_ewald_softcore_force_excl(r, softcore_rq, alpha),sep);
    printf("%le%c", mdg_ewald_softcore_potential_excl14(r, softcore_rq, alpha),sep);
    printf("%le%c", mdg_ewald_softcore_force_excl14(r, softcore_rq, alpha),sep);
    printf("%le%c", mdg_coulomb_cf_softcore_potential(r, softcore_rq),sep);
    printf("%le%c", mdg_coulomb_cf_softcore_force(r, softcore_rq),sep);
    printf("%le%c", mdg_coulomb_lf_softcore_potential(r, softcore_rq),sep);
    printf("%le%c", mdg_coulomb_lf_softcore_force(r, softcore_rq),sep);
    printf("%le%c", mdg_coulomb_gromacs_softcore_potential(r, softcore_rq/sqrt(2.0)),sep);
    printf("%le%c", mdg_coulomb_gromacs_softcore_force(r, softcore_rq/sqrt(2.0)),sep);
    printf("%le%c", mdg_lj_potential(r/vdw_sigma),sep);
    printf("%le%c", mdg_lj_force(r/vdw_sigma),sep);
    printf("\n");
  }    
}

#endif

#ifndef __NONBOND_FUNCTIONS_H
#define __NONBOND_FUNCTIONS_H

#ifdef __cplusplus
extern "C" {
#endif

  double mdg_ewald_potential_core(double r, double alpha);
  double mdg_ewald_potential(double r, double alpha);
  double mdg_ewald_force_core(double x, double alpha);
  double mdg_ewald_force(double r, double alpha);
  double mdg_ewald_potential_excl(double r, double alpha);
  double mdg_ewald_force_excl(double r, double alpha);
  double mdg_ewald_potential_excl14(double r, double alpha);
  double mdg_ewald_force_excl14(double r, double alpha);

  double mdg_zd_potential(double r, double rc, double alpha);
  double mdg_zd_force(double r, double rc, double alpha);

  double mdg_coulomb_softcore_potential(double r, double r0);
  double mdg_coulomb_softcore_force(double r, double r0);

  double mdg_ewald_softcore_potential(double r, double r0, double alpha);
  double mdg_ewald_softcore_force(double r, double r0, double alpha);
  double mdg_ewald_softcore_potential_excl(double r, double r0, double alpha);
  double mdg_ewald_softcore_force_excl(double r, double r0, double alpha);
  double mdg_ewald_softcore_potential_excl14(double r, double r0, double alpha);
  double mdg_ewald_softcore_force_excl14(double r, double r0, double alpha);

  double mfg_coulomb_cf_softcore_potential(double r, double r0);
  double mdg_coulomb_cf_softcore_force(double r, double r0);

  double mdg_coulomb_lf_softcore_potential(double r, double r0);
  double mdg_coulomb_lf_softcore_force(double r, double r0);

  double mdg_lf_softcore_potential(double r, double r0, double f0r0, double phi0);
  double mdg_lf_softcore_force(double r, double r0, double f0);

  double mdg_coulomb_gromacs_softcore_potential(double r, double sigma);
  double mdg_coulomb_gromacs_softcore_force(double r, double sigma);
  
  double mdg_lj_potential(double x);
  double mdg_lj_force(double x);

#ifdef __cplusplus
}
#endif


#endif

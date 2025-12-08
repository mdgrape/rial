#ifndef __SCALE_H
#define __SCALE_H

#ifdef __cplusplus
extern "C" {
#endif

/* scale.c */
typedef struct {
  double length;
  double charge;
  double vdw_sigma_max;
  double vdw_epsilon_max;
  double vdw_energy_factor;
  double coulomb_energy_factor; 
  int kernel_bp_force;
  int kernel_bp_potential;
  int norm_coulomb_potential;
  int norm_coulomb_force;
  int norm_vdw_potential;
  int norm_vdw_force;
  int bp_force;
  int bp_grid_potential;
} mdg_scale_conversion_constants_t;

extern mdg_scale_conversion_constants_t *MDG_SC;

/* scale.c */
uint32_t mdg_scale_r(double r);
uint32_t mdg_scale_rcut2(double rc);
uint32_t mdg_scale_rough_rcut2(double rc);
int32_t  mdg_scale_charge(double q);

double mdg_rev_scale_r(uint32_t s);
double mdg_rev_scale_r_ltrmax(uint32_t s,double rmax);
double mdg_rev_scale_charge(int32_t q);

uint32_t mdg_scale_vdw_factor(double gbc);
uint32_t mdg_scale_coulomb_factor(double gbc);
double mdg_rev_scale_vdw_factor(uint32_t s);		//Added (Komatsu)
double mdg_rev_scale_coulomb_factor(uint32_t s);	//Added (Komatsu)

double mdg_rev_scale_vdw_radius_arithmetic(uint32_t s);
double mdg_rev_scale_vdw_radius_geometric(uint32_t s);
double mdg_rev_scale_vdw_energy(uint32_t eps);

uint32_t mdg_scale_vdw_radius_arithmetic(double s);
uint32_t mdg_scale_vdw_radius_geometric(double s);
uint32_t mdg_scale_vdw_energy(double eps);

double mdg_scale_force(int32_t q);
double mdg_scale_force64(int64_t q);
double mdg_scale_potential(int64_t q);

double mdg_rev_scale_coulomb_exclusion(int32_t s); 
int32_t mdg_scale_coulomb_exclusion(double x);

void mdg_scale_coulomb_softcore(double r0, double f0, double phi0,
				uint32_t *r0i, uint32_t *f0r0i, uint32_t *k0i,
				int *ex_norm_force, int *ex_norm_potential);

void mdg_scale_vdw_softcore(double r0, double f0, double phi0,
			    uint32_t *r0i, uint32_t *f0r0i, uint32_t *k0i,
			    int *norm_force, int *norm_potential);

mdg_scale_conversion_constants_t * mdg_set_gromacs_conversion();

#ifdef __cplusplus
}
#endif

#endif

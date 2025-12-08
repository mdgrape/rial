#ifndef __FPMD_H
#define __FPMD_H

#include "gentable.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MAXPIPE 100

  extern int B_MDG_CHECK;
  extern int MDG_DEBUG_FLAG;

#ifdef DEBUG
#define DEBUG_PRINTF(...) { if (MDG_DEBUG_FLAG) printf(__VA_ARGS__); }
#else
#define DEBUG_PRINTF(...) 
#endif

  int mdg_set_check_flag(int flag);
  int mdg_set_debug_flag(int flag);
  void mdg_get_environment();
  
  uint32_t mdg_inv_sqrt(uint32_t x);
  double mdg_rev_inv_sqrt(uint32_t x); //for debug (Komatsu)
  uint32_t mdg_inv_sqrt3(uint32_t x);
  double mdg_rev_inv_sqrt3(uint32_t x); //for debug (Komatsu)
  uint32_t mdg_r2sqrt(uint32_t x);

  double mdg_rev_funcev_vdw_common(uint64_t x); //for debug (Komatsu)
  int32_t mdg_funcev_coulomb(uint32_t x, mdg_func_table *t);
  int32_t mdg_funcev_coulomb_array(uint32_t x, int64_t c[512][4], char *name);

  void mdg_conv_r2_fromDouble(double x, uint32_t *z1, uint32_t *z2);

  void mdg_rsquare(uint32_t rij[3], int rij_valid,
		   uint32_t *r2, uint32_t *r2p, int *r2_valid, 
		   uint64_t *r2_soft, int *r2_soft_gt1);
  void mdg_sub32(uint32_t ri, uint32_t rj,
		 uint32_t *z, int32_t *zs, int *inrange);

  uint64_t mdg_fmul(uint64_t x, uint64_t y, int xmw, int xew, int ymw, int yew, int zmw, int zew);
  uint32_t mdg_charge_mult(int32_t qi, int32_t qj);
  double mdg_rev_charge_mult(uint32_t qiqj); //for debug (Komatsu)

  uint64_t mdg_coulomb_prefactor(uint32_t qiqj_gbc, uint32_t rinv);
  double mdg_rev_coulomb_prefactor(uint64_t x); //for debug (Komatsu)

  uint64_t mdg_funcev_vdw(uint32_t x, int exp_max, uint32_t underflow_val, mdg_func_table *t);
  uint64_t mdg_funcev_vdw_array(uint32_t x, int exp_max, uint32_t underflow_val, 
				int64_t vdw_table[512][4], int vdw_exp_table[512],
				char *name);

  uint32_t mdg_vdw_eps_mult(uint32_t ei, uint32_t ej);

  uint32_t mdg_vdw_sigma_sum(uint32_t si, uint32_t sj);
  uint32_t mdg_vdw_sigma_inv(uint32_t sij);

  uint32_t mdg_vdw_sigma_mult(uint32_t si, uint32_t sj);
  double mdg_rev_vdw_sigma_mult(uint32_t x); //for debug (Komatsu)

  uint32_t mdg_vdw_eij_div_sij(uint32_t eij, uint32_t invsij, int pot_mode);
  uint32_t mdg_vdw_group_cutoff(uint32_t eij, uint32_t gbc);

  uint32_t mdg_vdw_rij_div_sij(uint32_t rij, uint32_t invsij);
  double mdg_rev_vdw_rij_div_sij(uint32_t x); //for debug (Komatsu)

  void mdg_vdw_kernel(uint32_t eij_sij, uint64_t res_fev_vdw, int norm_vdw, int round_64,
		      int64_t *z, int *zovf);

  void mdg_coulomb_kernel(uint64_t coulomb_core, int32_t fev_excluded,
			  int norm_coulomb, int round_64,
			  int64_t *z, int *zovf);

  void mdg_force_factor(int64_t ckernel, int ckernel_overflow,
			int64_t vkernel, int vkernel_overflow,
			int negate,
			int64_t *fkernel, int *fkernel_overflow);

  void mdg_force_mul(int32_t r, int64_t fkernel, int fkernel_overflow,
		     int32_t *f, int *f_ovf);

  uint32_t mdg_coulomb_group_cutoff(uint32_t x, uint32_t y);
  double mdg_rev_coulomb_group_cutoff(uint32_t x); //for debug (Komatsu)

  int32_t mdg_coulomb_exclusion(int32_t res_fev_coul, int32_t noexcl_factor,
				int32_t exclfull_factor, int32_t excl14_factor,
				int excl_full, int excl_14 );

  uint32_t mdg_coulomb_softcore( uint64_t r2, uint32_t k, int potential_mode);
  uint64_t mdg_vdw_softcore( uint32_t rs, uint32_t k0, int potential_mode );

  void mdg_add32_sat(int32_t x, int xovf, int32_t y, int yovf, int32_t *z, int *zovf);
  void mdg_sub32_sat(int32_t x, int xovf, int32_t y, int yovf, int32_t *z, int *zovf);
  void mdg_add64_sat(int64_t x, int xovf, int64_t y, int yovf, int64_t *z, int *zovf);
  void mdg_add_sat(int n, int64_t x, int xovf, int64_t y, int yovf, int64_t *z, int *zovf);

  int mdg_checkrange(int32_t ri[3], int32_t rj[3], int cutoff);

#ifdef __cplusplus
}
#endif

#endif

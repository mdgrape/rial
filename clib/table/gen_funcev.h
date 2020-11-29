#ifndef __GEN_FUNCEV_H
#define __GEN_FUNCEV_H

#include "gentable.h"

#ifdef __cplusplus
extern "C" {
#endif
/* gen_funcev.c */

#define MDG_FUNCEV_COULOMB_BP 26
#define MDG_FUNCEV_COULOMB_NORDER 3
#define MDG_FUNCEV_COULOMB_R2MAX 4.0
#define MDG_FUNCEV_COULOMB_ADRW 8
#define MDG_FUNCEV_COULOMB_EXPW 4
#define MDG_FUNCEV_COULOMB_DEPTH (1<<MDG_FUNCEV_COULOMB_ADRW)

  extern const int mdg_funcev_coulomb_table_width[4];
  extern const int mdg_funcev_coulomb_calc_width[4];
  extern const int mdg_funcev_coulomb_table_bit_div;
  
#define MDG_FUNCEV_VDW_BP 25
#define MDG_FUNCEV_VDW_NORDER 3
  // In the followintg setting,
  //   Definition range : [2^-3, 2^5) = [1/8, 32)
#define MDG_FUNCEV_VDW_RS2MAX_EXP 5
#define MDG_FUNCEV_VDW_RS2MAX ((double)(1<<MDG_FUNCEV_VDW_RS2MAX_EXP))
#define MDG_FUNCEV_VDW_ADRW 9
#define MDG_FUNCEV_VDW_EXPW 3
#define MDG_FUNCEV_VDW_DEPTH (1<<MDG_FUNCEV_COULOMB_ADRW)

  extern const int mdg_funcev_vdw_table_width[4];
  extern const int mdg_funcev_vdw_calc_width[4];
  extern const int mdg_funcev_vdw_table_bit_div;
  extern const int mdg_funcev_vdw_table_exp_width;
  
  void mdg_gentable_funcev_coulomb(mdg_func_table *t, double (*func)(double), const char *name);
  int mdg_pack_funcev_table(mdg_func_table *t, int n, const int *width, int packwidth, uint64_t *packed);
  void mdg_gentable_funcev_vdw(mdg_func_table *t, double (*func)(double), const char *name);

  int mdg_pack_coulomb_table(int64_t *p, uint64_t *packed);
  int mdg_pack_vdw_table(int64_t *p, int ex, uint64_t *packed);
  void mdg_unpack_coulomb_table(uint64_t *packed, int64_t *p);
  void mdg_unpack_vdw_table(uint64_t *packed, int64_t *p, int *ex);

  int mdg_pack32_coulomb_table(int64_t *p, uint32_t *packed);
  int mdg_pack32_vdw_table(int64_t *p, int ex, uint32_t *packed);

#ifdef __cplusplus
}
#endif

#endif

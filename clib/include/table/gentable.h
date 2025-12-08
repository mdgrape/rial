#ifndef __GENTABLE_H
#define __GENTABLE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
/* gentable.c */
typedef struct {
  int norder;
  int max_entries;
  int nentries;
  int nprec;
  int floating;
  double *xmin;
  double *w;
  double *c;
  double *c1;
  double *cp;
  int    *exp;
  int64_t *ci;
  int *cbit;
  int expbits;
  int exporg;
  int *sign;
  int *hsign;
  char *name;
  int allow_carryup;
  int exp_max;
} mdg_func_table;

  void mdg_table_init(mdg_func_table *t, int norder, int nmax, const char *name, int allow_carryup);
  void mdg_table_add_intervals(double (*f)(double), int n, double xmin, double w, mdg_func_table *t);
  void mdg_table_integer_conversion(int nprec, mdg_func_table *t);
  void mdg_table_integer_conversion_float( int nprec, mdg_func_table *t);
  double mdg_table_eval(double x, int k, mdg_func_table *t);
  int64_t mdg_table_eval_ii(int64_t dx, int dxprec, int k, mdg_func_table *t);
  double mdg_table_eval_i(double x, int dxprec, mdg_func_table *t);
  void mdg_table_plot(int n, int ylog, double (*f)(double), mdg_func_table *t);
  void mdg_table_output_verilog(char *name, int sign_mode, mdg_func_table *t, char *command);
  void mdg_table_output_c(char *name, mdg_func_table *t, char *command);
  void mdg_table_integer_conversion_coulomb(int one, mdg_func_table *t);
#ifdef __cplusplus
}
#endif

#endif

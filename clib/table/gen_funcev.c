/*
 Table Generation for function evaluator
*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include "defs.h"
#include "gentable.h"
#include "gen_funcev.h"
#include "fpmd.h"	//Added (Komatsu) for debug

//int MDG_DEBUG_FLAG;

const int mdg_funcev_coulomb_table_bit_div = 42;
const int mdg_funcev_coulomb_table_width[4] = { 27, 21, 18, 18 };
const int mdg_funcev_coulomb_calc_width[4] = { 28, 22, 19, 18 };

const int mdg_funcev_vdw_table_bit_div  = 50;
const int mdg_funcev_vdw_table_exp_width= 7;
const int mdg_funcev_vdw_table_width[4] = { 27, 27, 22, 17 };
const int mdg_funcev_vdw_calc_width[4]  = { 28, 28, 23, 17 };

#ifdef DEBUG
static int parity64(uint64_t x) {
  x ^= x>>32;
  x ^= x>>16;
  x ^= x>>8;
  x ^= x>>4;
  x ^= x>>2;
  x ^= x>>1;
  return x&1;
}
#endif

void
mdg_gentable_funcev_coulomb(mdg_func_table *t, double (*func)(double), const char *name) {
  int norder  = MDG_FUNCEV_COULOMB_NORDER;
  double xmax = MDG_FUNCEV_COULOMB_R2MAX;
  int adrw    = MDG_FUNCEV_COULOMB_ADRW;
  int nmax    = 1<<adrw;
  int expw    = MDG_FUNCEV_COULOMB_EXPW;
  int e_max   = 1<<expw;
  // xmax corresponds to exponent == 1<<expw
  // min exponent == 1 since 0 is used for fixed-point region
  double xmin = scalbn(xmax, -((1<<expw)-1));

  int manw = adrw-expw;

  mdg_table_init(t, norder, nmax, name, 1);

  int nblock = 1<<manw;
  if(MDG_DEBUG_FLAG)printf("INFO (%s): man_bits=%d nblock=%d\n", __func__, manw, nblock);
  int n=0;
  // Add interval from 0~xmin
  double w = scalbn(xmin, -manw);
  mdg_table_add_intervals(func, nblock, 0.0, w, t);
  if(MDG_DEBUG_FLAG)printf("INFO (%s): Block %d : xmin=%f width=%e xmax=%f\n",
	 __func__, n, 0.0, w, nblock*w);
  ++n;
  for(int ex=1; ex<e_max;++ex,++n) {
    double x = scalbn(xmin, ex-1);
    double w = scalbn(xmin, ex-1-manw);
    if(MDG_DEBUG_FLAG)printf("INFO (%s): Block %d : xmin=%f width=%e xmax=%f\n",
	   __func__, n, x, w, x+nblock*w);
    mdg_table_add_intervals(func, nblock, x, w, t);
  }

  mdg_table_integer_conversion_coulomb(MASK(MDG_FUNCEV_COULOMB_BP), t);
}

int
mdg_pack_funcev_table(mdg_func_table *t, int n, const int *width, int packwidth, uint64_t *packed) {
  int64_t *p=t->ci;
  p += (t->norder+1)*n;

  int w=0;
  int j=0;
  __uint128_t z=0;
  for(int i=t->norder; i>=0; --i) {
    __uint128_t c = p[i] & MASK64(width[i]);
    //    printf("%llx\n",c);
    z |= c<<w;
    w += width[i];
    //printf("%d %d %d\n", w, j, packwidth);
    if (w >= packwidth) {
      packed[j] = z & MASK64(packwidth);
      z >>= packwidth;
      w -= packwidth;
      ++j;
    }
  }
  if (t->floating) {
    __uint128_t c = t->exp[n] & MASK64(t->expbits);
    //    printf("%llx\n",c);
    z |= c<<w;
    w += t->expbits;
    //printf("%d %d %d\n", w, j, packwidth);
    if (w >= packwidth) {
      packed[j] = z & MASK64(packwidth);
      z >>= packwidth;
      w -= packwidth;
      ++j;
    }
  }
  if (w>0) {
    packed[j] = z;
    ++j;
  }
  return j;
}

void
mdg_gentable_funcev_vdw(mdg_func_table *t, double (*func)(double), const char *name) {
  int norder  = MDG_FUNCEV_VDW_NORDER;
  double xmax = MDG_FUNCEV_VDW_RS2MAX;
  int adrw    = MDG_FUNCEV_VDW_ADRW;
  int nmax    = 1<<adrw;
  int expw    = MDG_FUNCEV_VDW_EXPW;
  int e_max   = 1<<expw;
  // xmax corresponds to exponent == 1<<expw
  // min exponent == 1 since 0 is used for fixed-point region
  double xmin = scalbn(xmax, -(1<<expw));
  int exp_max = ilogb(xmax)-1;
  t->exp_max = exp_max;

  int manw = adrw-expw;

  mdg_table_init(t, norder, nmax, name, 1);

  int nblock = 1<<manw;
  if(MDG_DEBUG_FLAG)printf("INFO (%s): man_bits=%d nblock=%d\n", __func__, manw, nblock);
  int n=0;
  for(int ex=0; ex<e_max;++ex,++n) {
    double x = scalbn(xmin, ex);
    double w = scalbn(xmin, ex-manw);
    if(MDG_DEBUG_FLAG)printf("INFO (%s): Block %d : xmin=%f width=%e xmax=%f\n",
	   __func__, n, x, w, x+nblock*w);
    mdg_table_add_intervals(func, nblock, x, w, t);
  }

  mdg_table_integer_conversion_float(MDG_FUNCEV_VDW_BP, t);
}

int
mdg_pack_coulomb_table(int64_t *p, uint64_t *packed) {
  int w=0;
  int j=0;
  __uint128_t z=0;
  for(int i=3; i>=0; --i) {
    int width = mdg_funcev_coulomb_table_width[i];
    __uint128_t c = p[i] & MASK64(width);
    //    printf("%llx\n",c);
    z |= c<<w;
    w += width;
    //printf("%d %d %d\n", w, j, packwidth);
    if (w >= mdg_funcev_coulomb_table_bit_div) {
      packed[j] = z & MASK64(mdg_funcev_coulomb_table_bit_div);
      z >>= mdg_funcev_coulomb_table_bit_div;
      w -= mdg_funcev_coulomb_table_bit_div;
      ++j;
    }
  }
  if (w>0) {
    packed[j] = z;
    ++j;
  }
  return j;
}

void
mdg_unpack_coulomb_table(uint64_t *packed, int64_t *p) {
  __uint128_t z=packed[1] & MASK64(mdg_funcev_coulomb_table_bit_div);
  z <<= mdg_funcev_coulomb_table_bit_div;
  z += packed[0] & MASK64(mdg_funcev_coulomb_table_bit_div);
  for(int i=3; i>=0; --i) {
    int width = mdg_funcev_coulomb_table_width[i];
    p[i] = z & MASK64(width);
    if (i) p[i] = EXTENDLL(p[i], width);
    z >>=width;
  }
}

int
mdg_pack_vdw_table(int64_t *p, int ex, uint64_t *packed) {
  int width = mdg_funcev_vdw_table_exp_width;
  __uint128_t z=ex & MASK64(width);
  for(int i=0; i<4; ++i) {
    width = mdg_funcev_vdw_table_width[i];
    z = (z << width) + (p[i] & MASK64(width));
  }
  packed[0] = z & MASK64(mdg_funcev_vdw_table_bit_div);
  packed[1] = z >> mdg_funcev_vdw_table_bit_div;
  return 2;
}

void
mdg_unpack_vdw_table(uint64_t *packed, int64_t *p, int *ex) {
  __uint128_t z=packed[1] & MASK64(mdg_funcev_vdw_table_bit_div);
  z <<= mdg_funcev_vdw_table_bit_div;
  z += packed[0] & MASK64(mdg_funcev_vdw_table_bit_div);
  for(int i=3; i>=0; --i) {
    int width = mdg_funcev_vdw_table_width[i];
    p[i] = z & MASK64(width);
    if (i) p[i] = EXTENDLL(p[i], width);
    z >>=width;
  }
  *ex = z & MASK64(mdg_funcev_vdw_table_exp_width);
}

static int parity(uint32_t x) {
  x ^= x>>16;
  x ^= x>>8;
  x ^= x>>4;
  x ^= x>>2;
  x ^= x>>1;
  return x&1;
}

int
mdg_pack32_coulomb_table(int64_t *p, uint32_t *packed) {
  int w=0;
  const int w0=mdg_funcev_coulomb_table_bit_div/2;
  int j=0;
  uint64_t z=0;
  int par=0;
  for(int i=3; i>=0; --i) {
    int width = mdg_funcev_coulomb_table_width[i];
    uint64_t c = p[i] & MASK64(width);
    //    printf("%llx\n",c);
    z |= c<<w;
    w += width;
    //printf("%d %d %d\n", w, j, packwidth);
    if (w >= w0) {
      packed[j] = z & MASK64(w0);
      z >>= w0;
      w -= w0;
      par ^= parity(packed[j])<<(j/2);
      ++j;
    }
  }
  if (w>0) {
    packed[j] = z;
    par ^= parity(packed[j])<<(j/2);
    ++j;
  }
  // add parity bit
  for(int i=0;i<j;i+=2) {
    packed[i+1] += BIT(par,i/2)<<w0;
  }
  
  return j;
}

int
mdg_pack32_vdw_table(int64_t *p, int ex, uint32_t *packed) {
  int w=0;
  const int w0=mdg_funcev_vdw_table_bit_div/2;
  int j=0;
  uint64_t z=0;
  int par=0;
  for(int i=3; i>=0; --i) {
    int width = mdg_funcev_vdw_table_width[i];
    uint64_t c = p[i] & MASK64(width);
    //    printf("%llx\n",c);
    z |= c<<w;
    w += width;
    //printf("%d %d %d\n", w, j, packwidth);
    if (w >= w0) {
      packed[j] = z & MASK64(w0);
      z >>= w0;
      w -= w0;
      par ^= parity(packed[j])<<(j/2);
      ++j;
    }
  }
  int width = mdg_funcev_vdw_table_exp_width;
  uint64_t c=ex & MASK64(width);
  z |= c<<w;
  w += width;
  if (w>0) { // Always w < w0
    packed[j] = z;
    par ^= parity(packed[j])<<(j/2);
    ++j;
  }
  // add parity bit
  for(int i=0;i<j;i+=2) {
    packed[i+1] += BIT(par,i/2)<<w0;
  }  
  return j;
}

//  fpmd.c
//  
//  Simulator of arithmetic units
//  used in MD chip
//
//  Copyright (c) 1993, 1994 by Makoto Taiji/University of Tokyo
//  taiji@grape.c.u-tokyo.ac.jp
//  All rights reserved
//
//  Copyright (c) 2002, 2003 by Makoto Taiji/RIKEN GSC
//  taiji@gsc.riken.go.jp
//
//  Copyright (c) 2010 by Makoto Taiji/RIKEN ASI
//  taiji@riken.jp
//
//  Copyright (c) 2011 by Makoto Taiji/RIKEN QBiC
//  taiji@riken.jp
//
//  Copyright (c) 2016 by Makoto Taiji/RIKEN QBiC
//  taiji@riken.jp
//
//  (2017) bug fix by teruhisa.komatsu@riken.jp
//

/** @file fpmd.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include "defs.h"
#include "fpmd.h"
#include "gentable.h"
#include "gen_funcev.h"
//#include "md.h"

int B_MDG_CHECK=0;
int mdg_set_check_flag(int flag) {
  if (flag<0) return MDG_DEBUG_FLAG;
  else return MDG_DEBUG_FLAG=flag; 
}

int MDG_DEBUG_FLAG=0;
int mdg_set_debug_flag(int flag) {
  if (flag<0) return MDG_DEBUG_FLAG;
  else return MDG_DEBUG_FLAG=flag; 
}

void
mdg_get_environment() {
  if (getenv("B_MDG_CHECK")!=NULL) B_MDG_CHECK=1;
  if (getenv("MDG_DEBUG_FLAG")!=NULL) MDG_DEBUG_FLAG=atoi(getenv("MDG_DEBUG_FLAG"));
}

//
// Leading i position
// return 0-64
// return 0 if x==0
static int ilogbi(uint64_t x) {
  int n=64;
  while ((n>0)&&((x&0x8000000000000000ULL)==0)) {
    --n;
    x<<=1;
  }
  return n;
}

/**
   @brief unbiased rounding 
   @param[in] x   input value
   @param[in] pos rounding position
   @retval rounded result, shifted
*/
static uint64_t round_unbiased_u64( uint64_t x, int pos) {
  if (pos<=0) return x;
  int stickey=0;
  if (pos>1) stickey = (x & MASK(pos-1)) ? 1 : 0;
  int round=BIT(x,pos-1);
  int lsb = BIT(x,pos);
  x>>=pos;
  x += round & (stickey | lsb);
  return x;
}

/**
   @brief biased rounding 
   @param[in] x   input value
   @param[in] pos rounding position
   @retval rounded result, shifted
*/
static uint64_t round_biased_u64( uint64_t x, int pos) {
  if (pos<=0) return x;
  int round=BIT(x,pos-1);
  x>>=pos;
  x += round;
  return x;
}

/**
   @brief polynomial evaluation with fixed sign
   @param[in] adr   table address
   @param[in] d     fixed-point difference from center
   @param[in] dbit  to specify binary point
   @param[in] norder polynomial order
   @param[in] adrmax table entry max
   @param[in] sign  sign of coefficients 
   @param[in] cbit  width of coefficients
   @param[in] table pointer to table
   @param[in] name  name of function for debug
   @retval evaluation result in fixed point
*/
static uint32_t
polynomial_fixedsign32(int adr, int d, int dbit, int norder, int adrmax,
		       int *sign, int *cbit, int *table, char *name)
{
  if (adr>=adrmax) {
    printf("ERROR (%s:%s) :address out of range %d >= %d\n",
	   __func__, name, adr, adrmax);
    return 0;
  }
  int *c = table+adr*(norder+1);
  
  int64_t dmax=(1LL<<(dbit-1));
  int64_t dl=d-dmax;
  int64_t s=c[norder];

  for(int i=norder-1;i>=0;--i) {
    int64_t prod = s*dl;
    prod >>= dbit-1;
    s = c[i] + sign[i+1]*prod;
    if (llabs(s) >= (1ULL<<cbit[i])) {
      printf("WARNING (%s:%s) : carry up occurs at order=%d, adr=%d, d=%d\n",
	     __func__, name, i, adr, d);
    } else if (s<0) {
      printf("WARNING (%s:%s) : negative sum at order=%d, adr=%d, d=%d\n",
	     __func__, name, i, adr, d);
    }
  }
  s *= sign[0];

  return s;
}

static uint32_t
polynomial_withsign32(int adr, int d, int dbit, int norder, int adrmax,
		      int *cbit, int *table, const char *name, int ignore_overflow_at_c0)
{
  if (adr>=adrmax) {
    printf("ERROR (%s:%s) :address out of range %d >= %d\n",
	   __func__, name, adr, adrmax);
    return 0;
  }
  int *c = table+adr*(norder+1);
  
  int64_t dmax=(1LL<<(dbit-1));
  int64_t dl=d-dmax;
  int64_t s=c[norder];

  for(int i=norder-1;i>=0;--i) {
    int64_t prod = s*dl;
    //printf("prod %16llx ", prod);
    prod >>= dbit-1;
    //printf(" %16llx ", prod);
    s = c[i] + prod;
    //printf(" c[%d]=%8x s=%16llx cbit=%d\n", i, c[i], s, cbit[i]);
    if ( !(ignore_overflow_at_c0 && (i==0))
	 && ( (s >= (1LL<<cbit[i]))||(s < -(1LL<<cbit[i])) ) ) { 
      printf("WARNING (%s:%s) : carry up occurs at order=%d, adr=%d, d=%d\n",
	     __func__, name, i, adr, d);
    }
  }
  return s;
}

#include "rinv.h"
#include "rinv3.h"
#include "r2sqrt.h"

static const uint32_t sqrt2=47453133; 

double mdg_rev_inv_sqrt3(uint32_t x){ //for debug (Komatsu)
	return mdg_rev_inv_sqrt(x);
}
double mdg_rev_inv_sqrt(uint32_t x){ //for debug (Komatsu)
	const int mw=23;
	const int ew=7;
	const int bias=32;
	double z=TO_DOUBLE_U(x,mw,ew,bias);
	return z;
}

// Calculate 1/r from r^2
// input : 25-bit mantissa and 6-bit exponent, base=0x3e
// output: 23-bit mantissa and 7-bit exponent, base=32  =bias?
//
uint32_t mdg_inv_sqrt(uint32_t x)
{
  if (x>=(1U<<31)) {
    fprintf(stderr,"WARNING (%s) : illegal input %x >= 1<<31\n",
	    __func__, x);
  } 
  int ex  = (x>>25) & MASK(6);
  int adr = (x>>(25-7)) & MASK(7);
  int d   = x & MASK(25-7);
  //  printf("%s : %08x -> exp = %1x adr = %02x d = %05x\n",
  // __func__, x, exp, adr, d);
  uint32_t s = polynomial_withsign32(adr, d, 25-7, NORDER_rinv, ADRMAX_rinv, 
				     rinv_cbit, (int *)rinv_table, NAME_rinv, 1);
  uint64_t z = (ex&1) ? sqrt2 : (1<<25);
  z *= s;
  //fprintf(stderr,"ex=%d s=%8x z=%16llx\n", ex&1, s, z);

  int n = ilogbi(z);
  if ((n>52)||(n<50)) {
    fprintf(stderr,"WARNING (%s) : something wrong with leading one position %d , %" PRIx64 "\n",
	   __func__, n, z);
  }
  //  fprintf(stderr,"shift=%16llx\n", z>>(n-27));
  z = (z>>(n-25))&MASK(24);
  if (ex&1) --n;
  ++z; // rounding
  z >>= 1;
  if (BIT(z,23)) { ++n; z=0; }

  // input exponent : 0x3f -> [2,4)
  // input exponent : 0x3e -> [1,2) -> 1/r in (1/2,1]
  // n==52 -> exponent = 0, exponent offset = 32
  int zexp = 32+(0x1f-(ex>>1)) + (n-52);
  if (ex==0) {
    zexp = 0;
    z = 0;
  }
  if ((zexp<0)||(zexp>=128)) {
    fprintf(stderr,"WARNING (%s) : result exponent out of range %d for input exponent %d\n",
	    __func__, zexp, ex);
  }
  uint32_t res = (z&MASK(23)) + ( (zexp&MASK(7)) << 23);
  //printf("res = %08x\n", res);
  return res;
}

// Calculate 1/r^3 from r^2
// input : 25-bit mantissa and 6-bit exponent, base=0x3e
// output: 25-bit mantissa and 7-bit exponent, base=32
//
uint32_t mdg_inv_sqrt3(uint32_t x)
{
  if (x>=(1U<<31)) {
    fprintf(stderr,"WARNING (%s) : illegal input %x >= 1<<31\n",
	   __func__, x);
  } 
  int ex  = (x>>25) & MASK(6);
  int adr = (x>>(25-7)) & MASK(7);
  int d   = x & MASK(25-7);
  //  printf("%s : %08x -> exp = %1x adr = %02x d = %05x\n",
  // __func__, x, exp, adr, d);
  uint32_t s = polynomial_withsign32(adr, d, 25-7, NORDER_rinv3, ADRMAX_rinv3, 
				     rinv3_cbit, (int *)rinv3_table, NAME_rinv3, 1);
  uint64_t z = (ex&1) ? sqrt2 : (1<<25);
  z *= s;
  //fprintf(stderr,"ex=%d s=%8x z=%16llx\n", ex&1, s, z);
  
  int n = ilogbi(z);
  if ((n>52)||(n<50)) {
    fprintf(stderr,"WARNING (%s) : something wrong with leading one position %d , %" PRIx64 "\n",
	   __func__, n, z);
  }
  //  printf("shift=%16llx\n", z>>(n-27));
  z = (z>>(n-25))&MASK(24);
  if (ex&1) n-=2;
  ++z; // rounding
  z >>= 1;
  if (BIT(z,23)) { ++n; z=0; }

  // input exponent : 0x3f -> [2,4)
  // input exponent : 0x3e -> [1,2) -> 1/r in (1/2,1]
  // n==52 -> exponent = 0, exponent offset = 32
  int zexp = 32+(0x1f*3-(ex>>1)*3) + (n-52);
  if (ex==0) {
    zexp = 0;
    z = 0;
  }
  if ((zexp<0)||(zexp>=128)) {
    fprintf(stderr,"WARNING (%s) : result exponent out of range %d for input exponent %d\n",
	   __func__, zexp, ex);
  }
  uint32_t res = (z&MASK(23)) + ( (zexp&MASK(7)) << 23);
  //printf("res = %08x\n", res);
  return res;
}

/**
   @brief calc sqrt of r^2
   @param[in] x r^2 in floating point (5 bit exponent, 25 bit mantissa, exponent bias=0x1e)
   @retval evaluation result in floating point (6 bit exponent, 25 bit mantissa, exponent bias=0x20)
*/
uint32_t mdg_r2sqrt(uint32_t x)
{
  if (x>=(1U<<31)) {
    printf("WARNING (%s) : illegal input %x >= 1<<31\n",
	   __func__, x);
  } 
  int ex  = (x>>25) & MASK(6);
  int adr = (x>>(25-7)) & MASK(7);
  int d   = x & MASK(25-7);
  //  printf("%s : %08x -> exp = %1x adr = %02x d = %05x\n",
  // __func__, x, exp, adr, d);

  // r2sqrt table calculates sqrt(x)/2
  uint32_t s = polynomial_withsign32(adr, d, 25-7, NORDER_r2sqrt, ADRMAX_r2sqrt, 
				     r2sqrt_cbit, (int *)r2sqrt_table, NAME_r2sqrt, 1);
  uint64_t z = (ex&1) ? sqrt2 : (1<<25);
  z *= s;
  //  printf("ex=%d s=%8x z=%16llx ", ex&1, s, z);

  int n = ilogbi(z);
  if ((n>52)||(n<50)) {
    printf("WARNING (%s) : something wrong with leading one position %d , %" PRIx64 "\n",
	   __func__, n, z);
  }
  //  printf("shift=%16llx\n", z>>(n-27));
  z = (z>>(n-25))&MASK(24);
  ++z; // rounding
  z >>= 1;
  if (BIT(z,23)) { ++n; z=0; }

  // input exponent : 0x3f -> [2,4)
  // input exponent : 0x3e -> [1,2) -> 1/r in (1/2,1]
  // n==52 -> exponent = 0, exponent offset = 32
  // 1 added to multiply sqrt(x)/2 by 2 
  int zexp = 32+((ex>>1)-0x1f) + (n-52) + 1;
  if (ex==0) {
    zexp = 0;
    z = 0;
  }
  if ((zexp<0)||(zexp>=128)) {
    printf("WARNING (%s) : result exponent out of range %d for input exponent %d\n",
	   __func__, zexp, ex);
  }
  uint32_t res = (z&MASK(23)) + ( (zexp&MASK(7)) << 23);
  //printf("res = %08x\n", res);
  return res;
}

/**
   @brief polynomial function evaluator for funcev_coulomb/funcev_vdw
   @param[in] dx     fixed-point difference from center
   @param[in] dbit   to specify binary point
   @param[in] norder polynomial order
   @param[in] calc_width  width of intermediate results after addition with sign
   @param[in] table  [norder] table entry
   @param[in] name   function name (to display error)
   @param[in] x      input value (to display error)
*/
static int64_t
polynomial_funcev(int32_t dx, int dbit, int norder,
		  const int *calc_width, const int64_t *table, const char *name, double x)
{
  // Note: in C, a result of division is rounded toward 0.
  // That's why I use double for this calculation..

  double dxf = dx;
  double dxnorm = 1<<(dbit-1);
  double dxl;
  double dxw;

  double z = table[norder];
  for (int i=norder-1; i>=0; --i) {
    int width = calc_width[i+1];
    if (dbit > width) {
      dxw = 1<<(width-1);
      dxl = floor(dxf / (1<<(dbit-width)) );
    } else {
      dxw = dxnorm;
      dxl = dxf;
    }
#ifdef DEBUG_VERBOSE
    double z0 = z;
    double prod = dxl * z;
#endif
    double dz = floor(dxl * z / dxw);
    z = table[i] + dz;
#ifdef DEBUG_VERBOSE
    printf("%d: %llx %llx %d %llx %x %x %x\n", i, (int64_t) z, table[i], width, (int64_t) prod, (int32_t)dz, (int32_t) dxl, (uint32_t) dxw);
#endif
    int64_t cmax, cmin;
    cmax = (1<<(calc_width[i]-1))-1;
    cmin = -cmax-1;
    if ( (z>cmax) || (z<cmin) ) { 
      fprintf(stdout, "ERROR(%s): range error for %s(%f): order=%d, dx=%x(%f), z=%" PRIx64 "(%le), c=%" PRIx64 ", cmin=%" PRIx64 ", cmax=%" PRIx64 "\n",
	      __func__, name, x, i, dx, dxf/dxnorm, (int64_t) z,z, table[i], cmin, cmax);
    }
  }
  return (int64_t) z;
}

/**
   @brief function evaluator for coulomb, common routine
   @param[in] x r^2 in partial floating, 25-bit mantissa and 4-bit exponent
   @param[in] polynomial order (should be 3 now)
   @param[in] pointer to table entry
   @param[in] function name for debug
   @retval evaluation result in fixed point
*/
static int32_t
mdg_funcev_coulomb_common(uint32_t x, int norder, int64_t *c, char *name) {
  int32_t dx = (x & MASK(21)) - (1<<20);
  int exp = SLICE(x,25,4);
  double xf;
  if (exp==0) { // fixed point region
    xf = x & MASK(25);
    xf = scalbn(xf, -24 + 1 - 15);
  } else { // floating region
    xf = (x & MASK(25)) + (1<<25);
    // max 2^26-1 == 4-eps
    //     2^25 -> 2
    // exp==15 -> x in [2,4)
    xf = scalbn(xf, -24 + exp - 15);
  }
  
  int64_t z = polynomial_funcev(dx, 21, norder,
				mdg_funcev_coulomb_calc_width,
				c, name, xf);
  if (z<0) z=0;
  else if (z>=(1<<26)) z=MASK(26);
  return (int32_t) z;
}

/**
   @brief function evaluator for coulomb, use struct 
   @param[in] x r^2 in partial floating, 25-bit mantissa and 4-bit exponent
   @param[in] t pointer to table
   @retval evaluation result in fixed point
*/
int32_t
mdg_funcev_coulomb(uint32_t x, mdg_func_table *t) {
  int adr = SLICE(x,21,8);
  int64_t *c=t->ci+adr*(t->norder+1);
  //  printf("adr = %x, c = %llx, %llx, %llx, %llx\n", adr, c[0], c[1], c[2], c[3]);

  return mdg_funcev_coulomb_common(x, t->norder, c, t->name);
}

/**
   @brief function evaluator for coulomb, use array for table
   @param[in] x r^2 in partial floating, 25-bit mantissa and 4-bit exponent
   @param[in] c pointer to table 
   @param[in] function name
   @retval evaluation result in fixed point
*/
int32_t
mdg_funcev_coulomb_array(uint32_t x, int64_t c[512][4], char *name) {
  int adr = SLICE(x,21,8);
  //  printf("adr = %x, c = %llx, %llx, %llx, %llx\n", adr, c[adr][0], c[adr][1], c[adr][2], c[adr][3]);

  return mdg_funcev_coulomb_common(x, 3, c[adr], name);
}

/**
   @brief r2 format conversion from double (max 4.0)
   @param[in] x r^2 in partial floating, 25-bit mantissa and 4-bit exponent
   @param[out] z1 r^2 in floating, 25-bit mantissa and 6-bit exponent
   @param[out] z2 r^2 in partial floating, 25-bit mantissa and 4-bit exponent
*/
void
mdg_conv_r2_fromDouble(double x, uint32_t *z1, uint32_t *z2) {
  int ex = ilogb(x);
  // exponent: 1 ~ -62
  if ((ex>1)||(ex<-62)) {
    printf("WARNING(%s): exponent out of range %d\n",__func__, ex);
    *z1 = *z2 = 0;
  } else {
    // 53 bit leading 1
    uint64_t man = scalbn(x, 52-ex);
    *z1 = man>>27;
    // if exponent in 1..-13, z2 is floating
    // otherwise, fixed
    int ex2;
    if (ex>=-13) {
      *z2 = *z1;
      ex2 = ex+14;
    } else {
      *z2 = man>>(28 + (-ex-14));
      ex2 = 0;
    }
    *z1 &= MASK(25);
    *z2 &= MASK(25);

    *z1 += (ex+62)<<25;
    *z2 += ex2<<25;
  }
}

/**
   @brief calculate r^2 from xij 
   @param[in] rij rij in 28-bit fixed each
   @param[in] rij_valid 4-bit enable input
   @param[out] r2 r^2 in floating, 25-bit mantissa and 6-bit exponent
   @param[out] r2p r^2 in partial floating, 25-bit mantissa and 4-bit exponent
   @param[out] r2_valid enable output
   @param[out] r2_soft r2 for softcore, fixed 34-bit, range <1  
   @param[out] r2_soft_gt1 r2 for softcore is greater than 1
*/
void
mdg_rsquare(uint32_t rij[3], int rij_valid,
	    uint32_t *r2, uint32_t *r2p, int *r2_valid,
	    uint64_t *r2_soft, int *r2_soft_gt1
	    ) {
  *r2_valid = (rij_valid==0xf) ? 1 : 0;
  uint64_t sum = ((uint64_t)rij[0])*rij[0]
    + ((uint64_t)rij[1])*rij[1]
    + ((uint64_t)rij[2])*rij[2];
  if (sum>=(1ULL<<56)) {
    *r2_valid=0;
  }
  uint64_t x = sum & MASK64(56);
  int lop=0;
  while(x) { ++lop; x>>=1; }

  int ex0 = 0x3f - (56-lop);
  //printf("rij = %08x %08x %08x sum = %llx ex0=%d\n", rij[0],rij[1],rij[2],sum,ex0);
  uint64_t man0 = sum << (56-lop);
  int round = (man0>>(56-27)) &1;
  uint32_t man = (man0>>(56-26)) & MASK(25);
  man += round;

  int ex=ex0;
  if (man>=(1<<25)) {
    ex=ex0+1;
    man = 0;
    if (ex>=0x40) {
      *r2_valid=0;
      ex=0;
    }
  }
  *r2=(ex<<25) + man;
  
  int ex_p;
  uint32_t man_p;
  if (ex0>0x30) {
    ex_p = ex - 0x30;
    man_p = man;
  } else {
    ex_p = 0;
    uint64_t man_p0 = man0 >> (0x30-ex0+1);
    int round = (man_p0>>(56-27)) &1;
    man_p = (man_p0>>(56-26)) & MASK(25);
    man_p += round;
    if (man_p>=(1<<25)) {
      ex_p=1;
      man_p = 0;
    }
  }
  *r2p=(ex_p<<25) + man_p;

  // For softcore 
  // sum bp: at 54, take 24+7+1(round) = 32 bit : b53-b22
  *r2_soft = (sum>>22) & MASK64(32);
  *r2_soft_gt1 = (sum>>54)!=0;
}
 
/**
   @brief calculate rij from ri and rj (one component)
   @param[in] ri ri in 32-bit fixed
   @param[in] rj rj in 32-bit fixed
   @param[out] z   absolute rij in 28-bit fixed 
   @param[out] zs  rij with sign in 29-bit fixed
   @param[out] inrange rij < 2^28
*/
void
mdg_sub32(uint32_t ri, uint32_t rj,
	  uint32_t *z, int32_t *zs, int *inrange) {
  int32_t z0 = ri-rj;
  int32_t zabs = abs(z0);

  *inrange = (zabs < (1<<28)) ? 1 : 0;
  *z  = zabs & MASK(28);
  *zs = (z0 & MASK(28));
  if (z0<0) *zs += 1<<28;
}

/**
   @brief floating point multiplier 
   @param[in] x xman bit mantissa, xexp bit exponent, 1 bit sign
   @param[in] y yman bit mantissa, yexp bit exponent, 1 bit sign
   @param[in] xmw mantissa bit width of x
   @param[in] xew exponent bit width of x
   @param[in] ymw mantissa bit width of y
   @param[in] yew exponent bit width of y
   @param[in] zmw mantissa bit width of y
   @param[in] zew exponent bit width of y

   @retval x*y zman bit mantissa, zexp bit exponent, 1 bit sign
*/
uint64_t
mdg_fmul(uint64_t x, uint64_t y, int xmw, int xew, int ymw, int yew, int zmw, int zew) {
  __uint128_t xman = (x & MASK64(xmw)) + (1<<xmw);
  __uint128_t yman = (y & MASK64(ymw)) + (1<<ymw);

  int xexp = SLICE64(x, xmw, xew);
  int yexp = SLICE64(y, ymw, yew);

  int xsign = BIT(x, xmw+xew);
  int ysign = BIT(y, ymw+yew);
  
  __uint128_t zman = xman * yman;
  int gt2 = zman >> (xmw+ymw+2-1); // greater than 2
  if (gt2==0) zman <<= 1;
  uint64_t zman_sft = zman >> (xmw+ymw+2-zmw-1);
  int round = BIT(zman, xmw+ymw+2-zmw-2);
  zman_sft += round;
  int shift_after_round = BIT(zman_sft, zmw+1);
  zman_sft &=  MASK64(zmw);

  if (shift_after_round && gt2) {
    printf("WARNING(%s) : result becomes >= 4\n", __func__);
  }
  int zexp = 0;
  if (xexp && yexp) {
    zexp = xexp + yexp + gt2 + shift_after_round;
  } else {
    zexp = 0;
    zman_sft = 0;
  }

  if (zexp >= (1<<zew)) {	//(Komatsu)
    printf("ERROR(%s) : exponent overflow\n", __func__);
  }
  zexp &= MASK64(zew);
	       
  int sign = xsign ^ ysign;
  if (zexp==0) sign=0;
  uint64_t z = ( sign << (zmw+zew) ) + ( zexp << zmw ) + zman_sft;

  return z;
}

/**
   @brief calculate qiqj from qi and qj
   @param[in] qi qi in 32-bit fixed (26-bit fraction, 2's complement)
   @param[in] qj qj in 32-bit fixed (26-bit fraction, 2's complement)
   @retval qiqj 29-bit floating (23-bit mantissa, 5-bit exponent, 1-bit sign)
*/
uint32_t
mdg_charge_mult(int32_t qi, int32_t qj) {
  int64_t qiqj = ((int64_t)qi) * ((int64_t)qj);
  int sign = (qiqj<0) ? 1 : 0;
  uint64_t qiqj_abs = llabs(qiqj);
  //printf("%llx %llx ", qiqj, qiqj_abs);
  int lop=0;
  if (qiqj_abs==0) {
    lop=64;
  } else {
    while ((qiqj_abs & (1ULL<<63))==0) {
      ++lop;
      qiqj_abs<<=1;
    }
  }
  //printf("%d %llx\n", lop, qiqj_abs);
  uint32_t z;
  if (lop>=31) {
    z=0;
  } else {
    uint32_t ex = 31-lop;
    qiqj_abs >>= (64-25);
    uint32_t man = (qiqj_abs>>1) & MASK(23);
    man += qiqj_abs &1;
    if (BIT(man,23)) {
      man &= MASK(23);
      ++ex;
      if (ex>=32) {
	printf("ERROR (%s) : illegal exponent after rounding.\n",__func__);
	ex = 0;
      }
    }
    z= (sign<<28) + (ex<<23) + man;
  }
  return z;
}

/**
   @brief return qiqj in double for debug (Komatsu)
   @param[in] qiqj 29-bit floating (23-bit mantissa, 5-bit exponent, 1-bit sign) bias=20
   @retval qiqj in double
*/
double
mdg_rev_charge_mult(uint32_t qiqj){
  const int mw=23;
  const int ew=5;
  const int bias=20;
  double qiqj_d=TO_DOUBLE_U(qiqj,mw,ew,bias);
  if(BIT(qiqj,mw+ew))qiqj_d*=-1;
  return qiqj_d;
}

/**
   @brief return vdw_common for debug (Komatsu)
   @param[in] uint64_t 25-bit signed mantissa and 8-bit exponent, bias unknown
   @retval double
*/
double
mdg_rev_funcev_vdw_common(uint64_t x) {
  const int mw=25;
  const int ew=8;
  const int bias=20+16+32-1;
  int64_t c64 = x & MASK64(mw);
  c64 = EXTENDLL(c64,mw);
  int cex64 = SLICE64(x,mw,ew);
  double c64d = ldexp((double)c64,cex64-bias-mw+1);
  return c64d;
}

/**
   @brief function evaluator for van der Waals
   @param[in] x r^2/sigma^2 in floating, 25-bit mantissa and 7-bit exponent, bias?
   @param[in] exp_max max exponent
   @param[in] norder order of polynomial
   @param[in] c pointer to table
   @param[in] exponent value
   @param[in] function name
   @retval evaluation result in 25-bit signed mantissa and 8-bit exponent, bias? (not in fixedpoint)
*/
static uint64_t
mdg_funcev_vdw_common(uint32_t x,
		      int norder, int64_t *c, int ex_o, char *name) {
  int32_t dx = (x & MASK(19)) - (1<<18);
  //printf("%08x %08x\n", x, dx);

  // xf is only used for debug
  int ex_in = SLICE(x,25,7);
  double xf = (x & MASK(25)) + (1<<25);
  if (ex_in==0) xf=0.0;
  else xf = scalbn(xf, ex_in-26); // minimum is 1
  
  int64_t z = polynomial_funcev(dx, 19, norder,
				mdg_funcev_vdw_calc_width,
				c, name, xf);
  // normalization and rounding
  int round, stickey, lsb;
  int w=25; // signed mantissa width
  int nr = 28-w-1;
  if ((z>=(1<<28))||(z<-(1<<28))) {
    ++nr;
    ++ex_o;
  }
  round   = BIT(z,nr);
  stickey = (SLICE(z,0,nr)>0) ? 1 : 0;
  lsb     = BIT(z,nr+1);
  z >>= nr+1;
  int carry = round & (stickey | lsb);
  z += carry;

  if ((z>=(1<<(w-1)))||(z<-(1<<(w-1)))) {
    ++ex_o;
    z >>=1;
  }
  uint64_t res = (((uint64_t)ex_o & MASK64(8))<<w) + (z & MASK64(w));
  return res;
}

/**
   @brief function evaluator for van der Waals
   @param[in] x r^2/sigma^2 in floating, 25-bit mantissa and 7-bit exponent
   @param[in] exp_max max input exponent
   @param[in] underflow_val return value when input range is out of definition range
   @param[in] t pointer to table
   @retval evaluation result in floating point, 24bit mantissa unnormalized, 8bit exponent
	Inconsistent with mdg_funcev_vdw_common retval (25-bit signed mantissa, 8bit exponent, bias)? (Komatsu)
*/
uint64_t
mdg_funcev_vdw(uint32_t x, int exp_max, uint32_t underflow_val, mdg_func_table *t) {
  int ex_in = SLICE(x,25,7);
  int exp_diff = exp_max-ex_in;
  int overflow = exp_diff<0;
  int underflow = exp_diff>=8;

  if (overflow) return 0;
  if (underflow) return underflow_val;

  int ex_adr = 7- (exp_diff & MASK(3));
  int adr = (ex_adr << 6) + SLICE(x,19,6);
  int64_t *c=t->ci+adr*(t->norder+1);
  int ex_o = t->exp[adr];
  
  return mdg_funcev_vdw_common(x, t->norder, c, ex_o, t->name);
}

/**
   @brief function evaluator for van der Waals, using array 
   @param[in] x r^2/sigma^2 in floating, 25-bit mantissa and 7-bit exponent
   @param[in] exp_max max exponent
   @param[in] underflow_val return value when input range is out of definition range
   @param[in] vdW coefficient table
   @param[in] vdW exponent table
   @param[in] function name for debug
   @retval evaluation result in floating point, 24bit mantissa unnormalized, 8bit exponent, bias?
	Inconsistent with mdg_funcev_vdw_common retval (25-bit signed mantissa, 8bit exponent, bias)? (Komatsu)
*/
uint64_t
mdg_funcev_vdw_array(uint32_t x, int exp_max, uint32_t underflow_val, 
		     int64_t vdw_table[512][4], int vdw_exp_table[512],
		     char *name) {
  int ex_in = SLICE(x,25,7);
  int exp_diff = exp_max-ex_in;
  int overflow = exp_diff<0;
  int underflow = exp_diff>=8;

  if (overflow) return 0;
  if (underflow) return underflow_val;

  int ex_adr = 7- (exp_diff & MASK(3));
  int adr = (ex_adr << 6) + SLICE(x,19,6);
  int64_t *c=vdw_table[adr];
  int ex_o = vdw_exp_table[adr];

  return mdg_funcev_vdw_common(x, 3, c, ex_o, name);
}

/**
   @brief vdw epsilon multiplier
   @param[in] ei sqrt(eps_i) in floating, 23-bit mantissa and 5-bit exponent
   @param[in] ej the same as ei
   @retval evaluation result in floating, 23-bit mantissa and 6-bit exponent
*/
uint32_t
mdg_vdw_eps_mult(uint32_t ei, uint32_t ej) {
  uint32_t z = mdg_fmul((uint64_t)ei, (uint64_t)ej, 23, 5, 23, 5, 23, 6);

  return z;
}

/**
   @brief return coulomb (qiqj x gbc) x rinv in double for debug (Komatsu)
   @param[in] product in floating, 25-bit signed mantissa and 8-bit exponent,  bias=20+16+32-1
   @retval double
*/
double
mdg_rev_coulomb_prefactor(uint64_t x){
  const int mw=25;
  const int ew=8;
  const int bias=20+16+32-1;

  int32_t c32 = x & MASK(mw);
  c32 = EXTEND(c32,mw);
  int cex32 = SLICE(x,mw,ew);
  double c32d = ldexp((double)c32,cex32-bias-mw+1);

/*
	int64_t c64 = x & MASK64(mw);
	c64 = EXTENDLL(c64,mw);
	int cex64 = SLICE64(x,mw,ew);
	double c64d = ldexp((double)c64,cex64-bias-mw+1);

#define TO_DOUBLE_S32(x,m,n,bias)  (ldexp((double)(EXTEND((int32_t)((x)&MASK(m)),m)), SLICE(x,m,n)-(m)+1-(bias)));
#define TO_DOUBLE_S64(x,m,n,bias)  (ldexp((long double)(EXTENDLL((int64_t)((x)&MASK64(m)),m)), (int)(SLICE64(x,m,n))-(m)+1-(bias)));

	double z = TO_DOUBLE_S(x,mw,ew,bias);
	double z32 = TO_DOUBLE_S32(x,mw,ew,bias);
	double z64 = TO_DOUBLE_S64(x,mw,ew,bias);

#undef TO_DOUBLE_S32
#undef TO_DOUBLE_S64

	printf("%s: x=%lx z=%g z32=%g z64=%g c32d=%g c64d=%g\n",__func__,x,z,z32,z64,c32d,c64d); //in case negative, bug in TO_DOUBLE_S?
*/

  return c32d;
}
/**
   @brief coulomb (qiqj x gbc) x rinv paranoia_check for debug (Komatsu)
   @param[in] qiqj_gbc qiqjxgbc in floating, 23-bit mantissa and 6-bit exponent, 1-bit sign, bias=20+16
   @param[in] rinv rinv/rinv^3  in floating, 23-bit mantissa and 7-bit exponent, no sign, bias=32  (correct?)
   @retval product in floating, 25-bit signed mantissa and 8-bit exponent, bias=20+16+32-1?
*/
//  const int xmw = 23;	//24bit effective
//  const int ymw = 23;
//  const int zmw = 25;	//equivalent to 25bit integer = 24bit effective (1bit sign)
//  const int xew = 6;
//  const int yew = 7;
//  const int zew = 8;
//  const int expbiasdiff = 0;
uint64_t
mdg_coulomb_prefactor(uint32_t qiqj_gbc, uint32_t rinv) {
  uint64_t xman = (qiqj_gbc & MASK(23)) + (1<<23);
  uint64_t yman = (rinv & MASK(23)) + (1<<23);

  int xexp = SLICE(qiqj_gbc, 23, 6);
  int yexp = SLICE(rinv, 23, 7);

  int sign = BIT(qiqj_gbc, 29);
  
  uint64_t zman = xman * yman;
  int gt2 = zman >> (23+23+2-1); // greater than 2 (or equal)

  if(B_MDG_CHECK){
    double qiqj_gbc_double=mdg_rev_coulomb_group_cutoff(qiqj_gbc);
    double rinv_double=mdg_rev_inv_sqrt(rinv);
    printf("%s: qiqj_gbc=%g rinv(/rinv3)=%g\n",__func__,qiqj_gbc_double,rinv_double);
    printf("%s: zman=%lx(%ld) xman=%lx(%ld) yman=%lx(%ld)\n",__func__,zman,zman,xman,xman,yman,yman);
  }

  if (gt2==0) zman <<= 1;
  int64_t zman_sft = zman >> (23+23+2-23-1);

  if(B_MDG_CHECK)printf("%s: zman_sft=%lx\n",__func__,zman_sft);

  int round = BIT(zman, 23+23+2-23-2);
  zman_sft += round;
  int shift_after_round = 0;

  if ((sign==0) && round && BIT(zman_sft,25-1)) { //BIT(zman_sft,zmw-1) //bug fix (Komatsu)
    if(B_MDG_CHECK)printf("%s: zman_sft overflow adjusted\n",__func__);
    zman_sft >>= 1; shift_after_round++;
  }

  if(B_MDG_CHECK)if(round)printf("%s: (+round) zman_sft=%lx\n",__func__,zman_sft);

  //if (sign) shift_after_round = BIT(zman_sft, 23+1); //comment out (Komatsu)
  if (sign) zman_sft = -zman_sft;
  zman_sft &=  MASK64(25);

  if (shift_after_round && gt2) {
    printf("WARNING(%s) : result becomes >= 4\n", __func__);
  }
  uint64_t zexp = 0;
  if (xexp && yexp) {
    zexp = xexp + yexp + gt2 + shift_after_round;
  } else {
    zexp = 0;
    zman_sft = 0;
  }

  if (zexp >= (1<<8)) {
    printf("ERROR(%s) : exponent overflow\n", __func__);
  }
  zexp &= MASK64(8);
  uint64_t z = ( zexp << 25 ) + ( zman_sft & MASK(25) );

  if(B_MDG_CHECK)printf("%s: retval=%g\n",__func__,mdg_rev_coulomb_prefactor(z));

  return z;
}

/**
   @brief van der Waals combination rule for sigma - si+sj
   @param[in] si in fixed, 28-bit, no sign
   @param[in] sj in fixed, 28-bit, no sign
   @retval si+sj in floating, 23-bit mantissa and 5-bit exponent, no sign
*/
uint32_t
mdg_vdw_sigma_sum(uint32_t si, uint32_t sj) {
  si &= MASK(28);
  sj &= MASK(28);
  uint64_t sum = si+sj; // 29 bit

  int lop = ilogbi(sum); // 32..1,  0 for 0
  
  if (lop==0) return 0;
  
  sum <<= 29-lop; 
  sum = round_unbiased_u64(sum, 5); // 24b (except for carry after rounding)

  lop += 1; // max 30
  if (BIT(sum,24)) ++lop; 

  uint32_t z = (lop<<23) + (sum & MASK(23));
  return z;
}

#include "inv_table.h"
/**
   @brief van der Waals combination rule for sigma - inverse si+sj
   @param[in] sij  in floating, 23-bit mantissa and 5-bit exponent, no sign
   @retval 1.0/sij in floating, 23-bit mantissa and 5-bit exponent, no sign
*/
uint32_t
mdg_vdw_sigma_inv(uint32_t sij) {
  int ex_in = SLICE(sij,23,5);
  
  if (ex_in==0) return 0;
  
  uint64_t man = SLICE(sij,0,23)+(1<<23);

  uint64_t a = (MASK(7)<<17) & man; // 24b 
  uint64_t b = MASK(17) & man;
  
  int adr = (a>>17) & MASK(6);
  uint64_t a2inv = mdg4_inv_table[adr]; // 26b, 2^25=1
  uint64_t prod = a2inv * (a-b); // 50b
  prod = round_biased_u64(prod, 24); // 26b (but MSB should be always 0)
  //printf("adr=%d a=%llx b=%llx 1/a2=%llx a-b=%llx (a-b)/a^2=%llx ", adr, a, b, a2inv, a-b, prod);

  // b: 2^23
  uint64_t b2 = (b>>3) * (b>>3);  // b2: 2^40
  b2 = round_biased_u64(b2, 15);  // b2: 2^25
  b2 *= a2inv>>12;                // b2: 2^38
  b2 = round_biased_u64(b2, 13);  // b2: 2^25
  b2 += 1<<25;

  prod *= b2; //51b
  prod = round_biased_u64(prod, 25); //26b
  //printf("b2=%llx res=%llx\n", b2, prod);

  int ex = 31-ex_in;
  if (BIT(prod,24)) {
    ++ex;
    prod=0; // Correspond to 1
  }

  uint32_t z = (ex<<23) + (prod & MASK(23));
  return z;
}

/**
   @brief return vdw inverse sigma in double for debug (Komatsu)
   @param[in] product in floating, 23-bit mantissa and 6-bit exponent, no sign, bias?
   @return double
*/
double
mdg_rev_vdw_sigma_mult(uint32_t x) {
  const int mw=23;
  const int ew=6;
  const int bias=14;
  double z = TO_DOUBLE_U(x,mw,ew,bias);
  return z;
}

/**
   @brief vdw inverse sigma (geometric mode)
   @param[in] 1/si in floating, 23-bit mantissa and 5-bit exponent, no sign
   @param[in] 1/sj in floating, 23-bit mantissa and 5-bit exponent, no sign
   @retval product in floating, 23-bit mantissa and 6-bit exponent, no sign
*/
uint32_t
mdg_vdw_sigma_mult(uint32_t si, uint32_t sj) {
  uint32_t z = mdg_fmul((uint64_t)si, (uint64_t)sj, 23, 5, 23, 5, 23, 6);
  return z;
}

/**
   @brief vdw eij / sij^2
   @param[in] eij in floating, 23-bit mantissa and 6-bit exponent, no sign
   @param[in] 1/sij^2 in floating, 23-bit mantissa and 6-bit exponent, no sign
   @param[in] pot_mode skip multiplication of 1/sij^2 for potential (add 64 to exponent)
   @retval product in floating, 23-bit mantissa and 7-bit exponent, no sign
*/
uint32_t
mdg_vdw_eij_div_sij(uint32_t eij, uint32_t invsij, int pot_mode) {
  if (pot_mode) invsij = 14<<23; // exp bias set to 18
  uint32_t z = mdg_fmul((uint64_t)eij, (uint64_t)invsij, 23, 6, 23, 6, 23, 7);
  return z;
}

/**
   @brief return rij_div_sij for debug (Komatsu)
   @param[in] floating, 25-bit mantissa and 7-bit exponent, no sign, bias?
   @retval double
*/
double
mdg_rev_vdw_rij_div_sij(uint32_t x) {
  const int mw=25;
  const int ew=7;
  const int bias=62+14;
  double z = TO_DOUBLE_U(x,mw,ew,bias);
  return z;
}

/**
   @brief vdw rij^2 / sij^2
   @param[in] rij^2 in floating, 25-bit mantissa and 6-bit exponent, no sign
   @param[in] 1/sij^2 in floating, 23-bit mantissa and 6-bit exponent, no sign
   @retval product in floating, 25-bit mantissa and 7-bit exponent, no sign
*/
uint32_t
mdg_vdw_rij_div_sij(uint32_t rij, uint32_t invsij) {
  uint32_t z = mdg_fmul((uint64_t)rij, (uint64_t)invsij, 25, 6, 23, 6, 25, 7);
  return z;
}

/**
   @brief vdw group cutoff
   @param[in] eij in floating, 23-bit mantissa and 7-bit exponent, no sign
   @param[in] group cutoff in floating, 23-bit mantissa and 7-bit exponent, no sign
   @retval product in floating, 23-bit mantissa and 7-bit exponent, no sign
*/
uint32_t
mdg_vdw_group_cutoff(uint32_t eij, uint32_t gbc) {
  uint32_t z = mdg_fmul((uint64_t)eij, (uint64_t)gbc, 23, 7, 23, 7, 23, 8);
  return z;
}

/**
   @brief coulomb group cutoff
   @param[in] qiqj in floating, 23-bit mantissa and 5-bit exponent, 1 sign
   @param[in] group cutoff in floating, 23-bit mantissa and 5-bit exponent, 1 sign
   @retval product in floating, 23-bit mantissa and 6-bit exponent, 1 sign
*/
uint32_t
mdg_coulomb_group_cutoff(uint32_t x, uint32_t y) {
  uint32_t z = mdg_fmul((uint64_t)x, (uint64_t)y, 23, 5, 23, 5, 23, 6);
  return z;
}
/**
   @brief return coulomb group cutoff in double for debug (Komatsu)
   @param[in] product in floating, 23-bit mantissa and 6-bit exponent, 1 sign, bias=20+16
   @retval double value 
*/
double
mdg_rev_coulomb_group_cutoff(uint32_t x){
  const int mw=23;
  const int ew=6;
  const int bias=20+16;
  double z = TO_DOUBLE_U(x,mw,ew,bias);
  if(BIT(x,mw+ew))z*=-1;
  return z;
}


/**
   @brief vdw factor in total
   @param[in] eij/sij in floating, 23-bit mantissa and 8-bit exponent, no sign
   @param[in] result of function evaluator (sij^-7 - sij^-4 etc.) in floating, 25-bit signed mantissa and 8-bit exponent
   @param[in] normalization - max exponent, 9 bit
   @param[in] 1 for rounding at 64bit, 0 for 32bit
   @param[out] result van der Waals kernel 32/64bit
   @param[out] overflow flag
*/
void
mdg_vdw_kernel(uint32_t eij_sij, uint64_t res_fev_vdw,
	       int norm_vdw, int round_64,
	       int64_t *z, int *zovf) {
  int64_t m1 = MANTISSA(eij_sij,23);
  int e1 = SLICE(eij_sij,23,8);

  int64_t m2 = EXTENDLL(res_fev_vdw & MASK64(25), 25);
  int e2 = SLICE(res_fev_vdw, 25,8);

  int zero = (e1 == 0) || (e2 == 0);

  // max: (2^24-1) * (2^23-1) <  2^47
  // min: -(2^24-1) * 2^23    > -2^47
  int64_t m12 = m1*m2;
  const int64_t msb = 1LL<<48;
  if ((m12>=msb)||(m12<=-msb)) {
    printf("ERROR(%s): multiplied value %" PRIx64 " out of range - should never happen\n",
	   __func__, m12);
  }
  int e12 = e1+e2-norm_vdw;
  //printf("vkern m=%llx e1=%d e2=%d e=%d\n", m12, e1, e2, e12);
#ifdef DEBUG_VERBOSE
  printf("%llx & %llx = %llx\n",res_fev_vdw,MASK64(25),res_fev_vdw & MASK64(25));
  printf("((x)&(1LL<<((n)-1))) %llx\n",((res_fev_vdw & MASK64(25))&(1LL<<((25)-1))));
  printf("(((x)&MASK64(n))-(1LL<<(n))) = (((%llx)& %llx)-(1LL<<(25))) = ((%llx)-(%llx)) %llx\n", res_fev_vdw & MASK64(25),MASK64(25),(res_fev_vdw & MASK64(25))&MASK64(25), 1LL<<(25),(((res_fev_vdw & MASK64(25))&MASK64(25) )-(1LL<<(25))));
#endif
  int64_t z0;
  int ovf=0;
  int64_t round=0;
  if (zero ||(m12==0)) {
    z0=0; ovf=0;
  } else if (e12>=49) { // overflow
    z0=0; ovf=1;
  } else if (e12>0) { // shiftup
    __int128_t m128 = m12;
    m128 <<= e12; // less than 49
    if ((m128>=msb)||(m128<-msb)) {
      z0=0; ovf=1;
    } else {
      z0=m128<<15; ovf=0;
    }
  } else { // shiftdown
    ovf=0;
    if (-e12 > 64) {
      z0 = 0; round=0;
    } else {
      __int128_t m128 = m12;
      m128 <<= 16; // for rounding 
      m128 >>= -e12;
      // rounding 
      if ( round_64 && (m128 & 1)) {
	round = 1;
      }
      z0 = m128>>1;
    }
  }
  if ((!round_64) && BIT(z0, 31) ) {
    round = 1LL<<32;
  }
  if (round) {
    int sign = (z0<0) ? 1 : 0; 
    z0 += round;
    if ((sign==0)&&(z0<0)) {
      ovf = 1;
      z0 = 0;
    }
    //printf("%llx %d\n", z0, round);
  }
  // if (!round_64) z0 &= 0xffffffff00000000;
  *z = z0;
  *zovf=ovf;
}

/**
   @brief coulomb factor in total
   @param[in] coulomb_core in floating, 25-bit signed mantissa and 8-bit exponent
   @param[in] fev_excluded in fixed, 28-bit signed  
   @param[in] normalization - max exponent, 8 bit
   @param[in] 1 for rounding at 64bit, 0 for 32bit
   @param[out] result coulomb kernel 32/64bit
   @param[out] overflow flag
*/
void
mdg_coulomb_kernel(uint64_t coulomb_core, int32_t fev_excluded,
		   int norm_coulomb, int round_64,
		   int64_t *z, int *zovf) {

  int64_t m1 = EXTENDLL(coulomb_core & MASK64(25), 25);
  int e1 = SLICE(coulomb_core,25,8);

  int64_t m2 = EXTENDLL( fev_excluded & MASK64(28), 28);

  int zero = (e1 == 0) || (m2 == 0);

  // max: -2^24 * -2^27 = 2^51
  // min: -2^24 * (2^27-1) > -2^51
  int64_t m12 = m1*m2;
  const int64_t msb = 1LL<<51;
  if ((m12>msb)||(m12<=-msb)) {
    printf("ERROR(%s): multiplied value %" PRIx64 " out of range - should never happen\n",
	   __func__, m12);
  }
  // Since max multiplied result is 2^51, result is expressed by 53 bits.
  const int64_t msb2 = 1LL<<52;
  int e12 = e1-norm_coulomb;
  //printf("%llx %d\n", m12, e12);
  int64_t z0;
  int ovf;
  int64_t round=0;
  if (zero ||(m12==0)) {
    z0=0; ovf=0;
  } else if (e12>=52) { // overflow
    z0=0; ovf=1;
  } else if (e12>0) { // shiftup
    __int128_t m128 = m12;
    m128 <<= e12; // less than 52
    //printf("shiftup %llx %lx\n", m128, msb2);
    if ((m128>=msb2)||(m128<-msb2)) {
      z0=0; ovf=1;
    } else {
      z0=m128<<11; ovf=0;
    }
  } else { // shiftdown
    ovf=0;
    if (-e12 > 64) {
      z0 = 0; round=0;
    } else {
      __int128_t m128 = m12;
      m128 <<= 12; // for rounding 
      m128 >>= -e12;
      // rounding 
      if ( round_64 && (m128 & 1)) {
	round = 1;
      }
      z0 = m128>>1;
      //printf("shiftdown %llx %d\n", m128, e12);
    }
  }
  if ((!round_64) && BIT(z0, 31) ) {
    round = 1LL<<32;
  }
  if (round) {
    int sign = (z0<0) ? 1 : 0; 
    z0 += round;
    if ((sign==0)&&(z0<0)) {
      ovf = 1;
      z0 = 0;
    }
    //printf("%llx %d\n", z0, round);
  }
  // if (!round_64) z0 &= 0xffffffff00000000;
  *z = z0;
  *zovf=ovf;
}

/**
   @brief force factor in total
   @param[in] coulomb factor, 32-bit fixed
   @param[in] coulomb factor overflow flag
   @param[in] vdw factor, 32-bit fixed
   @param[in] vdw factor overflow flag
   @param[in] negate result
   @param[out] force factor, 33-bit fixed
   @param[out] force factor overflow flag
*/
void
mdg_force_factor(int64_t ckernel, int ckernel_overflow,
		 int64_t vkernel, int vkernel_overflow,
		 int negate,
		 int64_t *fkernel, int *fkernel_overflow) {

  //printf("ck %lx vk %lx", ckernel, vkernel);
  int ovf = ckernel_overflow || vkernel_overflow;
  int64_t sign_ovf=0;
  if (vkernel_overflow) sign_ovf = BIT(vkernel,31);
  else if (ckernel_overflow) sign_ovf = BIT(ckernel,31);
  ckernel = EXTENDLL(ckernel, 32);
  vkernel = EXTENDLL(vkernel, 32);
  //printf(" ck %lx vk %lx\n", ckernel, vkernel);
  int64_t sum = ckernel + vkernel;
  if (negate) sum = -sum;
  if (sum>=0x100000000L) {
    ovf = 1;
    sum &= MASK64(32);
  } else if (sum<-0x100000000L) {
    printf("ERROR(%s): something wrong; improbable overflow\n", __func__);
  } else if (ovf) {
    sum &= MASK64(32);
    sum += sign_ovf <<32;
    sum = EXTENDLL(sum, 33);
  }
  *fkernel = sum;
  *fkernel_overflow = ovf;
}

/**
   @brief force multilication
   @param[in] relative distance, 29-bit fixed
   @param[in] force factor, 33-bit fixed
   @param[in] force factor overflow flag
   @param[out] force, 32-bit fixed
   @param[out] force overflow
*/
void
mdg_force_mul(int32_t r, int64_t fkernel, int fkernel_overflow,
	      int32_t *f, int *f_ovf) {
  //int sign = BIT(r,28) ^ BIT(fkernel, 32);
  int64_t r64 = EXTENDLL((int64_t)r, 29);
  fkernel = EXTENDLL(fkernel, 33);
  int64_t prod = r64*fkernel;
  int stickey = prod & MASK64(28);
  int round = BIT(prod,28) && (BIT(prod,29) || stickey);
  prod >>= 29;
  if (round) ++prod;
  if ((prod>(1L<<31))||(prod<-(1L<<31))) {
    printf("ERROR(%s): something wrong.\n", __func__);
  }
  int ovf = fkernel_overflow || (prod>=(1L<<31)) || (prod<-(1L<<31));
  //*f = (sign<<31) + (prod & MASK(31));
  *f = prod & MASK64(32);
  *f_ovf = ovf;
}

/**
   @brief force accumulation
   @param[in] force vector, 32-bit fixed X 3
   @param[in] force overflow flag, 3 bit
   @param[out] force i, 32-bit fixed
   @param[out] force j, 32-bit fixed
   @param[out] overflow flag
*/
void
mdg_accum
(
 int64_t *fin,
 int fin_overflow,
 int64_t *fi,
 int64_t *fj,
 int *ovf
 )
{
  const int64_t max=1<<31;
  *ovf = fin_overflow;
  if (fin_overflow==0) {
    int64_t tmpi[3], tmpj[3];
    for(int i=0;i<3;++i) {
#ifdef DEBUG_VERBOSE
      int64_t tmp = EXTENDLL(fin[i],32);
#endif
      tmpi[i] = fi[i]+fin[i];
      tmpj[i] = fj[i]-fin[i];
      if ( (tmpi[i]>=max)||(tmpi[i]<-max)||
	   (tmpj[i]>=max)||(tmpj[i]<-max) ) *ovf=1;
    }
    if (*ovf==0) {
      for(int i=0;i<3;++i) {
	fi[i] = tmpi[i];
	fj[i] = tmpj[i];
      }
    }
  }
}

/**
   @brief coulomb exclusion
   @param[in] result of function evaluator, 26-bit fixed unsigned
   @param[in] noexcl_factor, 27-bit fixed
   @param[in] exclfull_factor, 27-bit fixed
   @param[in] excl14_factor, 27-bit fixed
   @param[in] full exclusion
   @param[in] 1-4 exclusion
   @retval excluded coulomb factor, 28-bit signed fixed
*/
int32_t
mdg_coulomb_exclusion
(
 int32_t res_fev_coul,
 int32_t noexcl_factor,
 int32_t exclfull_factor,
 int32_t excl14_factor,
 int excl_full,
 int excl_14
 )
{
  res_fev_coul &= MASK(26);

  int32_t excl_factor;
  if (excl_full)    excl_factor=exclfull_factor;
  else if (excl_14) excl_factor=excl14_factor;
  else              excl_factor=noexcl_factor;

  //printf("excl_factor=%x ",excl_factor);
  excl_factor = EXTEND(excl_factor,27);
  //printf("%x %x %x %x\n",excl_factor, noexcl_factor, excl_factor, excl14_factor);

  return res_fev_coul + excl_factor;
}

/**
   @brief r0^2 - r^2 for coulomb softcore potential
   @param[in] r^2, fixed 32bit, max=1 (lsb is used for rounding)
   @param[in] k: r0^2 + 2r0 phi(r0)/f(r0) 26bit + exp 3bit (exp==7 -> max=1)
   @retval r0^2 + 2r0 phi(r0)/f(r0) - r^2, 26bit
   // shift = 0 case -> minimum k, lsb = 1 x 2^-128 = 1/128, r0 ~ 1/30
*/
uint32_t
mdg_coulomb_softcore
(
 uint64_t r2,
 uint32_t k0,
 int potential_mode
 )
{
  // ex = 7: BP of k0 = 24, max k0 = 4 (-eps)
  // r2 max = 1, BP at 32
  //  ex = 7 -> k0 = 2 ~ 4 shift 8 
  //  ex = 6 -> k0 = 1 ~ 2 shift 7
  //  ex = 5 -> k0 = 0.5 ~ 1 shift 6
  //  ex = 4 -> k0 = 0.25 ~ 0.5 shift 5
  //  ex = 3 -> k0 = 0.125 ~ 0.25 shift 4
  //  shift amount -1 to incl. round bit
  if (potential_mode==0) return (1<<25); // This correspond to 1/2
  int ex = k0>>26;
  uint32_t k = k0 & MASK(26);
  // rounding for minus
  uint64_t x = r2>>ex;
  //printf("k0 = %x k = %x shift=%d r2=%lx r2s=%lx\n", k0, k, ex, r2, x);
  x=(x>>1)+(x&1);
  k -= x;
  return k & MASK(26);
}

/**
   @brief r0^2 - (r/s)^2 for van der Waals softcore potential
   @param[in] (r/s)^2, floating 25-bit mantissa, 7-bit exponent
   @param[in] k: r0^2 + 2r0 phi(r0)/f(r0), floating 25-bit mantissa, 7-bit exponent
   @retval r0^2 + 2r0 phi(r0)/f(r0) - r^2, 25bit signed mantissa, 8-bit exponent
   exponent bias: 76
*/
uint64_t
mdg_vdw_softcore
(
 uint32_t rs,
 uint32_t k0,
 int potential_mode
 )
{
  uint64_t bias = 76;
  if (potential_mode==0) return (bias<<25) + 0x800000; // This correspond to 1/2

  int32_t mank = MANTISSA(k0,25);
  int32_t manr = MANTISSA(rs,25);
  int expk = k0>>25;
  int expr = rs>>25;
  int ex_diff = expk - expr;
  // ex_diff should be positive if rs < rswitch

  // result should be positive
  int32_t man = mank - (manr>>ex_diff);
  if (man<0) {
    printf("ERROR(%s): negative result\n", __func__);
  }
  if (man>=(1<<25)) {
    int round = BIT(man,1) & (BIT(man,2) | BIT(man,0));
    man = (man>>2) + round;
    expk += 2;
  } else if (man>=(1<<24)) {
    int round = BIT(man,0) & BIT(man,1);
    man = (man>>1) + round;
    expk += 1;
  }
  if (man>=(1<<24)) {
    man >>=1;
    expk += 1;
  }
  if (expk>=256) {
    printf("ERROR(%s): exponent overflow %d\n", __func__, expk);
  }
  return ((uint64_t)expk<<25) + (man & MASK(25));
}

/**
   @brief addition with saturation, +inf + (-inf) = +inf
   @param[in] x 32bit
   @param[in] xovf 1bit
   @param[in] y 32bit
   @param[in] yovf 1bit
   @param[out] z 32bit
   @param[out] zovf 1bit
*/
void
mdg_add32_sat(int32_t x, int xovf, int32_t y, int yovf, int32_t *z, int *zovf) {
  int xsign = (x<0) ? 1 : 0;
  int ysign = (y<0) ? 1 : 0;
  int zsign = 0;
  *zovf=0;
  if (xovf && yovf) {
    zsign = xsign & ysign;
    *zovf = 1;
  } else if (xovf) {
    zsign = xsign;
    *zovf = 1;
  } else if (yovf) {
    zsign = ysign;
    *zovf = 1;
  } else {
    int64_t zl = x;
    zl += (int64_t)y;
    if (zl<-0x80000000L) {
      zsign =1;
      *zovf = 1;
    } else if (zl>=0x80000000L) { 
      zsign =0;
      *zovf = 1;
    } else {
      *z = zl;
    }
  }
  if (*zovf) {
    if (zsign) *z = -0x80000000;
    else *z = 0x7fffffff;
  }
}

/**
   @brief subtraction with saturation, +inf - (+inf)/ -inf - (-inf) = +inf
   @param[in] x 32bit
   @param[in] xovf 1bit
   @param[in] y 32bit
   @param[in] yovf 1bit
   @param[out] z = x-y 32bit
   @param[out] zovf 1bit
*/
void
mdg_sub32_sat(int32_t x, int xovf, int32_t y, int yovf, int32_t *z, int *zovf) {
  int xsign = (x<0) ? 1 : 0;
  int ysign = (y<0) ? 0 : 1;
  int zsign = 0;
  *zovf=0;
  if (xovf && yovf) {
    zsign = xsign & ysign;
    *zovf = 1;
  } else if (xovf) {
    zsign = xsign;
    *zovf = 1;
  } else if (yovf) {
    zsign = ysign;
    *zovf = 1;
  } else {
    int64_t zl = x;
    zl -= (int64_t)y;
    if (zl<-0x80000000L) {
      zsign =1;
      *zovf = 1;
    } else if (zl>=0x80000000L) { 
      zsign =0;
      *zovf = 1;
    } else {
      *z = zl;
    }
  }
  if (*zovf) {
    if (zsign) *z = -0x80000000;
    else *z = 0x7fffffff;
  }
}

/**
   @brief addition with saturation, +inf + (-inf) = +inf
   @param[in] x 64bit
   @param[in] xovf 1bit
   @param[in] y 64bit
   @param[in] yovf 1bit
   @param[out] z 64bit
   @param[out] zovf 1bit
*/
void
mdg_add64_sat(int64_t x, int xovf, int64_t y, int yovf, int64_t *z, int *zovf) {
  int xsign = (x<0) ? 1 : 0;
  int ysign = (y<0) ? 1 : 0;
  int zsign = 0;
  *zovf=0;
  if (xovf && yovf) {
    zsign = xsign & ysign;
    *zovf = 1;
  } else if (xovf) {
    zsign = xsign;
    *zovf = 1;
  } else if (yovf) {
    zsign = ysign;
    *zovf = 1;
  } else {
    int64_t zl = x + y;
    if ((x>0) && (y>0) && (zl<=0)) {
      zsign =0;
      *zovf = 1;
    } else if ((x<0) && (y<0) && (zl>=0)) {
      zsign =1;
      *zovf = 1;
    } else {
      *z = zl;
    }
  }
  if (*zovf) {
    if (zsign) *z = -0x8000000000000000L;
    else *z = 0x7fffffffffffffffL;
  }
}

/**
   @brief addition with saturation, overflow at n (<64) bit, +inf + (-inf) = +inf
   @param[in] x 64bit
   @param[in] xovf 1bit
   @param[in] y 64bit
   @param[in] yovf 1bit
   @param[out] z 64bit
   @param[out] zovf 1bit
*/
void
mdg_add_sat(int n, int64_t x, int xovf, int64_t y, int yovf, int64_t *z, int *zovf) {
  int xsign = (x<0) ? 1 : 0;
  int ysign = (y<0) ? 1 : 0;
  int zsign = 0;
  *zovf=0;
  if (xovf && yovf) {
    zsign = xsign & ysign;
    *zovf = 1;
  } else if (xovf) {
    zsign = xsign;
    *zovf = 1;
  } else if (yovf) {
    zsign = ysign;
    *zovf = 1;
  } else {
    int64_t zl = x + y;
    if (zl>= (1L<<(n-1))) {
      zsign =0;
      *zovf = 1;
    } else if (zl< -(1L<<(n-1))) {
      zsign =1;
      *zovf = 1;
    } else {
      *z = zl;
    }
  }
  if (*zovf) {
    if (zsign) *z = -(1L<<(n-1));
    else *z = MASK(n-1);
  }
}

/**
   @brief rough cutoff r^2 calculation
   @param[in] ri 32bit
   @param[in] rj 32bit
   @param[in] cutoff 6bit
   @param[out] inrange flag
*/
int
mdg_checkrange(int32_t ri[3], int32_t rj[3], int cutoff) {
  if (cutoff==0) return 1;
  int inrange = 1;
  int s = 0;
  for(int i=0; i<3; ++i) {
    int xi = SLICE(ri[i],24,8);
    int xj = SLICE(rj[i],24,8);
    int x = xi-xj;
    if (x<-128) x+=256;
    if (x>=128) x-=256;
    //printf("%d ", x);
    if ((x>=16)||(x<-16)) {
      inrange = 0;
    } else {
      if (x<0) x=-x;
      x *= x;
      x >>= 3;
      s += x;
    }
  }
  if (inrange) {
    if (s >= cutoff) inrange=0;
  }
  //printf("%x %x %d\n", s, cutoff, inrange);
  return inrange;
}

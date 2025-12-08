#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#include "defs.h"
#include "scale.h"

/** @file scale.c */

/*
  What we have to scale appropriately

  length
  charge
  vdw sigma
  vdw energy
  coulomb 1/4 pi epsilon
  group-based cutoff factor
  force
  potential

  For this, we specify
  Length scale
  Maximum charge
  Maximum vdw radius
  Maximum vdw energy

*/  

#define FD stdout

mdg_scale_conversion_constants_t *MDG_SC=0; //set 0 (Komatsu)

/**
   @brief Conversion of coordinates
   @param[in] r   coordinate
   @retval scaled coordinate
*/
// We assume r>=0.
// typical scale (sc->length) will be scaled to 2^27.
uint32_t
mdg_scale_r(double r)
{
  double r0 = ldexp(r/MDG_SC->length, 27);
  double max = ldexp(1.0, 32);
  if (r0 >=max) {
    fprintf(FD, "WARNING(%s) : input coordinate out of range.\n",__func__);
  }
  r0 = fmod(r0, max);
  if (r0<0.0) r0 += max;
  return (uint32_t) round(r0); //add round.(Komatsu)
}

/**
   @brief Reverse conversion of coordinates
   @param[in] s  scaled coordinate
   @retval coordinate in double
*/
double
mdg_rev_scale_r(uint32_t s) //added (Komatsu)
{
  double x = ldexp((double)(s), -27);
  return x * MDG_SC->length;
}
double
mdg_rev_scale_r_ltrmax(uint32_t s,double rmax) //added (Komatsu)
{
  double x = ldexp((double)(s), -27)*MDG_SC->length;
  double R = 32.0*MDG_SC->length;
  while(x>=rmax)x-=R;
  while(x<rmax-R)x+=R;
  return x;
}

/**
   @brief Conversion of cutoff^2 Rcut^2
   @param[in] rc  cutoff (not square rc)
   @retval scaled cutoff
*/
//   exp=0x3f -> r^2 = (1<<55)
//   exp=0x40 correspond to (2^28)^2 = (2*sc->length)^2
uint32_t
mdg_scale_rcut2(double rc)
{
  if (rc<=0.0) {
    fprintf(FD, "ERROR(%s): cutoff radius is negative or zero (%e).\n",
	    __func__, rc);
    return 0;
  }
  double x = rc/MDG_SC->length;
  x *= x;
  int ex = ilogb(x);
  uint32_t r_man = round(ldexp(x, -ex + 25)); //add round.(Komatsu) 
  // ex==2 -> exp=0x40
  int exp = 0x40 + ex - 2;
  if (exp<=0) {
    fprintf(FD, "WARNING(%s): cutoff radius %e is too small - set to zero.\n",
	    __func__, rc);
    return 0;
  }
  if (exp>=0x40) {
    fprintf(FD, "WARNING(%s): cutoff radius %e is too large - set to max.\n",
	    __func__, rc);
    return MASK31;
  }
  return (exp << 25) + ( r_man & MASK(25) );
} 

/**
   @brief Conversion of cutoff^2 Rcut^2
   @param[in] rc  cutoff (not square rc)
   @retval scaled rough cutoff
*/
uint32_t
mdg_scale_rough_rcut2(double rc)
{
  if (rc<=0.0) {
    fprintf(FD, "ERROR(%s): cutoff radius is negative or zero (%e).\n",
	    __func__, rc);
    return 0;
  }
  double x = rc/MDG_SC->length;
  x *= x;
  uint32_t cut = ceil(ldexp(x, 3));
  if (cut==0) {
    fprintf(FD, "WARNING(%s): cutoff radius %e is too small - set to 1.\n",
	    __func__, rc);
    return 1;
  }
  if (cut>=0x40) {
    fprintf(FD, "WARNING(%s): cutoff radius %e is too large - disable cutoff.\n",
	    __func__, rc);
    return 0;
  }
  return cut;
} 

/**
   @brief Conversion of charges
   @param[in] q charge
   @retval scaled charge, fixed point
*/
// Conversion of charges
// typical charge will be scaled to 1<<26
int32_t
mdg_scale_charge(double q)
{
  double q0 = ldexp(q/MDG_SC->charge, 26);
  int64_t qi = lround(q0);
  if ((qi>=0x80000000L)||(qi<=-0x80000000L)) {
    fprintf(FD, "ERROR(%s): charge %e exceeds max charge %f.\n",
	    __func__, q, MDG_SC->charge*(1<<(31-26)));
    if (q>0.0) return 0x7fffffff;
    else return 0x80000000;
  }
  return (int32_t)qi;	//OK without round.(Komatsu)
}

/**
   @brief Reverse conversion of charges
   @param[in] s scaled charge in fixed point
   @retval charge in double
*/
double
mdg_rev_scale_charge(int32_t s) //added (Komatsu)
{
  double x = ldexp((double)(s), -26);
  return x * MDG_SC->charge;
}

/**
   @brief Conversion of vdw radius for arithmetic mean
   @param[in] s: radius
   @retval scaled radius sigma, fixed point
*/
// 
uint32_t
mdg_scale_vdw_radius_arithmetic(double s)
{
  if (s<=0.0) {
    fprintf(FD, "ERROR(%s): sigma=%e <=0.0 is not allowed\n",
	    __func__, s);
    return 0;
  }
  s *= 0.5;
  double s0 = ldexp(s/MDG_SC->length, 23);
  uint32_t si = round(s0);	//add round.(Komatsu)
  if ( (s0>=ldexp(1.0, 28)) || (si>=(1<<28)) ) {
    fprintf(FD, "ERROR(%s): sigma=%e overflow.\n",
	    __func__, s);
    return MASK(28);
  }
  return si;
}

/**
   @brief Reverse conversion of vdw radius for arithmetic mean
   @param[in] s: radius
   @retval radius sigma,
*/
double
mdg_rev_scale_vdw_radius_arithmetic(uint32_t s)
{
  double x = ldexp(MASK(28)&s, -23);
  return x * 2.0 * MDG_SC->length;
}

#define VDW_R_GEO_BIAS 7

/**
   @brief Conversion of vdw sigma for geometric mean
   @param[in] sigma
   @retval scaled 1.0/sigma, floating point 23-bit mantissa, 5-bit exponent.
   exponent bias = 7 (modified from 9 by Komatsu)
   inverse sigma :  2^-8 /L <= 1/s < 2^22 /L
   range of sigma : typical length (SC.length) L -> 2^8 L ~ 2^-22 L
*/
uint32_t
mdg_scale_vdw_radius_geometric(double s)
{
  if (s<=0.0) {
    fprintf(FD, "ERROR(%s): sigma=%e <=0.0 is not allowed\n",
	    __func__, s);
    return 0;
  }
  double x = MDG_SC->length/s;
  int ex = ilogb(x);
  uint32_t s_man = round(ldexp(x, -ex + 23));	//add round.(Komatsu) 
  ex += VDW_R_GEO_BIAS; //9 (modified by Komatsu) 
  if (ex>31) {
    fprintf(FD, "ERROR(%s): sigma=%e too small - set to max.\n",
	    __func__, s);
    return MASK(28);
  }
  if (ex<=0) {
    fprintf(FD, "WARNING(%s): sigma=%e too large - set 1/s to zero.\n",
	    __func__, s);
    return 0;
  }
  uint32_t sinv = (ex<<23) + (MASK(23) & s_man);
  return sinv;
}

/**
   @brief Reverse conversion of vdw sigma for geometric mean
   @param[in] 1.0/sigma
   @retval sigma
*/
double
mdg_rev_scale_vdw_radius_geometric(uint32_t s)
{
  double sinv = TO_DOUBLE_U(s, 23, 5, VDW_R_GEO_BIAS); // bias 7 follow changeset:   104:c8f0bf8d1f1e
  return MDG_SC->length/sinv;
}

#define VDW_ENERGY_BIAS 27
/**
   @brief Conversion of vdw energy epsilon 
   @param[in] eps energy
   @retval scaled sqrt epsilon, floating point 23-bit mantissa, 5-bit exponent.
   range of sqrt epsilon : 2^5 > eps >= 2^-26 (kJ/mol)^1/2
   exponent bias = 27
*/
uint32_t
mdg_scale_vdw_energy(double eps)
{
  if (eps<0.0) {
    fprintf(FD, "ERROR(%s): negative eps=%e is not allowed, set to zero.\n",
	    __func__, eps);
    return 0;
  }
  if (eps==0.0) return 0;
  double x = sqrt(eps);
  int ex=ilogb(x);
  uint32_t z = round(ldexp(x, 23-ex));	//add round.(Komatsu)
  ex+=VDW_ENERGY_BIAS;
  if (ex>31) {
    fprintf(FD, "ERROR(%s): eps=%e (sqrt(eps)=%e) overflow.  Set to max.\n",
	    __func__, eps, x);
    return MASK(28);
  }
  if (ex<=0) {
    fprintf(FD, "WARNING(%s): eps=%e (sqrt(eps)=%e) underflow. Set to zero.\n",
	    __func__, eps, x);
    return 0;
  }
  z = (ex<<23) + (MASK(23) & z);
  return z;
}

/**
   @brief Reverse conversion of vdw energy epsilon 
   @param[in] eps energy
   @retval recovered sqrt epsilon
*/
double
mdg_rev_scale_vdw_energy(uint32_t eps)
{
  double e=TO_DOUBLE_U(eps, 23, 5, VDW_ENERGY_BIAS);
  return e*e;
}

/**
   @brief Reverse conversion of vdw group-based cutoff factor
   @param[in] scaled cutoff factor, floating point 23-bit mantissa, 7-bit exponent.
   exponent bias : 32
   @retval cutoff factor in double (divided by MDG_SC->vdw_energy_factor)
*///added (Komatsu)
double
mdg_rev_scale_vdw_factor(uint32_t s)
{
  double e=TO_DOUBLE_U(s, 23, 7, 32);
  return e/MDG_SC->vdw_energy_factor;
}

/**
   @brief Conversion of vdw group-based cutoff factor (also used for 1-4 exclusion)
   @param[in] factor cutoff factor, range 2^32 > factor >= 2^-31
   @retval scaled cutoff factor, floating point 23-bit mantissa, 7-bit exponent.
   exponent bias : 32
*/
uint32_t
mdg_scale_vdw_factor(double factor)
{
  if (factor<0.0) {
    fprintf(FD, "ERROR(%s): negative factor=%e is not allowed, set to zero.\n",
	    __func__,factor);
    return 0;
  }
  if (factor==0.0) return 0;

  factor *= MDG_SC->vdw_energy_factor;

  int ex = ilogb(factor);
  uint32_t z = round(ldexp(factor, -ex+23));
  if ( z>=(1<<24) ) { 
    z = 0;
    ++ex;
  } else {
    z &= MASK(23);
  }
  ex += 32;
  if (ex>127) {
    fprintf(FD, "ERROR(%s): factor=%e overflow.  Set to max.\n",
	    __func__, factor);
    return MASK(29);
  }
  if (ex<=0) {
    fprintf(FD, "WARNING(%s): factor=%e underflow. Set to zero.\n",
	    __func__, factor);
    return 0;
  }
  z += (ex<<23);
  return z;
}

/**
   @brief Reverse conversion of coulomb group-based cutoff factor
   @param[in] scaled cutoff factor, floating point 23-bit mantissa, 5-bit exponent, 1 sign.
   exponent bias : 16
   @retval cutoff factor in double (divided by MDG_SC->coulomb_energy_factor)
*///added (Komatsu)
double
mdg_rev_scale_coulomb_factor(uint32_t s)
{
  double e=TO_DOUBLE_U(s, 23, 5, 16);
  if((s>>28)&MASK(1))e=-e;
  return e/MDG_SC->coulomb_energy_factor;
}

/**
   @brief Conversion of coulomb group-based cutoff factor (also used for 1-4 exclusion)
   @param[in] factor cutoff factor, multiplied by 1/(4 pi eps0) in unit of kJ nm mol^-1 e^-2
   @retval scaled cutoff factor, floating point 23-bit mantissa, 5-bit exponent, 1 sign.
   exponent bias : 16
*/
uint32_t
mdg_scale_coulomb_factor(double factor)
{
  int sign = 0;
  double fac=factor;
  if (fac<0.0) {
    sign = 1<<28;
    fac = -fac;
  } else if (fac==0.0) 
    return 0;

  fac *= MDG_SC->coulomb_energy_factor;
  int ex = ilogb(fac);
  uint32_t z = round(ldexp(fac, -ex+23));
  if ( z>=(1<<24) ) { 
    z = 0;
    ++ex;
  } else {
    z &= MASK(23);
  }
  ex += 16;
  if (ex>31) {
    fprintf(FD, "ERROR(%s): factor=%e overflow.  Set to max.\n",
	    __func__, factor);
    return sign+MASK(28);
  }
  if (ex<=0) {
    fprintf(FD, "WARNING(%s): factor=%e underflow. Set to zero.\n",
	    __func__, factor);
    return 0;
  }
  z += (ex<<23);
  return z;
}

double
mdg_scale_potential(int64_t x)
{
  double z = x;
  z = ldexp(z, -MDG_SC->kernel_bp_potential);
  return z;
}

double
mdg_scale_force(int32_t x)
{
  double z = ldexp((double)x, -MDG_SC->bp_force);
  return z;
}

double
mdg_scale_force64(int64_t x)
{
  double z = ldexp((double)x, -MDG_SC->bp_force);
  return z;
}

/**
   @brief Reverse conversion of coulomb_exclusion_factor
   @param[in] s  scaled coordinate
   @retval coordinate in double
*/
double
mdg_rev_scale_coulomb_exclusion(int32_t s) //added (Komatsu) NOT FULLY CONFIRMED
{
  if(BIT(s,26))s|=0x1f<<27;
  double x = ldexp((double)(s), -26);
  return x;
}

int32_t
mdg_scale_coulomb_exclusion(double x)
{
  // set negated value
  if (x>1.0) {
    fprintf(FD, "ERROR(%s): excl=%e overflow. Set to 1.\n",
	    __func__, x);
    return -(1<<26);
  } else if (x<-1.0) {
    fprintf(FD, "ERROR(%s): excl=%e underflow. Set to -1.\n",
	    __func__, x);
    return MASK(26);
  }
  int32_t z = -round(ldexp(x,26));
  //printf("excl %f %x\n", x, z);
  if (z<-(1<<26)) z = -(1<<26);
  if (z>=(1<<26)) z = MASK(26);

  return z;
}

/**
   @brief Conversion of Coulomb softcore parameters
   @param[in] r0: switching point, should be less than typical length (1nm)
   @param[in] f0: magnitude of force at r0; Do not include 1/(4 pi eps0). For example, please set 1/r0^2 for direct
   @param[in] phi0: magnitude of potential at r0; Do not include 1/(4 pi eps0). For example, please set 1/r0 for direct
   @param[out] r0i: switching radius
   @param[out] f0r0i: f(r0)/r0, floating normalized to 1.0-2.0
   @param[out] k0i:   fixed, normalized <1.0, with shift factor
   @param[out] ex_norm_force: exponent shift factor for force
   @param[out] ex_norm_potential: exponent shift factor for potential
   @retval scaled cutoff factor, floating point 23-bit mantissa, 5-bit exponent, 1 sign.
   exponent bias : 16
*/
//Modified type uint64_t -> uint32_t for r0i (Komatsu)
void
mdg_scale_coulomb_softcore(double r0, double f0, double phi0,
			   uint32_t *r0i, uint32_t *f0r0i, uint32_t *k0i,
			   int *ex_norm_force, int *ex_norm_potential) 
{
  double r1 = r0 / MDG_SC->length;
  if (r1>=1.0) {
    fprintf(FD, "ERROR(%s): softcore switching point %e overflow.\n",
	    __func__, r0);
  }
  // The following may be wrong for non-standard scaling.
  double f0r0 = f0/r1;
  int exf = ilogb(f0r0);
  *f0r0i = mdg_scale_coulomb_factor(ldexp(f0r0, -exf));
  *ex_norm_force = exf + 1; // mult by 2
  double k0 = r1*r1 + 2.0*phi0/f0r0;
  //printf("r0 = %f f0 = %f f0/r0 = %f pot0 = %f r0^2 = %f 2phi/f0r0 =%f k0 = %f \n",
  //r1, f0, f0r0, phi0, r1*r1, 2.0*phi0/f0r0, k0);
  int ex = ilogb(k0); 
  // ex = 1..-6
  if (ex>=2) {
    fprintf(FD, "ERROR(%s): softcore potential constant k0=%e overflow; r0=%f f0=%f phi0=%f.\n",
	    __func__, k0, r0, f0, phi0);
    *k0i = MASK(29);
  } else {
    if (ex<-6) ex=-6;
    *ex_norm_potential = exf + ex;
    uint32_t z = round(ldexp(k0, 25-ex)); // 26 bit	//add round.(Komatsu)
    ex += 6;
    // ex = 0..7, binary point at 25 - (ex-6) = 31 - ex
    *k0i = (ex<<26) + z;
    //printf("k0=%f k0i=%x ex=%d z=%x\n", k0, *k0i, ex, z);
  }
  
  r1 *= r1;
  uint64_t r2 = round(ldexp(r1, 32));
  if (r2>0x100000000UL) {
    printf("ERROR(%s):switching radius %f overflow\n", __func__,  r0);
  } else {
    *r0i = r2;
  }
  //printf("%x %f\n", *r0i, r0);
}

/**
   @brief Conversion of vdw softcore parameters
   @param[in] r0: switching point, divided by vdw radius.
   @param[in] f0: magnitude of force at r0; Do not include 4eps_ij/sij^2. 
   @param[in] phi0: magnitude of potential at r0; Do not include 4eps_ij.
   @param[out] r0i: switching radius
   @param[out] f0r0i: f(r0)/r0, floating normalized to 1.0-2.0
   @param[out] k0i:   floating 25bit mantissa, 7bit exponent 
   @param[out] norm_force: exponent shift factor for force
   @param[out] norm_potential: exponent shift factor for potential
   @retval scaled cutoff factor, floating point 23-bit mantissa, 5-bit exponent, 1 sign.
   exponent bias : 16
*/
void
mdg_scale_vdw_softcore(double r0, double f0, double phi0,
		       uint32_t *r0i, uint32_t *f0r0i, uint32_t *k0i,
		       int *norm_force, int *norm_potential)
{
  if (r0>=1.0) {
    fprintf(FD, "ERROR(%s): softcore switching point %e overflow.\n",
	    __func__, r0);
  }
  // The following may be wrong for non-standard scaling.
  double f0r0 = f0/r0;
  int exf = ilogb(f0r0);
  *f0r0i = mdg_scale_vdw_factor(ldexp(f0r0, -exf));
  *norm_force = 237 - exf - MDG_SC->kernel_bp_force; 
  *norm_potential = 240 - exf - MDG_SC->kernel_bp_potential;
  double k0 = r0*r0 + 2.0*phi0/f0r0;
  int ex = ilogb(k0); 
  int man = round(ldexp(k0, 25-ex));	//add round.(Komatsu)
  man &= MASK(25);
  ex +=76;
  if ((ex<=0)||(ex>=128)) {
    fprintf(FD, "ERROR(%s): vdw softcore potential factor %e out of range.\n",
	    __func__, k0);
  }
  *k0i = (ex<<25)+man;

  double r2=r0*r0;
  ex = ilogb(r2);
  man = round(ldexp(r2, 25-ex));	//add round.(Komatsu)
  man &= MASK(25);
  ex +=76;
  if ((ex<=0)||(ex>=128)) {
    fprintf(FD, "ERROR(%s): vdw softcore switch radius (r/s) %e out of range.\n",
	    __func__, r0);
  }
  *r0i = (ex<<25)+man;
  //printf("%x %f\n", *r0i, r0);
}

mdg_scale_conversion_constants_t *
mdg_set_gromacs_conversion()
{
  if(MDG_SC!=0)return MDG_SC; //Skip if already assigned (Komatsu)
  mdg_scale_conversion_constants_t *sc;
  sc = (mdg_scale_conversion_constants_t *)malloc(sizeof(mdg_scale_conversion_constants_t));
  sc->length = 1.0; // nm
  sc->charge = 1.0; // e
  sc->vdw_energy_factor = 4.0; // kJ/mol; is this correct? if eps does not contain factor 4...

  // Dielectric constant of vaccum eps_0: 8.85418782 x 10^-12 F/m = m^-3 kg^-1 s^4 A^2 
  // e = 1.60217653 x 10^-19 C = s A
  // NA = 6.02214179 x 10^23 
  // 1/4 pi eps_0 * e^2 = 0.0230707725 x 10^-26 m^2 kg s^-2
  //                    = 2.30707725 x 10^-28 m^3 kg s^-2
  //                    = 2.30707725 x 10^-19 nm J
  //                    = 13.89354634 x 10^4 nm J/mol
  //                    = 138.9354634 nm kJ/mol

  sc->coulomb_energy_factor = 138.9354634; // nm kJ/mol; Energy at r=1nm

  sc->kernel_bp_force     = 16 + 32;
  sc->kernel_bp_potential = 32;

  sc->norm_coulomb_force     = 128 - sc->kernel_bp_force;
  sc->norm_coulomb_potential = 128 - sc->kernel_bp_potential;

  // should be subtracted exp_org of vdw table
  sc->norm_vdw_force     = 135 - sc->kernel_bp_force;
  sc->norm_vdw_potential = 135 - sc->kernel_bp_potential;

  sc->bp_force = sc->kernel_bp_force - 32 - 2;
  sc->bp_grid_potential = 20;

  MDG_SC = sc;
  return sc;
}

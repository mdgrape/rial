/*
    Random Number Generator
*/
#include "math.h"
#include "defs.h"
#include "mt.h"
#include "random.h"

uint64_t mdg_genrand64()
{
  return (((uint64_t)mdg_genrand())<<32) + mdg_genrand();
}

double
mdg_frnd() /* uniform random number from 0..1 */
{
  /*
  union { double f; longlong i; } r;
  r.i.h=(genrand()&0xfffff)+0x3ff00000;
  r.i.l=genrand();
  r.f=r.f-1.0;
  r.i.h=(genrand()&0xfffff)+(r.i.h & 0x3ff00000);
  r.i.l=genrand();
  return r.f;
  */
  uint64_t ri = mdg_genrand64();
  double r = ri;
  return scalbn(r,-64);
}

double
mdg_frnds() /* uniform random number from -1..1 */
{/*
  DOUBLE r;
  r.i=(mdg_genrand64()&0xfffffffffffffULL)+0x3ff00000;
  r.f=r.f-2.0;
  r.i=(mdg_genrand64()&0xfffffffffffffULL)+(r.i & 0x3ff0000000000000ULL);
  return r.f;
 */
  uint64_t ri = mdg_genrand64();
  double r = ri;
  return scalbn(r,-63)-1.0;
}





#include <stdio.h>
#include <stdint.h>

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static uint32_t mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */
static int pos=0; // for 1by1

/* initializes mt[N] with a seed */
void mdg_sgenrand(uint32_t seed)
{
    mt[0]= seed & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = (1664525UL * mt[mti-1] + 1UL);
        mt[mti] &= 0xffffffffUL;
    }
    pos=0;
}

uint32_t mdg_seed_get(int adr) {
  if ((adr<0)||(adr>=N)) {
    return 0;
  } 
  return mt[adr];
}

void mdg_seed_set(int adr, uint32_t seed) {
  if ((adr>=0)&&(adr<N)) {
    mt[adr]=seed;
  } 
}

/* generates a random number on [0,0xffffffff]-interval */
uint32_t mdg_genrand(void)
{
    uint32_t y;
    static uint32_t mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;
	
        if (mti == N+1)   /* if sgenrand() has not been called, */
	  mdg_sgenrand(5489UL); /* a default initial seed is used   */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
    //printf("%08x\n", mt[mti]);
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

// Generate one-by-one
uint32_t mdg_genrand_1by1(void)
{
    uint32_t y;
    static uint32_t mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if ((pos<0)||(pos>=N)) pos=0;
    int kk1=pos+1;
    if (kk1==N) kk1=0;
    int kkM=pos+M;
    if (kkM>=N) kkM-=N;
    
    y = (mt[pos]&UPPER_MASK)|(mt[kk1]&LOWER_MASK);
    y = mt[pos] = mt[kkM] ^ (y >> 1) ^ mag01[y & 0x1UL];
    //printf("%08x\n", mt[mti]);
    pos++;
    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
int32_t mdg_genrand31(void)
{
    return (int32_t)(mdg_genrand()>>1);
}

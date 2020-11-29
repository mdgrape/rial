#ifndef __DEFS_H
#define __DEFS_H

#include <stdint.h>

// For MacOS, use __BIG_ENDIAN__ / __LITTLE_ENDIAN__
//   (MacOS also has machine/endian.h, but value differs.)
// For Linux, use endian.h

#undef BYTE_ORDER
#if defined(__BIG_ENDIAN__)
#define BYTE_ORDER __BIG_ENDIAN
#elif defined(__LITTLE_ENDIAN__)
#define BYTE_ORDER __LITTLE_ENDIAN
#else
#include <endian.h>
#endif

typedef union { double f; uint64_t i; } DOUBLE;   ///< union to convert double <-> uint64_t
typedef union { float f; uint32_t i; } FLOAT; ///< union to convert float <-> uint32_t

#define MASK(n) ((0x1U<<(n)) -1)
#define MASK31 0x7fffffff
#define MASK64(n) ((0x1LLU<<(n)) -1)
#define BIT(n,b) (( (n)>>(b) )&1)
#define SLICE(n,pos,width) (( (n)>>(pos) )& MASK(width) )
#define SLICE64(n,pos,width) (( (n)>>(pos) )& MASK64(width) )
#define EXTEND(x,n) ( ((x)&(1<<((n)-1))) ? (((x)&MASK(n))-(1<<(n))) : x )
#define EXTENDLL(x,n) ( ((x)&(1LL<<((n)-1))) ? (((x)&MASK64(n))-(1LL<<(n))) : x )
#define MANTISSA(x,n) ( ((x)&MASK64(n)) + (1LL<<(n)) )
#define ROUND(x,n) ( ((x)>>(n)) + ((x>>((n)-1)) & 1) )
#define ROUNDUB(x,n) ( ((x)>>(n)) + ( BIT(x,n-1) && ( BIT(x,n) || ((x) & MASK64(n-1)) ) ) ? 1 : 0 )

#define TO_DOUBLE_U(x,m,n,bias)  (ldexp(MANTISSA(x,m), SLICE(x,m,n)-(m)-(bias)));
#define TO_DOUBLE_S(x,m,n,bias)  (ldexp((int64_t)EXTENDLL((x)&MASK64(m),m), SLICE(x,m,n)-(m)+1-(bias)));

#define SFTMASK(x,n,s) (((x)&MASK(n))<<(s))
#define PACK2(x0,w0,x1,w1) (SFTMASK(x0,w0,0)+SFTMASK(x1,w1,w0))
#define PACK3(x0,w0,x1,w1,x2,w2) (SFTMASK(x0,w0,0)+SFTMASK(x1,w1,w0)+SFTMASK(x2,w2,(w0)+(w1)))

//#define B_MDG_CHECK (getenv("B_MDG_CHECK")!=NULL)	//for debug (Komatsu)

#endif

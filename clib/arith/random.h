#ifndef __RANDOM_H
#define __RANDOM_H
#include <stdint.h>
/* random.c */
#ifdef __cplusplus
extern "C" {
#endif
  double mdg_frnd(void);
  double mdg_frnds(void);
  uint64_t mdg_genrand64();
#ifdef __cplusplus
}
#endif
#endif

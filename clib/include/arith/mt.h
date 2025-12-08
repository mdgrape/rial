#ifndef MT_H
#define MT_H
#include <stdint.h>

/* mt.c */
#ifdef __cplusplus
extern "C" {
#endif
  void mdg_sgenrand(uint32_t seed);
  uint32_t mdg_seed_get(int adr);
  void mdg_seed_set(int adr, uint32_t seed);
  uint32_t mdg_genrand(void);
  uint32_t mdg_genrand_1by1(void);
  int32_t mdg_genrand31(void);
#ifdef __cplusplus
}
#endif
#endif

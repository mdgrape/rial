#ifndef __MDG_THREEFRY_H
#define __MDG_THREEFRY_H

#include "Random123/threefry.h"

typedef struct {
  threefry4x32_ctr_t ctr;
  threefry4x32_key_t key;
} mdg_threefry4x32_t;

#ifdef __cplusplus
extern "C" {
#endif
  mdg_threefry4x32_t * mdg_threefry4x32_init(uint32_t * count, uint32_t * key);
  void mdg_threefry4x32(mdg_threefry4x32_t *p, uint32_t *r);
  void mdg_threefry4x32_free(mdg_threefry4x32_t *p);
#ifdef __cplusplus
}
#endif
#endif

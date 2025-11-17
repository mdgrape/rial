/*
    Random Number Generator using Random123
    for Scala interface
*/
#include <stdint.h>
#include <stdlib.h>
#include "mdg_threefry.h"
#include "Random123/threefry.h"

mdg_threefry4x32_t * mdg_threefry4x32_init(uint32_t * count, uint32_t * key) {
  mdg_threefry4x32_t *p = (mdg_threefry4x32_t *)malloc(sizeof(mdg_threefry4x32_t));
  for(int i=0;i<4;++i) {
    p->ctr.v[i] = count[i];
    p->key.v[i] = key[i];
  }
  return p;
}

void mdg_threefry4x32(mdg_threefry4x32_t *p, uint32_t *rand) {
  threefry4x32_ctr_t r = threefry4x32(p->ctr, p->key);
  // Copy result
  for(int i=0;i<4;++i) rand[i] = r.v[i];
  // Increment 128bit counter
  for(int i=0;i<4;++i) {
    ++p->ctr.v[i];
    if (p->ctr.v[i] != 0) break;
  }
}

void mdg_threefry4x32_free(mdg_threefry4x32_t *p) {
  free(p);
}

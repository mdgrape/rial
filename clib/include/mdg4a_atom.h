#ifndef __MDG4A_ATOM_H
#define __MDG4A_ATOM_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  union {
    uint32_t val[16];
    struct{ 
  uint32_t coord[3];
  uint32_t address;
  int32_t charge;
  uint32_t vdw_type;  // b[31:16] residue index, b[15:0] vdw_type
  uint32_t id;        // b[31:23] cellkey, b[22:0] atomID
  uint32_t exclude;
  int32_t force[3];
  int32_t inv_mass;
  int32_t vel[3];
  int32_t mass;
    };
  };
} mdg4a_atom_t;

#ifdef __cplusplus
}
#endif
#endif //__MDG4A_ATOM_H

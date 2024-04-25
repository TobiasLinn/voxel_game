#ifndef VOXEL_HPP
#define VOXEL_HPP

#include "shape.hpp"

namespace vxg {

// represent voxel as 32 bit integer
using voxel = u32;

// masks
static constexpr u32 MATERIAL_MASK = 0x00FFFFFF;
static constexpr u32 SHAPE_MASK    = 0xFF000000;
    
inline u32 get_material_index(voxel vx) {
    return (vx & MATERIAL_MASK);
}
inline voxel set_material_index(voxel vx, u32 mati) {
    return (vx & SHAPE_MASK) | mati;
}

inline u32 get_shape_index(voxel vx) {
    return (vx & SHAPE_MASK) >> 24;
}
inline voxel set_shape_index(voxel vx, u32 shi) {
    return (vx & MATERIAL_MASK) | (shi << 24);
}

inline voxel get_voxel(u32 shi, u32 mati) {
    return (shi << 24) | mati;
}

}

#endif

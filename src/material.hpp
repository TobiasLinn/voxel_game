#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include <unordered_map>

#include "util.hpp"

namespace vxg {

struct material {
    // color
    u32 col;
};

static constexpr u32 MAT_EMPTY = 0;
static constexpr u32 MAT_STONE = 1;

static const std::unordered_map<u32, material> materials = {
    { MAT_EMPTY, { 0x00000000 } },
    { MAT_STONE, { 0x00808080 } }
};

}

#endif

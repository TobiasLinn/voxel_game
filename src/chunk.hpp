#ifndef CHUNK_HPP
#define CHUNK_HPP

#include <algorithm>

#include "bounding_box.hpp"
#include "voxel.hpp"

namespace vxg {

template<class T, T EMPTY, i32 LS>
class chunk {
public:
    static constexpr i32 E  = EMPTY;
    static constexpr i32 S  = (1 <<      LS );
    static constexpr i32 S3 = (1 << (3 * LS));

private:
    // bounding box
    bounding_box<T,EMPTY,S> bb;

    // data
    alignas(16) T data[S3];

public:
    inline chunk() {
        std::fill(std::begin(data), std::end(data), EMPTY);
    }

    inline const bounding_box<T,EMPTY,S>& get_bounding_box() const {
        return bb;
    }

    inline T get_data(i32 x, i32 y, i32 z) const {
        return data[z * S * S + y * S + x];
    }
    inline T get_data(const i32 xyz[3]) const {
        return get_data(xyz[0], xyz[1], xyz[2]);
    }
    inline void set_data(i32 x, i32 y, i32 z, T d) {
        i32 xyz[3] = { x, y, z };
        set_data(xyz, d);
    }
    inline void set_data(const i32 xyz[3], T d) {
        T before = data[xyz[2] * S * S + xyz[1] * S + xyz[0]];
        data[xyz[2] * S * S + xyz[1] * S + xyz[0]] = d;

        // update bounding box
        bb.update(xyz, before, d, data);
    }
};

template<i32 LS>
using internal_chunk = chunk<i32, -1, LS>;

template<i32 LS>
using leaf_chunk = chunk<voxel, 0, LS>;

}

#endif

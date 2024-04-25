#ifndef BOUNDING_BOX_HPP
#define BOUNDING_BOX_HPP

#include "util.hpp"

namespace vxg {

template<class T, T EMPTY, i32 S>
struct bounding_box {
    alignas(16) i32 b0[4]; // minx, miny, minz, 0
    alignas(16) i32 b1[4]; // maxx, maxy, maxz, 0

    inline bounding_box()
        : b0 { S, S, S, 0 },
          b1 { 0, 0, 0, 0 } {
    }

    inline void get_m128i(m128i &min, m128i &max) const {
        min = load_m128i_aligned(reinterpret_cast<const m128i*>(b0));
        max = load_m128i_aligned(reinterpret_cast<const m128i*>(b1));
    }
    inline void update(const i32 xyz[3], T before, T after, const T data[S*S*S]) {
        if ((before == EMPTY) && (after != EMPTY)) {
            for (i32 dir = 0; dir < 3; ++dir) {
                if (xyz[dir] <  b0[dir]) b0[dir] = xyz[dir];
                if (xyz[dir] >= b1[dir]) b1[dir] = xyz[dir] + 1;
            }
        } else if ((before != EMPTY) && (after == EMPTY)) {
            i32 xyz2[3];

            for (i32 dir = 0; dir < 3; ++dir) {
                // get other 2 directions
                static constexpr i32 d1[3] = { 1, 0, 0 };
                static constexpr i32 d2[3] = { 2, 2, 1 };
                i32 dir1 = d1[dir];
                i32 dir2 = d2[dir];

                // update lower bound
                if (xyz[dir] == b0[dir]) {
                    for (i32 k = xyz[dir]; k < S; ++k) {
                        xyz2[dir] = k;
                        for (i32 u = b0[dir1]; u < b1[dir1]; ++u) {
                            xyz2[dir1] = u;
                            for (i32 v = b0[dir2]; v < b1[dir2]; ++v) {
                                xyz2[dir2] = v;
                                if (data[xyz2[2] * S * S + xyz2[1] * S + xyz2[0]] != EMPTY) goto end_k1;
                            }
                        }
                        ++b0[dir];
                    }
                end_k1:;
                }

                // update upper bound
                if (xyz[dir] == b1[dir] - 1) {
                    for (i32 k = xyz[dir]; k >= 0; --k) {
                        xyz2[dir] = k;
                        for (i32 u = b0[dir1]; u < b1[dir1]; ++u) {
                            xyz2[dir1] = u;
                            for (i32 v = b0[dir2]; v < b1[dir2]; ++v) {
                                xyz2[dir2] = v;
                                if (data[xyz2[2] * S * S + xyz2[1] * S + xyz2[0]] != EMPTY) goto end_k2;
                            }
                        }
                        --b1[dir];
                    }
                end_k2:;
                }
            }
        }
    }
};

}

#endif

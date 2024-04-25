#ifndef FAST_MATH_HPP
#define FAST_MATH_HPP

#include "util.hpp"

namespace vxg {

inline i32 fast_floor(f32 x) {
    m128 r = load_m128(x);
    r = floorss(r, r);
    return cvtss2si(r);
}

// fast cosine for x in [-pi/2; pi/2]
inline f32 fast_cos(f32 x) {
    static constexpr f32 a = f32(4.0 / (M_PI * M_PI * M_PI) * (4.0 / M_PI - 1.0));
    static constexpr f32 b = f32(1.0 / M_PI * (8.0 / M_PI - 1.0));

    f32 x2 = x * x;
    f32 x4 = x2 * x2;

    return a * x4 - b * x2 + 1;
}
inline m128 fast_cos(m128 x) {
    static constexpr f32 a = f32(4.0 / (M_PI * M_PI * M_PI) * (4.0 / M_PI - 1.0));
    static constexpr f32 b = f32(1.0 / M_PI * (8.0 / M_PI - 1.0));

    m128 aa = broadcast_m128(a);
    m128 bb = broadcast_m128(b);

    m128 x2 = mulps(x, x);
    m128 x4 = mulps(x2, x2);

    m128 c0 = broadcast_m128(1.0f);
    m128 c1 = mulps(bb, x2);
    m128 c2 = mulps(aa, x4);

    m128 c = subps(c0, c1);
    return addps(c, c2);
}

}

#endif

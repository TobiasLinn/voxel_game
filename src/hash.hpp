#ifndef HASH_HPP
#define HASH_HPP

#include "util.hpp"

namespace vxg {

// reversible integer hash with low bias
inline u32 hash(u32 x) {
    // prevent hash(0) == 0
    x ^= 0xE0AEDEC4u;

    // hash
    x ^= x >> 17;
    x *= 0xED5AD4BBu;
    x ^= x >> 11;
    x *= 0xAC4C1B51u;
    x ^= x >> 15;
    x *= 0x31848BABu;
    x ^= x >> 14;

    return x;
}
inline m128i hash(m128i x) {
    // constants
    m128i c = load_m128i(0xE0AEDEC4u, 0xED5AD4BBu, 0xAC4C1B51u, 0x31848BABu);

    // tmp
    m128i r;

    // prevent hash(0) == 0
    r = pshufd<0,0,0,0>(c);
    x = pxor(x, r);

    // hash
    r = psrld<17>(x);
    x = pxor(x, r);
    r = pshufd<1,1,1,1>(c);
    x = pmulld(x, r);
    r = psrld<11>(x);
    x = pxor(x, r);
    r = pshufd<2,2,2,2>(c);
    x = pmulld(x, r);
    r = psrld<15>(x);
    x = pxor(x, r);
    r = pshufd<3,3,3,3>(c);
    x = pmulld(x, r);
    r = psrld<14>(x);
    x = pxor(x, r);

    return x;
}

inline u32 hash(u32 seed, u32 x) {
    return hash(seed ^ x);
}
inline m128i hash(u32 seed, m128i x) {
    return hash(pxor(broadcast_m128i(seed), x));
}

}

#endif

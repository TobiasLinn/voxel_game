#ifndef NOISE_HPP
#define NOISE_HPP

#include "fast_math.hpp"
#include "hash.hpp"

namespace vxg {

namespace noise_impl {

// convert 32 bit integer to float in [-1.0f; +1.0f] (slight bias because of +/- 0.0)
inline f32 value(i32 x) {
    // get absolute and sign of x
    i32 a = x & 0x7FFFFFFF;
    i32 s = (x >> 31) | 1;

    return s * (f32(a) / f32(0x7FFFFFFF));
}
inline m128 value(m128i x) {
    // sign mask
    m128i msk = broadcast_m128i(0x7FFFFFFF);

    // get absolute and sign bit of x
    m128i a = pand(x, msk);
    m128i s = pandn(msk, x);

    // v in [-1.0f; 1.0f]
    m128 v = cvtdq2ps(a);
    m128 u = cvtdq2ps(msk);
    v = divps(v, u);
    v = orps(m128i_to_m128(s), v);

    return v;
}

inline f32 grad(i32 h, f32 x) {
    return value(h) * x;
}
inline m128 grad(m128i h, m128 x) {
    return mulps(value(h), x);
}

inline f32 grad(i32 h, f32 x, f32 y) {
    // take the two most significant bits of h as signs of cos and sin
    i32 sgn1 = ( h       >> 31) | 1;
    i32 sgn2 = ((h << 1) >> 31) | 1;
    h = h & 0x3FFFFFFF;

    // get angle between 0 and pi/2
    f32 a = (f32(h) / f32(0x3FFFFFFF)) * f32(M_PI_2);

    // get sine and cosine
    f32 c = sgn1 * fast_cos(a);
    f32 s = sgn2 * fast_cos(f32(M_PI_2) - a);

    // scalar product of (x,y) with vector on unit circle
    return c * x + s * y;
}
inline m128 grad(m128i h, m128 x, m128 y) {
    // constants
    m128i msk = broadcast_m128i(0x7FFFFFFF);
    m128 pi_2 = broadcast_m128(f32(M_PI_2));

    // take the two most significant bits of h as signs of cos and sin
    m128i sgn1 = pandn(msk, h);
    m128i sgn2 = pslld<1>(h);
    sgn2 = pandn(msk, sgn2);
    msk = psrld<1>(msk);
    h = pand(msk, h);

    // get angle between 0 and pi/2
    m128 a = cvtdq2ps(h);
    m128 b = cvtdq2ps(msk);
    a = divps(a, b);
    a = mulps(a, pi_2);

    // get sine and cosine
    m128 c = fast_cos(a);
    m128 s = subps(pi_2, a);
    s      = fast_cos(s);
    c      = orps(c, m128i_to_m128(sgn1));
    s      = orps(s, m128i_to_m128(sgn2));

    // scalar product of (x,y) with vector on unit circle
    c = mulps(c, x);
    s = mulps(s, y);

    return addps(c, s);
}

inline f32 grad(i32 h, f32 x, f32 y, f32 z) {
    // get second random number
    i32 h2 = hash(h);

    // take the two most significant bits of h as signs of cos and sin
    i32 sgn1 = ( h       >> 31) | 1;
    i32 sgn2 = ((h << 1) >> 31) | 1;
    h = h & 0x3FFFFFFF;

    // get phi angle between 0 and pi/2
    f32 phi = (f32(h) / f32(0x3FFFFFFF)) * f32(M_PI_2);

    // get sine and cosine
    f32 c = sgn1 * fast_cos(phi);
    f32 s = sgn2 * fast_cos(f32(M_PI_2) - phi);

    // get u in [-1; +1] and v = sqrt(1 - u²)
    f32 u = value(h2);
    f32 v = std::sqrt(1 - u * u);

    // scalar product of (x,y,z) with vector on unit sphere
    return c * v * x + s * v * y + u * z;
}
inline m128 grad(m128i h, m128 x, m128 y, m128 z) {
    // constants
    m128i msk = broadcast_m128i(0x7FFFFFFF);
    m128 pi_2 = broadcast_m128(f32(M_PI_2));

    // get second random number
    m128i h2 = hash(h);

    // take the two most significant bits of h as signs of cos and sin
    m128i sgn1 = pandn(msk, h);
    m128i sgn2 = pslld<1>(h);
    sgn2 = pandn(msk, sgn2);
    msk = psrld<1>(msk);
    h = pand(msk, h);

    // get phi angle between 0 and pi/2
    m128 phi = cvtdq2ps(h);
    m128 b = cvtdq2ps(msk);
    phi = divps(phi, b);
    phi = mulps(phi, pi_2);

    // get sine and cosine
    m128 c = fast_cos(phi);
    m128 s = subps(pi_2, phi);
    s      = fast_cos(s);
    c      = orps(c, m128i_to_m128(sgn1));
    s      = orps(s, m128i_to_m128(sgn2));

    // get u in [-1; +1] and v = sqrt(1 - u²)
    m128 u   = value(h2);
    m128 one = broadcast_m128(1.0f);
    m128 v   = mulps(u, u);
    v        = subps(one, v);
    v        = sqrtps(v);

    // scalar product of (x,y,z) with vector on unit sphere
    c = mulps(c, v);
    c = mulps(c, x);
    s = mulps(s, v);
    s = mulps(s, y);
    c = addps(c, s);
    s = mulps(u, z);

    return addps(c, s);
}

}

inline f32 value_noise(u32 seed, f32 x) {
    using namespace vxg::noise_impl;

    // corner coordinates (nearest integer values)
    i32 i0 = fast_floor(x);
    i32 i1 = i0 + 1;

    // hash corner coordinates
    i32 h0 = hash(seed, i0);
    i32 h1 = hash(seed, i1);

    // scale hash values to [-1.0f; +1.0f]
    f32 v0 = value(h0);
    f32 v1 = value(h1);

    // distance to left corner
    f32 xx = x - i0;

    // get weight
    f32 w = xx*xx*xx * (10 + 3*xx * (2*xx - 5));

    // interpolate
    return v0 + w * (v1 - v0);
}
inline m128 value_noise(u32 seed, m128 x) {
    using namespace vxg::noise_impl;

    // corner coordinates (nearest integer values)
    m128 one = broadcast_m128(1.0f);
    m128 i0f = floorps(x);
    m128 i1f = addps(i0f, one);
    m128i i0 = cvtps2dq(i0f);
    m128i i1 = cvtps2dq(i1f);

    // hash corner coordinates
    m128i h0 = hash(seed, i0);
    m128i h1 = hash(seed, i1);

    // scale hash values to [-1.0f; +1.0f]
    m128 v0 = value(h0);
    m128 v1 = value(h1);

    // distance to left corner
    m128 xx = subps(x, i0f);

    // get weight: use quintic S-curve w = dx0^3 * (2*5 + 3*dx0 * (2*dx0 - 5)
    m128 five   = broadcast_m128(5.0f);
    m128 ten    = addps(five, five);
    m128 twox   = addps(xx, xx);
    m128 threex = addps(twox, xx);
    m128 xcubed = mulps(xx, xx);
    xcubed      = mulps(xcubed, xx);
    m128 w      = subps(twox, five);
    w           = mulps(w, threex);
    w           = addps(w, ten);
    w           = mulps(w, xcubed);

    // interpolate: v = (1-w) * v0 + w * v1 = v0 + w * (v1 - v0)
    m128 v = subps(v1, v0);
    v = mulps(v, w);
    return addps(v0, v);
}

inline f32 value_noise(u32 seed, f32 x, f32 y) {
    using namespace vxg::noise_impl;

    // corner coordinates (nearest integer values)
    i32 i0 = fast_floor(x);
    i32 i1 = i0 + 1;
    i32 j0 = fast_floor(y);
    i32 j1 = j0 + 1;

    // hash corner coordinates
    i32 h0  = hash(seed, i0);
    i32 h1  = hash(seed, i1);
    i32 h00 = hash(seed, j0 + h0);
    i32 h01 = hash(seed, j0 + h1);
    i32 h10 = hash(seed, j1 + h0);
    i32 h11 = hash(seed, j1 + h1);

    // scale hash values to [-1.0f; +1.0f]
    f32 v00 = value(h00);
    f32 v01 = value(h01);
    f32 v10 = value(h10);
    f32 v11 = value(h11);

    // distance to left/lower corner
    f32 xx = x - i0;
    f32 yy = y - j0;

    // get weights
    f32 wx = xx*xx*xx * (10 + 3*xx * (2*xx - 5));
    f32 wy = yy*yy*yy * (10 + 3*yy * (2*yy - 5));

    // bilinear interpolation
    f32 vx0 = v00 + wx * (v01 - v00);
    f32 vx1 = v10 + wx * (v11 - v10);
    return vx0 + wy * (vx1 - vx0);
}
inline m128 value_noise(u32 seed, m128 x, m128 y) {
    using namespace vxg::noise_impl;

    // corner coordinates (nearest integer values)
    m128 one = broadcast_m128(1.0f);
    m128 i0f = floorps(x);
    m128 i1f = addps(i0f, one);
    m128 j0f = floorps(y);
    m128 j1f = addps(j0f, one);
    m128i i0 = cvtps2dq(i0f);
    m128i i1 = cvtps2dq(i1f);
    m128i j0 = cvtps2dq(j0f);
    m128i j1 = cvtps2dq(j1f);

    // hash corner coordinates
    m128i h0  = hash(seed,       i0     );
    m128i h1  = hash(seed,       i1     );
    m128i h00 = hash(seed, paddd(j0, h0));
    m128i h01 = hash(seed, paddd(j0, h1));
    m128i h10 = hash(seed, paddd(j1, h0));
    m128i h11 = hash(seed, paddd(j1, h1));

    // scale hash values to [-1.0f; +1.0f]
    m128 v00 = value(h00);
    m128 v01 = value(h01);
    m128 v10 = value(h10);
    m128 v11 = value(h11);

    // distance to left/lower corner
    m128 xx = subps(x, i0f);
    m128 yy = subps(y, j0f);

    // get weights
    m128 five   = broadcast_m128(5.0f);
    m128 ten    = addps(five, five);
    m128 twox   = addps(xx, xx);
    m128 twoy   = addps(yy, yy);
    m128 threex = addps(twox, xx);
    m128 threey = addps(twoy, yy);
    m128 xcubed = mulps(xx, xx);
    m128 ycubed = mulps(yy, yy);
    xcubed      = mulps(xcubed, xx);
    ycubed      = mulps(ycubed, yy);
    m128 wx     = subps(twox, five);
    m128 wy     = subps(twoy, five);
    wx          = mulps(wx, threex);
    wy          = mulps(wy, threey);
    wx          = addps(wx, ten);
    wy          = addps(wy, ten);
    wx          = mulps(wx, xcubed);
    wy          = mulps(wy, ycubed);

    // bilinear interpolation
    m128 vx0 = subps(v01, v00);
    m128 vx1 = subps(v11, v10);
    vx0      = mulps(vx0, wx);
    vx1      = mulps(vx1, wx);
    vx0      = addps(v00, vx0);
    vx1      = addps(v10, vx1);
    m128 v   = subps(vx1, vx0);
    v        = mulps(v, wy);
    return addps(vx0, v);
}

inline f32 value_noise(u32 seed, f32 x, f32 y, f32 z) {
    using namespace vxg::noise_impl;

    // corner coordinates (nearest integer values)
    i32 i0 = fast_floor(x);
    i32 i1 = i0 + 1;
    i32 j0 = fast_floor(y);
    i32 j1 = j0 + 1;
    i32 k0 = fast_floor(z);
    i32 k1 = k0 + 1;

    // hash corner coordinates
    i32 h0   = hash(seed, i0);
    i32 h1   = hash(seed, i1);
    i32 h00  = hash(seed, j0 + h0);
    i32 h01  = hash(seed, j0 + h1);
    i32 h10  = hash(seed, j1 + h0);
    i32 h11  = hash(seed, j1 + h1);
    i32 h000 = hash(seed, k0 + h00);
    i32 h001 = hash(seed, k0 + h01);
    i32 h010 = hash(seed, k0 + h10);
    i32 h011 = hash(seed, k0 + h11);
    i32 h100 = hash(seed, k1 + h00);
    i32 h101 = hash(seed, k1 + h01);
    i32 h110 = hash(seed, k1 + h10);
    i32 h111 = hash(seed, k1 + h11);

    // scale hash values to [-1.0f; +1.0f]
    f32 v000 = value(h000);
    f32 v001 = value(h001);
    f32 v010 = value(h010);
    f32 v011 = value(h011);
    f32 v100 = value(h100);
    f32 v101 = value(h101);
    f32 v110 = value(h110);
    f32 v111 = value(h111);

    // distance to left/lower/front corner
    f32 xx = x - i0;
    f32 yy = y - j0;
    f32 zz = z - k0;

    // get weights
    f32 wx = xx*xx*xx * (10 + 3*xx * (2*xx - 5));
    f32 wy = yy*yy*yy * (10 + 3*yy * (2*yy - 5));
    f32 wz = zz*zz*zz * (10 + 3*zz * (2*zz - 5));

    // trilinear interpolation
    f32 vx00 = v000 + wx * (v001 - v000);
    f32 vx01 = v010 + wx * (v011 - v010);
    f32 vx10 = v100 + wx * (v101 - v100);
    f32 vx11 = v110 + wx * (v111 - v110);
    f32 vxy0 = vx00 + wy * (vx01 - vx00);
    f32 vxy1 = vx10 + wy * (vx11 - vx10);
    return vxy0 + wz * (vxy1 - vxy0);
}
inline m128 value_noise(u32 seed, m128 x, m128 y, m128 z) {
    using namespace vxg::noise_impl;

    // corner coordinates (nearest integer values)
    m128 one = broadcast_m128(1.0f);
    m128 i0f = floorps(x);
    m128 i1f = addps(i0f, one);
    m128 j0f = floorps(y);
    m128 j1f = addps(j0f, one);
    m128 k0f = floorps(z);
    m128 k1f = addps(k0f, one);
    m128i i0 = cvtps2dq(i0f);
    m128i i1 = cvtps2dq(i1f);
    m128i j0 = cvtps2dq(j0f);
    m128i j1 = cvtps2dq(j1f);
    m128i k0 = cvtps2dq(k0f);
    m128i k1 = cvtps2dq(k1f);

    // hash corner coordinates
    m128i h0   = hash(seed,       i0      );
    m128i h1   = hash(seed,       i1      );
    m128i h00  = hash(seed, paddd(j0, h0 ));
    m128i h01  = hash(seed, paddd(j0, h1 ));
    m128i h10  = hash(seed, paddd(j1, h0 ));
    m128i h11  = hash(seed, paddd(j1, h1 ));
    m128i h000 = hash(seed, paddd(k0, h00));
    m128i h001 = hash(seed, paddd(k0, h01));
    m128i h010 = hash(seed, paddd(k0, h10));
    m128i h011 = hash(seed, paddd(k0, h11));
    m128i h100 = hash(seed, paddd(k1, h00));
    m128i h101 = hash(seed, paddd(k1, h01));
    m128i h110 = hash(seed, paddd(k1, h10));
    m128i h111 = hash(seed, paddd(k1, h11));

    // scale hash values to [-1.0f; +1.0f]
    m128 v000 = value(h000);
    m128 v001 = value(h001);
    m128 v010 = value(h010);
    m128 v011 = value(h011);
    m128 v100 = value(h100);
    m128 v101 = value(h101);
    m128 v110 = value(h110);
    m128 v111 = value(h111);

    // distance to left/lower/front corner
    m128 xx = subps(x, i0f);
    m128 yy = subps(y, j0f);
    m128 zz = subps(z, k0f);

    // get weights
    m128 five   = broadcast_m128(5.0f);
    m128 ten    = addps(five, five);
    m128 twox   = addps(xx, xx);
    m128 twoy   = addps(yy, yy);
    m128 twoz   = addps(zz, zz);
    m128 threex = addps(twox, xx);
    m128 threey = addps(twoy, yy);
    m128 threez = addps(twoz, zz);
    m128 xcubed = mulps(xx, xx);
    m128 ycubed = mulps(yy, yy);
    m128 zcubed = mulps(zz, zz);
    xcubed      = mulps(xcubed, xx);
    ycubed      = mulps(ycubed, yy);
    zcubed      = mulps(zcubed, zz);
    m128 wx     = subps(twox, five);
    m128 wy     = subps(twoy, five);
    m128 wz     = subps(twoz, five);
    wx          = mulps(wx, threex);
    wy          = mulps(wy, threey);
    wz          = mulps(wz, threez);
    wx          = addps(wx, ten);
    wy          = addps(wy, ten);
    wz          = addps(wz, ten);
    wx          = mulps(wx, xcubed);
    wy          = mulps(wy, ycubed);
    wz          = mulps(wz, zcubed);

    // trilinear interpolation
    m128 vx00 = subps(v001, v000);
    m128 vx01 = subps(v011, v010);
    m128 vx10 = subps(v101, v100);
    m128 vx11 = subps(v111, v110);
    vx00      = mulps(vx00, wx);
    vx01      = mulps(vx01, wx);
    vx10      = mulps(vx10, wx);
    vx11      = mulps(vx11, wx);
    vx00      = addps(vx00, v000);
    vx01      = addps(vx01, v010);
    vx10      = addps(vx10, v100);
    vx11      = addps(vx11, v110);
    m128 vxy0 = subps(vx01, vx00);
    m128 vxy1 = subps(vx11, vx10);
    vxy0      = mulps(vxy0, wy);
    vxy1      = mulps(vxy1, wy);
    vxy0      = addps(vxy0, vx00);
    vxy1      = addps(vxy1, vx10);
    m128 vxyz = subps(vxy1, vxy0);
    vxyz      = mulps(vxyz, wz);
    vxyz      = addps(vxyz, vxy0);

    return vxyz;
}

inline f32 simplex_noise(u32 seed, f32 x) {
    using namespace vxg::noise_impl;

    // corner coordinates (nearest integer values)
    i32 i0 = fast_floor(x);
    i32 i1 = i0 + 1;

    // hash corner coordinates
    i32 h0 = hash(seed, i0);
    i32 h1 = hash(seed, i1);

    // distances to corners
    f32 dx0 = x - i0;
    f32 dx1 = x - i1;

    // get both gradient scalar products
    f32 g0 = grad(h0, dx0);
    f32 g1 = grad(h1, dx1);

    // calculate the contribution from both corners
    f32 t0 = 1 - dx0 * dx0;
    f32 t1 = 1 - dx1 * dx1;
    t0 *= t0;
    t1 *= t1;
    f32 n0 = t0 * t0 * g0;
    f32 n1 = t1 * t1 * g1;

    // scale by (3/4)^(-4) to fit in [-1.0f; 1.0f]
    return 3.16049382716f * (n0 + n1);
}
inline m128 simplex_noise(u32 seed, m128 x) {
    using namespace vxg::noise_impl;

    // corner coordinates (nearest integer values)
    m128 one = broadcast_m128(1.0f);
    m128 i0f = floorps(x);
    m128 i1f = addps(i0f, one);
    m128i i0 = cvtps2dq(i0f);
    m128i i1 = cvtps2dq(i1f);

    // hash corner coordinates
    m128i h0 = hash(seed, i0);
    m128i h1 = hash(seed, i1);

    // distances to corners
    m128 dx0 = x - i0f;
    m128 dx1 = x - i1f;

    // get both gradient scalar products
    m128 g0 = grad(h0, dx0);
    m128 g1 = grad(h1, dx1);

    // calculate the contribution from both corners
    m128 t0 = mulps(dx0, dx0);
    m128 t1 = mulps(dx1, dx1);
    t0      = subps(one, t0);
    t1      = subps(one, t1);
    t0      = mulps(t0, t0);
    t1      = mulps(t1, t1);
    t0      = mulps(t0, t0);
    t1      = mulps(t1, t1);
    m128 n0 = mulps(t0, g0);
    m128 n1 = mulps(t1, g1);

    // scale by (3/4)^(-4) to fit in [-1.0f; 1.0f]
    m128 g = addps(n0, n1);
    m128 scale = broadcast_m128(3.16049382716f);
    return mulps(scale, g);
}

inline f32 simplex_noise(u32 seed, f32 x, f32 y) {
    using namespace vxg::noise_impl;

    // skewing/unskewing factors for 2D
    static constexpr f32 F2 = 0.366025403f;  // F2 = (sqrt(3) - 1) / 2
    static constexpr f32 G2 = 0.211324865f;  // G2 = (3 - sqrt(3)) / 6

    // skew input space to determine simplex cell
    f32 s  = (x + y) * F2;
    f32 xs = x + s;
    f32 ys = y + s;
    i32 i0 = fast_floor(xs);
    i32 j0 = fast_floor(ys);

    // unskew (i0,j0) back to get cell origin
    f32 t = (i0 + j0) * G2;
    f32 x0 = i0 - t;
    f32 y0 = j0 - t;

    // get distances from cell origin
    f32 dx0 = x - x0;
    f32 dy0 = y - y0;

    // get distances from third corner
    f32 dx2 = dx0 - 1 + 2 * G2;
    f32 dy2 = dy0 - 1 + 2 * G2;

    // hash values for j
    i32 hj0 = hash(seed, j0    );
    i32 hj1 = hash(seed, j0 + 1);

    // hash values for first and third corner
    i32 h0 = hash(seed, i0 + hj0    );
    i32 h2 = hash(seed, i0 + hj1 + 1);

    // determine middle corner
    f32 dx1, dy1;
    i32 h1;
    if (dx0 > dy0) { // lower triangle
        dx1 = dx0 - 1 + G2;
        dy1 = dy0     + G2;
        h1  = hash(seed, i0 + hj0 + 1);
    } else { // upper triangle
        dx1 = dx0     + G2;
        dy1 = dy0 - 1 + G2;
        h1  = hash(seed, i0 + hj1    );
    }

    // contributions from each corner
    f32 n0, n1, n2;

    // first corner
    f32 t0 = 0.5f - dx0*dx0 - dy0*dy0;
    if (t0 < 0) {
        n0 = 0;
    } else {
        t0 *= t0;
        n0 = t0 * t0 * grad(h0, dx0, dy0);
    }

    // second corner
    f32 t1 = 0.5f - dx1*dx1 - dy1*dy1;
    if (t1 < 0) {
        n1 = 0;
    } else {
        t1 *= t1;
        n1 = t1 * t1 * grad(h1, dx1, dy1);
    }

    // third corner
    f32 t2 = 0.5f - dx2*dx2 - dy2*dy2;
    if (t2 < 0) {
        n2 = 0;
    } else {
        t2 *= t2;
        n2 = t2 * t2 * grad(h2, dx2, dy2);
    }

    // add contributions and scale to interval [-1,1];
    return 98.0f * (n0 + n1 + n2);
}
inline m128 simplex_noise(u32 seed, m128 x, m128 y) {
    using namespace vxg::noise_impl;

    // skewing/unskewing factors for 2D
    static constexpr f32 F2 = 0.366025403f;  // F2 = (sqrt(3) - 1) / 2
    static constexpr f32 G2 = 0.211324865f;  // G2 = (3 - sqrt(3)) / 6
    m128 f = broadcast_m128(F2);
    m128 g = broadcast_m128(G2);

    // tmp
    m128 r, s, t;
    m128i u, v;
    m128i one = broadcast_m128i(1);
    m128 onef = cvtdq2ps(one);
    m128 half = broadcast_m128(0.5f);

    // skew input space to determine simplex cell
    s        = addps(x, y);
    s        = mulps(s, f);
    m128 xs  = addps(x, s);
    m128 ys  = addps(y, s);
    m128 i0f = floorps(xs);
    m128 j0f = floorps(ys);
    m128i i0 = cvtps2dq(i0f);
    m128i j0 = cvtps2dq(j0f);
    m128i i1 = paddd(i0, one);
    m128i j1 = paddd(j0, one);

    // unskew (i0,j0) back to get cell origin
    t       = addps(i0f, j0f);
    t       = mulps(t, g);
    m128 x0 = subps(i0f, t);
    m128 y0 = subps(j0f, t);

    // get distances from cell origin
    m128 dx0 = subps(x, x0);
    m128 dy0 = subps(y, y0);

    // get distances from third corner
    r        = subps(dx0, onef);
    s        = subps(dy0, onef);
    t        = addps(g, g);
    m128 dx2 = addps(r, t);
    m128 dy2 = addps(s, t);

    // hash values for j
    m128i hj0 = hash(seed, j0);
    m128i hj1 = hash(seed, j1);

    // hash values for first and third corner
    m128i h0 = paddd(i0, hj0);
    m128i h2 = paddd(i1, hj1);
    h0 = hash(seed, h0);
    h2 = hash(seed, h2);

    // determine middle corner
    m128 msk = cmpgtps(dx0, dy0);
    r = cvtdq2ps(m128_to_m128i(msk)); // r = (dx0 > dy0) ? -1.0f : 0.0f
    s = andnps(msk, onef);            // s = (dx0 > dy0) ?  0.0f : 1.0f
    m128 dx1 = addps(dx0, r);
    m128 dy1 = subps(dy0, s);
    dx1 = addps(dx1, g);
    dy1 = addps(dy1, g);
    u = pand(m128_to_m128i(msk), hj0);  // u = (dx0 > dy0) ? hj0 :   0
    v = pandn(m128_to_m128i(msk), hj1); // v = (dx0 > dy0) ?   0 : hj1
    u = por(u, v);                      // u = (dx0 > dy0) ? hj0 : hj1
    v = psrld<31>(m128_to_m128i(msk));  // v = (dx0 > dy0) ?   1 :   0
    u = paddd(u, v);                    // u = (dx0 > dy0) ? hj0 + 1 : hj1
    u = paddd(u, i0);                   // u = i0 + (dx0 > dy0) ? hj0 + 1 : hj1
    m128i h1 = hash(seed, u);

    // gradients
    m128 g0 = grad(h0, dx0, dy0);
    m128 g1 = grad(h1, dx1, dy1);
    m128 g2 = grad(h2, dx2, dy2);

    // first corner
    m128 t0, n0;
    r   = mulps(dx0, dx0); // r = dx0*dx0
    s   = mulps(dy0, dy0); // s = dy0*dy0
    t0  = subps(half, r);  // t0 = 0.5f - dx0*dx0
    t0  = subps(t0, s);    // t0 = 0.5f - dx0*dx0 - dy0*dy0
    msk = m128i_to_m128(psrad<31>(m128_to_m128i(t0))); // (t0 < 0)
    t0  = mulps(t0, t0);
    t0  = mulps(t0, t0);
    n0  = mulps(t0, g0);
    n0  = andnps(msk, n0); // n0 = (t0 < 0) ? 0 : ...

    // second corner
    m128 t1, n1;
    r   = mulps(dx1, dx1); // r = dx1*dx1
    s   = mulps(dy1, dy1); // s = dy1*dy1
    t1  = subps(half, r);  // t1 = 0.5f - dx1*dx1
    t1  = subps(t1, s);    // t1 = 0.5f - dx1*dx1 - dy1*dy1
    msk = m128i_to_m128(psrad<31>(m128_to_m128i(t1))); // (t1 < 0)
    t1  = mulps(t1, t1);
    t1  = mulps(t1, t1);
    n1  = mulps(t1, g1);
    n1  = andnps(msk, n1); // n1 = (t1 < 0) ? 0 : ...

    // third corner
    m128 t2, n2;
    r   = mulps(dx2, dx2); // r = dx2*dx2
    s   = mulps(dy2, dy2); // s = dy2*dy2
    t2  = subps(half, r);  // t2 = 0.5f - dx2*dx2
    t2  = subps(t2, s);    // t2 = 0.5f - dx2*dx2 - dy2*dy2
    msk = m128i_to_m128(psrad<31>(m128_to_m128i(t2))); // (t2 < 0)
    t2  = mulps(t2, t2);
    t2  = mulps(t2, t2);
    n2  = mulps(t2, g2);
    n2  = andnps(msk, n2); // n2 = (t2 < 0) ? 0 : ...

    // add contributions and scale to interval [-1,1];
    m128 n = addps(n0, n1);
    n = addps(n, n2);
    r = broadcast_m128(98.0f);
    return mulps(n, r);
}

inline f32 simplex_noise(u32 seed, f32 x, f32 y, f32 z) {
    using namespace vxg::noise_impl;

    // skewing/unskewing factors for 3D
    static constexpr f32 F3 = 1.0f / 3.0f;
    static constexpr f32 G3 = 1.0f / 6.0f;

    // skew input space to determine simplex cell
    f32 s = (x + y + z) * F3;
    f32 xs = x + s;
    f32 ys = y + s;
    f32 zs = z + s;
    i32 i0 = fast_floor(xs);
    i32 j0 = fast_floor(ys);
    i32 k0 = fast_floor(zs);

    // unskew (i0,j0,k0) back to get cell origin
    f32 t = (i0 + j0 + k0) * G3;
    f32 x0 = i0 - t;
    f32 y0 = j0 - t;
    f32 z0 = k0 - t;

    // get distances from cell origin
    f32 dx0 = x - x0;
    f32 dy0 = y - y0;
    f32 dz0 = z - z0;

    // get distances from 4-th corner
    f32 dx3 = dx0 - 1 + 3 * G3;
    f32 dy3 = dy0 - 1 + 3 * G3;
    f32 dz3 = dz0 - 1 + 3 * G3;

    // hash values for k
    i32 hk0 = hash(seed, k0    );
    i32 hk1 = hash(seed, k0 + 1);

    // hash values for j,k
    i32 hjk00 = hash(seed, j0     + hk0);
    i32 hjk11 = hash(seed, j0 + 1 + hk1);

    // hash values for first and fourth corner
    i32 h0 = hash(seed, i0     + hjk00);
    i32 h3 = hash(seed, i0 + 1 + hjk11);

    // determine second and third corner
    f32 dx1, dy1, dz1;
    f32 dx2, dy2, dz2;
    i32 h1, h2;
    if (dx0 >= dy0) {
        if (dy0 >= dz0) { // z y x
            dx1 = dx0 - 1 +     G3;
            dy1 = dy0     +     G3;
            dz1 = dz0     +     G3;
            dx2 = dx0 - 1 + 2 * G3;
            dy2 = dy0 - 1 + 2 * G3;
            dz2 = dz0     + 2 * G3;
            i32 hjk10 = hash(seed, j0 + 1 + hk0);
            h1 = hash(seed, i0 + 1 + hjk00);
            h2 = hash(seed, i0 + 1 + hjk10);
        } else if (dx0 >= dz0) { // y z x
            dx1 = dx0 - 1 +     G3;
            dy1 = dy0     +     G3;
            dz1 = dz0     +     G3;
            dx2 = dx0 - 1 + 2 * G3;
            dy2 = dy0     + 2 * G3;
            dz2 = dz0 - 1 + 2 * G3;
            i32 hjk01 = hash(seed, j0 + hk1);
            h1 = hash(seed, i0 + 1 + hjk00);
            h2 = hash(seed, i0 + 1 + hjk01);
        } else { // y x z
            dx1 = dx0     +     G3;
            dy1 = dy0     +     G3;
            dz1 = dz0 - 1 +     G3;
            dx2 = dx0 - 1 + 2 * G3;
            dy2 = dy0     + 2 * G3;
            dz2 = dz0 - 1 + 2 * G3;
            i32 hjk01 = hash(seed, j0 + hk1);
            h1 = hash(seed, i0     + hjk01);
            h2 = hash(seed, i0 + 1 + hjk01);
        }
    } else {
        if (dy0 < dz0) { // x y z
            dx1 = dx0     +     G3;
            dy1 = dy0     +     G3;
            dz1 = dz0 - 1 +     G3;
            dx2 = dx0     + 2 * G3;
            dy2 = dy0 - 1 + 2 * G3;
            dz2 = dz0 - 1 + 2 * G3;
            i32 hjk01 = hash(seed, j0 + hk1);
            h1 = hash(seed, i0 + hjk01);
            h2 = hash(seed, i0 + hjk11);
        } else if (dx0 < dz0) { // x z y
            dx1 = dx0     +     G3;
            dy1 = dy0 - 1 +     G3;
            dz1 = dz0     +     G3;
            dx2 = dx0     + 2 * G3;
            dy2 = dy0 - 1 + 2 * G3;
            dz2 = dz0 - 1 + 2 * G3;
            i32 hjk10 = hash(seed, j0 + 1 + hk0);
            h1 = hash(seed, i0 + hjk10);
            h2 = hash(seed, i0 + hjk11);
        } else { // z x y
            dx1 = dx0     +     G3;
            dy1 = dy0 - 1 +     G3;
            dz1 = dz0     +     G3;
            dx2 = dx0 - 1 + 2 * G3;
            dy2 = dy0 - 1 + 2 * G3;
            dz2 = dz0     + 2 * G3;
            i32 hjk10 = hash(seed, j0 + 1 + hk0);
            h1 = hash(seed, i0     + hjk10);
            h2 = hash(seed, i0 + 1 + hjk10);
        }
    }

    // contributions from each corner
    f32 n0, n1, n2, n3;

    // first corner
    f32 t0 = 0.6f - dx0*dx0 - dy0*dy0 - dz0*dz0;
    if (t0 < 0) {
        n0 = 0;
    } else {
        t0 *= t0;
        n0 = t0 * t0 * grad(h0, dx0, dy0, dz0);
    }

    // second corner
    f32 t1 = 0.6f - dx1*dx1 - dy1*dy1 - dz1*dz1;
    if (t1 < 0) {
        n1 = 0;
    } else {
        t1 *= t1;
        n1 = t1 * t1 * grad(h1, dx1, dy1, dz1);
    }

    // third corner
    f32 t2 = 0.6f - dx2*dx2 - dy2*dy2 - dz2*dz2;
    if (t2 < 0) {
        n2 = 0;
    } else {
        t2 *= t2;
        n2 = t2 * t2 * grad(h2, dx2, dy2, dz2);
    }

    // fourth corner
    f32 t3 = 0.6f - dx3*dx3 - dy3*dy3 - dz3*dz3;
    if (t3 < 0) {
        n3 = 0;
    } else {
        t3 *= t3;
        n3 = t3 * t3 * grad(h3, dx3, dy3, dz3);
    }

    // add contributions and scale to interval [-1,1];
    return 32 * (n0 + n1 + n2 + n3);
}
inline m128 simplex_noise(u32 seed, m128 x, m128 y, m128 z) {
    using namespace vxg::noise_impl;

    // skewing/unskewing factors for 3D
    static constexpr f32 F3 = 1.0f / 3.0f;
    static constexpr f32 G3 = 1.0f / 6.0f;
    m128 f = broadcast_m128(F3);
    m128 g = broadcast_m128(G3);

    // tmp
    m128 r, s, t;
    m128i u, v;
    m128i one = broadcast_m128i(1);
    m128 onef = cvtdq2ps(one);
    m128 p6   = broadcast_m128(0.6f);

    // skew input space to determine simplex cell
    s        = addps(x, y);
    s        = addps(s, z);
    s        = mulps(s, f);
    m128 xs  = addps(x, s);
    m128 ys  = addps(y, s);
    m128 zs  = addps(z, s);
    m128 i0f = floorps(xs);
    m128 j0f = floorps(ys);
    m128 k0f = floorps(zs);
    m128i i0 = cvtps2dq(i0f);
    m128i j0 = cvtps2dq(j0f);
    m128i k0 = cvtps2dq(k0f);
    m128i i1 = paddd(i0, one);
    m128i j1 = paddd(j0, one);
    m128i k1 = paddd(k0, one);

    // unskew (i0,j0,k0) back to get cell origin
    t       = addps(i0f, j0f);
    t       = addps(t, k0f);
    t       = mulps(t, g);
    m128 x0 = subps(i0f, t);
    m128 y0 = subps(j0f, t);
    m128 z0 = subps(k0f, t);

    // get distances from cell origin
    m128 dx0 = subps(x, x0);
    m128 dy0 = subps(y, y0);
    m128 dz0 = subps(z, z0);

    // get distances from 4-th corner
    m128 dx3 = subps(dx0, onef);
    m128 dy3 = subps(dy0, onef);
    m128 dz3 = subps(dz0, onef);
    r        = addps(g, g);
    r        = addps(r, g);
    dx3      = addps(dx3, r);
    dy3      = addps(dy3, r);
    dz3      = addps(dz3, r);

    // hash values for k
    m128i hk0 = hash(seed, k0);
    m128i hk1 = hash(seed, k1);

    // hash values for j,k
    m128i hjk00 = paddd(j0, hk0);
    m128i hjk11 = paddd(j1, hk1);
    hjk00 = hash(seed, hjk00);
    hjk11 = hash(seed, hjk11);

    // hash values for first and fourth corner
    m128i h0 = paddd(i0, hjk00);
    m128i h3 = paddd(i1, hjk11);
    h0 = hash(seed, h0);
    h3 = hash(seed, h3);

    // determine second and third corner
    m128 mskxy = cmpgeps(dx0, dy0);
    m128 mskxz = cmpgeps(dx0, dz0);
    m128 mskyz = cmpgeps(dy0, dz0);
    m128 dx1   = addps(dx0, g);
    m128 dy1   = addps(dy0, g);
    m128 dz1   = addps(dz0, g);
    m128 dx2   = addps(dx1, g);
    m128 dy2   = addps(dy1, g);
    m128 dz2   = addps(dz1, g);
    m128 xy0   = andnps(mskxy, onef);
    m128 xy1   = andps( mskxy, onef);
    m128 xz0   = andnps(mskxz, onef);
    m128 xz1   = andps( mskxz, onef);
    m128 yz0   = andnps(mskyz, onef);
    m128 yz1   = andps( mskyz, onef);
    m128 xx1   = andps(xy1, xz1); // xx1 = ( (dx0 >= dy0) &  (dx0 >= dz0)) & 1.0f
    m128 xx2   = orps( xy1, xz1); // xx2 = ( (dx0 >= dy0) |  (dx0 >= dz0)) & 1.0f
    m128 yy1   = andps(xy0, yz1); // yy1 = (~(dx0 >= dy0) &  (dy0 >= dz0)) & 1.0f
    m128 yy2   = orps( xy0, yz1); // yy2 = (~(dx0 >= dy0) |  (dy0 >= dz0)) & 1.0f
    m128 zz1   = andps(xz0, yz0); // zz1 = (~(dx0 >= dz0) & ~(dy0 >= dz0)) & 1.0f
    m128 zz2   = orps( xz0, yz0); // zz2 = (~(dx0 >= dz0) | ~(dy0 >= dz0)) & 1.0f
    dx1        = subps(dx1, xx1);
    dx2        = subps(dx2, xx2);
    dy1        = subps(dy1, yy1);
    dy2        = subps(dy2, yy2);
    dz1        = subps(dz1, zz1);
    dz2        = subps(dz2, zz2);

    // hjktmp = hash(seed, (dy0 >= dz0) ? (j1 + hk0) : (j0 + hk1))
    u = paddd(j1, hk0);
    v = paddd(j0, hk1);
    u = pand( m128_to_m128i(mskyz), u);
    v = pandn(m128_to_m128i(mskyz), v);
    m128i hjktmp = por(u, v);
    hjktmp = hash(seed, hjktmp);

    // h1 = hash(seed, ((dx0 >= dy0) & (dx0 >= dz0)) ? (i1 + hjk00) : (i0 + hjktmp))
    u = paddd(i1, hjk00);
    v = paddd(i0, hjktmp);
    r = andps(mskxy, mskxz);
    u = pand( m128_to_m128i(r), u);
    v = pandn(m128_to_m128i(r), v);
    m128i h1 = por(u, v);
    h1 = hash(seed, h1);

    // h2 = hash(seed, ((dx0 >= dy0) | (dx0 >= dz0)) ? (i1 + hjktmp) : (i0 + hjk11))
    u = paddd(i1, hjktmp);
    v = paddd(i0, hjk11);
    r = orps(mskxy, mskxz);
    u = pand( m128_to_m128i(r), u);
    v = pandn(m128_to_m128i(r), v);
    m128i h2 = por(u, v);
    h2 = hash(seed, h2);

    // gradients
    m128 g0 = grad(h0, dx0, dy0, dz0);
    m128 g1 = grad(h1, dx1, dy1, dz1);
    m128 g2 = grad(h2, dx2, dy2, dz2);
    m128 g3 = grad(h3, dx3, dy3, dz3);

    // first corner
    m128 t0, n0;
    r  = mulps(dx0, dx0); // r  = dx0*dx0
    s  = mulps(dy0, dy0); // s  = dy0*dy0
    t  = mulps(dz0, dz0); // t  = dz0*dz0
    t0 = subps(p6, r);    // t0 = 0.6f - dx0*dx0
    s  = addps(s, t);     // s  = dy0*dy0 + dz0*dz0
    t0 = subps(t0, s);    // t0 = 0.6f - dx0*dx0 - dy0*dy0 - dz0*dz0
    r  = m128i_to_m128(psrad<31>(m128_to_m128i(t0))); // r = (t0 < 0)
    t0 = mulps(t0, t0);
    t0 = mulps(t0, t0);
    n0 = mulps(t0, g0);
    n0 = andnps(r, n0); // n0 = (t0 < 0) ? 0 : ...

    // second corner
    m128 t1, n1;
    r  = mulps(dx1, dx1); // r  = dx1*dx1
    s  = mulps(dy1, dy1); // s  = dy1*dy1
    t  = mulps(dz1, dz1); // t  = dz1*dz1
    t1 = subps(p6, r);    // t1 = 0.6f - dx1*dx1
    s  = addps(s, t);     // s  = dy1*dy1 + dz1*dz1
    t1 = subps(t1, s);    // t1 = 0.6f - dx1*dx1 - dy1*dy1 - dz1*dz1
    r  = m128i_to_m128(psrad<31>(m128_to_m128i(t1))); // r = (t1 < 0)
    t1 = mulps(t1, t1);
    t1 = mulps(t1, t1);
    n1 = mulps(t1, g1);
    n1 = andnps(r, n1); // n1 = (t1 < 0) ? 0 : ...

    // third corner
    m128 t2, n2;
    r  = mulps(dx2, dx2); // r  = dx2*dx2
    s  = mulps(dy2, dy2); // s  = dy2*dy2
    t  = mulps(dz2, dz2); // t  = dz2*dz2
    t2 = subps(p6, r);    // t2 = 0.6f - dx2*dx2
    s  = addps(s, t);     // s  = dy2*dy2 + dz2*dz2
    t2 = subps(t2, s);    // t2 = 0.6f - dx2*dx2 - dy2*dy2 - dz2*dz2
    r  = m128i_to_m128(psrad<31>(m128_to_m128i(t2))); // r = (t2 < 0)
    t2 = mulps(t2, t2);
    t2 = mulps(t2, t2);
    n2 = mulps(t2, g2);
    n2 = andnps(r, n2); // n2 = (t2 < 0) ? 0 : ...

    // fourth corner
    m128 t3, n3;
    r  = mulps(dx3, dx3); // r  = dx3*dx3
    s  = mulps(dy3, dy3); // s  = dy3*dy3
    t  = mulps(dz3, dz3); // t  = dz3*dz3
    t3 = subps(p6, r);    // t3 = 0.6f - dx3*dx3
    s  = addps(s, t);     // s  = dy3*dy3 + dz3*dz3
    t3 = subps(t3, s);    // t3 = 0.6f - dx3*dx3 - dy3*dy3 - dz3*dz3
    r  = m128i_to_m128(psrad<31>(m128_to_m128i(t3))); // r = (t3 < 0)
    t3 = mulps(t3, t3);
    t3 = mulps(t3, t3);
    n3 = mulps(t3, g3);
    n3 = andnps(r, n3); // n3 = (t3 < 0) ? 0 : ...

    // add contributions and scale to interval [-1,1];
    n0 = addps(n0, n1);
    n2 = addps(n2, n3);
    m128 n = addps(n0, n2);
    r = broadcast_m128(32.0f);
    return mulps(n, r);
}

// value-simplex noise
inline f32 noise(u32 seed, f32 x) {
    // get hashed seed
    u32 seed2 = hash(seed);

    // combine value and gradient noise
    return 0.5f * (simplex_noise(seed, x) + value_noise(seed2, x));
}
inline m128 noise(u32 seed, m128 x) {
    // get hashed seed
    u32 seed2 = hash(seed);

    // combine value and gradient noise
    m128 g = simplex_noise(seed, x);
    m128 v = value_noise(seed2, x);
    g = addps(g, v);
    v = broadcast_m128(0.5f);
    return mulps(g, v);
}

inline f32 noise(u32 seed, f32 x, f32 y) {
    // get hashed seed
    u32 seed2 = hash(seed);

    // combine value and gradient noise
    return 0.5f * (simplex_noise(seed, x, y) + value_noise(seed2, x, y));
}
inline m128 noise(u32 seed, m128 x, m128 y) {
    // get hashed seed
    u32 seed2 = hash(seed);

    // combine value and gradient noise
    m128 g = simplex_noise(seed, x, y);
    m128 v = value_noise(seed2, x, y);
    g = addps(g, v);
    v = broadcast_m128(0.5f);
    return mulps(g, v);
}

inline f32 noise(u32 seed, f32 x, f32 y, f32 z) {
    // get hashed seed
    u32 seed2 = hash(seed);

    // combine value and gradient noise
    return 0.5f * (simplex_noise(seed, x, y, z) + value_noise(seed2, x, y, z));
}
inline m128 noise(u32 seed, m128 x, m128 y, m128 z) {
    // get hashed seed
    u32 seed2 = hash(seed);

    // combine value and gradient noise
    m128 g = simplex_noise(seed, x, y, z);
    m128 v = value_noise(seed2, x, y, z);
    g = addps(g, v);
    v = broadcast_m128(0.5f);
    return mulps(g, v);
}


// fractal noise
inline f32 fractal_noise(u32 seed, f32 x, i32 n, const f32 *w) {
    f32 r = 0.0f;

    for (i32 i = 0; i < n; ++i) {
        r += noise(seed, x) * w[i];
        x = x + x;
    }

    return r;
}
inline m128 fractal_noise(u32 seed, m128 x, i32 n, const f32 *w) {
    m128 r = zero_m128();

    for (i32 i = 0; i < n; ++i) {
        m128 s = noise(seed, x);
        s = mulps(s, broadcast_m128(w[i]));
        r = addps(r, s);
        x = addps(x, x);
    }

    return r;
}

inline f32 fractal_noise(u32 seed, f32 x, f32 y, i32 n, const f32 *w) {
    f32 r = 0.0f;

    for (i32 i = 0; i < n; ++i) {
        r += noise(seed, x, y) * w[i];
        x = x + x;
        y = y + y;
    }

    return r;
}
inline m128 fractal_noise(u32 seed, m128 x, m128 y, i32 n, const f32 *w) {
    m128 r = zero_m128();

    for (i32 i = 0; i < n; ++i) {
        m128 s = noise(seed, x, y);
        s = mulps(s, broadcast_m128(w[i]));
        r = addps(r, s);
        x = addps(x, x);
        y = addps(y, y);
    }

    return r;
}

inline f32 fractal_noise(u32 seed, f32 x, f32 y, f32 z, i32 n, const f32 *w) {
    f32 r = 0.0f;

    for (i32 i = 0; i < n; ++i) {
        r += noise(seed, x, y, z) * w[i];
        x = x + x;
        y = y + y;
        z = z + z;
    }

    return r;
}
inline m128 fractal_noise(u32 seed, m128 x, m128 y, m128 z, i32 n, const f32 *w) {
    m128 r = zero_m128();

    for (i32 i = 0; i < n; ++i) {
        m128 s = noise(seed, x, y, z);
        s = mulps(s, broadcast_m128(w[i]));
        r = addps(r, s);
        x = addps(x, x);
        y = addps(y, y);
        z = addps(z, z);
    }

    return r;
}

}

#endif

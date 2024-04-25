#ifndef QUAD_HPP
#define QUAD_HPP

#include <algorithm>
#include <utility>

#include "frustum.hpp"
#include "tree.hpp"
#include "var_to_template.hpp"

namespace vxg {

template<i32 QSIZE>
class quad {
    // number of pixels already set
    i32 npx;

    // signs of corner ray directions
    i32 sgn[3];

    // main axis
    i32 ax;

    // origin
    m128 org0;

    // corner rays (dx, dy, dz)
    m128 d[3];

    // inverse of corner rays (invx, invy, invz)
    m128 inv[3];

    // 2 rows of 4 rays each at the top left (ddx, ddy, ddz)
    m256 dd[3];
    m256 next_x[2];
    m256 next_y[3];

public:
    // pixel data
    alignas(16) u32 px[QSIZE * QSIZE];

    // z-buffer
    alignas(16) f32 zb[QSIZE * QSIZE];

    inline quad(const frustum &fr, u32 wx0, u32 wy0) {
        // tmp
        m128 r, s, t;
        m256 u;

        // reset number of pixels set
        npx = 0;

        // clean up pixels and z-buffer
        std::fill(std::begin(px), std::end(px), 0);
        std::fill(std::begin(zb), std::end(zb), 3e38f);

        // origin
        org0 = fr.org;

        // shift ray to (wx0, wy0)
        r = broadcast_m128((f32)wx0); // r  = wx0
        r = mulps(r, fr.sx);          // r  = wx0 * sx
        s = broadcast_m128((f32)wy0); // s  = wy0
        s = mulps(s, fr.sy);          // s  = wy0 * sy
        r = addps(r, s);              // r  = wx0 * sx + wy0 * sy
        m128 d0 = addps(r, fr.d0);    // d0 = wx0 * sx + wy0 * sy + fr.d0

        // unpack ray
        m128 dx = shufps<0, 0, 0, 0>(d0, d0); // dx = d0x + wx0 * sxx + wy0 * syx
        m128 dy = shufps<1, 1, 1, 1>(d0, d0); // dy = d0y + wx0 * sxy + wy0 * syy
        m128 dz = shufps<2, 2, 2, 2>(d0, d0); // dz = d0z             + wy0 * syz

        // init corner rays
        d[0] = addps(dx, fr.ax); // b.dx = d0x + wx0 * sxx + wy0 * syx + ax
        d[1] = addps(dy, fr.ay); // b.dy = d0y + wx0 * sxy + wy0 * syy + ay
        d[2] = addps(dz, fr.az); // b.dz = d0z             + wy0 * syz + az

        // get inverse of corner rays
        r = broadcast_m128(1.0f);
        inv[0] = divps(r, d[0]);
        inv[1] = divps(r, d[1]);
        inv[2] = divps(r, d[2]);

        // init line of rays
        r = addps(dx, fr.bx);
        s = shufps<0,0,0,0>(fr.sy, fr.sy);
        s = addps(s, r);
        dd[0] = vperm2f128<0,2>(m128_to_m256(r), m128_to_m256(s));
        r = addps(dy, fr.by);
        s = shufps<1,1,1,1>(fr.sy, fr.sy);
        s = addps(s, r);
        dd[1] = vperm2f128<0,2>(m128_to_m256(r), m128_to_m256(s));
        r = dz;
        s = shufps<2,2,2,2>(fr.sy, fr.sy);
        s = addps(s, r);
        dd[2] = vperm2f128<0,2>(m128_to_m256(r), m128_to_m256(s));

        // go to next 4 columns (move in x direction)
        next_x[0] = m128_to_m256(fr.cx);
        next_x[0] = vperm2f128<0,0>(next_x[0], next_x[0]);
        next_x[1] = m128_to_m256(fr.cy);
        next_x[1] = vperm2f128<0,0>(next_x[1], next_x[1]);

        // go to next 2 rows (move in y direction)
        r = addps(fr.sy, fr.sy); // r = 2 * sy
        u = m128_to_m256(r);
        u = vperm2f128<0,0>(u, u);
        next_y[0] = vshufps<0,0,0,0>(u, u);
        next_y[1] = vshufps<1,1,1,1>(u, u);
        next_y[2] = vshufps<2,2,2,2>(u, u);

        // set direction signs
        for (i32 dir = 0; dir < 3; ++dir) {
            i32 msk = movmskps(d[dir]);
            if (msk == 0) {
                sgn[dir] = 1;
            } else if (msk == 0x0F) {
                sgn[dir] = -1;
            } else {
                sgn[dir] = 0;
            }
        }

        // get absolute value of d0
        r = zero_m128();          // r = 0
        r = subps(r, d0);         // r = -d0
        r = andps(r, d0);         // r = abs(d0)
        s = shufps<1,1,1,1>(r,r); // s = d0y, ...
        t = movhlps(r, r);        // t = d0z, ...

        // set main axis
        if (ucomigtss(r, s)) {
            if (ucomigtss(r, t)) ax = 0; // x-axis
            else ax = 2; // z-axis
        } else {
            if (ucomigtss(s, t)) ax = 1; // y-axis
            else ax = 2; // z-axis
        }
    }

    template<i32 ... LS>
    inline void render(const tree<LS...> &tr) {
        // convert axis parameters to template values
        var_to_template<render_callback, 4>::template impl<void>(
            ax, sgn[0], sgn[1], sgn[2], this, tr
        );
    }

private:
    template<class ... TPARAMS>
    struct render_callback {
        template<class R, class ... ARGS>
        static inline R impl(quad *q, ARGS&& ... args) {
            return q->render<TPARAMS::value...>(std::forward<ARGS>(args)...);
        }
    };
    template<i32 AX, i32 SX, i32 SY, i32 SZ, i32 ... LS>
    inline void render(const tree<LS...> &tr) {
        static constexpr i32 LSa[] = { 0, LS... }; // LSa = LS[LVL-1]

        // choose k, u, v axes
        static constexpr i32 IK = AX;
        static constexpr i32 IU = (AX == 0) ? (SZ == 0 ? 1 : 2) :
                                  (AX == 1) ? (SZ == 0 ? 0 : 2) :
                                              (SY == 0 ? 0 : 1);
        static constexpr i32 IV = (AX == 0) ? (SZ == 0 ? 2 : 1) :
                                  (AX == 1) ? (SZ == 0 ? 2 : 0) :
                                              (SY == 0 ? 1 : 0);
        static constexpr i32 SK = (IK == 0) ? SX :
                                  (IK == 1) ? SY :
                                              SZ;
        static constexpr i32 SU = (IU == 0) ? SX :
                                  (IU == 1) ? SY :
                                              SZ;
        static constexpr i32 SV = (IV == 0) ? SX :
                                  (IV == 1) ? SY :
                                              SZ;
        static constexpr i32 NK = SK;
        static constexpr i32 NU = (SU != 0) ? SU : 1;
        static constexpr i32 NV = (SV != 0) ? SV : 1;

        // enable front to back rendering if directions are unambiguous
        static constexpr bool F2B = (SU != 0) && (SV != 0);

        // tmp
        m128i a, b, c;

        // scale origin by 2^{- sum_{i=0}^{sizeof...(LS)-2} LS[i]}; (division is slow)
        a = m128_to_m128i(org0);
        static constexpr i32 sum_LS = (... + LS) - LSa[sizeof...(LS)];
        b = load_m128i(sum_LS<<23, sum_LS<<23, sum_LS<<23, 0); // shift to float exponent
        c = zero_m128i();
        c = pcmpeqd(a, c); // check for zero org
        b = pandn(c, b);   // do not scale zero
        a = psubd(a, b);   // scale org by (neg.) power of 2
        m128 org = m128i_to_m128(a);

        // render root chunk
        render<F2B,IK,IU,IV,NK,NU,NV,SU,std::integral_constant<i32,sizeof...(LS)-1>,LS...>(tr, tr.get_root(), org, std::integral_constant<i32,sizeof...(LS)-1>());
    }

    // render chunk
    template<bool F2B, i32 IK, i32 IU, i32 IV, i32 NK, i32 NU, i32 NV, i32 SU, class LVL_TYPE, i32 ... LS>
    inline void render(const tree<LS...> &tr, i32 chptr, m128 org, LVL_TYPE) {
        static constexpr i32 LSa[] = { 0, LS... }; // LSa = LS[LVL-1]
        static constexpr i32 LVL = LVL_TYPE::value;

        // tmp
        m128i a, b;
        m128 r;

        // get chunk
        auto &ch = tr.get_chunk(chptr, LVL_TYPE());

        // get bounding box of chunk
        m128i min, max;
        ch.get_bounding_box().get_m128i(min, max);

        // check if beam hits chunk at all
        if (!beam_box_intersection(org, min, max)) {
            return;
        }

        // get limits for k
        i32 k0, k1;
        get_k_bounds<IK,NK>(org, min, max, k0, k1);

        // init quad frustum
        m128 dfr, qfr;
        m128i qfr0;
        init_quad_frustum<IK,IU,IV,NK,NU,NV>(org, ch.get_bounding_box().b0, ch.get_bounding_box().b1, k0, dfr, qfr0, qfr);

        // scale origin by 2^{LS[LVL-1]}
        if (LVL > 0) {
            a   = load_m128i(1<<LSa[LVL], 1<<LSa[LVL], 1<<LSa[LVL], 1);
            r   = cvtdq2ps(a);
            org = mulps(org, r);
        }

        // move through k slices
        for (i32 k = k0; k != k1; k += NK) {
            // get uv bounds
            i32 u0, u1, v0, v1;
            get_uv_bounds<NU,NV>(qfr0, qfr, u0, u1, v0, v1);

            // check children in slice
            for (i32 u = u0; u != u1; u += NU) {
                for (i32 v = v0; v != v1; v += NV) {
                    // get coordinates
                    i32 x, y, z;
                    get_xyz<IK,IU,IV>(k, u, v, x, y, z);

                    // load child (either chunk pointer or voxel)
                    auto child = ch.get_data(x, y, z);

                    // check whether child is empty
                    if (child == ch.E) {
                        continue;
                    }

                    // shift origin by - 2^{LS[LVL-1]} * (x, y, z)
                    a = load_m128i(x, y, z, 0);
                    if (LVL > 0) {
                        a = pslld<LSa[LVL]>(a);
                    }
                    r = cvtdq2ps(a);
                    m128 org1 = subps(org, r);

                    // render child
                    render<F2B,IK,IU,IV,NK,NU,NV,SU>(tr, child, org1, std::integral_constant<i32,LVL-1>());
                    if (F2B && (npx == QSIZE * QSIZE)) return;
                }
                if ((SU != 0) && (npx == QSIZE * QSIZE)) return;
            }
            if (npx == QSIZE * QSIZE) return;

            // update quad frustum
            qfr = addps(qfr, dfr);
        }
    }

    // render voxel
    template<bool F2B, i32 IK, i32 IU, i32 IV, i32 NK, i32 NU, i32 NV, i32 SU, i32 ... LS>
    inline void render(const tree<LS...>&, voxel vx, m128 org, std::integral_constant<i32,-1>) {
        // get voxel shape
        const shape &sh = shapes[get_shape_index(vx)];

        // render triangles and rectangles
        for (i32 i = 0; i < sh.tr.size(); ++i) {
            render_face<3,F2B>(org, triangles[sh.tr[i]], 0xFF1F1F1F + ((1+i) << 21) + ((1+i) << 13) + ((1+i) << 5));
            if (F2B && (npx == QSIZE * QSIZE)) return;
        }
        for (i32 i = 0; i < sh.rc.size(); ++i) {
            render_face<4,F2B>(org, rectangles[sh.rc[i]], 0xFF1F1F1F + ((sh.tr.size()+1+i) << 21) + ((sh.tr.size()+1+i) << 13) + ((sh.tr.size()+1+i) << 5));
            if (F2B && (npx == QSIZE * QSIZE)) return;
        }
    }

    // test whether beam can hit axis aligned box (with some false positives)
    bool beam_box_intersection(m128 org, m128i min, m128i max) const {
        // tmp
        m128 r, s;

        // subtract origin from planes
        m128 minf = cvtdq2ps(min);
        m128 maxf = cvtdq2ps(max);
        minf = subps(minf, org);
        maxf = subps(maxf, org);

        // get intersection params
        r = shufps<0,0,0,0>(minf, minf);
        s = shufps<0,0,0,0>(maxf, maxf);
        r = mulps(r, inv[0]);
        s = mulps(s, inv[0]);
        m128 t0x = minps(r, s);
        m128 t1x = maxps(r, s);
        r = shufps<1,1,1,1>(minf, minf);
        s = shufps<1,1,1,1>(maxf, maxf);
        r = mulps(r, inv[1]);
        s = mulps(s, inv[1]);
        m128 t0y = minps(r, s);
        m128 t1y = maxps(r, s);
        r = shufps<2,2,2,2>(minf, minf);
        s = shufps<2,2,2,2>(maxf, maxf);
        r = mulps(r, inv[2]);
        s = mulps(s, inv[2]);
        m128 t0z = minps(r, s);
        m128 t1z = maxps(r, s);

        // get max intersection params (t1 must be > 0)
        m128 t1 = minps(t1x, t1y);
        t1 = minps(t1, t1z);
        if (movmskps(t1) == 0xF) return false;

        // get min intersection params
        m128 t0 = maxps(t0x, t0y);
        t0 = maxps(t0, t0z);

        // check whether t0 is the same for all corner rays (this creates false positives)
        r = cmpeqps(t0, t0x);
        i32 msk = movmskps(r);
        if ((msk != 0) && (msk != 0xF)) return true;
        r = cmpeqps(t0, t0y);
        msk = movmskps(r);
        if ((msk != 0) && (msk != 0xF)) return true;
        r = cmpeqps(t0, t0z);
        msk = movmskps(r);
        if ((msk != 0) && (msk != 0xF)) return true;

        // intersection if t0 <= t1
        r = cmpleps(t0, t1);
        s = m128i_to_m128(psrad<31>(m128_to_m128i(t1)));
        r = andnps(s, r);

        return movmskps(r) != 0;
    }

    template<i32 IK, i32 NK>
    inline void get_k_bounds(m128 org, m128i min, m128i max, i32 &k0, i32 &k1) {
        // get floor(org_k) as integer
        m128 r = (IK == 0) ? org : (IK == 1) ? shufps<1,0,2,3>(org, org) : movhlps(org, org);
        r      = floorss(r, r);
        i32 k2 = cvtss2si(r);

        // get chunk bounds
        i32 kk0, kk1;
        if (IK == 0) {
            kk0 = m128i_to_i32(min);
            kk1 = m128i_to_i32(max);
        } else if (IK == 1) {
            kk0 = pextrd<1>(min);
            kk1 = pextrd<1>(max);
        } else {
            kk0 = pextrd<2>(min);
            kk1 = pextrd<2>(max);
        }

        if (NK > 0) {
            k0 = std::max(k2, kk0);
            k1 = std::max(k2, kk1);
        } else {
            k0 = std::min(k2, kk1 - 1);
            k1 = std::min(k2, kk0 - 1);
        }
    }

    template<i32 IK, i32 IU, i32 IV, i32 NK, i32 NU, i32 NV>
    inline void init_quad_frustum(m128 org, const i32 b0[3], const i32 b1[3], i32 k0, m128 &dfr, m128i &qfr0, m128 &qfr) {
        // tmp
        m128 r, s, t;
        m128i l;

        // get |dk|
        m128 adk = d[IK];
        if (NK < 0) {
            r = zero_m128();     // r   = 0
            adk = subps(r, adk); // adk = - dk
        }

        // get 1 / |dk|
        r = broadcast_m128(1.0f);   // r = 1
        r = divps(r, adk);          // r = 1 / |dk|

        // get dudk, dvdk
        s = mulps(r, d[IU]);          // s = dudk0,dudk1,dudk2,dudk3
        t = mulps(r, d[IV]);          // t = dvdk0,dvdk1,dvdk2,dvdk3
        r = movlhps(s, t);            // r = dudk0,dudk1,dvdk0,dvdk1
        s = movhlps(t, s);            // s = dudk2,dudk3,dvdk2,dvdk3
        t = minps(r, s);              // t = min(u0,u2),min(u1,u3),min(v0,v2),min(v1,v3)
        r = maxps(r, s);              // r = max(u0,u2),max(u1,u3),max(v0,v2),max(v1,v3)
        s = shufps<1,0,3,2>(t, t);    // s = min(u1,u3),min(u0,u2),min(v1,v3),min(v0,v2)
        t = minps(s, t);              // t = min(u0,u1,u2,u3),min(u0,u1,u2,u3),min(v0,v1,v2,v3),min(v0,v1,v2,v3)
        s = shufps<1,0,3,2>(r, r);    // s = max(u1,u3),max(u0,u2),max(v1,v3),max(v0,v2)
        r = maxps(r, s);              // r = max(u0,u1,u2,u3),max(u0,u1,u2,u3),max(v0,v1,v2,v3),max(v0,v1,v2,v3)
        dfr = shufps<0,2,0,2>(t, r);  // dudkmin,dvdkmin,dudkmax,dvdkmax

        // setup qfr0
        qfr0 = load_m128i(b0[IU], b0[IV], b1[IU], b1[IV]);

        // get k coordinate of origin
        m128 org_k = shufps<IK,IK,IK,IK>(org, org);

        // get distance to initial plane
        l = broadcast_m128i(k0); // l = k0
        r = cvtdq2ps(l);         // r = k0
        r = subps(r, org_k);     // r = k0 - org_k

        // get sign(dk) * dfr
        if (NK > 0) {
            s = dfr;
        } else {
            s = zero_m128();
            s = subps(s, dfr);
        }

        // get quad frustum
        r = mulps(r, s);                   // r = sgn(dk) * (k0     - org_k) * (dudkmin,dvdkmin,dudkmax,dvdkmax) = fr0
        s = addps(r, s);                   // s = sgn(dk) * (k0 + 1 - org_k) * (dudkmin,dvdkmin,dudkmax,dvdkmax) = fr1
        t = minps(r, s);                   // t = min(fr0, fr1)
        r = maxps(r, s);                   // r = max(fr0, fr1)
        r = shufps<0,1,2,3>(t, r);         // r = minu,minv,maxu,maxv
        s = shufps<IU,IV,IU,IV>(org, org); // s = org_u,org_v,org_u,org_v
        qfr = addps(r, s);                 // uA,vA,uB,vB
    }

    template<i32 NU, i32 NV>
    inline void get_uv_bounds(m128i qfr0, m128 qfr, i32 &u0, i32 &u1, i32 &v0, i32 &v1) {
        // tmp
        m128 r;
        m128i l, m;

        r = floorps(qfr);        // r = floor(uA,vA,uB,vB)
        l = cvtps2dq(r);         // l = floor(uA,vA,uB,vB)
        m = broadcast_m128i(1);  // m = 1
        m = paddd(m, l);         // m = floor(uA+1,vA+1,uB+1,vB+1)

        // merge with chunk bounds
        m = pminsd(m, qfr0);
        l = pmaxsd(l, qfr0);

        if (NU > 0) {
            u0 = m128i_to_i32(l);
            u1 = pextrd<2>(m);
            u0 = std::min(u0, u1); // u0 <= u1
        } else {
            u1 = m128i_to_i32(l) - 1;
            u0 = pextrd<2>(m) - 1;
            u0 = std::max(u0, u1); // u0 >= u1
        }
        if (NV > 0) {
            v0 = pextrd<1>(l);
            v1 = pextrd<3>(m);
            v0 = std::min(v0, v1); // v0 <= v1
        } else {
            v1 = pextrd<1>(l) - 1;
            v0 = pextrd<3>(m) - 1;
            v0 = std::max(v0, v1); // v0 >= v1
        }
    }

    template<i32 IK, i32 IU, i32 IV>
    inline void get_xyz(i32 k, i32 u, i32 v, i32 &x, i32 &y, i32&z) {
        x = (IK == 0) ? k :
            (IU == 0) ? u :
                        v;
        y = (IK == 1) ? k :
            (IU == 1) ? u :
                        v;
        z = (IK == 2) ? k :
            (IU == 2) ? u :
                        v;
    }

    // render triangle/rectangle, update number of pixels set
    template<i32 N, bool F2B>
    inline void render_face(m128 org, const face<N> &f, u32 col) {
        static_assert((N == 3) || (N == 4), "N must be 3 (triangle) or 4 (rectangle!");

        // temp variables
        m128 l, m;
        m256 r, s, t;
        m128i u;

        // load face plane equation
        u = load_m128i(*(reinterpret_cast<const i32*>(f.n)));
        u = pmovsxbd(u); // sign extend 4 bytes to 4 dwords (SSE 4.1)
        m128 n = cvtdq2ps(u);

        // precompute n . org + b = [nx, ny, nz, b] . [ox, oy, oz, 1]
        m128 n_o_b_128 = dpps<0xFF>(n, org); // dot product (SSE 4.1)

        // check if org is on backside of face
        if (m128_to_f32(n_o_b_128) < 0) return;

        // calculate n . d for corner rays
        l = shufps<0,0,0,0>(n, n);
        l = mulps(d[0], l);
        m = shufps<1,1,1,1>(n, n);
        m = mulps(d[1], m);
        l = addps(l, m);
        m = shufps<2,2,2,2>(n, n);
        m = mulps(d[2], m);
        m128 n_d_128 = addps(l, m);

        // init ray mask: rays must point towards face (n . d < 0)
        u8 msk = movmskps(n_d_128);
        if (msk == 0) return;

        // check edges for corner rays
        for (i32 k = 0; k < N; ++k) {
            u8 mski = check_edge(org, *(const i32*)(f.nn[k]), d[0], d[1], d[2], n_d_128, n_o_b_128);
            // if all corner rays miss on same side => beam does not hit face
            if (mski == 0) return;
            msk &= mski;
        }

        // if all corner rays hit face, no need to check edges for rest of rays
        bool check_rays = (msk != 0xF);

        // expand n . org + b
        r = m128_to_m256(n_o_b_128);
        m256 n_o_b = vperm2f128<0,0>(r,r);

        // unpack face normal vector
        r = m128_to_m256(n);
        r = vperm2f128<0,0>(r,r);
        m256 nx = vshufps<0,0,0,0>(r, r);
        m256 ny = vshufps<1,1,1,1>(r, r);
        m256 nz = vshufps<2,2,2,2>(r, r);

        // get 2 lines of 4 rays each in the top left corner
        m256 dx0 = dd[0];
        m256 dy0 = dd[1];
        m256 dz0 = dd[2];

        // loop over rows of pixels
        i32 i = 0;
        goto for_i_0;
        for (; i < QSIZE; i += 2) {
            // go to next 2 rows
            dx0 = vaddps(dx0, next_y[0]);
            dy0 = vaddps(dy0, next_y[1]);
            dz0 = vaddps(dz0, next_y[2]);

        for_i_0:
            // start with leftmost rays in this row
            m256 dx = dx0;
            m256 dy = dy0;
            m256 dz = dz0;

            // loop over columns (2 rows with 4 rays each at a time)
            i32 j = 0;
            goto for_j_0;
            for (; j < QSIZE; j += 4) {
                // go to next 4 columns (move in x direction, z does not change here)
                dx = vaddps(dx, next_x[0]);
                dy = vaddps(dy, next_x[1]);

            for_j_0:
                // init ray mask for 8 rays
                msk = 0xFF;

                // load old pixel data and get old pixel mask
                l = load_m128_aligned(reinterpret_cast<const f32*>(px +  i     * QSIZE + j));
                m = load_m128_aligned(reinterpret_cast<const f32*>(px + (i + 1)* QSIZE + j));
                i32 msk0 = movmskps(l) | (movmskps(m) << 4);

                if (F2B) {
                    // update ray mask (front-to-back: if a pixel is already filled, do not update it)
                    msk &= (~msk0);
                    if (msk == 0) continue;
                }

                // calculate n . d
                r = vmulps(dx, nx);
                s = vmulps(dy, ny);
                t = vmulps(dz, nz);
                r = vaddps(r, s);
                m256 n_d = vaddps(r, t);

                // check rays if necessary
                if (check_rays) {
                    // update ray mask: rays must point towards face (n . d < 0)
                    msk &= vmovmskps(n_d);
                    if (msk == 0) continue;

                    // check edges
                    for (i32 k = 0; k < N; ++k) {
                        msk &= check_edge(org, *(const i32*)(f.nn[k]), dx, dy, dz, n_d, n_o_b);
                        if (msk == 0) goto for_j_continue;
                    }
                }

                // get length of d vector = sqrt(d . d)
                r = vmulps(dx, dx);
                s = vmulps(dy, dy);
                t = vmulps(dz, dz);
                r = vaddps(r, s);
                r = vaddps(r, t);
                r = vsqrtps(r);

                // get distance to plane (zb = - (n . o + b) / (n . d) * sqrt(d . d))
                t = zero_m256();
                s = vshufps<0,0,0,0>(n_o_b, n_o_b);
                s = vsubps(t, s);
                s = vdivps(s, n_d);
                r = vmulps(r, s);

                if (!F2B) {
                    // load old z buffer values
                    s = m128_to_m256(load_m128_aligned(zb +  i      * QSIZE + j));
                    t = m128_to_m256(load_m128_aligned(zb + (i + 1) * QSIZE + j));
                    s = vperm2f128<0,2>(s, t);

                    // get z buffer mask: zb_new < zb_old
                    s = vcmpps<_CMP_LT_OQ>(r, s);

                    // update ray mask
                    msk &= vmovmskps(s);
                    if (msk == 0) continue;
                }

                if (msk == 0xFF) {
                    // fill all 8 pixels
                    for (i32 k = 0; k < 4; ++k) {
                        px[ i      * QSIZE + j + k] = col;
                        px[(i + 1) * QSIZE + j + k] = col;
                    }
                    store_m128_aligned(zb +  i      * QSIZE + j, m256_to_m128(r));
                    store_m128_aligned(zb + (i + 1) * QSIZE + j, vextractf128<1>(r));
                } else {
                    // test each ray individually
                    l = m256_to_m128(r);
                    for (i32 k = 0; k < 2; ++k) {
                        if (msk & (1 << (4 * k))) {
                            px[(i + k) * QSIZE + j    ] = col;
                            zb[(i + k) * QSIZE + j    ] = m128_to_f32(l);
                        }
                        if (msk & (1 << (4 * k + 1))) {
                            px[(i + k) * QSIZE + j + 1] = col;
                            zb[(i + k) * QSIZE + j + 1] = m128_to_f32(shufps<1,1,1,1>(l,l));
                        }
                        if (msk & (1 << (4 * k + 2))) {
                            px[(i + k) * QSIZE + j + 2] = col;
                            zb[(i + k) * QSIZE + j + 2] = m128_to_f32(movhlps(l,l));
                        }
                        if (msk & (1 << (4 * k + 3))) {
                            px[(i + k) * QSIZE + j + 3] = col;
                            zb[(i + k) * QSIZE + j + 3] = m128_to_f32(shufps<3,3,3,3>(l,l));
                        }
                        l = vextractf128<1>(r);
                    }
                }

                // update number of pixels set
                npx += __builtin_popcount(msk);
                if (F2B && (npx == QSIZE * QSIZE)) return;

                // do not count pixels that were previously set
                if (!F2B) {
                    npx -= __builtin_popcount(msk & msk0);
                }

            for_j_continue:;
            }
        }
    }

    inline i32 check_edge(m128 org, i32 ed, m256 dx, m256 dy, m256 dz, m256 n_d, m256 n_o_b) const {
        // tmp
        m256 r, s, t;
        m128 l;
        m128i u;

        // expand edge plane equation from the 4 bytes in ed
        u = load_m128i(ed);
        u = pmovsxbd(u); // sign extend 4 bytes to 4 dwords (SSE 4.1)
        m128 nn = cvtdq2ps(u);

        // compute nn . org + bb = [nnx, nny, nnz, bb] . [ox, oy, oz, 1]
        l = dpps<0xFF>(nn, org); // dot product (SSE 4.1)
        r = m128_to_m256(l);
        m256 nn_o_bb = vperm2f128<0,0>(r,r); // broadcast to all 8 values

        // compute nn . d for all eight rays
        t = m128_to_m256(nn);       // t = [nnx, nny, nnz,  bb,   0,   0,   0,   0]
        t = vperm2f128<0,0>(t, t);  // t = [nnx, nny, nnz,  bb, nnx, nny, nnz,  bb]
        r = vshufps<0,0,0,0>(t, t); // r = [nnx, nnx, nnx, nnx, nnx, nnx, nnx, nnx]
        s = vshufps<1,1,1,1>(t, t); // s = [nny, nny, nny, nny, nny, nny, nny, nny]
        t = vshufps<2,2,2,2>(t, t); // t = [nnz, nnz, nnz, nnz, nnz, nnz, nnz, nnz]
        r = vmulps(r, dx);
        s = vmulps(s, dy);
        t = vmulps(t, dz);
        r = vaddps(r, s);
        m256 nn_d = vaddps(r, t);

        // get mask for (nn . org + bb) * (n . d) < (n . org + b) * (nn . d)
        r = vmulps(nn_o_bb, n_d);
        s = vmulps(n_o_b, nn_d);
        r = vcmpps<_CMP_LT_OQ>(r, s);

        // return mask
        return vmovmskps(r);
    }

    inline i32 check_edge(m128 org, i32 ed, m128 dx, m128 dy, m128 dz, m128 n_d, m128 n_o_b) const {
        // tmp
        m128 r, s, t;
        m128i u;

        // expand edge plane equation from the 4 bytes in ed
        u = load_m128i(ed);
        u = pmovsxbd(u); // sign extend 4 bytes to 4 dwords (SSE 4.1)
        m128 nn = cvtdq2ps(u);

        // compute nn . org + bb = [nnx, nny, nnz, bb] . [ox, oy, oz, 1]
        m128 nn_o_bb = dpps<0xFF>(nn, org); // dot product (SSE 4.1)

        // compute nn . d for all four rays
        r = shufps<0, 0, 0, 0>(nn, nn); // r = nnx
        s = shufps<1, 1, 1, 1>(nn, nn); // s = nny
        t = shufps<2, 2, 2, 2>(nn, nn); // t = nnz
        r = mulps(r, dx);
        s = mulps(s, dy);
        t = mulps(t, dz);
        r = addps(r, s);
        m128 nn_d = addps(r, t);

        // get mask for (nn . org + bb) * (n . d) < (n . org + b) * (nn . d)
        r = mulps(nn_o_bb, n_d);
        s = mulps(n_o_b, nn_d);
        r = cmpltps(r, s);

        // return mask
        return movmskps(r);
    }

};

}

#endif

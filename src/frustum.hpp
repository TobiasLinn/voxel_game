#ifndef FRUSTUM_HPP
#define FRUSTUM_HPP

#include "camera.hpp"

namespace vxg {

struct frustum {
    m128 org; //    ox,    oy,    oz,           1

    m128 d0;  //   d0x,   d0y,   d0z,           0

    m128 sx;  //   sxx,   sxy,     0,           0
    m128 sy;  //   syx,   syy,   syz,           0

    m128 ax;  //     0, 3*sxx, 3*syx, 3*(sxx+syx)
    m128 ay;  //     0, 3*sxy, 3*syy, 3*(sxy+syy)
    m128 az;  //     0,     0, 3*syz, 3*(    syz)

    m128 bx;  //     0,   sxx, 2*sxx, 3*(sxx    )
    m128 by;  //     0,   sxy, 2*sxy, 3*(sxy    )

    m128 cx;  // 4*sxx, 4*sxx, 4*sxx, 4*(sxx    )
    m128 cy;  // 4*sxy, 4*sxy, 4*sxy, 4*(sxy    )

    // update values depending on camera position and angles
    inline void update(const camera & cam, i32 resx, i32 resy, i32 qsize) {
        const m128 zeros = broadcast_m128(0.0f);
        const m128 halfs = broadcast_m128(0.5f);
        const m128 ones  = broadcast_m128(1.0f);

        m128 r, s, sca;
        f32 fov, inv_resx_2;

        // load position and angles from camera
        fov = cam.get_fov();
        org = load_m128(cam.get_x(), cam.get_y(), cam.get_z(), 1.0f);
        sca = load_m128(cam.get_sin_theta(), cam.get_cos_theta(), cam.get_sin_phi(), cam.get_cos_phi());

        // sx, sy
        inv_resx_2 = 1.0f / ((f32)(resx >> 1));
        r  = broadcast_m128(inv_resx_2);     // r  = 2 / RESX
        r  = mulps(r, sca);                  // r  = 2 / RESX * ( sin_th, cos_th,  sin_ph, cos_ph)
        r  = shufps<0, 1, 3, 2>(r, r);       // r  = 2 / RESX * ( sin_th, cos_th,  cos_ph, sin_ph)
        r  = addsubps(zeros, r);             // r  = 2 / RESX * (-sin_th, cos_th, -cos_ph, sin_ph)
        sx = shufps<3, 2, 0, 0>(r, zeros);   // sx = 2 / RESX * (sin_ph, -cos_ph, 0, 0)
        
        s  = shufps<1, 1, 1, 1>(r, r);       // s  = 2 / RESX * cos_th
        r  = shufps<0, 0, 0, 0>(r, zeros);   // r  = 2 / RESX * (-sin_th, -sin_th, 0, 0)
        s  = mulps(s, sca);                  // s  = 2 / RESX * cos_th * (sin_th, cos_th, sin_ph, cos_ph)
        sy = shufps<3, 2, 0, 2>(s, r);       // sy = 2 / RESX * (cos_th * cos_ph, cos_th * sin_ph, -sin_th, 0)

        // h
        r  = load_m128(fov,fov,0.0f,0.0f);   // r  = (fov, fov, 0, 0)
        r  = mulps(r, sca);                  // r  = fov * (sin_th, cos_th, 0, 0)
        r  = shufps<0, 0, 1, 2>(r, r);       // r  = fov * (sin_th, sin_th, cos_th, 0)
        s  = shufps<3, 2, 0, 0>(sca, ones);  // s  = (cos_ph, sin_ph, 1.0f, 1.0f)
        d0 = mulps(r, s);                    // d0 = fov * (sin_th * cos_ph, sin_th * sin_ph, cos_th, 0) = h

        // h - (RESX - 1) / 2 * sx
        r  = broadcast_m128((f32)(resx - 1));
        r  = mulps(r, halfs);                // r = (RESX - 1) / 2
        r  = mulps(r, sx);                   // r  = (RESX - 1) / 2 * sx
        d0 = subps(d0, r);                   // d0 = h - (RESX - 1) / 2 * sx

        // h - (RESX - 1) / 2 * sx - (RESY - 1) / 2 * sy
        r  = broadcast_m128((f32)(resy - 1));
        r  = mulps(r, halfs);                // r  = (RESY - 1) / 2
        r  = mulps(r, sy);                   // r  = (RESY - 1) / 2 * sy
        d0 = subps(d0, r);                   // d0 = h - (RESX - 1) / 2 * sx - (RESY - 1) / 2 * sy

        // ax, ay, az => create quad corner rays
        r  = shufps<0, 2, 1, 2>(sx, sx);     // r  = (    sxx,   0,     sxy,       0)
        s  = shufps<0, 0, 1, 1>(sy, sy);     // s  = (    syx, syx,     syy,     syy)
        s  = addps(s, r);                    // s  = (sxx+syx, syx, sxy+syy,     syy)
        ax = shufps<1, 0, 1, 0>(r, s);       // ax = (      0, sxx,     syx, sxx+syx)
        ay = shufps<3, 2, 3, 2>(r, s);       // ay = (      0, sxy,     syy, sxy+syy)
        az = shufps<0, 0, 2, 2>(zeros, sy);  // az = (      0,   0,     syz,     syz)
        r  = broadcast_m128((f32)(qsize-1)); // r  = qsize - 1
        ax = mulps(ax, r);                   // ax = (      0, sxx,     syx, sxx+syx) * (qsize - 1)
        ay = mulps(ay, r);                   // ay = (      0, sxy,     syy, sxy+syy) * (qsize - 1)
        az = mulps(az, r);                   // az = (      0,   0,     syz,     syz) * (qsize - 1)

        // bx, by => create line of rays
        r  = shufps<2, 0, 2, 0>(sx, sx);     // r  = (      0, sxx,       0,     sxx)
        s  = shufps<2, 2, 0, 0>(sx, sx);     // s  = (      0,   0,     sxx,     sxx)
        r  = addps(r, s);                    // r  = (      0, sxx,     sxx,   2*sxx)
        bx = addps(r, s);                    // bx = (      0, sxx,   2*sxx,   3*sxx)
        r  = shufps<2, 1, 2, 1>(sx, sx);     // r  = (      0, sxy,       0,     sxy)
        s  = shufps<2, 2, 1, 1>(sx, sx);     // s  = (      0,   0,     sxy,     sxy)
        r  = addps(r, s);                    // r  = (      0, sxy,     sxy,   2*sxy)
        by = addps(r, s);                    // by = (      0, sxy,   2*sxy,   3*sxy)

        // cx, cy => go to next 4 rays
        r  = shufps<0, 0, 0, 0>(sx, sx);     // r  = (  sxx,   sxx,   sxx,   sxx)
        s  = broadcast_m128(4.0f);           // s  = 4
        cx = mulps(r, s);                    // cx = (4*sxx, 4*sxx, 4*sxx, 4*sxx)
        r  = shufps<1, 1, 1, 1>(sx, sx);     // r  = (  sxy,   sxy,   sxy,   sxy)
        cy = mulps(r, s);                    // cy = (4*sxy, 4*sxy, 4*sxy, 4*sxy)
    }

};

}

#endif

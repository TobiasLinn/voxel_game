#ifndef HEIGHT_MAP_HPP
#define HEIGHT_MAP_HPP

#include "noise.hpp"

namespace vxg {

class height_map {
    i32 nx;
    i32 ny;

    std::vector<f32> w;
    std::vector<i32> h;

public:
    inline height_map(i32 nnx, i32 nny, i32 noct, f32 foct)
        : nx(nnx), ny(nny), w(noct), h(nnx*nny) {

        // generate octave weights
        w[0] = 1.0f;
        f32 ws = w[0];
        for (i32 i = 1; i < noct; ++i) {
            w[i] = foct * w[i-1];
            ws += w[i];
        }
        for (i32 i = 0; i < noct; ++i) {
            w[i] /= ws;
        }
    }

    inline void generate(u32 seed, f32 x0, f32 x1, f32 y0, f32 y1, f32 z0, f32 z1) {
        std::vector<f32> x(nx), y(ny);
        linspace(x0, x1, x);
        linspace(y0, y1, y);

        for (i32 j = 0; j < ny; ++j) {
            for (i32 i = 0; i < nx; i = i + 4) {
                m128 xx = load_m128(x.data() + i);
                m128 yy = broadcast_m128(y[j]);
                m128 hh = z0 + 0.5f * (z1 - z0) * (fractal_noise(seed, xx, yy, w.size(), w.data()) + 1.0f);
                m128i hhi = cvttps2dq(floorps(hh));
                if (i + 4 > nx) {
                    for (i32 k = i; k < nx; ++k) {
                        h[j * nx + k] = hhi[k-i];
                    }
                } else {
                    store_m128i(h.data() + j * nx + i, hhi);
                }
            }
        }
    }

    inline i32 get(i32 ix, i32 iy) const {
        return h[iy * nx + ix];
    }
};

}

#endif

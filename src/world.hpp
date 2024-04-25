#ifndef WORLD_HPP
#define WORLD_HPP

#include <ctime>

#include "height_map.hpp"
#include "tree.hpp"

namespace vxg {

class world {
public:
    // world seed
    u32 seed;

    // min, max
    f32 zmin = 0.0f;
    f32 zmax = 256.0f;

    // voxel data in tree format
    tree<4, 6> tr;


    inline world()
        : world(hash(std::time(nullptr))) {
    }
    inline world(u32 s)
        : seed(s) {
        generate();
    }

private:
    inline void generate() {
        height_map hm(17, 17, 8, 0.35f);
        f32 fs = 1.0f / 256.0f;
        std::get<0>(tr.ichunks).resize(1);
        for (i32 y = 0; y < 64; ++y) {
            for (i32 x = 0; x < 64; ++x) {
                // generate heightmap
                hm.generate(seed, x*16*fs, ((x+1)*16)*fs, y*16*fs, ((y+1)*16)*fs, zmin, zmax);

                for (i32 z = 0; z < 64; ++z) {
                    leaf_chunk<4> ch;
                    bool empty;
                    generate(hm, z*16, ch, empty);

                    if (empty) {
                        break;
                    }

                    tr.lchunks.push_back(ch);
                    std::get<0>(tr.ichunks)[0].set_data(x, y, z, tr.lchunks.size() - 1);
                }
            }
        }
    }

    template<i32 LS>
    inline void generate(const height_map &hm, i32 z0, leaf_chunk<LS> &ch, bool &empty) {
        empty = true;

        // set voxels
        for (i32 y = 0; y < ch.S; ++y) {
            for (i32 x = 0; x < ch.S; ++x) {
                // get min and max of 4 height values (integers)
                i32 hmin = 0x7FFFFFFF;
                i32 hmax = 0x80000000;
                for (i32 yy = 0; yy < 2; ++yy) {
                    for (i32 xx = 0; xx < 2; ++xx) {
                        i32 h = hm.get(x + xx, y + yy) - z0;
                        if (h < hmin) hmin = h;
                        if (h > hmax) hmax = h;
                    }
                }

                // limit to chunk size
                if (hmin < 0) hmin = 0;
                if (hmin > ch.S) hmin = ch.S;
                if (hmax < 0) hmax = 0;
                if (hmax > ch.S) hmax = ch.S;

                // set full voxels up to hmin
                if (hmin > 0) {
                    empty = false;
                    for (i32 z = 0; z < hmin; ++z) {
                        ch.set_data(x, y, z, get_voxel(0xFF, 0x00FFFFFF));
                    }
                }

                // set partial voxels
                if (hmax > hmin) {
                    for (i32 z = hmin; z < hmax; ++z) {
                        // get shape
                        i32 sh = 0;
                        i32 bit = 1;
                        for (i32 zz = 0; zz < 2; ++zz) {
                            for (i32 yy = 0; yy < 2; ++yy) {
                                for (i32 xx = 0; xx < 2; ++xx) {
                                    i32 h = hm.get(x + xx, y + yy) - z0;

                                    if (z + zz <= h) {
                                        sh |= bit;
                                    }
                                    bit = bit << 1;
                                }
                            }
                        }

                        // smooth terrain
                        if (sh == 0x1F) {
                            sh = 0x17;
                        } else if (sh == 0x2F) {
                            sh = 0x2B;
                        } else if (sh == 0x4F) {
                            sh = 0x4D;
                        } else if (sh == 0x8F) {
                            sh = 0x8E;
                        }

                        // set voxel
                        if ((sh != 0) && (shapes[sh].tr.size() + shapes[sh].rc.size() >= 4)) {
                            empty = false;
                            ch.set_data(x, y, z, get_voxel(sh, 0x00FFFFFF));
                        }
                    }
                }
            }
        }
    }
};

}

#endif

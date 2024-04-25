#ifndef RENDERER_HPP
#define RENDERER_HPP

#include <string>
#include <thread>
#include <vector>
#include <SDL2/SDL.h>

#include "camera.hpp"
#include "frustum.hpp"
#include "quad.hpp"
#include "world.hpp"

namespace vxg {

class renderer {
    // resolution
    u32 resx;
    u32 resy;

    // quad size (must be divisible by 4)
    static constexpr i32 QSIZE = 4;

    // pixels and z buffer
    std::vector<u32> px;
    std::vector<f32> zb;    

    // SDL related objects
    SDL_Renderer *sdl_renderer = nullptr;
    SDL_Texture *sdl_texture = nullptr;
    
    // frustum
    frustum fr;

public:
    // constructor
    inline renderer() {
    }
    inline renderer(SDL_Window *wnd, u32 rx, u32 ry, const camera &cam) :
        resx(rx),
        resy(ry),
        px(rx * ry),
        zb(rx * ry) {

        // create SDL renderer
        sdl_renderer = SDL_CreateRenderer(wnd, -1, 0);//SDL_RENDERER_PRESENTVSYNC);

        // set logical resolution, texture will be scaled to fit window resolution
        SDL_RenderSetLogicalSize(sdl_renderer, resx, resy);

        // setup SDL texture to allow software rendering
        sdl_texture = SDL_CreateTexture(sdl_renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, resx, resy);
        
        // update frustum
        update_frustum(cam);
    }

    // delete copy constructor to prevent double usage of SDL
    inline renderer(const renderer&) = delete;

    // destructor
    inline ~renderer() {
        // destroy SDL related objects
        if (sdl_texture) {
            SDL_DestroyTexture(sdl_texture);
            sdl_texture = nullptr;
        }
        if (sdl_renderer) {
            SDL_DestroyRenderer(sdl_renderer);
            sdl_renderer = nullptr;
        }
    }

    inline void render(const world &wrld) {
        // render quads in parallel
        #pragma omp parallel for schedule(dynamic) collapse(2)
        for (u32 wy0 = 0; wy0 <= resy - QSIZE; wy0 += QSIZE) {
            for (u32 wx0 = 0; wx0 <= resx - QSIZE; wx0 += QSIZE) {
                // render quad
                quad<QSIZE> q(fr, wx0, wy0);
                q.render(wrld.tr);

                // copy quad pixel data to px
                for (u32 i = 0; i < QSIZE; ++i) {
                    for (u32 j = 0; j < QSIZE; ++j) {
                        px[(wy0 + i) * resx + wx0 + j] = q.px[i * QSIZE + j];
                    }
                }
            }
        }

        // upload pixel data to gpu
        SDL_UpdateTexture(sdl_texture, nullptr, px.data(), resx * 4);

        // clear screen
        SDL_RenderClear(sdl_renderer);

        // render texture to framebuffer, scale it to fit window resolution
        SDL_RenderCopy(sdl_renderer, sdl_texture, 0, 0);

        // show new framebuffer
        SDL_RenderPresent(sdl_renderer);
    }

    inline void update_frustum(const camera &cam) {
        fr.update(cam, resx, resy, QSIZE);
    }

};

}

#endif

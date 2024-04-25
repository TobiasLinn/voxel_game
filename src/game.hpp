#ifndef GAME_HPP
#define GAME_HPP

#include <memory>

#include <SDL2/SDL.h>
#include <omp.h>

#include "noise.hpp"
#include "renderer.hpp"

namespace vxg {

class game {
    // run flag
    bool running = false;

    // window resolution
    u32 wndresx = 1600;
    u32 wndresy = 960;

    // length of each frame in ms
    u32 framelength = 16;

    // move and turn speed
    f64 movespeed = 0.1;
    f64 fastspeed = 0.3;
    f64 turnspeed = 0.025;

    // mouse turning
    bool turnmouse = false;

    // render threads
    i32 nthreads = 8;

    // fps related
    u32 frames0 = 0;
    u32 frames = 0;
    u32 fps = 0;

    // window and context
    SDL_Window *wnd = nullptr;

    // camera
    camera cam;

    // renderer
    renderer rnd;

    // world
    world wrld;

public:
    // constructor
    inline game()
        : wrld(12342134) {
        // init sdl
        if (SDL_Init(SDL_INIT_EVERYTHING) < 0) {
            throw "game::game(): Failed to init SDL!";
        }

        // create window
        wnd = SDL_CreateWindow("", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, wndresx, wndresy, SDL_WINDOW_OPENGL);
        if (!wnd) {
            throw "game::game(): Unable to create window!";
        }

        // init camera
        cam = camera(-0.5f, -0.5f, 128.0f, 0.75f * M_PI, 0.25f * M_PI, 0.5f * M_PI);

        // init renderer in place
        new(&rnd) renderer(wnd, wndresx, wndresy, cam);

        // set number of render threads
        omp_set_num_threads(nthreads);
    }

    // destructor
    inline ~game() {
        // delete sdl_window
        SDL_DestroyWindow(wnd);

        // quit sdl
        SDL_Quit();
    }

    // run the game
    inline void run() {
        running = true;

        SDL_Event event;
        u32 ticks = SDL_GetTicks();

        if (turnmouse) {
            SDL_SetRelativeMouseMode(SDL_TRUE);
        }

        // game loop
        while (running) {
            // handle events
            while (SDL_PollEvent(&event)) {
                switch (event.type) {
                case SDL_KEYDOWN:
                    process_key_event(event);
                    break;
                case SDL_QUIT:
                    return;
                }
            }

            // game update
            for (int a = 0; (ticks + framelength < SDL_GetTicks()) && (a < 10); ++a) {
                ticks += framelength;
                update();
            }

            // render
            rnd.render(wrld);

            // calculate frames per second
            update_fps();
        }
    }

private:
    // handle keyboard
    inline void process_key_event(const SDL_Event &e) {
        switch(e.key.keysym.sym) {
        case SDLK_ESCAPE:
            running = false;
            break;
        case SDLK_i:
            turnmouse = !turnmouse;
            SDL_SetRelativeMouseMode(turnmouse ? SDL_TRUE : SDL_FALSE);
            break;
        }
    }

    // update game state
    inline void update() {
        auto * keystate = SDL_GetKeyboardState(0);

        // moving
        f32 speed = movespeed;
        if (keystate[SDL_SCANCODE_LSHIFT]) {
            speed = fastspeed;
        }
        if (keystate[SDL_SCANCODE_W]) {
            cam.move_fb(speed);
        }
        if (keystate[SDL_SCANCODE_A]) {
            cam.move_lr(speed);
        }
        if (keystate[SDL_SCANCODE_S]) {
            cam.move_fb(-speed);
        }
        if (keystate[SDL_SCANCODE_D]) {
            cam.move_lr(-speed);
        }
        if (keystate[SDL_SCANCODE_SPACE]) {
            cam.move_ud(speed);
        }
        if (keystate[SDL_SCANCODE_LCTRL]) {
            cam.move_ud(-speed);
        }

        // turning
        if (turnmouse) {
            int dx = 0;
            int dy = 0;
            SDL_GetRelativeMouseState(&dx, &dy);

            if (dx != 0) {
                cam.turn_lr(-dx * turnspeed * 0.1f);
            }
            if (dy != 0) {
                cam.turn_ud(dy * turnspeed * 0.1f);
            }
        } else {
            if (keystate[SDL_SCANCODE_LEFT]) {
                cam.turn_lr(turnspeed);
            }
            if (keystate[SDL_SCANCODE_RIGHT]) {
                cam.turn_lr(-turnspeed);
            }
            if (keystate[SDL_SCANCODE_UP]) {
                cam.turn_ud(-turnspeed);
            }
            if (keystate[SDL_SCANCODE_DOWN]) {
                cam.turn_ud(turnspeed);
            }
        }

        // update camera
        bool m = cam.was_moved();
        bool t = cam.was_turned();
        if (m) {
            cam.finish_move();
        }
        if (t) {
            cam.finish_turn();
        }
        if (m || t) {
            rnd.update_frustum(cam);
        }
    }

    // update fps counter
    inline void update_fps() {
        // update frame counter
        ++frames;

        // test if 1 second has passed
        if (frames0 + 1000 <= SDL_GetTicks()) {
            // set fps
            fps = frames;

            // display info in window title
            std::string s = "FPS = " + std::to_string(fps);
            SDL_SetWindowTitle(wnd, s.c_str());

            // reset frame counter
            frames0 = SDL_GetTicks();
            frames = 0;
        }
    }
};

}

#endif

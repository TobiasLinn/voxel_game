#ifndef CAMERA_HPP
#define CAMERA_HPP

#include "util.hpp"

namespace vxg {

class camera {
    // position
    f32 x;
    f32 y;
    f32 z;

    // direction angles
    f32 th;
    f32 ph;

    // field of view
    f32 fov;

    // sine/cosine of angles
    f32 sin_th;
    f32 cos_th;
    f32 sin_ph;
    f32 cos_ph;

    // flags to indicate whether camera was moved or turned
    bool moved;
    bool turned;

public:

    // constructor
    inline camera() {
    }
    inline camera(f32 x, f32 y, f32 z, f32 theta, f32 phi, f32 fov_ang) {
        set_pos(x, y, z);
        set_theta(theta);
        set_phi(phi);
        set_fov(fov_ang);
        
        finish_move();
        finish_turn();
    }

    // xyz position
    inline f32 get_x() const {
        return x;
    }
    inline void set_x(f32 x) {
        this->x = x;
        moved = true;
    }
    inline f32 get_y() const {
        return y;
    }
    inline void set_y(f32 y) {
        this->y = y;
        moved = true;
    }
    inline f32 get_z() const {
        return z;
    }
    inline void set_z(f32 z) {
        this->z = z;
        moved = true;
    }
    inline void get_pos(f32 &x, f32 &y, f32 &z) const {
        x = this->x;
        y = this->y;
        z = this->z;
    }
    inline void set_pos(f32 x, f32 y, f32 z) {
        this->x = x;
        this->y = y;
        this->z = z;
        moved = true;
    }

    // theta angle
    inline f32 get_theta() const {
        return th;
    }
    inline void set_theta(f32 th) {
        if (th < 0) {
            th = 0;
        } else if (th > PI) {
            th = PI;
        }
        this->th = th;
        turned = true;
    }
    inline f32 get_sin_theta() const {
        return sin_th;
    }
    inline f32 get_cos_theta() const {
        return cos_th;
    }

    // phi angle
    inline f32 get_phi() const {
        return ph;
    }
    inline void set_phi(f32 ph) {
        this->ph = ph - TWO_PI * std::floor(ph / TWO_PI);
        turned = true;
    }
    inline f32 get_sin_phi() const {
        return sin_ph;
    }
    inline f32 get_cos_phi() const {
        return cos_ph;
    }

    // field of view
    inline f32 get_fov() const {
        return fov;
    }
    inline void set_fov(f32 fov_ang) {
        fov = 1.0 / std::tan(0.5 * fov_ang);
        turned = true;
    }

    // move
    inline void move_fb(f32 d) {
        x += d * sin_th * cos_ph;
        y += d * sin_th * sin_ph;
        z += d * cos_th;
        moved = true;
    }
    inline void move_lr(f32 d) {
        x -= d * sin_ph;
        y += d * cos_ph;
        moved = true;
    }
    inline void move_ud(f32 d) {
        z += d;
        moved = true;
    }
    inline bool was_moved() const {
        return moved;
    }
    inline void finish_move() {
        moved = false;
    }

    // turn
    inline void turn_lr(f32 a) {
        set_phi(ph + a);
        turned = true;
    }
    inline void turn_ud(f32 a) {
        set_theta(th + a);
        turned = true;
    }
    inline bool was_turned() const {
        return turned;
    }
    inline void finish_turn() {
        sin_th = std::sin(th);
        cos_th = std::cos(th);
        sin_ph = std::sin(ph);
        cos_ph = std::cos(ph);

        turned = false;
    }
};

}

#endif

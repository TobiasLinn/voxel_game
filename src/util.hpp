#ifndef UTIL_HPP
#define UTIL_HPP

#include <cmath>
#include <cstdint>
#include <vector>

#include "simd/simd.hpp"

namespace vxg {

using namespace simd;

// unsigned integers
using u8  = std::uint8_t;
using u16 = std::uint16_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

// signed integers
using i8  = std::int8_t;
using i16 = std::int16_t;
using i32 = std::int32_t;
using i64 = std::int64_t;

// floating point numbers
using f32 = float;
using f64 = double;

// constants
static constexpr f64 PI     = 3.14159265358979323846;
static constexpr f64 TWO_PI = 6.28318530717958647693;

// generate linear spaced values
template<class T>
inline void linspace(T x0, T x1, std::vector<T> &x) {
    i32 n = x.size();
    for (i32 i = 0; i < n; ++i) {
        x[i] = x0 + ((x1 - x0) * i) / (n - 1);
    }
}

}

#endif

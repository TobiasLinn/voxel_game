#ifndef SIMD_SIMD_HPP
#define SIMD_SIMD_HPP

#include "mmx.hpp"
#include "sse.hpp"
#include "sse2.hpp"
#include "sse3.hpp"
#include "ssse3.hpp"
#include "sse4_1.hpp"
#include "sse4_2.hpp"

#ifdef __AVX__
#include "avx.hpp"
#endif

#ifdef __AVX2__
#include "avx2.hpp"
#include "fma.hpp"
#endif

#endif

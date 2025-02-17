# voxel_game
Small CPU-based Voxel Raycasting Demo written in header-only C++ using Multithreading and SSE/AVX intrinsics for low-level optimization.

The main rendering algorithm is taken from <https://dl.acm.org/doi/10.1145/1179352.1141913> and implemented in [quad.hpp](src/quad.hpp).

Voxels are not simply cubes, but can be cut along a fixed number of predefined planes, resulting in many possible shapes.

![screenshot](screenshot.png)

#ifndef RAYTRACE_MANDELBROT_GUARD
#define RAYTRACE_MANDELBROT_GUARD

#include "raytrace/quaternion.h"
// #include "raytrace/gradient_traceable.h"
#include "raytrace/ray_traceable.h"

const real MANDELBROT_BAILOUT = 1e5;
const int MANDELBROT_APPROACH_ITERATIONS = 1024;
const real MANDELBROT_APPROACH_ALPHA = 0.001;
const real MANDELBROT_APPROACH_EPSILON = 1e-6;

// class Mandelbrot: public GradientTraceable {
class Mandelbrot: public RayTraceable {
public:
    int left_order {1};
    int right_order {1};
    int max_iter {10};
    color pigment {0, 1, 1, 1};
    quaternion twist {0, 0, 0, 1};
    quaternion displacement {0, 0, 1, 0};
    real threshold = 0;
    real s_distance(quaternion);
    // quaternion normal(quaternion);
    // quaternion gradient(quaternion);
    color get_color(quaternion, quaternion);
    std::tuple<quaternion, quaternion> trace(quaternion, quaternion);
    // std::tuple<quaternion, quaternion, real> trace_S3(quaternion, quaternion);
};

#endif

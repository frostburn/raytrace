#ifndef RAYTRACE_MANDELBULB_GUARD
#define RAYTRACE_MANDELBULB_GUARD

#include "raytrace/quaternion.h"
// #include "raytrace/gradient_traceable.h"
#include "raytrace/ray_traceable.h"

const real MANDELBULB_BAILOUT = 1e5;
const int MANDELBULB_APPROACH_ITERATIONS = 1024;
const real MANDELBULB_APPROACH_ALPHA = 0.01;
const real MANDELBULB_APPROACH_EPSILON = 1e-1;

// class Mandelbulb: public GradientTraceable {
class Mandelbulb: public RayTraceable {
public:
    int order {6};
    int max_iter {5};
    color pigment {0, 1, 1, 1};
    real threshold {0};
    real s_distance(quaternion);
    // quaternion normal(quaternion);
    // quaternion gradient(quaternion);
    color get_color(quaternion, quaternion);
    std::tuple<quaternion, quaternion> trace(quaternion, quaternion);
    // std::tuple<quaternion, quaternion, real> trace_S3(quaternion, quaternion);
};

#endif

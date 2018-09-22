#ifndef RAYTRACE_SPHERE_GUARD
#define RAYTRACE_SPHERE_GUARD

#include "raytrace/quaternion.h"
#include "raytrace/gradient_traceable.h"

const int SPHERE_APPROACH_ITERATIONS = 128;
const real SPHERE_EPSILON = 1e-12;

class Sphere: public GradientTraceable {
public:
    color pigment {0, 1, 1, 1};
    real s_distance(quaternion);
    quaternion normal(quaternion);
    quaternion gradient(quaternion);
    color get_color(quaternion, quaternion);
    std::tuple<quaternion, quaternion> trace(quaternion, quaternion);
    std::tuple<quaternion, quaternion, real> trace_S3(quaternion, quaternion);
};

#endif

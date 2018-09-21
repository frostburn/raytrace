#ifndef RAYTRACE_SPHERE_GUARD
#define RAYTRACE_SPHERE_GUARD

#include "raytrace/quaternion.h"
#include "raytrace/gradient_traceable.h"

const int SPHERE_NEWTON_WARMUP = 32;
const int SPHERE_NEWTON_ITERATIONS = 16;
const int SPHERE_NEWTON_DESCENDS = 8;
const int SPHERE_NEWTON_SECOND_ROOT = 8;

class Sphere: public GradientTraceable {
public:
    real s_distance(quaternion);
    quaternion normal(quaternion);
    quaternion gradient(quaternion);
    color get_color(quaternion, quaternion);
    std::tuple<quaternion, quaternion> trace(quaternion, quaternion);
    std::tuple<quaternion, quaternion, real> trace_S3(quaternion, quaternion);
};

#endif

#ifndef RAYTRACE_SPHERE_GUARD
#define RAYTRACE_SPHERE_GUARD

#include "raytrace/quaternion.h"
#include "raytrace/gradient_traceable.h"

const int SPHERE_NEWTON_WARMUP = 512;
const int SPHERE_NEWTON_ITERATIONS = 32;

class Sphere: public GradientTraceable {
private:
    real newton(real, real, real, real, real, quaternion, quaternion, quaternion);
public:
    real s_distance(quaternion);
    quaternion normal(quaternion);
    quaternion gradient(quaternion);
    color get_color(quaternion, quaternion);
    std::tuple<quaternion, quaternion> trace(quaternion, quaternion);
    std::tuple<quaternion, quaternion, real> trace_S3(quaternion, quaternion);
};

#endif

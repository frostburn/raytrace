#ifndef RAYTRACE_SPHERE_GUARD
#define RAYTRACE_SPHERE_GUARD

#include "raytrace/quaternion.h"
#include "raytrace/gradient_traceable.h"

const int SPHERE_RAY_MARCH_DIVISIONS = 128;
const int SPHERE_BIJECT_ITERATIONS = 8;
const int SPHERE_NEWTON_ITERATIONS = 0;

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

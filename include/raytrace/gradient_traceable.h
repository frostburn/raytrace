#ifndef RAYTRACE_GRADIENT_TRACEABLE_H_GUARD
#define RAYTRACE_GRADIENT_TRACEABLE_H_GUARD

#include "raytrace/quaternion.h"
#include "raytrace/ray_traceable.h"

const real NEWTON_NUDGE = 0.01;
const int NEWTON_ITERATIONS = 10;
const real NEWTON_EPSILON = 1e-5;

class GradientTraceable: public RayTraceable {
public:
    virtual quaternion gradient(quaternion) = 0;
    virtual std::tuple<quaternion, quaternion, real> trace_S3(quaternion, quaternion);
};

#endif

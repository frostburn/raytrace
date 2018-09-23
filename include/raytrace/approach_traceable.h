#ifndef RAYTRACE_APPROACH_TRACEABLE_H_GUARD
#define RAYTRACE_APPROACH_TRACEABLE_H_GUARD

#include "raytrace/quaternion.h"
#include "raytrace/ray_traceable.h"

const int APPROACH_ITERATIONS = 128;
const real APPROACH_EPSILON = 1e-12;

class ApproachTraceable: public RayTraceable {
public:
    virtual std::tuple<quaternion, quaternion> trace(quaternion, quaternion);
    virtual std::tuple<quaternion, quaternion, real> trace_S3(quaternion, quaternion);
};

#endif

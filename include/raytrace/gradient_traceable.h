#ifndef RAYTRACE_GRADIENT_TRACEABLE_H_GUARD
#define RAYTRACE_GRADIENT_TRACEABLE_H_GUARD

#include "raytrace/quaternion.h"
#include "raytrace/ray_traceable.h"

class GradientTraceable: public RayTraceable {
public:
    virtual quaternion gradient(quaternion) = 0;
    virtual std::tuple<quaternion, quaternion, real> trace_S3(quaternion, quaternion);
};

#endif

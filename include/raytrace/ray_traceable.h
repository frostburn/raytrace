#ifndef RAYTRACE_RAY_TRACEABLE_H_GUARD
#define RAYTRACE_RAY_TRACEABLE_H_GUARD

#include "raytrace/quaternion.h"

class RayTraceable {
public:
    quaternion scale {1, 1, 1, 1};
    quaternion left_transform {1, 0, 0, 0};
    quaternion right_transform {1, 0, 0, 0};
    quaternion location {0, 0, 0, 0};
    bool reflective {false};
    virtual real s_distance(quaternion) = 0;
    virtual quaternion normal(quaternion) = 0;
    virtual color get_color(quaternion, quaternion) = 0;
    virtual std::tuple<quaternion, quaternion> trace(quaternion, quaternion);
    virtual std::tuple<quaternion, quaternion, real> trace_S3(quaternion, quaternion);
    quaternion project(quaternion);
    quaternion inverse_project(quaternion);
};

#endif

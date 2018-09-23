#ifndef RAYTRACE_HYPERTORUS_GUARD
#define RAYTRACE_HYPERTORUS_GUARD

#include "raytrace/quaternion.h"
#include "raytrace/approach_traceable.h"

class HyperTorus: public ApproachTraceable {
public:
    real minor_radius {0.2};
    color pigment_a {0, 1, 1, 1};
    color pigment_b {0, 0.2, 0.2, 0.2};
    real s_distance(quaternion);
    quaternion normal(quaternion);
    quaternion gradient(quaternion);
    color get_color(quaternion, quaternion);
    // std::tuple<quaternion, quaternion> trace(quaternion, quaternion);
    // std::tuple<quaternion, quaternion, real> trace_S3(quaternion, quaternion);
};

#endif

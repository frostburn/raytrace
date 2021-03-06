#ifndef RAYTRACE_CLIFFORD_TORUS_GUARD
#define RAYTRACE_CLIFFORD_TORUS_GUARD

#include "raytrace/quaternion.h"
#include "raytrace/approach_traceable.h"

class CliffordTorus: public ApproachTraceable {
public:
    color pigment_a {0, 1, 1, 1};
    color pigment_b {0, 0.2, 0.2, 0.2};
    color pigment_c {0, 0.5, 0.5, 0.5};
    real s_distance(quaternion);
    quaternion normal(quaternion);
    quaternion gradient(quaternion);
    color get_color(quaternion, quaternion);
};

#endif

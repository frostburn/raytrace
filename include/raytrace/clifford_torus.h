#ifndef RAYTRACE_CLIFFORD_TORUS_GUARD
#define RAYTRACE_CLIFFORD_TORUS_GUARD

#include "raytrace/quaternion.h"
#include "raytrace/gradient_traceable.h"

class CliffordTorus: public GradientTraceable {
public:
    real s_distance(quaternion);
    quaternion normal(quaternion);
    quaternion gradient(quaternion);
    color get_color(quaternion, quaternion);
};

#endif

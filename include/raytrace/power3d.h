#ifndef RAYTRACE_POWER3D_GUARD
#define RAYTRACE_POWER3D_GUARD

#include "raytrace/quaternion.h"

quaternion pow3d(const quaternion& v, const int& n);
// std::tuple<quaternion, quaternion, quaternion, quaternion> pow3d_grad(const quaternion&, const int&);
quaternion* pow3d_grad(const quaternion&, const int&);

#endif

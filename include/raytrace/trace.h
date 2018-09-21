#ifndef RAYTRACE_TRACE_GUARD
#define RAYTRACE_TRACE_GUARD

#include "raytrace/quaternion.h"
#include "raytrace/ray_traceable.h"

color raytrace(quaternion source, quaternion target, int depth, const std::vector<std::shared_ptr<RayTraceable>>& objects);

color raytrace_S3(quaternion source, quaternion target, int depth, const std::vector<std::shared_ptr<RayTraceable>>& objects);

#endif

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <vector>
#include <memory>
#include <random>
#include <cassert>

#include "raytrace/quaternion.h"
#include "raytrace/ray_traceable.h"
#include "raytrace/trace.h"

const real RAYTRACE_NUDGE = 0.01;

color raytrace(quaternion source, quaternion target, int depth, const std::vector<std::shared_ptr<RayTraceable>>& objects) {
    if (depth <= 0) {
        return (color){0, 0.05, 0, 0};
    }
    std::vector<std::shared_ptr<RayTraceable>>::const_iterator obj;
    color result = {0, 0.05, 0, 0};
    quaternion closest_surface = infinity();

    for (obj = objects.begin(); obj != objects.end(); ++obj) {
        auto [surface, ray] = (*obj)->trace(source, target);
        if (norm(surface - source) >= norm(closest_surface - source)) {
            continue;
        }
        closest_surface = surface;
        surface = (*obj)->inverse_project(surface - (*obj)->location);
        ray = (*obj)->inverse_project(ray);
        result = (*obj)->get_color(surface, ray / norm(ray));
        if ((*obj)->reflective) {
            quaternion normal = (*obj)->normal(surface);
            ray = ray - 2 * dot(ray, normal) * normal;
            // surface = (*obj)->project(surface) + (*obj)->location;
            ray = (*obj)->project(ray);
            result = 0.1 * result + 0.9 * raytrace(
                closest_surface + RAYTRACE_NUDGE * ray,
                closest_surface + ray,
                depth - 1,
                objects
            );
        }
    }
    return result;
}

color raytrace_S3(quaternion source, quaternion target, int depth, const std::vector<std::shared_ptr<RayTraceable>>& objects) {
    if (depth <= 0) {
        return (color){0, 0.05, 0, 0};
    }
    std::vector<std::shared_ptr<RayTraceable>>::const_iterator obj;
    color result = {0, 0.05, 0, 0};
    quaternion closest_surface = infinity();
    real closest_distance = std::numeric_limits<real>::infinity();

    for (obj = objects.begin(); obj != objects.end(); ++obj) {
        auto [surface, ray, distance] = (*obj)->trace_S3(source, target);
        if (distance >= closest_distance) {
            continue;
        }
        closest_surface = surface;
        closest_distance = distance;
        surface = (*obj)->inverse_project(surface - (*obj)->location);
        ray = (*obj)->inverse_project(ray);
        result = (*obj)->get_color(surface, ray / norm(ray));
        if ((*obj)->reflective) {
            quaternion normal = (*obj)->normal(surface);
            ray = ray - 2 * dot(ray, normal) * normal;
            // surface = (*obj)->project(surface) + (*obj)->location;
            ray = (*obj)->project(ray);
            source = closest_surface + RAYTRACE_NUDGE * ray;
            target = source + ray;
            source = source / norm(source);
            target = target / norm(target);
            result = 0.1 * result + 0.9 * raytrace_S3(
                source,
                target,
                depth - 1,
                objects
            );
        }
    }
    return result;
}

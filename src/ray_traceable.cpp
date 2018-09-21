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

const real RAY_MARCH_ACCURACY = 0.01;
const int RAY_MARCH_ITERATIONS = 10000;
const int BISECT_ITERATIONS = 100;

quaternion RayTraceable::project(quaternion location) {
    location = (quaternion){
        location.r * this->scale.r,
        location.x * this->scale.x,
        location.y * this->scale.y,
        location.z * this->scale.z,
    };
    return this->left_transform * location * this->right_transform;
}

quaternion RayTraceable::inverse_project(quaternion location) {
    location = (1.0 / this->left_transform) * location / this->right_transform;
    return (quaternion){
        location.r / this->scale.r,
        location.x / this->scale.x,
        location.y / this->scale.y,
        location.z / this->scale.z,
    };
}

std::tuple<quaternion, quaternion> RayTraceable::trace(quaternion source, quaternion target) {
    quaternion direction = target - source;
    direction = direction / norm(direction);
    quaternion dt = direction * RAY_MARCH_ACCURACY;

    target = source + dt;

    source = this->inverse_project(source - this->location);
    target = this->inverse_project(target - this->location);

    real a = this->s_distance(source);
    real b = this->s_distance(target);
    // Ray march
    for (int i = 0; i < RAY_MARCH_ITERATIONS; ++i) {
        if (a*b <= 0) {
            break;
        }
        source = target;
        target = target + dt;
        a = b;
        b = this->s_distance(target);
    }
    if (a*b > 0) {
        return {infinity(), direction};
    }
    // Ray bisect
    for (int i = 0; i < BISECT_ITERATIONS; ++i) {
        quaternion half_way = 0.5 * (source + target);
        real c = this->s_distance(half_way);
        if (a*c >= 0) {
            a = c;
            source = half_way;
        } else if (b*c >= 0) {
            b = c;
            target = half_way;
        } else {
            break;
        }
    }

    return {this->project(0.5 * (source + target)) + this->location, direction};
}

std::tuple<quaternion, quaternion, real> RayTraceable::trace_S3(quaternion source, quaternion target) {
    target = cross_align(source, target);
    // Ray geodesic: source * cos(t) + target * sin(t)

    quaternion sp = this->inverse_project(source);
    quaternion tp = this->inverse_project(target);
    quaternion c = this->inverse_project(-this->location);

    // Ray march
    quaternion v = sp + c;
    real a = this->s_distance(v);
    real b;
    real t = 0;
    real dt = 2 * M_PI / (real) RAY_MARCH_ITERATIONS;
    for (int i = 0; i < RAY_MARCH_ITERATIONS; ++i) {
        t += dt;
        v = sp*cos(t) + tp*sin(t) +  c;
        b = this->s_distance(v);
        if (a*b <= 0) {
            break;
        }
        a = b;
    }
    if (a*b > 0) {
        return {infinity(), infinity(), std::numeric_limits<real>::infinity()};
    }
    // Ray bisect
    for (int i = 0; i < BISECT_ITERATIONS; ++i) {
        dt *= 0.5;
        real half_way = t + dt;
        v = sp*cos(t) + tp*sin(t) +  c;
        real c = this->s_distance(v);
        if (a*c >= 0) {
            a = c;
            t = half_way;
        } else if (b*c >= 0) {
            b = c;
            t = half_way;
        } else {
            break;
        }
    }

    return {source*cos(t) + target*sin(t), target*cos(t) - source*sin(t), t};
}

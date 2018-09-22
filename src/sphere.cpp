#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <vector>
#include <memory>
#include <random>
#include <cassert>

#include "raytrace/quaternion.h"
#include "raytrace/sphere.h"

real Sphere::s_distance(quaternion location) {
    return norm(location) - 1;
}

quaternion Sphere::normal(quaternion location) {
    return location / norm(location);
}

quaternion Sphere::gradient(quaternion location) {
    return this->normal(location);
}

color Sphere::get_color(quaternion location, quaternion ray) {
    quaternion n = this->normal(location);
    real angle = -dot(n, ray);
    if (angle < 0.1) {
        angle = 0.1 - 0.1 * fabs(angle);
    }
    return (color) {0, fabs(angle), 0, 0};
}

std::tuple<quaternion, quaternion> Sphere::trace(quaternion source, quaternion target) {
    // Solve for t
    // |v|^2 = 1
    // v = source + t * (target - source)
    // ~ |u + t * w|^2 = 1
    // (u + t * w)*(u + t * w)' = 1
    // |u|^2 + t * (u*w' + w*u') + t^2*|w|^2 = 1
    // |w|^2*t^2 + 2*dot(u, w) * t + |u|^2 - 1 = 0
    // t = (-dot(u, w) +- sqrt(dot(u, w)^2+|w|^2*(1-|u|^2))) / |w|^2

    quaternion direction = target - source;
    direction = direction / norm(direction);

    source = this->inverse_project(source - this->location);
    target = this->inverse_project(target - this->location);

    quaternion u = source;
    quaternion w = target - source;
    w = w / norm(w);
    real b = dot(u, w);
    real d = b*b + 1.0 - dot(u, u);
    if (d < 0) {
        return {infinity(), direction};
    }
    d = sqrt(d);
    real t = -b + d;
    real t1 = -b - d;
    if (t1 < t && t1 > 0) {
        t = t1;
    }
    if (t < 0) {
        return {infinity(), direction};
    }
    return {this->project(u + t * w) + this->location, direction};
}

std::tuple<quaternion, quaternion, real> Sphere::trace_S3(quaternion source, quaternion target) {
    target = cross_align(source, target);

    quaternion sp = this->inverse_project(source);
    quaternion tp = this->inverse_project(target);
    quaternion c = this->inverse_project(-this->location);

    real s2 = dot(sp, sp);
    real st = 2*dot(sp, tp);
    real sc = 2*dot(sp, c);
    real t2 = dot(tp, tp);
    real tc = 2*dot(tp, c);
    real c2 = dot(c, c);
    // real f = c2 + cos_t * (s2 * cos_t + sc + st * sin_t) + sin_t * (tc + t2 * sin_t) - 1;
    // real df/dt = -sin_t*(sc + st * sin_t) + cos_t * (st * cos_t + tc - 2 * (s2 - t2) * sin_t);

    // real estimated_min_distance = c2 - s2 - fabs(sc) - fabs(st) - fabs(tc) - t2;
    // if (estimated_min_distance > 1) {
    //     return {infinity(), infinity(), std::numeric_limits<real>::infinity()};
    // }

    real t = 0;
    real cos_t = 1;
    real sin_t = 0;
    real squared_distance_to_center = c2 + s2 + sc;
    for (int i = 0; i < SPHERE_APPROACH_ITERATIONS; ++i) {
        real orbital_velocity = norm(-sin_t * sp + cos_t * tp);
        t += (sqrt(squared_distance_to_center) - 1) / orbital_velocity;
        cos_t = cos(t);
        sin_t = sin(t);
        squared_distance_to_center = c2 + cos_t * (s2 * cos_t + sc + st * sin_t) + sin_t * (tc + t2 * sin_t);
        if (squared_distance_to_center < 1 + SPHERE_EPSILON) {
            break;
        }
        if (t > 2*M_PI) {
            return {infinity(), infinity(), std::numeric_limits<real>::infinity()};
        }
    }
    if (squared_distance_to_center > 1 + SPHERE_EPSILON) {
        return {infinity(), infinity(), std::numeric_limits<real>::infinity()};
    }

    return {source*cos_t + target*sin_t, target*cos_t - source*sin_t, t};
}

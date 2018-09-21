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
    return RayTraceable::trace(source, target);
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

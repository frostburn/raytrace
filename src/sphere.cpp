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
    // Ray geodesic: source * cos(t) + target * sin(t)

    quaternion sp = this->inverse_project(source);
    quaternion tp = this->inverse_project(target);
    quaternion c = this->inverse_project(-this->location);

    // Signed distance estimate: s(source*cos(t) + target*sin(t))
    // ds / dt = dot(-source*sin(t) + target*cos(t), grad s)

    quaternion v;
    real t = 0;
    real cos_t = 1;
    real sin_t = 0;

    real t_closest = 0;
    real f_closest = std::numeric_limits<real>::infinity();
    for (int i = 0; i < SPHERE_NEWTON_WARMUP; ++i) {
        t = 2 * M_PI / (real) SPHERE_NEWTON_WARMUP;
        cos_t = cos(t);
        sin_t = sin(t);
        v = sp*cos_t + tp*sin_t + c;
        real f = this->s_distance(v);
        if (fabs(f) < f_closest) {
            f_closest = f;
            t_closest = t;
        }
    }
    t = t_closest;

    real min_gradient = std::numeric_limits<real>::infinity();
    real t_min = 0;

    // Newton's method
    for (int i = 0; i < SPHERE_NEWTON_ITERATIONS; ++i) {
        cos_t = cos(t);
        sin_t = sin(t);
        v = sp*cos_t + tp*sin_t + c;
        quaternion grad = this->gradient(v);
        real dt = dot(tp*cos_t - sp*sin_t, grad);
        if (fabs(dt) < min_gradient) {
            min_gradient = fabs(dt);
            t_min = t;
        }
        if (fabs(dt) < NEWTON_EPSILON) {
            dt = NEWTON_NUDGE;
        } else {
            real f = this->s_distance(v);
            dt = -f / dt;
        }
        if (dt > 0.5 * M_PI) {
            dt = 0.5 * M_PI;
        } else if (dt < -0.5 * M_PI) {
            dt = -0.5 * M_PI;
        }
        t += dt;
    }
    if (this->s_distance(v) > NEWTON_EPSILON) {
        return {infinity(), infinity(), std::numeric_limits<real>::infinity()};
    }

    // Gradient descent to find the midpoint
    for (int i = 0; i < SPHERE_NEWTON_DESCENDS; ++i) {
        cos_t = cos(t_min);
        sin_t = sin(t_min);
        v = sp*cos_t + tp*sin_t + c;
        quaternion grad = this->gradient(v);
        real dt = dot(tp*cos_t - sp*sin_t, grad);
        t_min -= 0.1 * dt;
    }
    // Relect the root around the midpoint to find another root.
    real t1 = 2*t_min - t;

    // A few more Newton's iterations to focus on the second root
    for (int i = 0; i < SPHERE_NEWTON_SECOND_ROOT; ++i) {
        cos_t = cos(t1);
        sin_t = sin(t1);
        v = sp*cos_t + tp*sin_t + c;
        quaternion grad = this->gradient(v);
        real dt = dot(tp*cos_t - sp*sin_t, grad);
        real f = this->s_distance(v);
        dt = -f / dt;
        t1 += dt;
    }

    // Snap to [0, tau] range
    t -= 2*M_PI * floor(0.5*t/M_PI);
    t1 -= 2*M_PI * floor(0.5*t1/M_PI);

    if (t1 < t) {
        t = t1;
    }

    cos_t = cos(t);
    sin_t = sin(t);
    return {source*cos_t + target*sin_t, target*cos_t - source*sin_t, t};
}

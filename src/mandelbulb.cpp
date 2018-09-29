#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <vector>
#include <memory>
#include <random>
#include <cassert>

#include "raytrace/quaternion.h"
#include "raytrace/power3d.h"
#include "raytrace/mandelbulb.h"

real Mandelbulb::s_distance(quaternion location) {
    quaternion z = location;
    int i;
    real magnitude;
    for (i = 0; i < this->max_iter; ++i) {
        z = pow3d(z, this->order) + location;
        // std::cerr << i << "," << z << std::endl;
        magnitude = z.x*z.x + z.y*z.y + z.z*z.z;
        if (magnitude > MANDELBULB_BAILOUT) {
            break;
        }
    }
    if (i >= this->max_iter - 1) {
        return -this->threshold;
    }
    real magnitude_norm = log(log(MANDELBULB_BAILOUT) * this->order) - log(log(MANDELBULB_BAILOUT));
    real v = (log(log(magnitude)) - log(log(MANDELBULB_BAILOUT))) / magnitude_norm;
    return v - i  + this->max_iter - 2 - this->threshold;
}

color Mandelbulb::get_color(quaternion location, quaternion ray) {
    // return this->pigment;
    quaternion n = this->normal(location);
    real angle = -dot(n, ray);
    if (angle < 0) {
        return (color) {-1, -angle, 0, 0};
    }
    return angle * this->pigment;
}

std::tuple<quaternion, quaternion> Mandelbulb::trace(quaternion source, quaternion target) {
    quaternion direction = target - source;
    direction = direction / norm(direction);
    source = this->inverse_project(source - this->location);
    quaternion dt = this->inverse_project(direction);
    dt = dt / norm(dt);

    real distance;
    for (int i = 0; i < MANDELBULB_APPROACH_ITERATIONS; ++i) {
        distance = this->s_distance(source);
        source = source + dt * distance * MANDELBULB_APPROACH_ALPHA;
        if (fabs(distance) <= MANDELBULB_APPROACH_EPSILON) {
            break;
        }
    }

    if (distance > MANDELBULB_APPROACH_EPSILON) {
        return {infinity(), direction};
    }

    return {this->project(source) + this->location, direction};
}

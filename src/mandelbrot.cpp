#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <vector>
#include <memory>
#include <random>
#include <cassert>

#include "raytrace/quaternion.h"
#include "raytrace/mandelbrot.h"

real Mandelbrot::s_distance(quaternion location) {
    quaternion z = location;
    real magnitude;
    int i;
    for (i = 0; i < this->max_iter; ++i) {
        quaternion zl = z;
        quaternion zr = z;
        for (int j = 1; j < this->left_order; ++j) {
            zl = zl * z;
        }
        for (int j = 1; j < this->right_order; ++j) {
            zr = zr * z;
        }
        z = zl * this->twist * zr + this->displacement;
        magnitude = z.r*z.r + z.x*z.x + z.y*z.y + z.z*z.z;
        if (magnitude > MANDELBROT_BAILOUT) {
            break;
        }
    }
    if (i >= this->max_iter - 1) {
        return -this->threshold;
    }
    real magnitude_norm = log(log(MANDELBROT_BAILOUT) * (this->left_order + this->right_order)) - log(log(MANDELBROT_BAILOUT));
    real v = (log(log(magnitude)) - log(log(MANDELBROT_BAILOUT))) / magnitude_norm;
    return v - i  + this->max_iter - 2 - this->threshold;
}

color Mandelbrot::get_color(quaternion location, quaternion ray) {
    quaternion n = this->normal(location);
    real angle = -dot(n, ray);
    if (angle < 0) {
        return (color) {-1, -angle, 0, 0};
    }
    return angle * this->pigment;
}

std::tuple<quaternion, quaternion> Mandelbrot::trace(quaternion source, quaternion target) {
    quaternion direction = target - source;
    direction = direction / norm(direction);
    source = this->inverse_project(source - this->location);
    quaternion dt = this->inverse_project(direction);
    dt = dt / norm(dt);

    real distance;
    for (int i = 0; i < MANDELBROT_APPROACH_ITERATIONS; ++i) {
        distance = this->s_distance(source);
        source = source + dt * distance * MANDELBROT_APPROACH_ALPHA;
        if (fabs(distance) < MANDELBROT_APPROACH_EPSILON) {
            break;
        }
    }

    if (distance > MANDELBROT_APPROACH_EPSILON) {
        return {infinity(), direction};
    }

    return {this->project(source) + this->location, direction};
}

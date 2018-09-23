#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include <cmath>

#include "raytrace/quaternion.h"
#include "raytrace/approach_traceable.h"

std::tuple<quaternion, quaternion> ApproachTraceable::trace(quaternion source, quaternion target) {
    quaternion direction = target - source;
    direction = direction / norm(direction);
    source = this->inverse_project(source - this->location);
    quaternion dt = this->inverse_project(direction);
    dt = dt / norm(dt);

    real distance;
    for (int i = 0; i < APPROACH_ITERATIONS; ++i) {
        distance = this->s_distance(source);
        source = source + dt * distance;
        if (fabs(distance) < APPROACH_EPSILON) {
            break;
        }
    }

    if (distance > APPROACH_EPSILON) {
        return {infinity(), direction};
    }

    return {this->project(source) + this->location, direction};
}

std::tuple<quaternion, quaternion, real> ApproachTraceable::trace_S3(quaternion source, quaternion target) {
    target = cross_align(source, target);

    quaternion sp = this->inverse_project(source);
    quaternion tp = this->inverse_project(target);
    quaternion c = this->inverse_project(-this->location);

    real t = 0;
    real cos_t = 1;
    real sin_t = 0;
    real distance = this->s_distance(sp);

    for (int i = 0; i < APPROACH_ITERATIONS; ++i) {
        real orbital_velocity = norm(-sin_t * sp + cos_t * tp);
        t += distance / orbital_velocity;
        if (t > 2*M_PI) {
            return {infinity(), infinity(), std::numeric_limits<real>::infinity()};
        }
        cos_t = cos(t);
        sin_t = sin(t);
        if (fabs(distance) < APPROACH_EPSILON) {
            break;
        }
        distance = this->s_distance(cos_t * sp + sin_t * tp + c);
    }
    if (distance > APPROACH_EPSILON) {
        return {infinity(), infinity(), std::numeric_limits<real>::infinity()};
    }

    return {source*cos_t + target*sin_t, target*cos_t - source*sin_t, t};
}

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

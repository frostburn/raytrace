#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <vector>
#include <memory>
#include <random>
#include <cassert>

#include "raytrace/quaternion.h"
#include "raytrace/gradient_traceable.h"

std::tuple<quaternion, quaternion, real> GradientTraceable::trace_S3(quaternion source, quaternion target) {
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

    // Newton's method
    for (int i = 0; i < NEWTON_ITERATIONS; ++i) {
        cos_t = cos(t);
        sin_t = sin(t);
        // v = inverse_project(source*cost_t + target*sin_t - this->location)
        v = sp*cos_t + tp*sin_t + c;
        quaternion grad = this->gradient(v);
        real dt = dot(tp*cos_t - sp*sin_t, grad);
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
        if (t + dt < 0) {
            dt = NEWTON_NUDGE;
        } else if (t + dt > 2 * M_PI) {
            dt = -NEWTON_NUDGE;
        }
        t += dt;
    }
    if (this->s_distance(v) > NEWTON_EPSILON) {
        return {infinity(), infinity(), std::numeric_limits<real>::infinity()};
    }

    return {source*cos_t + target*sin_t, target*cos_t - source*sin_t, t};
}

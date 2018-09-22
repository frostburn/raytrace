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

real Sphere::newton(real t, real cos_t, real sin_t, real f, real dt, quaternion sp, quaternion tp, quaternion c) {
    real initial_dt = dt;
    int i = 0;
    while (true) {
        if (fabs(dt) < NEWTON_EPSILON) {
            dt = NEWTON_NUDGE;
        } else {
            dt = -f / dt;
        }
        t += dt;
        i += 1;
        if (i >= SPHERE_NEWTON_ITERATIONS) {
            break;
        }
        // while (true) {
            cos_t = cos(t);
            sin_t = sin(t);
            quaternion v = sp*cos_t + tp*sin_t + c;
            f = this->s_distance(v);
            quaternion grad = this->gradient(v);
            dt = dot(tp*cos_t - sp*sin_t, grad);
        //     if (dt * initial_dt >= 0) {
        //         break;
        //     } else {
        //         if (initial_dt > 0) {
        //             t += NEWTON_NUDGE;
        //         } else {
        //             t -= NEWTON_NUDGE;
        //         }
        //     }
        // }
    }
    if (fabs(f) > NEWTON_EPSILON) {
        return std::numeric_limits<real>::infinity();
    }
    return t;
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


    real f_min_estimate = s2 + t2 + c2 + fabs(sc) + fabs(st) + fabs(tc) + 1;
    if (f_min_estimate < 0) {
        return {infinity(), infinity(), std::numeric_limits<real>::infinity()};
    }

    quaternion v;
    real t = 0;
    real cos_t = 1;
    real sin_t = 0;

    real a;
    real b = c2 + s2 + sc - 1;
    real f_min = b;
    real t_min = 0;
    for (int i = 1; i < SPHERE_RAY_MARCH_DIVISIONS; ++i) {
        t = i * 2*M_PI / (real) SPHERE_RAY_MARCH_DIVISIONS;
        sin_t = sin(t);
        cos_t = cos(t);
        a = b;
        b = c2 + cos_t * (s2 * cos_t + sc + st * sin_t) + sin_t * (tc + t2 * sin_t) - 1;
        if (b < f_min) {
            f_min = b;
            t_min = t;
        }
        if(a*b <= 0) {
            break;
        }
    }
    real t0 = t - 2*M_PI / (real) SPHERE_RAY_MARCH_DIVISIONS;
    if (a*b <= 0) {
        for (int i = 0; i < SPHERE_BIJECT_ITERATIONS; ++i) {
            real half_way = 0.5 * (t0 + t);
            real c = c2 + cos_t * (s2 * cos_t + sc + st * sin_t) + sin_t * (tc + t2 * sin_t) - 1;
            if (c < f_min) {
                f_min = c;
            }
            if (a*c >= 0) {
                a = c;
                t0 = half_way;
            } else if (b*c >= 0) {
                b = c;
                t = half_way;
            } else {
                break;
            }
        }
        t = 0.5 * (t + t0);
    } else {
        t = t_min;
    }

    for (int i = 0; i < SPHERE_NEWTON_ITERATIONS; ++i) {
        real f = c2 + cos_t * (s2 * cos_t + sc + st * sin_t) + sin_t * (tc + t2 * sin_t) - 1;
        real dt = -sin_t*(sc + st * sin_t) + cos_t * (st * cos_t + tc - 2 * (s2 - t2) * sin_t);
        if (f < f_min) {
            f_min = f;
        }
        t -= f/dt;
    }

    if (f_min > 0) {
        return {infinity(), infinity(), std::numeric_limits<real>::infinity()};
    }

    cos_t = cos(t);
    sin_t = sin(t);
    return {source*cos_t + target*sin_t, target*cos_t - source*sin_t, t};
}

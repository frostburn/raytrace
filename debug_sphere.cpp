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
#include "raytrace/gradient_traceable.h"
#include "raytrace/sphere.h"
#include "raytrace/clifford_torus.h"
#include "raytrace/trace.h"

using namespace std;

const int SPHERE_NEWTON_ITERATIONS = 10;

real sphere_newton(Sphere *sphere, real t, quaternion sp, quaternion tp, quaternion c) {
    // TODO: Reuse f and grad from warmup
    real cos_t = 1;
    real sin_t = 0;
    // Newton's method
    cout << "xy = [";
    for (int i = 0; i < SPHERE_NEWTON_ITERATIONS; ++i) {
        cos_t = cos(t);
        sin_t = sin(t);
        quaternion v = sp*cos_t + tp*sin_t + c;
        quaternion grad = sphere->gradient(v);
        real f = sphere->s_distance(v);
        cout << "[" << t << "," << f << "]";
        if (i < SPHERE_NEWTON_ITERATIONS - 1) {
            cout << ",";
        }
        real dt = dot(tp*cos_t - sp*sin_t, grad);
        if (fabs(dt) < NEWTON_EPSILON) {
            dt = NEWTON_NUDGE;
        } else {
            dt = -f / dt;
        }
        if (dt > 0.5 * M_PI) {
            dt = 0.5 * M_PI;
        } else if (dt < -0.5 * M_PI) {
            dt = -0.5 * M_PI;
        }
        t += dt;
    }
    cout << "]" << endl;
    return t;
}



int main()
{
    default_random_engine generator;
    generator.seed(random_device()());
    normal_distribution<real> distribution(0, 1);

    quaternion source = {
        distribution(generator),
        distribution(generator),
        distribution(generator),
        distribution(generator)
    };
    quaternion target = {
        distribution(generator),
        distribution(generator),
        distribution(generator),
        distribution(generator)
    };
    source = source / norm(source);
    target = target / norm(target);

    Sphere *sphere = new Sphere;
    sphere->location = {
        distribution(generator),
        distribution(generator),
        distribution(generator),
        distribution(generator)
    };
    sphere->scale = {
        1 + 0.2 * distribution(generator),
        1 + 0.2 * distribution(generator),
        1 + 0.2 * distribution(generator),
        1 + 0.2 * distribution(generator)
    };
    target = cross_align(source, target);

    quaternion sp = sphere->inverse_project(source);
    quaternion tp = sphere->inverse_project(target);
    quaternion c = sphere->inverse_project(-sphere->location);

    real s2 = dot(sp, sp);
    real st = 2*dot(sp, tp);
    real sc = 2*dot(sp, c);
    real t2 = dot(tp, tp);
    real tc = 2*dot(tp, c);
    real c2 = dot(c, c);

    quaternion v;
    real t = 0;
    real cos_t = 1;
    real sin_t = 0;

    real t_plus = 0;
    real t_minus = 0;
    real f_plus = std::numeric_limits<real>::infinity();
    real f_minus = std::numeric_limits<real>::infinity();
    cout << "f = [";
    for (int i = 0; i < 100; ++i) {
        t = i* 2*M_PI / (real) 100;
        cos_t = cos(t);
        sin_t = sin(t);
        v = sp*cos_t + tp*sin_t + c;
        real f = sphere->s_distance(v);
        f = (f+1) * (f+1) - 1;
        // cout << f;
        // real ff = sin_t*sin_t*t2 + sin_t*cos_t*st + sin_t*tc + cos_t*cos_t*s2 + cos_t*sc + c2 - 1;
        real ff = cos_t*cos_t*s2 + cos_t*sin_t*st + cos_t*sc + sin_t*sin_t*t2 + sin_t*tc + c2 - 1;
        real ffs = -1 + c2 + cos_t * (s2 * cos_t + sc + st * sin_t) + sin_t * (tc + t2 * sin_t);
        // real ff = dot(v, v) - 1;
        // cout << (ff - ffs);
        f = fabs(f);
        // quaternion grad = sphere->gradient(v);
        quaternion grad = 2*v;
        real dt = dot(tp*cos_t - sp*sin_t, grad);
        real dtt = -sin_t*(sc + st * sin_t) + cos_t * (st * cos_t + tc - 2 * (s2 - t2) * sin_t);
        cout << (dt - dtt);
        if (dt >= 0 && f < f_plus) {
            f_plus = f;
            t_plus = t;
        }
        if (dt < 0 && f < f_minus) {
            f_minus = f;
            t_minus = t;
        }
        if (i < 100 - 1) {
            cout << ",";
        }
    }
    cout << "]" << endl;

    sphere_newton(sphere, t_minus, sp, tp, c);
    sphere_newton(sphere, t_plus, sp, tp, c);

    return EXIT_SUCCESS;
}


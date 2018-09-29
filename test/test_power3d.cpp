#include <cassert>
#include <cmath>
#include <limits>
#include <random>

#include "raytrace/power3d.h"

const real EPSILON = 1e-12;

void test_identity() {
    std::default_random_engine generator;
    generator.seed(std::random_device()());
    std::normal_distribution<real> distribution(0.0, 1);
    quaternion v = {
        0,
        distribution(generator),
        distribution(generator),
        distribution(generator)
    };
    quaternion u = pow3d(v, 1);
    assert(fabs(u.x - v.x) < EPSILON);
    assert(fabs(u.y - v.y) < EPSILON);
    assert(fabs(u.z - v.z) < EPSILON);
}

void test_zero() {
    quaternion v = {0, 0, 0, 0};
    for (int i = 1; i < 10; ++i) {
        quaternion u = pow3d(v, i);
        assert(fabs(u.x) < EPSILON);
        assert(fabs(u.y) < EPSILON);
        assert(fabs(u.z) < EPSILON);
    }
}

void test_known() {
    quaternion v = {0, 1, 0, 0};
    quaternion u = pow3d(v, 2);
    assert(u.x == 1 && u.y == 0 && u.z == 0);
}

void test_on_axis() {
    std::default_random_engine generator;
    generator.seed(std::random_device()());
    std::normal_distribution<real> distribution(0.0, 1);
    quaternion v;
    for (int i = 0; i < 10; ++i) {
        v = {0, distribution(generator), 0, 0};
        for (int n = 1; n < 10; ++n) {
            assert(fabs(pow3d(v, n).x - pow(v.x, n)) < EPSILON);
        }
        v = {0, 0, distribution(generator), 0};
        for (int n = 1; n < 10; ++n) {
            if (n % 4 == 0) {
                assert(fabs(pow3d(v, n).x - pow(v.y, n)) < EPSILON);
            } else if(n % 4 == 1) {
                assert(fabs(pow3d(v, n).y - pow(v.y, n)) < EPSILON);
            } else if(n % 4 == 2) {
                assert(fabs(pow3d(v, n).x + pow(v.y, n)) < EPSILON);
            } else if(n % 4 == 3) {
                assert(fabs(pow3d(v, n).y + pow(v.y, n)) < EPSILON);
            }
        }
    }
}

void test_gradient() {
    std::default_random_engine generator;
    generator.seed(std::random_device()());
    std::normal_distribution<real> distribution(0.0, 1);
    quaternion v;
    quaternion h;
    quaternion d;
    for (int i = 0; i < 10; ++i) {
        v = {0, distribution(generator), distribution(generator), distribution(generator)};
        for (int n = 1; n < 10; ++n) {
            // auto [a, x, y, z] = pow3d_grad(v, n);
            quaternion *a = pow3d_grad(v, n);

            h = {0, EPSILON, 0, 0};
            d = (pow3d(v + h, n) - pow3d(v - h, n)) * 0.5 / EPSILON;
            assert(norm(a[1] - d) < norm(d) * 1e-2);

            h = {0, 0, EPSILON, 0};
            d = (pow3d(v + h, n) - pow3d(v - h, n)) * 0.5 / EPSILON;
            assert(norm(a[2] - d) < norm(d) * 1e-2);

            h = {0, 0, 0, EPSILON};
            d = (pow3d(v + h, n) - pow3d(v - h, n)) * 0.5 / EPSILON;
            assert(norm(a[3] - d) < norm(d) * 1e-2);
        }
    }
}

int main() {
    test_identity();
    test_zero();
    test_known();
    test_on_axis();
    test_gradient();
    return 0;
}

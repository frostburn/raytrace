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

const real ACCURACY = 0.01;

void test_trace() {
    quaternion source = {0, 0, 0, -2};
    quaternion target = {0, 0, 0, 0};
    Sphere *sphere = new Sphere;
    auto [surface, direction] = sphere->trace(source, target);
    assert(norm(surface - (quaternion){0, 0, 0, -1}) < ACCURACY);
}

void test_scale() {
    quaternion source = {0, 0, 3, 0};
    quaternion target = {0, 0, 0, 0};
    Sphere *sphere = new Sphere;
    sphere->scale = sphere->scale * 0.5;
    auto [surface, direction] = sphere->trace(source, target);
    assert(norm(surface - (quaternion){0, 0, 0.5, 0}) < ACCURACY);
}

void test_location() {
    quaternion source = {0, 1, 0.2, 0};
    quaternion target = {0, 0, 0, 0.1};
    Sphere *sphere = new Sphere;
    sphere->location = {0, 0, 0.2, 0.1};
    sphere->scale = sphere->scale * 0.5;
    auto [surface, direction] = sphere->trace(source, target);
    assert(norm(surface - (quaternion){0, 0.486928, 0.0973857, 0.0513072}) < ACCURACY);
    real radius = norm(surface - sphere->location);
    assert(radius > 0.5 - ACCURACY);
    assert(radius < 0.5 + ACCURACY);
}

void test_trace_S3() {
    default_random_engine generator;
    generator.seed(random_device()());
    normal_distribution<real> distribution(0.0, 1);

    quaternion source = {-1, 0, 0, 0};
    quaternion target = {0, 1, 0, 0};
    Sphere *sphere = new Sphere;
    sphere->location = (quaternion){1, 0, 0, 0};
    auto [surface, direction, distance] = sphere->trace_S3(source, target);
    real radius = norm(surface - sphere->location);
    assert(radius < 1 + ACCURACY);
    assert(radius > 1 - ACCURACY);
    assert(distance < M_PI);
    assert(direction.r > 0);
    assert(direction.x < 0);
}

void test_trace_S3_random() {
    default_random_engine generator;
    generator.seed(random_device()());
    normal_distribution<real> distribution(0.0, 1);

    Sphere *sphere = new Sphere;
    quaternion scale = sphere->scale;
    for (int i = 0; i < 100; ++i) {
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
        sphere->location = {
            distribution(generator),
            distribution(generator),
            distribution(generator),
            distribution(generator)
        };
        source = source / norm(source);
        target = target / norm(target);
        sphere->scale = scale * fabs(distribution(generator));

        auto [surface, direction, distance] = sphere->trace_S3(source, target);
        if (norm(target - sphere->location) < sphere->scale.r - ACCURACY) {
            quaternion v = sphere->inverse_project(target) + sphere->inverse_project(-sphere->location);
            assert(norm(v) < 1 + ACCURACY);

            assert(distance <= 2*M_PI);
            real radius = norm(surface - sphere->location);
            assert(radius < sphere->scale.r + ACCURACY);
            assert(radius > sphere->scale.r - ACCURACY);
        }
    }
};

int main() {
    test_trace();
    test_scale();
    test_location();
    test_trace_S3();
    // test_trace_S3_random();  // TODO: FIX
    return 0;
}

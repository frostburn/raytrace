#include <cassert>
#include <cmath>
#include <limits>
#include <random>

#include "raytrace/quaternion.h"

void test_great_circle() {
    std::default_random_engine generator;
    generator.seed(std::random_device()());
    std::normal_distribution<real> distribution(0.0, 1);
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

    quaternion omega = cross_align(source, target);

    real min_distance = std::numeric_limits<real>::infinity();
    real closest_theta = 0;
    for (int i = 0; i < 100; ++i) {
        real theta = M_PI * i * 0.02;
        quaternion v = cos(theta) * source + sin(theta) * omega;
        if (norm(v - target) < min_distance) {
            min_distance = norm(v - target);
            closest_theta = theta;
        }
    }
    assert(closest_theta < M_PI);
    assert(min_distance < 0.1);
}

int main() {
    test_great_circle();
    return 0;
}

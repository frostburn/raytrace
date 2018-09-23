#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <vector>
#include <memory>
#include <random>
#include <cassert>

#include "raytrace/quaternion.h"
#include "raytrace/clifford_torus.h"

real CliffordTorus::s_distance(quaternion location) {
    real rx = sqrt(location.r * location.r + location.x * location.x);
    real yz = sqrt(location.y * location.y + location.z * location.z);
    if (rx > yz) {
        return rx - 1;
    }
    return yz - 1;
}

quaternion CliffordTorus::normal(quaternion location) {
    real rx = sqrt(location.r * location.r + location.x * location.x);
    real yz = sqrt(location.y * location.y + location.z * location.z);
    if (rx > yz) {
        return (quaternion){location.r / rx, location.x / rx, 0, 0};
    }
    return (quaternion){0, 0, location.y / yz, location.z / yz};
}

quaternion CliffordTorus::gradient(quaternion location) {
    return this->normal(location);
}

color CliffordTorus::get_color(quaternion location, quaternion ray) {
    quaternion n = this->normal(location);
    real angle = -dot(n, ray);
    if (angle < 0) {
        return (color) {-1, -angle, 0, 0};
    }
    n = {location.y, location.z, 0, 0};
    n = n*n;
    n = n*n;
    real r2 = n.r*n.r;
    real x2 = n.x*n.x;
    real y = 1 - r2 - x2;
    return angle * (r2 * this->pigment_a + x2 * this->pigment_b + y * this->pigment_c);
}

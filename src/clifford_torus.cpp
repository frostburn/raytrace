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
    if (angle < 0.1) {
        angle = 0.1 - 0.1 * fabs(angle);
    }
    real pattern = cos(10*location.r)*cos(10*location.x)*cos(10*location.y)*cos(10*location.z);
    return (color) {0, (0.5 + 0.4 * pattern) * angle, 0, 0};
}

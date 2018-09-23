#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <vector>
#include <memory>
#include <random>
#include <cassert>

#include "raytrace/quaternion.h"
#include "raytrace/hyper_torus.h"

real HyperTorus::s_distance(quaternion location) {
    real rxy = 1 - sqrt(location.r * location.r + location.x * location.x + location.y * location.y);
    return sqrt(rxy*rxy + location.z*location.z) - this->minor_radius;
}

quaternion HyperTorus::gradient(quaternion location) {
    real rxy = sqrt(location.r * location.r + location.x * location.x + location.y * location.y);
    real r = sqrt((1-rxy)*(1-rxy) + location.z*location.z);
    quaternion n = {location.r, location.x, location.y, 0};
    n = n * (rxy - 1) / rxy;
    n.z = location.z;
    return n / r;
}

quaternion HyperTorus::normal(quaternion location) {
    quaternion n = this->gradient(location);
    return n / norm(n);
}

color HyperTorus::get_color(quaternion location, quaternion ray) {
    quaternion n = this->normal(location);
    real angle = -dot(n, ray);
    if (angle < 0) {
        return (color) {-1, -angle, 0, 0};
    }
    n = {location.x, location.y, 0, 0};
    n = n*n;
    n = n*n;
    n = n*n;
    if (n.r < 0) {
        return angle * this->pigment_b;
    }
    return angle * this->pigment_a;
}

real HyperTorus2::s_distance(quaternion location) {
    real xy = 1 - sqrt(location.x * location.x + location.y * location.y);
    return sqrt(xy*xy + location.z*location.z + location.r*location.r) - this->minor_radius;
}

quaternion HyperTorus2::gradient(quaternion location) {
    real xy = sqrt(location.x * location.x + location.y * location.y);
    real r = sqrt((1-xy)*(1-xy) + location.z*location.z + location.r*location.r);
    quaternion n = {0, location.x, location.y, 0};
    n = n * (xy - 1) / xy;
    n.r = location.r;
    n.z = location.z;
    return n / r;
}

quaternion HyperTorus2::normal(quaternion location) {
    quaternion n = this->gradient(location);
    return n / norm(n);
}

color HyperTorus2::get_color(quaternion location, quaternion ray) {
    quaternion n = this->normal(location);
    real angle = -dot(n, ray);
    if (angle < 0) {
        return (color) {-1, -angle, 0, 0};
    }
    n = {location.x, location.y, 0, 0};
    n = n*n;
    n = n*n;
    n = n*n;
    if (n.r < 0) {
        return angle * this->pigment_b;
    }
    return angle * this->pigment_a;
}

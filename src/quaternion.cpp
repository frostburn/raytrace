#include <iostream>
#include <limits>
#include <cmath>

#include "raytrace/quaternion.h"

const real EPSILON = 1e-12;

std::ostream& operator<<(std::ostream& os, const quaternion& v)
{
    return os << "(quaternion){" <<  v.r << ", " << v.x << ", " << v.y << ", " << v.z << "}"; 
}

real norm(quaternion v) {
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z + v.r*v.r);
}

quaternion conjugate(quaternion v) {
    return (quaternion){v.r, -v.x, -v.y, -v.z};
}

quaternion operator -(const quaternion& v) {
    return (quaternion){-v.r, -v.x, -v.y, -v.z};
}

quaternion operator *(const quaternion& a, const real b) {
    return (quaternion){a.r * b, a.x * b, a.y * b, a.z * b};
}

quaternion operator *(const real b, const quaternion& a) {
    return (quaternion){a.r * b, a.x * b, a.y * b, a.z * b};
}

quaternion operator /(const quaternion& a, const real b) {
    return a * (1.0 / b);
}

quaternion operator /(const real a, const quaternion& b) {
    real r = b.x*b.x + b.y*b.y + b.z*b.z + b.r*b.r;
    return (a / r) * conjugate(b);
}

quaternion operator +(const quaternion& a, const quaternion& b) {
    return (quaternion){a.r + b.r, a.x + b.x, a.y + b.y, a.z + b.z};
}

quaternion operator +(const real& a, const quaternion& b) {
    return (quaternion){a + b.r, b.x, b.y, b.z};
}

quaternion operator -(const quaternion& a, const quaternion& b) {
    return (quaternion){a.r - b.r, a.x - b.x, a.y - b.y, a.z - b.z};
}

quaternion operator *(const quaternion& a, const quaternion& b) {
    return (quaternion){
        a.r*b.r - a.x*b.x - a.y*b.y - a.z*b.z,
        a.r*b.x + a.x*b.r + a.y*b.z - a.z*b.y,
        a.r*b.y - a.x*b.z + a.y*b.r + a.z*b.x,
        a.r*b.z + a.x*b.y - a.y*b.x + a.z*b.r
    };
}

quaternion operator /(const quaternion& a, const quaternion& b) {
    real r = b.x*b.x + b.y*b.y + b.z*b.z + b.r*b.r;
    return a * conjugate(b) / r;
}

quaternion infinity () {
    real inf = std::numeric_limits<real>::infinity();
    return (quaternion){inf, inf, inf, inf};
}

real dot(const quaternion& a, const quaternion& b) {
    return a.r * b.r + a.x * b.x + a.y * b.y + a.z * b.z;
}

quaternion project(const quaternion& a, const quaternion& b) {
    return dot(a, b) * a;
}

quaternion cross_align(const quaternion& a, const quaternion& b) {
    // Assumes |a| = 1
    quaternion c = b;
    real angle = dot(a, b);
    if (angle > EPSILON) {
        c = b / angle - a;
    } else if (angle < -EPSILON) {
        c = a - b / angle;
    }
    return c / norm(c);
}

color clip_color(const color& a) {
    color pixel = a;
    if (pixel.r > 1) {
        pixel.r = 1;
    }
    if (pixel.r < 0) {
        pixel.r = 0;
    }
    if (pixel.x > 1) {
        pixel.x = 1;
    }
    if (pixel.x < 0) {
        pixel.x = 0;
    }
    if (pixel.y > 1) {
        pixel.y = 1;
    }
    if (pixel.y < 0) {
        pixel.y = 0;
    }
    if (pixel.z > 1) {
        pixel.z = 1;
    }
    if (pixel.z < 0) {
        pixel.z = 0;
    }
    return pixel;
}
